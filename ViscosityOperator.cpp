#include <ViscosityOperator.h>
#include <GradientCalculatorFD3.h>
#include <Vector5D.h>
#include <Utils.h>
#include <bits/stdc++.h>
using std::unique_ptr;

extern int verbose;
extern int INACTIVE_MATERIAL_ID;


//--------------------------------------------------------------------------

ViscosityOperator::ViscosityOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, 
                                     EquationsData &iod_eqs_, vector<VarFcnBase*>& varFcn_,
                                     SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_,
                                     InterpolatorBase &interpolator, GradientCalculatorBase &grad_)
                  : iod_eqs(iod_eqs_), coordinates(coordinates_), delta_xyz(delta_xyz_),
                    varFcn(varFcn_),
                    V_i_minus_half(comm_, &(dm_all_.ghosted1_3dof)),
                    V_j_minus_half(comm_, &(dm_all_.ghosted1_3dof)),
                    V_k_minus_half(comm_, &(dm_all_.ghosted1_3dof)),
                    dVdx_i_minus_half(comm_, &(dm_all_.ghosted1_3dof)),
                    dVdx_j_minus_half(comm_, &(dm_all_.ghosted1_3dof)),
                    dVdx_k_minus_half(comm_, &(dm_all_.ghosted1_3dof)),
                    dVdy_i_minus_half(comm_, &(dm_all_.ghosted1_3dof)),
                    dVdy_j_minus_half(comm_, &(dm_all_.ghosted1_3dof)),
                    dVdy_k_minus_half(comm_, &(dm_all_.ghosted1_3dof)),
                    dVdz_i_minus_half(comm_, &(dm_all_.ghosted1_3dof)),
                    dVdz_j_minus_half(comm_, &(dm_all_.ghosted1_3dof)),
                    dVdz_k_minus_half(comm_, &(dm_all_.ghosted1_3dof)),
                    interpolator(interpolator), grad(grad_),
                    cylindrical_symmetry(false),
                    dudx(NULL), dvdx(NULL), dudy(NULL), dvdy(NULL), 
                    Lam(NULL), Mu(NULL), dLamdx(NULL), dLamdy(NULL),
                    grad_minus(NULL), grad_plus(NULL)
{

  // Get i0, j0, etc.
  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);

  // Set internal variables to 0 (including ghost layer)
  V_i_minus_half.SetConstantValue(0.0, true);
  V_j_minus_half.SetConstantValue(0.0, true);
  V_k_minus_half.SetConstantValue(0.0, true);
  dVdx_i_minus_half.SetConstantValue(0.0, true);
  dVdx_j_minus_half.SetConstantValue(0.0, true);
  dVdx_k_minus_half.SetConstantValue(0.0, true);
  dVdy_i_minus_half.SetConstantValue(0.0, true);
  dVdy_j_minus_half.SetConstantValue(0.0, true);
  dVdy_k_minus_half.SetConstantValue(0.0, true);
  dVdz_i_minus_half.SetConstantValue(0.0, true);
  dVdz_j_minus_half.SetConstantValue(0.0, true);
  dVdz_k_minus_half.SetConstantValue(0.0, true);

  // Create a ViscoFcn for each material.
  // For problems involving multiple materials, it could happen that the user specified viscosity
  // model for only some of the materials, but not all. In this case, we still create a dummy 
  // (i.e. base) ViscoFcn for those "inviscid" materials. If all the materials are inviscid,
  // This constructor should not be called in the first place.
  //
  for(auto it = iod_eqs.materials.dataMap.begin(); it != iod_eqs.materials.dataMap.end(); it++)
    visFcn.push_back(NULL); //allocate space

  for(auto it = iod_eqs.materials.dataMap.begin(); it != iod_eqs.materials.dataMap.end(); it++) {
    int matid = it->first;
    if(matid < 0 || matid >= (int)visFcn.size()) {
      print_error("*** Error: Detected error in the specification of material indices (id = %d).\n", matid);
      exit_mpi();
    }
    switch (it->second->viscosity.type) {
      case ViscosityModelData::NONE :
        visFcn[matid] = new ViscoFcnBase(*varFcn[matid]);
        break;
      case ViscosityModelData::CONSTANT :
        visFcn[matid] = new ViscoFcnConstant(it->second->viscosity, *varFcn[matid]);
        break;
      case ViscosityModelData::SUTHERLAND :
        visFcn[matid] = new ViscoFcnSutherland(it->second->viscosity, *varFcn[matid]);
        break;
      case ViscosityModelData::ARTIFICIAL_RODIONOV :
        visFcn[matid] = new ViscoFcnRodionov(it->second->viscosity, *varFcn[matid]);
        break;
    }
  }


  // Initialize variables for cylindrical symmetry
  if(iod.mesh.type == MeshData::CYLINDRICAL) {

    cylindrical_symmetry = true;

    dudx = new SpaceVariable3D(comm_, &(dm_all_.ghosted1_1dof));
    dvdx = new SpaceVariable3D(comm_, &(dm_all_.ghosted1_1dof));
    dudy = new SpaceVariable3D(comm_, &(dm_all_.ghosted1_1dof));
    dvdy = new SpaceVariable3D(comm_, &(dm_all_.ghosted1_1dof));
    Lam = new SpaceVariable3D(comm_, &(dm_all_.ghosted1_1dof));
    Mu  = new SpaceVariable3D(comm_, &(dm_all_.ghosted1_1dof));
    dLamdx = new SpaceVariable3D(comm_, &(dm_all_.ghosted1_1dof));
    dLamdy = new SpaceVariable3D(comm_, &(dm_all_.ghosted1_1dof));

    double NX,NY,NZ;
    dudx->GetGlobalSize(&NX,&NY,&NZ);
    if(NZ != 1)
      print_warning("Warning: Enforcing cylindrical symmetry, but mesh has %d layers "
                    "in z-direction.\n", NZ); 

    grad_minus = new GradientCalculatorFD3(comm_, dm_all_, coordinates, delta_xyz, -1);
    grad_plus  = new GradientCalculatorFD3(comm_, dm_all_, coordinates, delta_xyz,  1);

  }
  else if(iod.mesh.type == MeshData::SPHERICAL) {
    print_error("*** Error: Currently, ViscosityOperator does not support spherical symmetry.\n");
    exit(-1);
  }

}

//--------------------------------------------------------------------------

ViscosityOperator::~ViscosityOperator()
{
  for(int i=0; i<(int)visFcn.size(); i++)
    delete visFcn[i];

  if(dudx) delete dudx;
  if(dvdx) delete dvdx;
  if(dudy) delete dudy;
  if(dvdy) delete dvdy;
  if(Lam)  delete Lam;
  if(Mu)   delete Mu;
  if(dLamdx) delete dLamdx;
  if(dLamdy) delete dLamdy;

  if(grad_minus) delete grad_minus;
  if(grad_plus) delete grad_plus;
}

//--------------------------------------------------------------------------

void
ViscosityOperator::Destroy()
{

  V_i_minus_half.Destroy();
  V_j_minus_half.Destroy();
  V_k_minus_half.Destroy();

  dVdx_i_minus_half.Destroy();
  dVdx_j_minus_half.Destroy();
  dVdx_k_minus_half.Destroy();
  dVdy_i_minus_half.Destroy();
  dVdy_j_minus_half.Destroy();
  dVdy_k_minus_half.Destroy();
  dVdz_i_minus_half.Destroy();
  dVdz_j_minus_half.Destroy();
  dVdz_k_minus_half.Destroy();


  // members for cylindrical symmetry

  if(dudx) dudx->Destroy();
  if(dvdx) dvdx->Destroy();
  if(dudy) dudy->Destroy();
  if(dvdy) dvdy->Destroy();

  if(Lam) Lam->Destroy();
  if(Mu)  Mu->Destroy();
  if(dLamdx) dLamdx->Destroy();
  if(dLamdy) dLamdy->Destroy();

  if(grad_minus)  grad_minus->Destroy();
  if(grad_plus)   grad_plus->Destroy();

}

//--------------------------------------------------------------------------
// Add diffusion fluxes on the left hand side of the N-S equations
void ViscosityOperator::AddDiffusionFluxes(SpaceVariable3D &V, SpaceVariable3D &ID, 
                                           vector<unique_ptr<EmbeddedBoundaryDataSet> > *EBDS,
                                           SpaceVariable3D &R)
{

  if(EBDS) 
    print_warning("Warning: ViscosityOperator::AddDiffusionFluxes: Not able to account for embedded surfaces.\n");


  std::vector<int> i123{1,2,3}, i012{0,1,2}; 
  //1. Calculate the x, y, and z velocities at cell interfaces by interpolation 
  interpolator.InterpolateAtCellInterfaces(0/*x-dir*/, V, i123, V_i_minus_half, i012);
  interpolator.InterpolateAtCellInterfaces(1/*y-dir*/, V, i123, V_j_minus_half, i012);
  interpolator.InterpolateAtCellInterfaces(2/*z-dir*/, V, i123, V_k_minus_half, i012);

  //2. Calculate velocity derivatives at cell interfaces
  grad.CalculateFirstDerivativeAtCellInterfaces(0, 0, V, i123, dVdx_i_minus_half, i012);
  grad.CalculateFirstDerivativeAtCellInterfaces(0, 1, V, i123, dVdx_j_minus_half, i012);
  grad.CalculateFirstDerivativeAtCellInterfaces(0, 2, V, i123, dVdx_k_minus_half, i012);
  grad.CalculateFirstDerivativeAtCellInterfaces(1, 0, V, i123, dVdy_i_minus_half, i012);
  grad.CalculateFirstDerivativeAtCellInterfaces(1, 1, V, i123, dVdy_j_minus_half, i012);
  grad.CalculateFirstDerivativeAtCellInterfaces(1, 2, V, i123, dVdy_k_minus_half, i012);
  grad.CalculateFirstDerivativeAtCellInterfaces(2, 0, V, i123, dVdz_i_minus_half, i012);
  grad.CalculateFirstDerivativeAtCellInterfaces(2, 1, V, i123, dVdz_j_minus_half, i012);
  grad.CalculateFirstDerivativeAtCellInterfaces(2, 2, V, i123, dVdz_k_minus_half, i012);

  //3. Loop through cell interfaces and calculate viscous fluxes
  Vec3D*** vi = (Vec3D***)V_i_minus_half.GetDataPointer();
  Vec3D*** vj = (Vec3D***)V_j_minus_half.GetDataPointer();
  Vec3D*** vk = (Vec3D***)V_k_minus_half.GetDataPointer();
  Vec3D*** dvdx_i = (Vec3D***)dVdx_i_minus_half.GetDataPointer();
  Vec3D*** dvdx_j = (Vec3D***)dVdx_j_minus_half.GetDataPointer();
  Vec3D*** dvdx_k = (Vec3D***)dVdx_k_minus_half.GetDataPointer();
  Vec3D*** dvdy_i = (Vec3D***)dVdy_i_minus_half.GetDataPointer();
  Vec3D*** dvdy_j = (Vec3D***)dVdy_j_minus_half.GetDataPointer();
  Vec3D*** dvdy_k = (Vec3D***)dVdy_k_minus_half.GetDataPointer();
  Vec3D*** dvdz_i = (Vec3D***)dVdz_i_minus_half.GetDataPointer();
  Vec3D*** dvdz_j = (Vec3D***)dVdz_j_minus_half.GetDataPointer();
  Vec3D*** dvdz_k = (Vec3D***)dVdz_k_minus_half.GetDataPointer();

  Vec3D*** dxyz = (Vec3D***)delta_xyz.GetDataPointer();
  Vec5D*** res  = (Vec5D***)R.GetDataPointer();
  double*** id  = (double***)ID.GetDataPointer();

  int myid = 0;
  double dx = 0.0, dy = 0.0, dz = 0.0;
  Vec5D flux;
  for(int k=k0; k<kkmax; k++)
    for(int j=j0; j<jjmax; j++)
      for(int i=i0; i<iimax; i++) {

        myid = id[k][j][i];
        dx   = dxyz[k][j][i][0];
        dy   = dxyz[k][j][i][1];
        dz   = dxyz[k][j][i][2];
         
        //*****************************************
        //calculate flux function F_{i-1/2,j,k}
        //*****************************************
        if(k!=kkmax-1 && j!=jjmax-1) {        
          visFcn[myid]->EvaluateViscousFluxFunction_F(flux, dvdx_i[k][j][i], dvdy_i[k][j][i],
                                                        dvdz_i[k][j][i], vi[k][j][i], &dx);
          flux *= dy*dz;
          res[k][j][i]   += flux;
          res[k][j][i-1] -= flux;
        }

        //*****************************************
        //calculate flux function G_{i,j-1/2,k}
        //*****************************************
        if(i!=iimax-1 && k!=kkmax-1) {        
          visFcn[myid]->EvaluateViscousFluxFunction_G(flux, dvdx_j[k][j][i], dvdy_j[k][j][i],
                                                        dvdz_j[k][j][i], vj[k][j][i], &dy);
          flux *= dx*dz;
          res[k][j][i]   += flux;
          res[k][j-1][i] -= flux;
        }

        //*****************************************
        //calculate flux function H_{i,j,k-1/2}
        //*****************************************
        if(i!=iimax-1 && j!=jjmax-1) {        
          visFcn[myid]->EvaluateViscousFluxFunction_H(flux, dvdx_k[k][j][i], dvdy_k[k][j][i],
                                                        dvdz_k[k][j][i], vk[k][j][i], &dz);
          flux *= dx*dy;
          res[k][j][i]   += flux;
          res[k-1][j][i] -= flux;
        }

      }


  if(cylindrical_symmetry) //Directly use id, dxyz, and res (local vars) 
    AddCylindricalSymmetryTerms(V, id, dxyz, EBDS, res);


  V_i_minus_half.RestoreDataPointerToLocalVector();
  V_j_minus_half.RestoreDataPointerToLocalVector();
  V_k_minus_half.RestoreDataPointerToLocalVector();
  dVdx_i_minus_half.RestoreDataPointerToLocalVector();
  dVdx_j_minus_half.RestoreDataPointerToLocalVector();
  dVdx_k_minus_half.RestoreDataPointerToLocalVector();
  dVdy_i_minus_half.RestoreDataPointerToLocalVector();
  dVdy_j_minus_half.RestoreDataPointerToLocalVector();
  dVdy_k_minus_half.RestoreDataPointerToLocalVector();
  dVdz_i_minus_half.RestoreDataPointerToLocalVector();
  dVdz_j_minus_half.RestoreDataPointerToLocalVector();
  dVdz_k_minus_half.RestoreDataPointerToLocalVector();
  delta_xyz.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();

  R.RestoreDataPointerToLocalVector(); //NOTE: although R has been updated, there is no need of
                                       //      cross-subdomain communications. So, no need to
                                       //      update the global vec.
}

//--------------------------------------------------------------------------

void
ViscosityOperator::AddCylindricalSymmetryTerms(SpaceVariable3D &V, double*** id, Vec3D*** dxyz,
                       vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS,
                       Vec5D*** res)
{

  assert(kmax-k0==1); //should be single layer in z-direction

  Vec5D*** v = (Vec5D***)V.GetDataPointer();


  // Step 1: Calculate lambda and mu (only in domain interior)
  double*** lam = Lam.GetDataPointer();
  double*** mu  = Lam.GetDataPointer();
  int myid;
  for(int k=k0; k<kmax; k++) //NOTE: not going to calculate d/dz. No need to copy ghosts!
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++) {

        myid = id[k][j][i];

        if(myid == INACTIVE_MATERIAL_ID) {
          mu[k][j][i]  = 0.0;
          lam[k][j][i] = 0.0;
          continue;
        } 

        switch(visFcn[myid]->type) {
          case ViscoFcnBase::NONE :
            mu[k][j][i]  = 0.0;
            lam[k][j][i] = 0.0;
            break;
          case ViscoFcnBase::CONSTANT :
          case ViscoFcnBase::SUTHERLAND :
            mu[k][j][i]  = visFcn[myid]->GetMu(V);
            lam[k][j][i] = visFcn[myid]->GetLambda(V);
            break;
          case ViscoFcnBase::ARTIFICIAL_RODIONOV :
            double div = CalculateLocalDiv2D(v, id, coords, i, j, k);
            double h   = 0.5*(dxyz[k][j][i][0] + dxyz[k][j][i][1]);
            mu[k][j][i]  = visFcn[myid]->GetMu(V);
            lam[k][j][i] = visFcn[myid]->GetLambda(V);
            break;
          default:
            fprintf(stderr,"\033[0;31m*** Error: ViscosityOperator detected unsupported viscosity model.\033[0m\n");
            exit(-1);
            break;
        }

      }


 
  // Step 2: Calculate dudx (i.e. dw/dz in KW notes) and dudy (i.e. dw/dr in KW notes)
  vector<int> ind0{0};
  double*** s = scalarG2.GetDataPointer();
  for(int k=k0; k<kmax; k++) //NOTE: not going to calculate d/dz. No need to copy ghosts!
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++)
        s[k][j][i] = v[k][j][i][1]; //extracting the axial velocity
  scalarG2.RestoreDataPointerAndInsert(); //need to exchange

  grad_minus->CalculateFirstDerivativeAtNodes(0/*x*/, scalarG2, ind0, dudxl, ind0);
  grad_plus->CalculateFirstDerivativeAtNodes(0/*x*/, scalarG2, ind0, dudxr, ind0);
  grad_minus->CalculateFirstDerivativeAtNodes(1/*y*/, scalarG2, ind0, dudyl, ind0);
  grad_plus->CalculateFirstDerivativeAtNodes(1/*y*/, scalarG2, ind0, dudyr, ind0);

  // dudx --> dw/dz in KW's notes
  double*** uxl = dudxl.GetDataPointer();
  double*** uxr = dudxr.GetDataPointer();
  // dudy --> dw/dr in KW's notes
  double*** uyl = dudyl.GetDataPointer();
  double*** uyr = dudyr.GetDataPointer();
  
  double a,b,c,d;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        double r = coords[k][j][i][1];
        assert(r>0);

        a = (i-2>=-1) ? phil[k][j][i] : (phi[k][j][i]-phi[k][j][i-1])/(coords[k][j][i][0]-coords[k][j][i-1][0]);
        b = (i+2<=NX) ? phir[k][j][i] : (phi[k][j][i+1]-phi[k][j][i])/(coords[k][j][i+1][0]-coords[k][j][i][0]);
        c = (j-2>=-1) ? phib[k][j][i] : (phi[k][j][i]-phi[k][j-1][i])/(coords[k][j][i][1]-coords[k][j-1][i][1]);
        d = (j+2<=NY) ? phit[k][j][i] : (phi[k][j+1][i]-phi[k][j][i])/(coords[k][j+1][i][1]-coords[k][j][i][1]);

I AM HERE
        res[k][j][i] = (v[k][j][i][1]>=0 ? a*v[k][j][i][1] : b*v[k][j][i][1])
                     + (v[k][j][i][2]>=0 ? c*v[k][j][i][2] : d*v[k][j][i][2])
                     + (v[k][j][i][3]>=0 ? e*v[k][j][i][3] : f*v[k][j][i][3]);

        res[k][j][i] *= -1.0; //moves the residual to the RHS

      }

}

//--------------------------------------------------------------------------
// calculates dudx + dvdy using a first order method.
double
ViscosityOperator::CalculateLocalDiv2D(Vec5D*** v, double*** id, Vec3D*** coords, int i, int j, int k)
{
  assert(id[k][j][i] != INACTIVE_MATERIAL_ID);

  double div(0.0);

  int il = (i-1>=ii0  && id[k][j][i-1] != INACTIVE_MATERIAL_ID) ? i-1 : i;
  int ir = (i+1<iimax && id[k][j][i+1] != INACTIVE_MATERIAL_ID) ? i+1 : i;
  if(il!=ir) //otherwise, i-1 and i+1 are both unavailable ==> set dudx = 0
    div += (v[k][j][ir][1] - v[k][j][il][1])/(coords[k][j][ir][0] - coords[k][j][il][0]);
    
  int jb = (j-1>=jj0  && id[k][j-1][i] != INACTIVE_MATERIAL_ID) ? j-1 : j;
  int jt = (j+1<jjmax && id[k][j+1][i] != INACTIVE_MATERIAL_ID) ? j+1 : j;
  if(jb!=jt)
    div += (v[k][jt][i][2] - v[k][jb][i][2])/(coords[k][jt][i][1] - coords[k][jb][i][1]);

  return div;
}

//--------------------------------------------------------------------------


