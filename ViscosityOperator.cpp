#include <ViscosityOperator.h>
#include <Vector5D.h>
#include <Utils.h>
#include <bits/stdc++.h>
using std::unique_ptr;
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
                    Lam(NULL), Mu(NULL), dLamdx(NULL), dLamdy(NULL)
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

  if(dudx) dudx->Destroy();
  if(dvdx) dvdx->Destroy();
  if(dudy) dudy->Destroy();
  if(dvdy) dvdy->Destroy();

  if(Lam) Lam->Destroy();
  if(Mu)  Mu->Destroy();
  if(dLamdx) dLamdx->Destroy();
  if(dLamdy) dLamdy->Destroy();

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

  XXX TODO I AM HERE XXX

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






