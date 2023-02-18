/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <ViscosityOperator.h>
#include <GradientCalculatorFD3.h>
#include <Vector5D.h>
#include <Utils.h>
#include <bits/stdc++.h>
using std::unique_ptr;

extern int verbose;
extern int INACTIVE_MATERIAL_ID;


//--------------------------------------------------------------------------

ViscosityOperator::ViscosityOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_,
                                     vector<VarFcnBase*>& varFcn_, GlobalMeshInfo &global_mesh_,
                                     SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_,
                                     SpaceVariable3D &volume_,
                                     InterpolatorBase &interpolator, GradientCalculatorBase &grad_,
                                     bool with_embedded_boundary)
                  : iod_eqs(iod_.eqs), coordinates(coordinates_), delta_xyz(delta_xyz_),
                    volume(volume_), varFcn(varFcn_), global_mesh(global_mesh_),
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
                    interpolator(interpolator), grad(grad_), gfo(NULL), Vgf(NULL),
                    cylindrical_symmetry(false),
                    DDXm(NULL), DDXp(NULL), DDYm(NULL), DDYp(NULL), 
                    Lam(NULL), Mu(NULL), scalarG2(NULL),
                    grad_minus(NULL), grad_plus(NULL)
{

  // Get i0, j0, etc.
  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);
  coordinates.GetGlobalSize(&NX,&NY,&NZ);

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

  // Initialize Vgf if ebo is not null.
  if(with_embedded_boundary) {
    gfo = new GhostFluidOperator(comm_, dm_all_, global_mesh);
    Vgf = new SpaceVariable3D(comm_, &(dm_all_.ghosted1_5dof));
  }

  // Initialize variables for cylindrical symmetry
  if(iod_.mesh.type == MeshData::CYLINDRICAL) {

    cylindrical_symmetry = true;

    DDXm = new SpaceVariable3D(comm_, &(dm_all_.ghosted1_1dof));
    DDXp = new SpaceVariable3D(comm_, &(dm_all_.ghosted1_1dof));
    DDYm = new SpaceVariable3D(comm_, &(dm_all_.ghosted1_1dof));
    DDYp = new SpaceVariable3D(comm_, &(dm_all_.ghosted1_1dof));

    Lam = new SpaceVariable3D(comm_, &(dm_all_.ghosted1_1dof));
    Mu  = new SpaceVariable3D(comm_, &(dm_all_.ghosted1_1dof));

    scalarG2 = new SpaceVariable3D(comm_, &(dm_all_.ghosted2_1dof)); //w/ 2 ghost layers

    if(NZ != 1) {
      print_error("*** Error: Enforcing cylindrical symmetry, but mesh has %d layers "
                    "in z-direction.\n", NZ); 
      exit_mpi();
    }

    grad_minus = new GradientCalculatorFD3(comm_, dm_all_, coordinates, delta_xyz, -1);
    grad_plus  = new GradientCalculatorFD3(comm_, dm_all_, coordinates, delta_xyz,  1);

  }
  else if(iod_.mesh.type == MeshData::SPHERICAL) {
    print_error("*** Error: Currently, ViscosityOperator does not support spherical symmetry.\n");
    exit(-1);
  }

}

//--------------------------------------------------------------------------

ViscosityOperator::~ViscosityOperator()
{
  for(int i=0; i<(int)visFcn.size(); i++)
    delete visFcn[i];

  if(gfo) delete gfo;
  if(Vgf) delete Vgf;

  if(DDXm) delete DDXm;
  if(DDXp) delete DDXp;
  if(DDYm) delete DDYm;
  if(DDYp) delete DDYp;
  if(Lam)  delete Lam;
  if(Mu)   delete Mu;

  if(scalarG2) delete scalarG2;

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

  if(gfo) gfo->Destroy();
  if(Vgf) Vgf->Destroy();

  if(DDXm) DDXm->Destroy();
  if(DDXp) DDXp->Destroy();
  if(DDYm) DDYm->Destroy();
  if(DDYp) DDYp->Destroy();

  if(Lam) Lam->Destroy();
  if(Mu)  Mu->Destroy();

  if(scalarG2) scalarG2->Destroy();

  if(grad_minus)  grad_minus->Destroy();
  if(grad_plus)   grad_plus->Destroy();

}

//--------------------------------------------------------------------------
// Add diffusion fluxes on the left hand side of the N-S equations
void ViscosityOperator::AddDiffusionFluxes(SpaceVariable3D &V, SpaceVariable3D &ID, 
                                           vector<unique_ptr<EmbeddedBoundaryDataSet> > *EBDS,
                                           SpaceVariable3D &R)
{

  SpaceVariable3D *VV = &V;

  if(EBDS && EBDS->size()>0) {
    assert(gfo); 
    int numGhostPopulated = gfo->PopulateGhostNodesForViscosityOperator(V, ID, EBDS, *Vgf);
    if(numGhostPopulated>0) // In this case, Vgf should have been filled, and it is the one to use below.
      VV = Vgf;

    //print_warning("Warning: ViscosityOperator::AddDiffusionFluxes: Not able to account for embedded surfaces.\n");
  }


  std::vector<int> i123{1,2,3}, i012{0,1,2}; 
  //1. Calculate the x, y, and z velocities at cell interfaces by interpolation 
  interpolator.InterpolateAtCellInterfaces(0/*x-dir*/, *VV, i123, V_i_minus_half, i012);
  interpolator.InterpolateAtCellInterfaces(1/*y-dir*/, *VV, i123, V_j_minus_half, i012);
  interpolator.InterpolateAtCellInterfaces(2/*z-dir*/, *VV, i123, V_k_minus_half, i012);

  //2. Calculate velocity derivatives at cell interfaces
  grad.CalculateFirstDerivativeAtCellInterfaces(0, 0, *VV, i123, dVdx_i_minus_half, i012);
  grad.CalculateFirstDerivativeAtCellInterfaces(0, 1, *VV, i123, dVdx_j_minus_half, i012);
  grad.CalculateFirstDerivativeAtCellInterfaces(0, 2, *VV, i123, dVdx_k_minus_half, i012);
  grad.CalculateFirstDerivativeAtCellInterfaces(1, 0, *VV, i123, dVdy_i_minus_half, i012);
  grad.CalculateFirstDerivativeAtCellInterfaces(1, 1, *VV, i123, dVdy_j_minus_half, i012);
  grad.CalculateFirstDerivativeAtCellInterfaces(1, 2, *VV, i123, dVdy_k_minus_half, i012);
  grad.CalculateFirstDerivativeAtCellInterfaces(2, 0, *VV, i123, dVdz_i_minus_half, i012);
  grad.CalculateFirstDerivativeAtCellInterfaces(2, 1, *VV, i123, dVdz_j_minus_half, i012);
  grad.CalculateFirstDerivativeAtCellInterfaces(2, 2, *VV, i123, dVdz_k_minus_half, i012);

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



  if(cylindrical_symmetry) //Directly use id, dxyz, and res (local vars) 
    AddCylindricalSymmetryTerms(*VV, id, dxyz, EBDS, res);


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

  assert(cylindrical_symmetry);
  assert(kmax-k0==1); //should be single layer in z-direction

  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
  double*** vol   = volume.GetDataPointer();

  Vec5D*** v = (Vec5D***)V.GetDataPointer();

  // Step 1: Calculate lambda and mu (for z-direction, only in domain interior)
  double*** lam = Lam->GetDataPointer();
  double*** mu  = Mu->GetDataPointer();
  int myid;
  double div, h;
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
            mu[k][j][i]  = visFcn[myid]->GetMu();
            lam[k][j][i] = visFcn[myid]->GetLambda();
            break;
          case ViscoFcnBase::SUTHERLAND :
            mu[k][j][i]  = visFcn[myid]->GetMu(v[k][j][i]);
            lam[k][j][i] = visFcn[myid]->GetLambda(v[k][j][i]);
            break;
          case ViscoFcnBase::ARTIFICIAL_RODIONOV :
            div = CalculateLocalDiv2D(v, id, coords, i, j, k);
            h   = 0.5*(dxyz[k][j][i][0] + dxyz[k][j][i][1]);
            mu[k][j][i]  = visFcn[myid]->GetMu(v[k][j][i], div, h);
            lam[k][j][i] = visFcn[myid]->GetLambda(v[k][j][i], div, h);
            break;
          default:
            fprintf(stderr,"\033[0;31m*** Error: ViscosityOperator detected unsupported viscosity model.\033[0m\n");
            exit(-1);
            break;
        }

      }

  // exchange data, so internal ghost nodes (i.e. owned by neighbors) get correct mu and lam.
  Lam->RestoreDataPointerAndInsert();
  Mu->RestoreDataPointerAndInsert();


 
  // Step 2: Add the terms that involve dudx (i.e. dw/dz in KW notes), dudy (i.e. dw/dr);
  //         Apply an upwinding 3rd order FD.

  lam = Lam->GetDataPointer();
  mu  = Mu->GetDataPointer();

  vector<int> ind0{0};
  double*** s = scalarG2->GetDataPointer();
  for(int k=k0; k<kmax; k++) //NOTE: not going to calculate d/dz. No need to copy ghosts!
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++)
        s[k][j][i] = v[k][j][i][1]; //extracting the axial velocity
  scalarG2->RestoreDataPointerAndInsert(); //need to exchange

  grad_minus->CalculateFirstDerivativeAtNodes(0/*x*/, *scalarG2, ind0, *DDXm, ind0); //dudxl
  grad_plus->CalculateFirstDerivativeAtNodes(0/*x*/, *scalarG2, ind0, *DDXp, ind0); //dudxr
  grad_minus->CalculateFirstDerivativeAtNodes(1/*y*/, *scalarG2, ind0, *DDYm, ind0); //dudyb
  grad_plus->CalculateFirstDerivativeAtNodes(1/*y*/, *scalarG2, ind0, *DDYp, ind0); //dudyt

  // dudx --> dw/dz in KW's notes
  double*** dwdzl = DDXm->GetDataPointer();
  double*** dwdzr = DDXp->GetDataPointer();
  // dudy --> dw/dr in KW's notes
  double*** dwdrb = DDYm->GetDataPointer();
  double*** dwdrt = DDYp->GetDataPointer();

  double r,w,ur,cv,c;
  double mu_over_r;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        if(id[k][j][i] == INACTIVE_MATERIAL_ID)
          continue;

        ur = v[k][j][i][2]; // radial velocity 
        w  = v[k][j][i][1]; // axial velocity
        cv = vol[k][j][i]; // cell volume
        r = coords[k][j][i][1]; // radial coord
        assert(r>0);

        //TODO: double-check! I think grad_minus/plus already takes care of i=0/NX-1
        //      So, no need of a,b,c,d (cf. LevelSetOperator)

        // (mu/r)*dw/dr ~ axial momentum eq
        mu_over_r = mu[k][j][i]/r;
        res[k][j][i][1] -= mu_over_r>0 ?
                           cv*mu_over_r*dwdrb[k][j][i] : cv*mu_over_r*dwdrt[k][j][i];
        
        // (mu*w/r)*dw/dr + (2*lam*ur/r)*dw/dz ~ energy eq.
        c = mu_over_r*w;
        res[k][j][i][4] -= c>0 ? cv*c*dwdrb[k][j][i] : cv*c*dwdrt[k][j][i];

        c = 2.0*lam[k][j][i]*ur/r;
        res[k][j][i][4] -= c>0 ? cv*c*dwdzl[k][j][i] : cv*c*dwdzr[k][j][i];
      }

  DDXm->RestoreDataPointerToLocalVector();
  DDXp->RestoreDataPointerToLocalVector();
  DDYm->RestoreDataPointerToLocalVector();
  DDYp->RestoreDataPointerToLocalVector();



 
  // Step 3: Add the terms that involve dvdx (i.e. dur/dz in KW notes), dvdy (i.e. dur/dr);
  //         Apply an upwinding 3rd order FD.

  s = scalarG2->GetDataPointer();
  for(int k=k0; k<kmax; k++) //NOTE: not going to calculate d/dz. No need to copy ghosts!
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++)
        s[k][j][i] = v[k][j][i][2]; //extracting the radial velocity
  scalarG2->RestoreDataPointerAndInsert(); //need to exchange

  grad_minus->CalculateFirstDerivativeAtNodes(0/*x*/, *scalarG2, ind0, *DDXm, ind0); //dvdxl
  grad_plus->CalculateFirstDerivativeAtNodes(0/*x*/, *scalarG2, ind0, *DDXp, ind0); //dvdxr
  grad_minus->CalculateFirstDerivativeAtNodes(1/*y*/, *scalarG2, ind0, *DDYm, ind0); //dvdyb
  grad_plus->CalculateFirstDerivativeAtNodes(1/*y*/, *scalarG2, ind0, *DDYp, ind0); //dvdyt

  // dvdx --> dur/dz in KW's notes
  double*** durdzl = DDXm->GetDataPointer();
  double*** durdzr = DDXp->GetDataPointer();
  // dvdy --> dur/dr in KW's notes
  double*** durdrb = DDYm->GetDataPointer();
  double*** durdrt = DDYp->GetDataPointer();

  double lam_mu_over_r;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        if(id[k][j][i] == INACTIVE_MATERIAL_ID)
          continue;

        ur = v[k][j][i][2]; // radial velocity 
        w  = v[k][j][i][1]; // axial velocity
        cv = vol[k][j][i]; // cell volume
        r = coords[k][j][i][1]; // radial coord
        assert(r>0);

        // (lam+mu)/r*dur/dz ~ axial momentum eq
        lam_mu_over_r = (lam[k][j][i] + mu[k][j][i])/r;
        res[k][j][i][1] -= lam_mu_over_r>0 ?
                           cv*lam_mu_over_r*durdzl[k][j][i] : cv*lam_mu_over_r*durdzr[k][j][i];
        
        // (lam+2*mu)/r*dur/dr ~ radial momentum eq
        c = (lam[k][j][i] + 2.0*mu[k][j][i])/r;
        res[k][j][i][2] -= c>0 ? cv*c*durdrb[k][j][i] : cv*c*durdrt[k][j][i];

        // (3*lam+2*mu)*ur/r*durdr + (lam+mu)*w/r*durdz ~ energy eq.
        c = (3.0*lam[k][j][i] + 2.0*mu[k][j][i])*ur/r;
        res[k][j][i][4] -= c>0 ? cv*c*durdrb[k][j][i] : cv*c*durdrt[k][j][i];

        c = lam_mu_over_r*w;
        res[k][j][i][4] -= c>0 ? cv*c*durdzl[k][j][i] : cv*c*durdzr[k][j][i];

      }

  DDXm->RestoreDataPointerToLocalVector();
  DDXp->RestoreDataPointerToLocalVector();
  DDYm->RestoreDataPointerToLocalVector();
  DDYp->RestoreDataPointerToLocalVector();




 
  // Step 4: Add the terms that involve dlamdx (i.e. dlam/dz in KW notes), dlamdy (dlam/dr);
  //         Apply an upwinding 3rd order FD.
  //         Also add the source term (on the LHS of the N-S equations)

  s = scalarG2->GetDataPointer();
  for(int k=k0; k<kmax; k++) //NOTE: not going to calculate d/dz. No need to copy ghosts!
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++)
        s[k][j][i] = lam[k][j][i]; 
  scalarG2->RestoreDataPointerAndInsert(); //need to exchange

  grad_minus->CalculateFirstDerivativeAtNodes(0/*x*/, *scalarG2, ind0, *DDXm, ind0); //dlamdxl
  grad_plus->CalculateFirstDerivativeAtNodes(0/*x*/, *scalarG2, ind0, *DDXp, ind0); //dlamdxr
  grad_minus->CalculateFirstDerivativeAtNodes(1/*y*/, *scalarG2, ind0, *DDYm, ind0); //dlamdyb
  grad_plus->CalculateFirstDerivativeAtNodes(1/*y*/, *scalarG2, ind0, *DDYp, ind0); //dlamdyt

  // dvdx --> dlam/dz in KW's notes
  double*** dlamdzl = DDXm->GetDataPointer();
  double*** dlamdzr = DDXp->GetDataPointer();
  // dvdy --> dlam/dr in KW's notes
  double*** dlamdrb = DDYm->GetDataPointer();
  double*** dlamdrt = DDYp->GetDataPointer();

  double ur_over_r;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        if(id[k][j][i] == INACTIVE_MATERIAL_ID)
          continue;

        ur = v[k][j][i][2]; // radial velocity 
        w  = v[k][j][i][1]; // axial velocity
        cv = vol[k][j][i]; // cell volume
        r = coords[k][j][i][1]; // radial coord
        assert(r>0);

        // ur/r*dlam/dz ~ axial momentum eq
        ur_over_r = ur/r;
        res[k][j][i][1] -= ur_over_r>0 ?
                           cv*ur_over_r*dlamdzl[k][j][i] : cv*ur_over_r*dlamdzr[k][j][i];
        
        // ur/r*dlam/dr - (lam+2mu)*ur/(r*r) ~ radial momentum eq
        res[k][j][i][2] -= ur_over_r>0 ? 
                           cv*ur_over_r*dlamdrb[k][j][i] : cv*ur_over_r*dlamdrt[k][j][i];

        res[k][j][i][2] += cv*(lam[k][j][i] + 2.0*mu[k][j][i])*ur/(r*r); //source term

        // (u_r*u_r/r)*dlam/dr + (ur*w/r)*dlam/dz ~ energy eq.
        c = ur_over_r*ur;
        res[k][j][i][4] -= c>0 ? cv*c*dlamdrb[k][j][i] : cv*c*dlamdrt[k][j][i];

        c = ur_over_r*w;
        res[k][j][i][4] -= c>0 ? cv*c*dlamdzl[k][j][i] : cv*c*dlamdzr[k][j][i];

      }

  DDXm->RestoreDataPointerToLocalVector();
  DDXp->RestoreDataPointerToLocalVector();
  DDYm->RestoreDataPointerToLocalVector();
  DDYp->RestoreDataPointerToLocalVector();


  Lam->RestoreDataPointerToLocalVector();
  Mu->RestoreDataPointerToLocalVector();

  V.RestoreDataPointerToLocalVector();
  volume.RestoreDataPointerToLocalVector();
  coordinates.RestoreDataPointerToLocalVector();

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


