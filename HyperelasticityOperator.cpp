/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<HyperelasticityOperator.h>
#include<EmbeddedBoundaryDataSet.h>
#include<linear_algebra.h>
#include<Vector5D.h>

//------------------------------------------------------------

HyperelasticityOperator::HyperelasticityOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, 
                             IoData &iod_, vector<VarFcnBase*> &varFcn_, SpaceVariable3D &coordinates_,
                             SpaceVariable3D &delta_xyz_, GlobalMeshInfo &global_mesh_,
                             InterpolatorBase &interpolator_, GradientCalculatorBase &grad_,
                             std::vector<GhostPoint> &ghost_nodes_inner_,
                             std::vector<GhostPoint> &ghost_nodes_outer_)
                       : comm(comm_), iod(iod_), varFcn(varFcn_), global_mesh(global_mesh_),
                         interpolator(interpolator_), grad(grad_),
                         refmap(comm_, dm_all_, iod_, coordinates_, delta_xyz_,
                                global_mesh_, ghost_nodes_inner_, ghost_nodes_outer_),
                         F(comm_, &(dm_all_.ghosted1_9dof)),
                         DetF(comm_, &(dm_all_.ghosted1_1dof)),
                         Var1(comm_, &(dm_all_.ghosted1_3dof)),
                         Var2(comm_, &(dm_all_.ghosted1_3dof)),
                         Var3(comm_, &(dm_all_.ghosted1_3dof)),
                         V_i_minus_half(comm_, &(dm_all_.ghosted1_3dof)),
                         V_j_minus_half(comm_, &(dm_all_.ghosted1_3dof)),
                         V_k_minus_half(comm_, &(dm_all_.ghosted1_3dof)),
                         dXidx_i_minus_half(comm_, &(dm_all_.ghosted1_3dof)),
                         dXidx_j_minus_half(comm_, &(dm_all_.ghosted1_3dof)),
                         dXidx_k_minus_half(comm_, &(dm_all_.ghosted1_3dof)),
                         dXidy_i_minus_half(comm_, &(dm_all_.ghosted1_3dof)),
                         dXidy_j_minus_half(comm_, &(dm_all_.ghosted1_3dof)),
                         dXidy_k_minus_half(comm_, &(dm_all_.ghosted1_3dof)),
                         dXidz_i_minus_half(comm_, &(dm_all_.ghosted1_3dof)),
                         dXidz_j_minus_half(comm_, &(dm_all_.ghosted1_3dof)),
                         dXidz_k_minus_half(comm_, &(dm_all_.ghosted1_3dof))
{
  coordinates_.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates_.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);


  // Create a HyperelasticitFcn for each material.
  // For problems involving multiple materials, it could happen that the user specified
  // hyperelasticity model for only some of the materials, but not all. In this case, we
  // still create a dummy (i.e. base) function for those unspecified materials.

  for(auto it = iod.eqs.materials.dataMap.begin(); it != iod.eqs.materials.dataMap.end(); it++)
    hyperFcn.push_back(NULL); //allocate space

  for(auto it = iod.eqs.materials.dataMap.begin(); it != iod.eqs.materials.dataMap.end(); it++) {
    int matid = it->first;
    if(matid < 0 || matid >= (int)hyperFcn.size()) {
      print_error("*** Error: Detected error in the specification of material indices (id = %d).\n", matid);
      exit_mpi();
    }
    switch (it->second->hyperelasticity.type) {
      case HyperelasticityModelData::NONE :
        hyperFcn[matid] = new HyperelasticityFcnBase(*varFcn[matid]);
        break;
      case HyperelasticityModelData::SAINTVENANT_KIRCHHOFF :
        hyperFcn[matid] = new HyperelasticityFcnSaintVenantKirchhoff(it->second->hyperelasticity,
                                                                     *varFcn[matid]);
        break;
      case HyperelasticityModelData::MODIFIED_SAINTVENANT_KIRCHHOFF :
        hyperFcn[matid] = new HyperelasticityFcnModifiedSaintVenantKirchhoff(it->second->hyperelasticity,
                                                                             *varFcn[matid]);
        break;
      case HyperelasticityModelData::NEO_HOOKEAN :
        hyperFcn[matid] = new HyperelasticityFcnNeoHookean(it->second->hyperelasticity, *varFcn[matid]);
        break;
      case HyperelasticityModelData::MOONEY_RIVLIN :
        hyperFcn[matid] = new HyperelasticityFcnMooneyRivlin(it->second->hyperelasticity, *varFcn[matid]);
        break;
    }
  }

}

//------------------------------------------------------------

HyperelasticityOperator::~HyperelasticityOperator()
{
  for(int i=0; i<(int)hyperFcn.size(); i++)
    delete hyperFcn[i];
}

//------------------------------------------------------------

void
HyperelasticityOperator::Destroy()
{
  refmap.Destroy();
  F.Destroy();
  DetF.Destroy();
  Var1.Destroy();
  Var2.Destroy();
  Var3.Destroy();
  
  V_i_minus_half.Destroy();
  V_j_minus_half.Destroy();
  V_k_minus_half.Destroy();
  dXidx_i_minus_half.Destroy();
  dXidx_j_minus_half.Destroy();
  dXidx_k_minus_half.Destroy();
  dXidy_i_minus_half.Destroy();
  dXidy_j_minus_half.Destroy();
  dXidy_k_minus_half.Destroy();
  dXidz_i_minus_half.Destroy();
  dXidz_j_minus_half.Destroy();
  dXidz_k_minus_half.Destroy();
}

//------------------------------------------------------------

void
HyperelasticityOperator::InitializeReferenceMap(SpaceVariable3D &Xi)
{
  refmap.SetInitialCondition(Xi);
}

//------------------------------------------------------------

void
HyperelasticityOperator::ApplyBoundaryConditionsToReferenceMap(SpaceVariable3D &Xi)
{
  refmap.ApplyBoundaryConditions(Xi);
}

//------------------------------------------------------------

void
HyperelasticityOperator::ComputeReferenceMapResidual(SpaceVariable3D &V, 
                                                     SpaceVariable3D &Xi, SpaceVariable3D &R)
{
  refmap.ComputeResidual(V,Xi,R);
}

//------------------------------------------------------------

void
HyperelasticityOperator::ComputeDeformationGradientAtNodes(SpaceVariable3D &Xi)
{

  // NOTE: This function computes deformation gradient (F) ONLY WITHIN THE PHYSICAL DOMAIN,
  //       NOT in the ghost boundary layer.

  // ------------------------------------
  // Step 1: Calculate the Jacobian of Xi
  // ------------------------------------
  std::vector<int> deriv_dofs{0,1,2};    
  // dxi/dx
  grad.CalculateFirstDerivativeAtNodes(0 /*dx*/, Xi, deriv_dofs, Var1, deriv_dofs);
  // dxi/dy
  grad.CalculateFirstDerivativeAtNodes(1 /*dy*/, Xi, deriv_dofs, Var2, deriv_dofs);
  // dxi/dz
  grad.CalculateFirstDerivativeAtNodes(2 /*dy*/, Xi, deriv_dofs, Var2, deriv_dofs);


  // ------------------------------------
  // Step 2: Calculate the deformation gradient: F = inv(grad.Xi) and its determinant
  // ------------------------------------
  double*** f    = F.GetDataPointer(); //column-first (aka. column-major)
  double*** detf = DetF.GetDataPointer();
  Vec3D*** dXidx = (Vec3D***)Var1.GetDataPointer();  
  Vec3D*** dXidy = (Vec3D***)Var2.GetDataPointer();  
  Vec3D*** dXidz = (Vec3D***)Var3.GetDataPointer();  

  double gradxi[9]; //local values of Jacobian of Xi
  bool invertible = false;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        for(int dim=0; dim<3; dim++)
          gradxi[dim]   = dXidx[k][j][i][dim];
        for(int dim=0; dim<3; dim++)
          gradxi[3+dim] = dXidy[k][j][i][dim];
        for(int dim=0; dim<3; dim++)
          gradxi[6+dim] = dXidz[k][j][i][dim];

        invertible = MathTools::LinearAlgebra::
                     CalculateMatrixInverseAndDeterminant3x3(gradxi, &f[k][j][i], &detf[k][j][i]);
        if(!invertible)
          fprintf(stdout,"\033[0;35mWarning: Jacobian of ref. map at (%d,%d,%d) is not invertible."
                         " determinant = %e.\033[0m\n", i,j,k, detf[k][j][i]);

        detf[k][j][i] = 1.0/detf[k][j][i]; // we want the determinant of F, not grad-xi!

        if(detf[k][j][i]<=0.0)
          fprintf(stdout,"\033[0;35mWarning: Determinant of deformation gradient at (%d,%d,%d)"
                         " is negative (%e).\033[0m\n", i,j,k, detf[k][j][i]);
      }

  F.RestoreDataPointerAndInsert();
  DetF.RestoreDataPointerAndInsert();
  Var1.RestoreDataPointerToLocalVector();
  Var2.RestoreDataPointerToLocalVector();
  Var3.RestoreDataPointerToLocalVector();

}

//------------------------------------------------------------

void
HyperelasticityOperator::AddHyperelasticityFluxes(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &Xi,
                                                  vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS,
                                                  SpaceVariable3D &R)
{

  if(EBDS)
    print_warning("Warning: AddHyperelasticityFluxes: Not able to account for embedded surfaces.\n");

  std::vector<int> i123{1,2,3}, i012{0,1,2};

  //1. Calculate the x, y, and z velocities at cell interfaces by interpolation
  interpolator.InterpolateAtCellInterfaces(0/*x-dir*/, V, i123, V_i_minus_half, i012);
  interpolator.InterpolateAtCellInterfaces(1/*y-dir*/, V, i123, V_j_minus_half, i012);
  interpolator.InterpolateAtCellInterfaces(2/*z-dir*/, V, i123, V_k_minus_half, i012); 

  //2. Calculate Xi derivatives at cell interfaces
  grad.CalculateFirstDerivativeAtCellInterfaces(0, 0, Xi, i012, dXidx_i_minus_half, i012);
  grad.CalculateFirstDerivativeAtCellInterfaces(0, 1, Xi, i012, dXidx_j_minus_half, i012);
  grad.CalculateFirstDerivativeAtCellInterfaces(0, 2, Xi, i012, dXidx_k_minus_half, i012);
  grad.CalculateFirstDerivativeAtCellInterfaces(1, 0, Xi, i012, dXidy_i_minus_half, i012);
  grad.CalculateFirstDerivativeAtCellInterfaces(1, 1, Xi, i012, dXidy_j_minus_half, i012);
  grad.CalculateFirstDerivativeAtCellInterfaces(1, 2, Xi, i012, dXidy_k_minus_half, i012);
  grad.CalculateFirstDerivativeAtCellInterfaces(2, 0, Xi, i012, dXidz_i_minus_half, i012);
  grad.CalculateFirstDerivativeAtCellInterfaces(2, 1, Xi, i012, dXidz_j_minus_half, i012);
  grad.CalculateFirstDerivativeAtCellInterfaces(2, 2, Xi, i012, dXidz_k_minus_half, i012);

  //3. Loop through cell interfaces and calculate fluxes
  Vec3D*** dxidx_i = (Vec3D***)dXidx_i_minus_half.GetDataPointer();
  Vec3D*** dxidx_j = (Vec3D***)dXidx_j_minus_half.GetDataPointer();
  Vec3D*** dxidx_k = (Vec3D***)dXidx_k_minus_half.GetDataPointer();
  Vec3D*** dxidy_i = (Vec3D***)dXidy_i_minus_half.GetDataPointer();
  Vec3D*** dxidy_j = (Vec3D***)dXidy_j_minus_half.GetDataPointer();
  Vec3D*** dxidy_k = (Vec3D***)dXidy_k_minus_half.GetDataPointer();
  Vec3D*** dxidz_i = (Vec3D***)dXidz_i_minus_half.GetDataPointer();
  Vec3D*** dxidz_j = (Vec3D***)dXidz_j_minus_half.GetDataPointer();
  Vec3D*** dxidz_k = (Vec3D***)dXidz_k_minus_half.GetDataPointer();

  Vec5D*** v    = (Vec5D***)V.GetDataPointer();
  Vec5D*** res  = (Vec5D***)R.GetDataPointer();
  double*** id  = (double***)ID.GetDataPointer();

  int myid = 0;
  double gradxi[9], f[9]; //nabla xi and deformation gradient
  double dx = 0.0, dy = 0.0, dz = 0.0, detf;
  Vec5D flux;
  bool invertible;
  for(int k=k0; k<kkmax; k++)
    for(int j=j0; j<jjmax; j++)
      for(int i=i0; i<iimax; i++) {

        myid = id[k][j][i];
        dx   = global_mesh.dx_glob[i];
        dy   = global_mesh.dy_glob[j];
        dz   = global_mesh.dz_glob[k];

        //*****************************************
        //calculate flux function F_{i-1/2,j,k}
        //*****************************************
        if(k!=kkmax-1 && j!=jjmax-1) {

          for(int dim=0; dim<3; dim++)
            gradxi[dim]   = dxidx_i[k][j][i][dim];
          for(int dim=0; dim<3; dim++)
            gradxi[3+dim] = dxidy_i[k][j][i][dim];
          for(int dim=0; dim<3; dim++)
            gradxi[6+dim] = dxidz_i[k][j][i][dim];

          invertible = MathTools::LinearAlgebra::
                       CalculateMatrixInverseAndDeterminant3x3(gradxi, f, &detf);
          if(!invertible)
            fprintf(stdout,"\033[0;35mWarning: Jacobian of ref. map at (%d-1/2,%d,%d) is not invertible."
                           " determinant = %e.\033[0m\n", i,j,k, detf);

          hyperFcn[myid]->EvaluateHyperelasticFluxFunction_F(flux, f, v[k][j][i], true);//TODO: multi-material
          flux *= dy*dz;
          res[k][j][i]   += flux;
          res[k][j][i-1] -= flux;
        }


        //*****************************************
        //calculate flux function G_{i,j-1/2,k}
        //*****************************************
        if(i!=iimax-1 && k!=kkmax-1) {

          for(int dim=0; dim<3; dim++)
            gradxi[dim]   = dxidx_j[k][j][i][dim];
          for(int dim=0; dim<3; dim++)
            gradxi[3+dim] = dxidy_j[k][j][i][dim];
          for(int dim=0; dim<3; dim++)
            gradxi[6+dim] = dxidz_j[k][j][i][dim];

          invertible = MathTools::LinearAlgebra::
                       CalculateMatrixInverseAndDeterminant3x3(gradxi, f, &detf);
          if(!invertible)
            fprintf(stdout,"\033[0;35mWarning: Jacobian of ref. map at (%d,%d-1/2,%d) is not invertible."
                           " determinant = %e.\033[0m\n", i,j,k, detf);

          hyperFcn[myid]->EvaluateHyperelasticFluxFunction_G(flux, f, v[k][j][i], true);//TODO: multi-material
          flux *= dx*dz;
          res[k][j][i]   += flux;
          res[k][j-1][i] -= flux;
        }

        //*****************************************
        //calculate flux function H_{i,j,k-1/2}
        //*****************************************
        if(i!=iimax-1 && j!=jjmax-1) {

          for(int dim=0; dim<3; dim++)
            gradxi[dim]   = dxidx_k[k][j][i][dim];
          for(int dim=0; dim<3; dim++)
            gradxi[3+dim] = dxidy_k[k][j][i][dim];
          for(int dim=0; dim<3; dim++)
            gradxi[6+dim] = dxidz_k[k][j][i][dim];

          invertible = MathTools::LinearAlgebra::
                       CalculateMatrixInverseAndDeterminant3x3(gradxi, f, &detf);
          if(!invertible)
            fprintf(stdout,"\033[0;35mWarning: Jacobian of ref. map at (%d,%d,%d-1/2) is not invertible."
                           " determinant = %e.\033[0m\n", i,j,k, detf);

          hyperFcn[myid]->EvaluateHyperelasticFluxFunction_H(flux, f, v[k][j][i], true);//TODO: multi-material!
          flux *= dx*dy;
          res[k][j][i]   += flux;
          res[k-1][j][i] -= flux;
        }

      }

  V_i_minus_half.RestoreDataPointerToLocalVector();
  V_j_minus_half.RestoreDataPointerToLocalVector();
  V_k_minus_half.RestoreDataPointerToLocalVector();
  dXidx_i_minus_half.RestoreDataPointerToLocalVector();
  dXidx_j_minus_half.RestoreDataPointerToLocalVector();
  dXidx_k_minus_half.RestoreDataPointerToLocalVector();
  dXidy_i_minus_half.RestoreDataPointerToLocalVector();
  dXidy_j_minus_half.RestoreDataPointerToLocalVector();
  dXidy_k_minus_half.RestoreDataPointerToLocalVector();
  dXidz_i_minus_half.RestoreDataPointerToLocalVector();
  dXidz_j_minus_half.RestoreDataPointerToLocalVector();
  dXidz_k_minus_half.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();
  V.RestoreDataPointerToLocalVector();

  R.RestoreDataPointerToLocalVector(); //NOTE: although R has been updated, there is no need of
                                       //      cross-subdomain communications. So, no need to
                                       //      update the global vec.
}

//------------------------------------------------------------


//------------------------------------------------------------

