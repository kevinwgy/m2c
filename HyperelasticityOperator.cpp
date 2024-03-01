/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<HyperelasticityOperator.h>
#include<EmbeddedBoundaryDataSet.h>
#include<linear_algebra.h>
#include<trilinear_interpolation.h>
#include<Vector5D.h>

extern int INACTIVE_MATERIAL_ID;

//------------------------------------------------------------

HyperelasticityOperator::HyperelasticityOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, 
                             IoData &iod_, vector<VarFcnBase*> &varFcn_, SpaceVariable3D &coordinates_,
                             SpaceVariable3D &delta_xyz_, SpaceVariable3D &volume_,
                             GlobalMeshInfo &global_mesh_,
                             InterpolatorBase &interpolator_, GradientCalculatorBase &grad_,
                             std::vector<GhostPoint> &ghost_nodes_inner_,
                             std::vector<GhostPoint> &ghost_nodes_outer_)
                       : comm(comm_), iod(iod_), varFcn(varFcn_), global_mesh(global_mesh_),
                         coordinates(coordinates_), volume(volume_),
                         interpolator(interpolator_), grad(grad_),
                         refmap(comm_, dm_all_, iod_, coordinates_, delta_xyz_,
                                global_mesh_, ghost_nodes_inner_, ghost_nodes_outer_),
                         F(comm_, &(dm_all_.ghosted1_9dof)),
                         J(comm_, &(dm_all_.ghosted1_1dof)),
                         Var1(comm_, &(dm_all_.ghosted1_3dof)),
                         Var2(comm_, &(dm_all_.ghosted1_3dof)),
                         Var3(comm_, &(dm_all_.ghosted1_3dof)),
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

  cylindrical_symmetry = false;
  if(iod.mesh.type == MeshData::CYLINDRICAL) {
    cylindrical_symmetry = true;
    if(!global_mesh.two_dimensional_xy) {
      print_error("*** Error: HyperelasticityOperator detected incorrect mesh, "
                  "with cylindrical symmetry requested.\n");
      exit_mpi();
    }
  } else if(iod.mesh.type != MeshData::THREEDIMENSIONAL) {
    print_error("*** Error: HyperelasticityOperator cannot handle user-specified domain/mesh type.\n");
    exit_mpi();
  }


  // Create a HyperelasticityFcn for each material.
  // For problems involving multiple materials, it could happen that the user specified
  // hyperelasticity model for only some of the materials, but not all. In this case, we
  // still create a dummy (i.e. base) function for those unspecified materials.

  for(auto it = iod.eqs.materials.dataMap.begin(); it != iod.eqs.materials.dataMap.end(); it++) {
    hyperFcn.push_back(NULL); //allocate space
    deviator_only.push_back(true); //default is true
  }

  for(auto it = iod.eqs.materials.dataMap.begin(); it != iod.eqs.materials.dataMap.end(); it++) {
    int matid = it->first;
    if(matid < 0 || matid >= (int)hyperFcn.size()) {
      print_error("*** Error: Detected error in the specification of material indices (id = %d).\n", matid);
      exit_mpi();
    }

    deviator_only[matid] = (it->second->hyperelasticity.stress_option 
                            == HyperelasticityModelData::DEVIATOR_ONLY);

    switch (it->second->hyperelasticity.type) {
      case HyperelasticityModelData::NONE :
        if(cylindrical_symmetry)
          hyperFcn[matid] = new HyperelasticityFcnBase2DCyl(*varFcn[matid]);
        else
          hyperFcn[matid] = new HyperelasticityFcnBase(*varFcn[matid]);                       
        break;
      case HyperelasticityModelData::SAINTVENANT_KIRCHHOFF :
        if(cylindrical_symmetry)
          hyperFcn[matid] = new HyperelasticityFcnSaintVenantKirchhoff2DCyl
                                    (it->second->hyperelasticity, *varFcn[matid]);
        else
          hyperFcn[matid] = new HyperelasticityFcnSaintVenantKirchhoff
                                    (it->second->hyperelasticity, *varFcn[matid]);
        break;
      case HyperelasticityModelData::MODIFIED_SAINTVENANT_KIRCHHOFF :
        if(cylindrical_symmetry)
          hyperFcn[matid] = new HyperelasticityFcnModifiedSaintVenantKirchhoff2DCyl
                                    (it->second->hyperelasticity, *varFcn[matid]);
        else
          hyperFcn[matid] = new HyperelasticityFcnModifiedSaintVenantKirchhoff
                                    (it->second->hyperelasticity, *varFcn[matid]);
        break;
      case HyperelasticityModelData::NEO_HOOKEAN :
        if(cylindrical_symmetry)
          hyperFcn[matid] = new HyperelasticityFcnNeoHookean2DCyl
                                    (it->second->hyperelasticity, *varFcn[matid]);
        else
          hyperFcn[matid] = new HyperelasticityFcnNeoHookean
                                    (it->second->hyperelasticity, *varFcn[matid]);
        break;
      case HyperelasticityModelData::MOONEY_RIVLIN :
        if(cylindrical_symmetry)
          hyperFcn[matid] = new HyperelasticityFcnMooneyRivlin2DCyl
                                    (it->second->hyperelasticity, *varFcn[matid]);
        else
          hyperFcn[matid] = new HyperelasticityFcnMooneyRivlin
                                    (it->second->hyperelasticity, *varFcn[matid]);
        break;
      default :
        print_error("*** Error: HyperelasticityOperator detected unknown model type.\n");
        exit_mpi();
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
  J.Destroy();
  Var1.Destroy();
  Var2.Destroy();
  Var3.Destroy();
  
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
HyperelasticityOperator::ComputeReferenceMapResidual(SpaceVariable3D &V, SpaceVariable3D &Xi,
                                                     SpaceVariable3D &R, [[maybe_unused]] double time)
{

#ifdef HYPERELASTICITY_TEST
  PrescribeVelocityForTesting(V, time);
#endif

  refmap.ComputeResidual(V,Xi,R);
}

//------------------------------------------------------------

void
HyperelasticityOperator::ComputeDeformationGradientAtNodes(SpaceVariable3D &Xi)
{
  if(cylindrical_symmetry)
    ComputeDeformGradAtNodes2DCylindrical(Xi);
  else
    ComputeDeformGradAtNodes3D(Xi);
}

//------------------------------------------------------------

void
HyperelasticityOperator::ComputeDeformGradAtNodes3D(SpaceVariable3D &Xi)
{

  // NOTE: This function computes deformation gradient (F) ONLY WITHIN THE PHYSICAL DOMAIN,
  //       NOT in the ghost boundary layer.

  // ------------------------------------
  // Step 1: Calculate the Jacobian of Xi
  // ------------------------------------
  std::vector<int> i012{0,1,2};    
  // dxi/dx
  grad.CalculateFirstDerivativeAtNodes(0 /*dx*/, Xi, i012, Var1, i012);
  // dxi/dy
  grad.CalculateFirstDerivativeAtNodes(1 /*dy*/, Xi, i012, Var2, i012);
  // dxi/dz
  grad.CalculateFirstDerivativeAtNodes(2 /*dz*/, Xi, i012, Var3, i012);


  // ------------------------------------
  // Step 2: Calculate the deformation gradient: F = inv(grad.Xi) and its determinant
  // ------------------------------------
  double*** f    = F.GetDataPointer(); //column-first (aka. column-major)
  double*** Jloc = J.GetDataPointer();
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
                     CalculateMatrixInverseAndDeterminant3x3(gradxi, &f[k][j][9*i], &Jloc[k][j][i]);
        if(!invertible)
          fprintf(stdout,"\033[0;35mWarning: Jacobian of ref. map at (%d,%d,%d) is not invertible."
                         " determinant = %e.\033[0m\n", i,j,k, Jloc[k][j][i]);

        Jloc[k][j][i] = 1.0/Jloc[k][j][i]; // we want the determinant of F, not grad-xi!

        if(Jloc[k][j][i]<=0.0)
          fprintf(stdout,"\033[0;35mWarning: Determinant of deformation gradient at (%d,%d,%d)"
                         " is negative (%e).\033[0m\n", i,j,k, Jloc[k][j][i]);
      }

  F.RestoreDataPointerAndInsert();
  J.RestoreDataPointerAndInsert();
  Var1.RestoreDataPointerToLocalVector();
  Var2.RestoreDataPointerToLocalVector();
  Var3.RestoreDataPointerToLocalVector();

}

//------------------------------------------------------------

void
HyperelasticityOperator::ComputeDeformGradAtNodes2DCylindrical(SpaceVariable3D &Xi)
{

  // NOTE: This function computes deformation gradient (F) ONLY WITHIN THE PHYSICAL DOMAIN,
  //       NOT in the ghost boundary layer.

  // ------------------------------------
  // Step 1: Calculate the Jacobian of Xi
  // ------------------------------------
  std::vector<int> i01{0,1};    
  // dxi/dx (i.e., dxi/dz in cylindrical coords)
  grad.CalculateFirstDerivativeAtNodes(0 /*dx*/, Xi, i01, Var1, i01);
  // dxi/dy (i.e., dxi/dr in cylindrical coords)
  grad.CalculateFirstDerivativeAtNodes(1 /*dy*/, Xi, i01, Var2, i01);


  // ------------------------------------
  // Step 2: Calculate the deformation gradient: F = inv(grad.Xi) and its determinant
  // ------------------------------------
  double*** f    = F.GetDataPointer(); //column-first (aka. column-major)
  double*** Jloc = J.GetDataPointer();
  Vec3D*** xi    = (Vec3D***)Xi.GetDataPointer();
  Vec3D*** dXidx = (Vec3D***)Var1.GetDataPointer();  
  Vec3D*** dXidy = (Vec3D***)Var2.GetDataPointer();  

  double gradxi[4]; //local values of Jacobian of Xi
  double f2[4], Jloc2, r0;
  bool invertible = false;
  for(int k=k0; k<kmax; k++) {
    for(int j=j0; j<jmax; j++) {
      double r =  global_mesh.GetY(j);
      for(int i=i0; i<imax; i++) {

        for(int dim=0; dim<2; dim++)
          gradxi[dim]   = dXidx[k][j][i][dim];
        for(int dim=0; dim<2; dim++)
          gradxi[2+dim] = dXidy[k][j][i][dim];

        invertible = MathTools::LinearAlgebra::
                     CalculateMatrixInverseAndDeterminant2x2(gradxi, f2, &Jloc2);
        if(!invertible)
          fprintf(stdout,"\033[0;35mWarning: Jacobian of ref. map at (%d,%d,%d) is not invertible."
                         " determinant = %e.\033[0m\n", i,j,k, Jloc2);

        // Note: f = [dz/dZ  0  dz/dR;  0  r/R  0;  dr/dZ  0  dr/dR]; //"x = z", "y = r"!
        r0 = xi[k][j][i][1];
        double* floc(&f[k][j][9*i]);
        floc[0] = f2[0];   floc[3] = 0.0;    floc[6] = f2[2];
        floc[1] = 0.0;     floc[4] = r/r0;   floc[7] = 0.0;
        floc[2] = f2[1];   floc[5] = 0.0;    floc[8] = f2[3];

        Jloc[k][j][i] = 1.0/Jloc2*r/r0; // we want the determinant of F, not grad-xi!

        if(Jloc[k][j][i]<=0.0)
          fprintf(stdout,"\033[0;35mWarning: Determinant of deformation gradient at (%d,%d,%d)"
                         " is negative (%e).\033[0m\n", i,j,k, Jloc[k][j][i]);
      }
    }
  }

  F.RestoreDataPointerAndInsert();
  J.RestoreDataPointerAndInsert();
  Xi.RestoreDataPointerToLocalVector();
  Var1.RestoreDataPointerToLocalVector();
  Var2.RestoreDataPointerToLocalVector();

}

//------------------------------------------------------------

void
HyperelasticityOperator::ComputePrincipalStresses(SpaceVariable3D &Xi, SpaceVariable3D &V,
                                                  SpaceVariable3D &ID, SpaceVariable3D &PS)
{
  if(cylindrical_symmetry)
    ComputePrincipalStresses2DCylindrical(Xi,V,ID,PS);
  else
    ComputePrincipalStresses3D(Xi,V,ID,PS);
}

//------------------------------------------------------------

void
HyperelasticityOperator::ComputePrincipalStresses3D(SpaceVariable3D &Xi, SpaceVariable3D &V,
                                                    SpaceVariable3D &ID, SpaceVariable3D &PS)
{
  ComputeDeformGradAtNodes3D(Xi); //fill F

  double*** f  = F.GetDataPointer();
  double*** id = ID.GetDataPointer();
  Vec5D*** v   = (Vec5D***)V.GetDataPointer();
  Vec3D*** ps  = (Vec3D***)PS.GetDataPointer();

  bool success;
  int myid;
  double sigma6[6], sigma9[9]; //Cauchy stress tensor, column-first, symmetric
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        myid = id[k][j][i];
        if(myid == INACTIVE_MATERIAL_ID) {
          ps[k][j][i] = 0.0;
          continue;
        }
        assert(myid>=0 && myid<(int)hyperFcn.size());

        double* floc(&f[k][j][9*i]);
        hyperFcn[myid]->GetCauchyStressTensor(floc, v[k][j][i], sigma6);

        // Get a general 3x3 matrix
        sigma9[0] = sigma6[0];  sigma9[3] = sigma6[1];  sigma9[6] = sigma6[2];
        sigma9[1] = sigma6[1];  sigma9[4] = sigma6[3];  sigma9[7] = sigma6[4];
        sigma9[2] = sigma6[2];  sigma9[5] = sigma6[4];  sigma9[8] = sigma6[5];

        success = MathTools::LinearAlgebra::
                  CalculateEigenSymmetricMatrix3x3(sigma9, ps[k][j][i]);
        if(!success) {
          fprintf(stdout,"\033[0;31m*** Error: Unable to calculate matrix "
                         "eigenvalues in HyperelasticityOperator.\n");
          exit(-1);
        }

        std::swap(ps[k][j][i][0], ps[k][j][i][2]); // should be in descending order

/*
        if(k==75) {
          double r = sqrt(global_mesh.GetX(i)*global_mesh.GetX(i)+global_mesh.GetY(j)*global_mesh.GetY(j));
          if(fabs(r-0.205)<1e-4)         
            fprintf(stdout,"At (x,y,z) = (%e,%e,%e), r = %e, phi = %e, F = [%e %e %e; %e %e %e; %e %e %e;]. "
                    "sigma6 = %e %e %e %e %e %e. ps = %e %e %e.\n",
                    global_mesh.GetX(i), global_mesh.GetY(j), global_mesh.GetZ(k), r,
                    atan(global_mesh.GetY(j)/global_mesh.GetX(i)), floc[0], floc[3], floc[6],
                    floc[1], floc[4], floc[7], floc[2], floc[5], floc[8],
                    sigma6[0], sigma6[1], sigma6[2], sigma6[3], sigma6[4], sigma6[5],
                    ps[k][j][i][0], ps[k][j][i][1], ps[k][j][i][2]);
        }
*/
      }

  F.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();
  V.RestoreDataPointerToLocalVector();
  PS.RestoreDataPointerAndInsert();
}

//------------------------------------------------------------

void
HyperelasticityOperator::ComputePrincipalStresses2DCylindrical(SpaceVariable3D &Xi, SpaceVariable3D &V,
                                                               SpaceVariable3D &ID, SpaceVariable3D &PS)
{
  ComputeDeformGradAtNodes2DCylindrical(Xi); //fills F

  double*** f  = F.GetDataPointer();
  double*** id = ID.GetDataPointer();
  Vec5D*** v   = (Vec5D***)V.GetDataPointer();
  Vec3D*** ps  = (Vec3D***)PS.GetDataPointer();

  bool success;
  int myid = 0;
  double sigma2d[3], sigma9[9] = {0.0}, sigma_phiphi; //"sigma_2D" and \sigma_{\phi\phi}
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        myid = id[k][j][i];
        if(myid == INACTIVE_MATERIAL_ID) {
          ps[k][j][i] = 0.0;
          continue;
        }
        assert(myid>=0 && myid<(int)hyperFcn.size());

        // Get sigma_2D and sigma_phiphi
        double* floc(&f[k][j][9*i]);
        dynamic_cast<HyperelasticityFcnBase2DCyl*>
            (hyperFcn[myid])->GetCauchyStressTensor(floc, v[k][j][i], sigma2d, sigma_phiphi);

        // Get a general 3x3 matrix
        sigma9[0] = sigma2d[0];                             sigma9[6] = sigma2d[1];
                                 sigma9[4] = sigma_phiphi; 
        sigma9[2] = sigma2d[1];                             sigma9[8] = sigma2d[2];

        success = MathTools::LinearAlgebra::
                  CalculateEigenSymmetricMatrix3x3(sigma9, ps[k][j][i]);
        if(!success) {
          fprintf(stdout,"\033[0;31m*** Error: Unable to calculate matrix "
                         "eigenvalues in HyperelasticityOperator.\n");
          exit(-1);
        }

        std::swap(ps[k][j][i][0], ps[k][j][i][2]); // should be in descending order
/*
        if(i==75 && j==20)
          fprintf(stdout,"At (z,r) = (%e,%e), F = [%e %e %e; %e %e %e; %e %e %e;]. sigma2d = [%e %e; %e %e]"
                  ". sigma_phiphi = %e. ps = %e %e %e.\n",
                  global_mesh.GetX(i), global_mesh.GetY(j), floc[0], floc[3], floc[6],
                  floc[1], floc[4], floc[7], floc[2], floc[5], floc[8],
                  sigma2d[0], sigma2d[1], sigma2d[1], sigma2d[2], sigma_phiphi,
                  ps[k][j][i][0], ps[k][j][i][1], ps[k][j][i][2]);
*/

      }

  F.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();
  V.RestoreDataPointerToLocalVector();
  PS.RestoreDataPointerAndInsert();
}

//------------------------------------------------------------

void
HyperelasticityOperator::ComputePrincipalStressesAtProbes(SpaceVariable3D &Xi, SpaceVariable3D &ID,
           std::vector<Int3> &ijk, std::vector<std::pair<int, std::array<bool,8> > > &ijk_valid,
           std::vector<Vec3D> &trilinear_coords, Vec5D*** v, vector<Vec3D> &sol)
{

  int nProbes = ijk.size();
  if(nProbes == 0) {
    sol.clear();
    return; //no probes
  }

  ComputeDeformationGradientAtNodes(Xi); // fill "F"

  double*** id = ID.GetDataPointer();
  double*** f  = F.GetDataPointer();
  
  sol.assign(nProbes, 0.0);
  int dim = 5;
  for(int iProbe=0; iProbe<nProbes; iProbe++) {
    int i = ijk[iProbe][0], j = ijk[iProbe][1], k = ijk[iProbe][2];

    if(!coordinates.IsHere(i,j,k,true))
      continue;

    //NOTE: Unlike other space variables, "Xi" is only populated within the physical domain.
    //      Therefore, we should calculate stresses only at nodes within the physical domain

    //this probe node is in the current subdomain (including ghost)
    // c000
    Vec3D c000 = 0.0;
    vector<bool> valid(8,false);
    int valid_count = 0;
    if(ijk_valid[iProbe].second[0] && coordinates.IsHereOrInternalGhost(i,j,k)) {
      valid[0] = true;
      valid_count++;
      c000 = ComputePrincipalStressesAtPoint(&f[k][j][i*9], v[k][j][i*dim], id[k][j][i]);
    }
    // c100
    Vec3D c100 = 0.0;
    if(ijk_valid[iProbe].second[1] && coordinates.IsHereOrInternalGhost(i+1,j,k)) {
      valid[1] = true;
      valid_count++;
      c100 = ComputePrincipalStressesAtPoint(&f[k][j][(i+1)*9], v[k][j][(i+1)*dim], id[k][j][i+1]);
    }
    // c010
    Vec3D c010 = 0.0;
    if(ijk_valid[iProbe].second[2] && coordinates.IsHereOrInternalGhost(i,j+1,k)) {
      valid[2] = true;
      valid_count++;
      c010 = ComputePrincipalStressesAtPoint(&f[k][j+1][i*9], v[k][j+1][i*dim], id[k][j+1][i]);
    }
    // c110
    Vec3D c110 = 0.0;
    if(ijk_valid[iProbe].second[3] && coordinates.IsHereOrInternalGhost(i+1,j+1,k)) {
      valid[3] = true;
      valid_count++;
      c110 = ComputePrincipalStressesAtPoint(&f[k][j+1][(i+1)*9], v[k][j+1][(i+1)*dim], id[k][j+1][i+1]);
    }
    // c001
    Vec3D c001 = 0.0;
    if(ijk_valid[iProbe].second[4] && coordinates.IsHereOrInternalGhost(i,j,k+1)) {
      valid[4] = true;
      valid_count++;
      c001 = ComputePrincipalStressesAtPoint(&f[k+1][j][i*9], v[k+1][j][i*dim], id[k+1][j][i]);
    }
    // c101
    Vec3D c101 = 0.0;
    if(ijk_valid[iProbe].second[5] && coordinates.IsHereOrInternalGhost(i+1,j,k+1)) {
      valid[5] = true;
      valid_count++;
      c101 = ComputePrincipalStressesAtPoint(&f[k+1][j][(i+1)*9], v[k+1][j][(i+1)*dim], id[k+1][j][i+1]);
    }
    // c011
    Vec3D c011 = 0.0;
    if(ijk_valid[iProbe].second[6] && coordinates.IsHereOrInternalGhost(i,j+1,k+1)) {
      valid[6] = true;
      valid_count++;
      c011 = ComputePrincipalStressesAtPoint(&f[k+1][j+1][i*9], v[k+1][j+1][i*dim], id[k+1][j+1][i]);
    }
    // c111
    Vec3D c111 = 0.0;
    if(ijk_valid[iProbe].second[7] && coordinates.IsHereOrInternalGhost(i+1,j+1,k+1)) {
      valid[7] = true;
      valid_count++;
      c111 = ComputePrincipalStressesAtPoint(&f[k+1][j+1][(i+1)*9], v[k+1][j+1][(i+1)*dim],
                                             id[k+1][j+1][i+1]);
    }

    if(valid_count<8) {//fill invalid slots with average value
      Vec3D c_avg = (c000+c100+c010+c110+c001+c101+c011+c111)/valid_count;
      if(!valid[0])  c000 = c_avg;
      if(!valid[1])  c100 = c_avg;
      if(!valid[2])  c010 = c_avg;
      if(!valid[3])  c110 = c_avg;
      if(!valid[4])  c001 = c_avg;
      if(!valid[5])  c101 = c_avg;
      if(!valid[6])  c011 = c_avg;
      if(!valid[7])  c111 = c_avg;
    }

    sol[iProbe] = MathTools::trilinear_interpolation(trilinear_coords[iProbe], c000, c100, c010, c110,
                                                     c001, c101, c011, c111);
  }
  MPI_Allreduce(MPI_IN_PLACE, (double*)sol.data(), 3*sol.size(), MPI_DOUBLE, MPI_SUM, comm);

  F.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();

}

//------------------------------------------------------------

Vec3D
HyperelasticityOperator::ComputePrincipalStressesAtPoint(double *f, Vec5D &v, int id)
{
  if(cylindrical_symmetry)
    return ComputePrincipalStressesAtPoint2DCylindrical(f, v, id);

  return ComputePrincipalStressesAtPoint3D(f, v, id);
}

//------------------------------------------------------------

Vec3D
HyperelasticityOperator::ComputePrincipalStressesAtPoint3D(double *f, Vec5D &v, int id)
{
  Vec3D ps(0.0);

  if(id == INACTIVE_MATERIAL_ID)
    return ps;

  assert(id>=0 && id<(int)hyperFcn.size());

  double sigma6[6], sigma9[9]; //Cauchy stress tensor, column-first, symmetric
  hyperFcn[id]->GetCauchyStressTensor(f, v, sigma6);

  // Get a general 3x3 matrix
  sigma9[0] = sigma6[0];  sigma9[3] = sigma6[1];  sigma9[6] = sigma6[2];
  sigma9[1] = sigma6[1];  sigma9[4] = sigma6[3];  sigma9[7] = sigma6[4];
  sigma9[2] = sigma6[2];  sigma9[5] = sigma6[4];  sigma9[8] = sigma6[5];

  bool success = MathTools::LinearAlgebra::CalculateEigenSymmetricMatrix3x3(sigma9, ps);
  if(!success) {
    fprintf(stdout,"\033[0;31m*** Error: Unable to calculate matrix "
                   "eigenvalues in HyperelasticityOperator (1P3D).\n");
    exit(-1);
  }

  std::swap(ps[0], ps[2]); // should be in descending order

  return ps;
}

//------------------------------------------------------------

Vec3D
HyperelasticityOperator::ComputePrincipalStressesAtPoint2DCylindrical(double *f, Vec5D &v, int id)
{
  Vec3D ps(0.0);

  if(id == INACTIVE_MATERIAL_ID)
    return ps;

  assert(id>=0 && id<(int)hyperFcn.size());

  // Get sigma_2D and sigma_phiphi
  double sigma2d[3], sigma9[9] = {0.0}, sigma_phiphi; //"sigma_2D" and \sigma_{\phi\phi}
  dynamic_cast<HyperelasticityFcnBase2DCyl*>
      (hyperFcn[id])->GetCauchyStressTensor(f, v, sigma2d, sigma_phiphi);

  // Get a general 3x3 matrix
  sigma9[0] = sigma2d[0];                             sigma9[6] = sigma2d[1];
                           sigma9[4] = sigma_phiphi; 
  sigma9[2] = sigma2d[1];                             sigma9[8] = sigma2d[2];

  bool success = MathTools::LinearAlgebra::CalculateEigenSymmetricMatrix3x3(sigma9, ps);
  if(!success) {
    fprintf(stdout,"\033[0;31m*** Error: Unable to calculate matrix "
                   "eigenvalues in HyperelasticityOperator (1P2DCyl).\n");
    exit(-1);
  }

  std::swap(ps[0], ps[2]); // should be in descending order

  return ps;
}

//------------------------------------------------------------
//Note: Fluxes are added on the left-hand-side of the Navier-Stokes equations
void
HyperelasticityOperator::AddHyperelasticityFluxes(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &Xi,
                                                  vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS,
                                                  SpaceVariable3D &R)
{
  if(cylindrical_symmetry) {
    AddFluxes2DCylindrical(V,ID,Xi,EBDS,R);
    AddCylindricalSourceTerms(V,ID,Xi,EBDS,R);
  } else
    AddFluxes3D(V,ID,Xi,EBDS,R);
}

//------------------------------------------------------------
//Note: Fluxes are added on the left-hand-side of the Navier-Stokes equations
void
HyperelasticityOperator::AddFluxes3D(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &Xi,
                                     vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS,
                                     SpaceVariable3D &R)
{
  if(EBDS)
    print_warning("Warning: AddHyperelasticityFluxes: Not able to account for embedded surfaces.\n");

  std::vector<int> i123{1,2,3}, i012{0,1,2};

  //1. Calculate Xi derivatives at cell interfaces
  grad.CalculateFirstDerivativeAtCellInterfaces(0, 0, Xi, i012, dXidx_i_minus_half, i012);
  grad.CalculateFirstDerivativeAtCellInterfaces(0, 1, Xi, i012, dXidx_j_minus_half, i012);
  grad.CalculateFirstDerivativeAtCellInterfaces(0, 2, Xi, i012, dXidx_k_minus_half, i012);
  grad.CalculateFirstDerivativeAtCellInterfaces(1, 0, Xi, i012, dXidy_i_minus_half, i012);
  grad.CalculateFirstDerivativeAtCellInterfaces(1, 1, Xi, i012, dXidy_j_minus_half, i012);
  grad.CalculateFirstDerivativeAtCellInterfaces(1, 2, Xi, i012, dXidy_k_minus_half, i012);
  grad.CalculateFirstDerivativeAtCellInterfaces(2, 0, Xi, i012, dXidz_i_minus_half, i012);
  grad.CalculateFirstDerivativeAtCellInterfaces(2, 1, Xi, i012, dXidz_j_minus_half, i012);
  grad.CalculateFirstDerivativeAtCellInterfaces(2, 2, Xi, i012, dXidz_k_minus_half, i012);

  //2. Loop through cell interfaces and calculate fluxes
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
  double dx = 0.0, dy = 0.0, dz = 0.0;
  Vec5D flux;
  bool invertible;
  for(int k=k0; k<kkmax; k++) {
    dz = global_mesh.GetDz(k);
    for(int j=j0; j<jjmax; j++) {
      dy = global_mesh.GetDy(j);
      for(int i=i0; i<iimax; i++) {

        dx   = global_mesh.GetDx(i);
        myid = id[k][j][i];

        if(myid == INACTIVE_MATERIAL_ID)
          continue;
        assert(myid>=0 && myid<(int)hyperFcn.size());

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
                       CalculateMatrixInverseAndDeterminant3x3(gradxi, f);
          if(!invertible)
            fprintf(stdout,"\033[0;35mWarning: Jacobian of ref. map at (%d-1/2,%d,%d)"
                           " is not invertible.\033[0m\n", i,j,k);

          hyperFcn[myid]->EvaluateHyperelasticFluxFunction_F(flux, f, v[k][j][i],
                                                             deviator_only[myid]);//TODO: multi-material
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
                       CalculateMatrixInverseAndDeterminant3x3(gradxi, f);
          if(!invertible)
            fprintf(stdout,"\033[0;35mWarning: Jacobian of ref. map at (%d,%d-1/2,%d)"
                           " is not invertible.\033[0m\n", i,j,k);

          hyperFcn[myid]->EvaluateHyperelasticFluxFunction_G(flux, f, v[k][j][i],
                                                             deviator_only[myid]);//TODO: multi-material
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
                       CalculateMatrixInverseAndDeterminant3x3(gradxi, f);
          if(!invertible)
            fprintf(stdout,"\033[0;35mWarning: Jacobian of ref. map at (%d,%d,%d-1/2)"
                           " is not invertible.\033[0m\n", i,j,k);

          hyperFcn[myid]->EvaluateHyperelasticFluxFunction_H(flux, f, v[k][j][i],
                                                             deviator_only[myid]);//TODO: multi-material
          flux *= dx*dy;
          res[k][j][i]   += flux;
          res[k-1][j][i] -= flux;
        }
      }
    }
  }

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
//Note: Fluxes are added on the left-hand-side of the Navier-Stokes equations
void
HyperelasticityOperator::AddFluxes2DCylindrical(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &Xi,
                                                vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS,
                                                SpaceVariable3D &R)
{
  if(EBDS)
    print_warning("Warning: AddHyperelasticityFluxes: Not able to account for embedded surfaces.\n");

  std::vector<int> i01{0,1};

  //1. Calculate Xi derivatives at cell interfaces
  grad.CalculateFirstDerivativeAtCellInterfaces(0, 0, Xi, i01, dXidx_i_minus_half, i01);
  grad.CalculateFirstDerivativeAtCellInterfaces(0, 1, Xi, i01, dXidx_j_minus_half, i01);
  grad.CalculateFirstDerivativeAtCellInterfaces(1, 0, Xi, i01, dXidy_i_minus_half, i01);
  grad.CalculateFirstDerivativeAtCellInterfaces(1, 1, Xi, i01, dXidy_j_minus_half, i01);

  //2. Loop through cell interfaces and calculate fluxes
  Vec3D*** dxidx_i = (Vec3D***)dXidx_i_minus_half.GetDataPointer();
  Vec3D*** dxidx_j = (Vec3D***)dXidx_j_minus_half.GetDataPointer();
  Vec3D*** dxidy_i = (Vec3D***)dXidy_i_minus_half.GetDataPointer();
  Vec3D*** dxidy_j = (Vec3D***)dXidy_j_minus_half.GetDataPointer();

  Vec3D*** xi  = (Vec3D***)Xi.GetDataPointer();
  Vec5D*** v   = (Vec5D***)V.GetDataPointer();
  Vec5D*** res = (Vec5D***)R.GetDataPointer();
  double*** id = (double***)ID.GetDataPointer();

  int myid = 0;
  double gradxi[4], f2[4], f[9] = {0.0}; //nabla xi and deformation gradient
  double dx = 0.0, dy = 0.0, dz = 0.0, r0;
  Vec5D flux;
  bool invertible;
  for(int k=k0; k<kkmax; k++) {
    dz = global_mesh.GetDz(k);
    for(int j=j0; j<jjmax; j++) {
      dy = global_mesh.GetDy(j);
      for(int i=i0; i<iimax; i++) {

        dx   = global_mesh.GetDx(i);
        myid = id[k][j][i];

        if(myid == INACTIVE_MATERIAL_ID)
          continue;
        assert(myid>=0 && myid<(int)hyperFcn.size());

        //*****************************************
        //calculate flux function F_{i-1/2,j,k}
        //*****************************************
        if(k!=kkmax-1 && j!=jjmax-1) {

          for(int dim=0; dim<2; dim++)
            gradxi[dim]   = dxidx_i[k][j][i][dim];
          for(int dim=0; dim<2; dim++)
            gradxi[2+dim] = dxidy_i[k][j][i][dim];

          invertible = MathTools::LinearAlgebra::
                       CalculateMatrixInverseAndDeterminant2x2(gradxi, f2);
          if(!invertible)
            fprintf(stdout,"\033[0;35mWarning: Jacobian of ref. map at (%d-1/2,%d,%d)"
                           " is not invertible.\033[0m\n", i,j,k);

          // Note: f = [dz/dZ  0  dz/dR;  0  r/R  0;  dr/dZ  0  dr/dR]; //"x = z", "y = r"!
          double r =  global_mesh.GetY(j);
          r0 = (global_mesh.GetDx(i-1)*xi[k][j][i][1] + dx*xi[k][j][i-1][1])
             / (global_mesh.GetDx(i-1) + dx);
          f[0] = f2[0];                    f[6] = f2[2];
                           f[4] = r/r0;    
          f[2] = f2[1];                    f[8] = f2[3];

          hyperFcn[myid]->EvaluateHyperelasticFluxFunction_F(flux, f, v[k][j][i],
                                                             deviator_only[myid]);//TODO: multi-material
          flux *= dy*dz;
          res[k][j][i]   += flux;
          res[k][j][i-1] -= flux;

        }


        //*****************************************
        //calculate flux function G_{i,j-1/2,k}
        //*****************************************
        if(i!=iimax-1 && k!=kkmax-1) {

          for(int dim=0; dim<2; dim++)
            gradxi[dim]   = dxidx_j[k][j][i][dim];
          for(int dim=0; dim<2; dim++)
            gradxi[2+dim] = dxidy_j[k][j][i][dim];

          invertible = MathTools::LinearAlgebra::
                       CalculateMatrixInverseAndDeterminant2x2(gradxi, f2);
          if(!invertible)
            fprintf(stdout,"\033[0;35mWarning: Jacobian of ref. map at (%d,%d-1/2,%d)"
                           " is not invertible.\033[0m\n", i,j,k);

          // Note: f = [dz/dZ  0  dz/dR;  0  r/R  0;  dr/dZ  0  dr/dR]; //"x = z", "y = r"!
          double r =  global_mesh.GetY(j) - 0.5*dy;
          r0 = (global_mesh.GetDy(j-1)*xi[k][j][i][1] + dy*xi[k][j-1][i][1])
             / (global_mesh.GetDy(j-1) + dy);
          f[0] = f2[0];                    f[6] = f2[2];
                           f[4] = r/r0;    
          f[2] = f2[1];                    f[8] = f2[3];

          hyperFcn[myid]->EvaluateHyperelasticFluxFunction_G(flux, f, v[k][j][i],
                                                             deviator_only[myid]);//TODO: multi-material
          flux *= dx*dz;
          res[k][j][i]   += flux;
          res[k][j-1][i] -= flux;
        }

      }
    }
  }

  dXidx_i_minus_half.RestoreDataPointerToLocalVector();
  dXidx_j_minus_half.RestoreDataPointerToLocalVector();
  dXidy_i_minus_half.RestoreDataPointerToLocalVector();
  dXidy_j_minus_half.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();
  V.RestoreDataPointerToLocalVector();
  Xi.RestoreDataPointerToLocalVector();

  R.RestoreDataPointerToLocalVector(); //NOTE: although R has been updated, there is no need of
                                       //      cross-subdomain communications. So, no need to
                                       //      update the global vec.
}

//------------------------------------------------------------
// Note: Source terms are placed on the left-hand-side of the Navier-Stokes equations
void
HyperelasticityOperator::AddCylindricalSourceTerms(SpaceVariable3D &V, SpaceVariable3D &ID,
                                                   SpaceVariable3D &Xi,
                                                   vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS,
                                                   SpaceVariable3D &R)
{
  if(EBDS)
    print_warning("Warning: AddHyperelasticityFluxes: Not able to account for embedded surfaces.\n");


  // NOTE: This function computes deformation gradient (F) ONLY WITHIN THE PHYSICAL DOMAIN,
  //       NOT in the ghost boundary layer.

  // ------------------------------------
  // Step 1: Calculate the Jacobian of Xi
  // ------------------------------------
  std::vector<int> i01{0,1};    
  // dxi/dx (i.e., dxi/dz in cylindrical coords)
  grad.CalculateFirstDerivativeAtNodes(0 /*dx*/, Xi, i01, Var1, i01);
  // dxi/dy (i.e., dxi/dr in cylindrical coords)
  grad.CalculateFirstDerivativeAtNodes(1 /*dy*/, Xi, i01, Var2, i01);


  // ------------------------------------
  // Step 2: Calculate the deformation gradient: F = inv(grad.Xi), its determinant.
  //         Then, evaluate the source terms due to cylindrical symmetry.
  // ------------------------------------

  Vec3D*** xi    = (Vec3D***)Xi.GetDataPointer();
  Vec3D*** dXidx = (Vec3D***)Var1.GetDataPointer();  
  Vec3D*** dXidy = (Vec3D***)Var2.GetDataPointer();  
  double*** id   = (double***)ID.GetDataPointer();
  Vec5D*** v     = (Vec5D***)V.GetDataPointer();
  Vec5D*** res   = (Vec5D***)R.GetDataPointer();
  double*** cv   = (double***)volume.GetDataPointer();

  int myid = 0;
  double gradxi[4]; //local values of Jacobian of Xi
  double f2[4], f[9] = {0.0}; //deform. grad
  double sigma[3], sigma_phiphi; //"sigma_2D" and \sigma_{\phi\phi}
  double Jloc2, r0, mycv;
  bool invertible = false;
  for(int k=k0; k<kmax; k++) {
    for(int j=j0; j<jmax; j++) {
      double r =  global_mesh.GetY(j);
      for(int i=i0; i<imax; i++) {

        myid = id[k][j][i];

        if(myid == INACTIVE_MATERIAL_ID)
          continue;
        assert(myid>=0 && myid<(int)hyperFcn.size());

        for(int dim=0; dim<2; dim++)
          gradxi[dim]   = dXidx[k][j][i][dim];
        for(int dim=0; dim<2; dim++)
          gradxi[2+dim] = dXidy[k][j][i][dim];

        invertible = MathTools::LinearAlgebra::
                     CalculateMatrixInverseAndDeterminant2x2(gradxi, f2, &Jloc2);
        if(!invertible || Jloc2<=0.0)
          fprintf(stdout,"\033[0;35mWarning: Jacobian of ref. map at (%d,%d,%d) is invalid."
                         " determinant = %e.\033[0m\n", i,j,k, Jloc2);
 
        // Note: f = [dz/dZ  0  dz/dR;  0  r/R  0;  dr/dZ  0  dr/dR]; //"x = z", "y = r"!
        r0 = xi[k][j][i][1];
        f[0] = f2[0];                f[6] = f2[2];
                       f[4] = r/r0;  
        f[2] = f2[1];                f[8] = f2[3];

        // Get sigma_2D and sigma_phiphi
        dynamic_cast<HyperelasticityFcnBase2DCyl*>
            (hyperFcn[myid])->GetCauchyStressTensor(f, v[k][j][i], sigma, sigma_phiphi);
        if(deviator_only[myid]) {                                                   
          double p = -1.0/3.0*(sigma[0] + sigma[2] + sigma_phiphi); //hydrostatic pressure
          sigma[0] += p;
          sigma[2] += p;
          sigma_phiphi += p;                 
        }
        
        mycv = cv[k][j][i];
        res[k][j][i][1] -= mycv/r*sigma[1]; //1/r*sigma_rz*cv
        res[k][j][i][2] -= mycv/r*(sigma[2] - sigma_phiphi); //1/r*(sigma_rr-sigma_phiphi)*cv
        res[k][j][i][4] -= mycv/r*(sigma[2]*v[k][j][i][2] + sigma[1]*v[k][j][i][1]);

      }
    }
  }

  Xi.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();
  V.RestoreDataPointerToLocalVector();
  volume.RestoreDataPointerToLocalVector();
  Var1.RestoreDataPointerToLocalVector();
  Var2.RestoreDataPointerToLocalVector();

  R.RestoreDataPointerToLocalVector(); //NOTE: although R has been updated, there is no need of
                                       //      cross-subdomain communications. So, no need to
                                       //      update the global vec.
}


//------------------------------------------------------------

void
HyperelasticityOperator::PrescribeVelocityForTesting([[maybe_unused]] SpaceVariable3D &V,
                                                     [[maybe_unused]] double time)
{

#if HYPERELASTICITY_TEST == 1
  double pi = acos(-1.0);
  double dmax = 0.5, Rmax = 1.0;
  Vec5D*** v = (Vec5D***)V.GetDataPointer(); 

  if(cylindrical_symmetry) {
    for(int k=k0; k<kmax; k++) {
      for(int j=j0; j<jmax; j++) {
        double r = global_mesh.GetY(j);
        for(int i=i0; i<imax; i++) {
          v[k][j][i][1] = r<=Rmax ? dmax*sin(pi*r/(2.0*Rmax)) : 0.0;
          v[k][j][i][2] = 0.0; //r-velocity
          v[k][j][i][3] = 0.0; //unused
        }
      }
    }
  }
  else { //true 3D
    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++) {
          double r = sqrt(pow(global_mesh.GetX(i),2) + pow(global_mesh.GetY(j),2));
          v[k][j][i][1] = 0.0; //x-velocity
          v[k][j][i][2] = 0.0; //y-velocity
          v[k][j][i][3] = r<=Rmax ? dmax*sin(pi*r/(2.0*Rmax)) : 0.0; //z-velocity
        }
  }

  V.RestoreDataPointerAndInsert();

#endif

}

//------------------------------------------------------------


//------------------------------------------------------------

