#include<HyperelasticityOperator.h>
#include<linear_algebra.h>

//------------------------------------------------------------

HyperelasticityOperator::HyperelasticityOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, 
                             IoData &iod_, SpaceVariable3D &coordinates_,
                             SpaceVariable3D &delta_xyz_, GlobalMeshInfo &global_mesh_,
                             InterpolatorBase &interpolator_, GradientCalculatorBase &grad_,
                             std::vector<GhostPoint> &ghost_nodes_inner_,
                             std::vector<GhostPoint> &ghost_nodes_outer_)
                       : comm(comm_), iod(iod_), global_mesh(global_mesh_),
                         interpolator(interpolator_), grad(grad_),
                         refmap(comm_, dm_all_, iod_, coordinates_, delta_xyz_,
                                global_mesh_, ghost_nodes_inner_, ghost_nodes_outer_),
                         F(comm_, &(dm_all_.ghosted1_9dof)),
                         DetF(comm_, &(dm_all_.ghosted1_1dof)),
                         Var1(comm_, &(dm_all_.ghosted1_3dof)),
                         Var2(comm_, &(dm_all_.ghosted1_3dof)),
                         Var3(comm_, &(dm_all_.ghosted1_3dof))
                         
{
  coordinates_.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates_.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);
}

//------------------------------------------------------------

HyperelasticityOperator::~HyperelasticityOperator()
{ }

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
HyperelasticityOperator::ComputeDeformationGradient(SpaceVariable3D &Xi)
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
          fprintf(stderr,"\033[0;35mWarning: Jacobian of ref. map at (%d,%d,%d) is not invertible."
                         " determinant = %e.\033[0m\n", i,j,k, detf[k][j][i]);

        detf[k][j][i] = 1.0/detf[k][j][i]; // we want the determinant of F, not grad-xi!

        if(detf[k][j][i]<=0.0)
          fprintf(stderr,"\033[0;35mWarning: Determinant of deformation gradient at (%d,%d,%d)"
                         " is negative (%e).\033[0m\n", i,j,k, detf[k][j][i]);
      }

  F.RestoreDataPointerAndInsert();
  DetF.RestoreDataPointerAndInsert();
  Var1.RestoreDataPointerToLocalVector();
  Var2.RestoreDataPointerToLocalVector();
  Var3.RestoreDataPointerToLocalVector();

}

//------------------------------------------------------------

