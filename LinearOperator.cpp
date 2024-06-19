/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<LinearOperator.h>
#include<petscksp.h> //for computing singular values & condition number
#include<cassert>
#include<cfloat>

//-----------------------------------------------------

LinearOperator::LinearOperator(MPI_Comm &comm_, DM &dm_)
              : comm(comm_)
{
  DMClone(dm_, &dm);
  DMSetMatType(dm, MATAIJ);
  DMSetMatrixPreallocateOnly(dm, PETSC_TRUE);
  DMCreateMatrix(dm, &A);

  // -------------------------------------------------------
  // Get info about the domain decomposition
  DMDAGetInfo(dm, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &dof, NULL, NULL, NULL, NULL, NULL);
 
  int nx(0), ny(0), nz(0);
  DMDAGetCorners(dm, &i0, &j0, &k0, &nx, &ny, &nz);
  imax = i0 + nx;
  jmax = j0 + ny;
  kmax = k0 + nz;

  int ghost_nx(0), ghost_ny(0), ghost_nz(0);
  DMDAGetGhostCorners(dm, &ii0, &jj0, &kk0, &ghost_nx, &ghost_ny, &ghost_nz);
  iimax = ii0 + ghost_nx;
  jjmax = jj0 + ghost_ny;
  kkmax = kk0 + ghost_nz;
  // -------------------------------------------------------
}

//-----------------------------------------------------

LinearOperator::~LinearOperator()
{ }

//-----------------------------------------------------

void
LinearOperator::Destroy()
{
  MatDestroy(&A);
  DMDestroy(&dm);
}

//-----------------------------------------------------

void
LinearOperator::SetLinearOperator(vector<RowEntries>& row_entries)
{

  //MatZeroEntries(A); //!< zeros all entries but retains the previous nonzero structure (KW: NOT WORKING?)
  MatDestroy(&A);
  DMCreateMatrix(dm, &A);

  for(auto&& entries : row_entries)
    MatSetValuesStencil(A, 1, &entries.row, entries.cols.size(), entries.cols.data(),
                        entries.vals.data(), ADD_VALUES);

  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

}

//-----------------------------------------------------

void
LinearOperator::ApplyLinearOperator(SpaceVariable3D &x, SpaceVariable3D &y)
{
  Vec &xx(x.GetRefToGlobalVec());
  Vec &yy(y.GetRefToGlobalVec());
  assert(&xx != &yy); //cannot be the same vector

  MatMult(A, xx, yy); 

  y.SyncLocalToGlobal(); //update the localVec of y 
}

//-----------------------------------------------------

void
LinearOperator::ApplyLinearOperatorAndAdd(SpaceVariable3D &x, SpaceVariable3D &b,
                                          SpaceVariable3D &y)
{
  Vec &xx(x.GetRefToGlobalVec());
  Vec &bb(b.GetRefToGlobalVec());
  Vec &yy(y.GetRefToGlobalVec());
  assert(&xx != &yy); //cannot be the same vector

  MatMultAdd(A, xx, bb, yy); 

  y.SyncLocalToGlobal(); //update the localVec of y 
}

//-----------------------------------------------------

void
LinearOperator::CalculateExtremeSingularValues(double &lambda_max, double &lambda_min)
{
  KSP ksp_tmp; //create a temporary Krylov-subspace based solver
  KSPCreate(comm, &ksp_tmp);
  KSPSetComputeSingularValues(ksp_tmp, PETSC_TRUE);
  PC pc;
  KSPGetPC(ksp_tmp, &pc);
  PCSetType(pc, PCNONE);
  KSPGMRESSetRestart(ksp_tmp, 1000000); //disable restart
  KSPSetOperators(ksp_tmp, A, A);
  KSPSetTolerances(ksp_tmp, 1.0e-8, 1.0e-8, PETSC_DEFAULT, 3000);
  KSPSetUp(ksp_tmp);

  //singular values are estimated while solving a linear system...
  Vec x,b;
  DMCreateGlobalVector(dm, &x);
  DMCreateGlobalVector(dm, &b);
  VecSet(x, 1.0);
  VecSet(b, 2.0);
  KSPSolve(ksp_tmp, b, x);

  KSPComputeExtremeSingularValues(ksp_tmp, &lambda_max, &lambda_min); //min can be quite inaccurate

  //cleanup
  VecDestroy(&x);
  VecDestroy(&b);
  KSPDestroy(&ksp_tmp);
}

//-----------------------------------------------------

double
LinearOperator::EstimateConditionNumber()
{
  double lambda_max(1.0), lambda_min(1.0);
  CalculateExtremeSingularValues(lambda_max, lambda_min);
  if(lambda_min == 0.0)
    return DBL_MAX;
  return lambda_max/lambda_min;
}

//-----------------------------------------------------

double
LinearOperator::CalculateMatrixOneNorm()
{
  double norm(0.0);
  MatNorm(A, NORM_1, &norm);
  return norm;
}

//-----------------------------------------------------

double
LinearOperator::CalculateMatrixTwoNorm()
{
  double lambda_max(0.0), lambda_min(0.0);
  CalculateExtremeSingularValues(lambda_max, lambda_min);
  return lambda_max; //matrix 2-norm is just max singular value
}

//-----------------------------------------------------

double
LinearOperator::CalculateMatrixInfNorm()
{
  double norm(0.0);
  MatNorm(A, NORM_INFINITY, &norm);
  return norm;
}

//-----------------------------------------------------

double
LinearOperator::CalculateMatrixFrobeniusNorm()
{
  double norm(0.0);
  MatNorm(A, NORM_FROBENIUS, &norm);
  return norm;
}

//-----------------------------------------------------

bool
LinearOperator::IsSymmetric(double tol)
{
  assert(tol>=0.0);

  PetscBool flg(PETSC_FALSE);
  MatIsSymmetric(A, tol, &flg);

  return flg==PETSC_TRUE;
}

//-----------------------------------------------------

void
LinearOperator::SetOutputVariableName(const char *name)
{
  PetscObjectSetName((PetscObject)A, name);
}

//-----------------------------------------------------

void
LinearOperator::WriteToMatlabFile(const char *filename, const char *varname)
{
  if(varname)
    SetOutputVariableName(varname);

  PetscViewer viewer;
  int code = PetscViewerASCIIOpen(comm, filename, &viewer);
  if(code) {
    print_error("*** Error: Cannot open file '%s' for output. (code: %d)\n", filename, code);
    exit_mpi();
  }

  PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
  MatView(A, viewer);

  PetscViewerDestroy(&viewer);
  MPI_Barrier(comm);
}


//-----------------------------------------------------
