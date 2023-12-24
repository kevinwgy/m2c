/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<LinearOperator.h>
#include<cassert>

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

  MatZeroEntries(A); //!< zeros all entries but retains the previous nonzero structure

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

double
LinearOperator::CalculateMatrixOneNorm()
{
  double norm(0.0);
  MatNorm(A, NORM_1, &norm);
  return norm;
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
