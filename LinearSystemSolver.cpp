/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<LinearSystemSolver.h>
#include<cassert>
#include<iostream>

//-----------------------------------------------------

LinearSystemSolver::LinearSystemSolver(MPI_Comm &comm_, DM &dm_, LinearSolverData &lin_input)
                  : LinearOperator(comm_, dm_)
{

  KSPCreate(comm, &ksp);
  KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); //!< initial guess is passed to KSPSolve

  SetTolerances(lin_input);

  if(lin_input.ksp == LinearSolverData::KSP_DEFAULT) {
    /* nothing to do */
  } else if(lin_input.ksp == LinearSolverData::GMRES) {
    KSPSetType(ksp, KSPGMRES);
  } else if(lin_input.ksp == LinearSolverData::FLEXIBLE_GMRES) {
    KSPSetType(ksp, KSPFGMRES);
  } else {
    print_error("*** Error: Detected unknown PETSc KSP type.\n");
    exit_mpi();
  }

  if(lin_input.pc == LinearSolverData::PC_DEFAULT) {
    /* nothing to do*/
  } 
  else {

    PC pc;
    KSPGetPC(ksp, &pc);
    //PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);

    if(lin_input.pc == LinearSolverData::PC_NONE)
      PCSetType(pc, PCNONE);
    else if(lin_input.pc == LinearSolverData::JACOBI)
      PCSetType(pc, PCJACOBI);
    else if(lin_input.pc == LinearSolverData::INCOMPLETE_LU)
      PCSetType(pc, PCILU);
    else if(lin_input.pc == LinearSolverData::INCOMPLETE_CHOLESKY)
      PCSetType(pc, PCICC);
    else if(lin_input.pc == LinearSolverData::MG)
      PCSetType(pc, PCMG);
    else { 
      print_error("*** Error: Detected unknown PETSc KSP preconditioner type.\n");
      exit_mpi();
    }
  }

  if(strcmp(lin_input.options_file, ""))
    PetscOptionsInsert(NULL, NULL, NULL, lin_input.options_file);

  KSPSetFromOptions(ksp); //overrides any options specified above


  PC pc;
  KSPGetPC(ksp, &pc);
  PCType pctype;
  PCGetType(pc, &pctype);
  std::cout << "Precondition: " << pctype << std::endl;
}

//-----------------------------------------------------

LinearSystemSolver::~LinearSystemSolver()
{ }

//-----------------------------------------------------

void
LinearSystemSolver::Destroy()
{
  KSPDestroy(&ksp);
  LinearOperator::Destroy();
}

//-----------------------------------------------------

void
LinearSystemSolver::SetTolerances(LinearSolverData &lin_input) 
{
  double relative_error = lin_input.rtol;
  double absolute_error = lin_input.abstol;
  double divergence_tol = lin_input.dtol;
  int    max_iterations = lin_input.maxits;

  KSPSetTolerances(ksp,
                   relative_error>0 ? relative_error : PETSC_DEFAULT,
                   absolute_error>0 ? absolute_error : PETSC_DEFAULT,
                   divergence_tol>0 ? divergence_tol : PETSC_DEFAULT,
                   max_iterations>0 ? max_iterations : PETSC_DEFAULT);
}

//-----------------------------------------------------

void
LinearSystemSolver::GetTolerances(double *rtol, double *abstol, double *dtol, int *maxits)
{
  KSPGetTolerances(ksp, rtol, abstol, dtol, maxits);
}

//-----------------------------------------------------

void
LinearSystemSolver::SetLinearOperator(vector<RowEntries>& row_entries)
{
  LinearOperator::SetLinearOperator(row_entries); //build A
  KSPSetOperators(ksp, A, A);
}

//-----------------------------------------------------

int
LinearSystemSolver::Solve(SpaceVariable3D &b, SpaceVariable3D &x)
{
  // --------------------------------------------------
  // Sanity checks
  int dof_ = b.NumDOF();
  assert(dof_ == dof);
  dof_ = x.NumDOF();
  assert(dof_ == dof);

  int i0_, j0_, k0_, imax_, jmax_, kmax_;
  b.GetCornerIndices(&i0_, &j0_, &k0_, &imax_, &jmax_, &kmax_);
  assert(i0_==i0 && j0_==j0 && k0_==k0 && imax_==imax && jmax_==jmax && kmax_==kmax);
  x.GetCornerIndices(&i0_, &j0_, &k0_, &imax_, &jmax_, &kmax_);
  assert(i0_==i0 && j0_==j0 && k0_==k0 && imax_==imax && jmax_==jmax && kmax_==kmax);

  b.GetGhostedCornerIndices(&i0_, &j0_, &k0_, &imax_, &jmax_, &kmax_);
  assert(i0_==ii0 && j0_==jj0 && k0_==kk0 && imax_==iimax && jmax_==jjmax && kmax_==kkmax);
  x.GetGhostedCornerIndices(&i0_, &j0_, &k0_, &imax_, &jmax_, &kmax_);
  assert(i0_==ii0 && j0_==jj0 && k0_==kk0 && imax_==iimax && jmax_==jjmax && kmax_==kkmax);
  // ---------------------------------------------------
  

  Vec &bb(b.GetRefToGlobalVec());
  Vec &xx(x.GetRefToGlobalVec());

  PetscErrorCode error_code = KSPSolve(ksp, bb, xx);

  return (int)error_code; //0 means no error
}


//-----------------------------------------------------

//-----------------------------------------------------
