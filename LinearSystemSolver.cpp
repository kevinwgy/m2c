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


  rnorm_history.resize(1000); //1000 entries should be more than enough
  KSPSetResidualHistory(rnorm_history.data(), rnorm_history.size(), true); //reset for each Solve

/*
  PC pc;
  KSPGetPC(ksp, &pc);
  PCType pctype;
  PCGetType(pc, &pctype);
  std::cout << "Precondition: " << pctype << std::endl;
*/
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

void
LinearSystemSolver::UsePreviousPreconditioner(bool reuse_or_not)
{
  KSPSetReusePreconditioner(ksp, reuse_or_not);
}

//-----------------------------------------------------

bool
LinearSystemSolver::Solve(SpaceVariable3D &b, SpaceVariable3D &x,
                          ConvergenceReason *reason, int *numIts, std::vector<double> *rnorm)
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
  

  // ---------------------------------------------------
  // Solve!
  Vec &bb(b.GetRefToGlobalVec());
  Vec &xx(x.GetRefToGlobalVec());
  KSPSolve(ksp, bb, xx);
  // ---------------------------------------------------


  KSPConvergedReason ksp_code; //positive if convergence; negative if diverged
  KSPGetConvergedReason(ksp, &ksp_code);
  bool success = ksp_code>0;

  if(reason) { //user requested convergence/divergence reason
    ConvergenceReason reason(NONE);
    if(ksp_code == KSP_CONVERGED_RTOL)
      *reason = CONVERGED_REL_TOL;
    else if(ksp_code == KSP_CONVERGED_ATOL)
      *reason = CONVERGED_ABS_TOL;
    else if((int)ksp_code>0)
      *reason = CONVERGED_OTHER;
    else if(ksp_code == KSP_DIVERGED_ITS)
      *reason = DIVERGED_ITS;
    else if(ksp_code == KSP_DIVERGED_DTOL)
      *reason = DIVERGED_DTOL;
    else
      *reason = DIVERGED_OTHER;
  }


  if(numIts) //user requested output of number of iterations
    KSPGetIterationNumber(ksp, numIts);

  if(rnorm) {//user requested output of residual norm history
    int nEntries(0);
    KSPGetResidualHistory(ksp, NULL, &nEntries);
    rnorm->resize(nEntries, 0);
    for(int i=0; i<nEntries) //copy data instead of passing rnorm_history (safer)
      (*rnorm)[i] = rnorm_history[i];
  }

  return success;
}


//-----------------------------------------------------

//-----------------------------------------------------