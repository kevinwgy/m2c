/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<LinearSystemSolver.h>
#include<cassert>
#include<iostream>

//-----------------------------------------------------

LinearSystemSolver::LinearSystemSolver(MPI_Comm &comm_, DM &dm_, LinearSolverData &lin_input,
                                       const char *equation_name_)
                  : LinearOperator(comm_, dm_), log_filename(lin_input.logfile)
{

  KSPCreate(comm, &ksp);
  KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); //!< initial guess is passed to KSPSolve

  SetTolerances(lin_input);

  if(lin_input.ksp == LinearSolverData::PETSC_KSP_DEFAULT) {
    /* nothing to do */
  } else if(lin_input.ksp == LinearSolverData::FLEXIBLE_GMRES) {
    KSPSetType(ksp, KSPFGMRES);
  } else if(lin_input.ksp == LinearSolverData::STAB_BI_CG) {
    KSPSetType(ksp, KSPBCGSL);
  } else if(lin_input.ksp == LinearSolverData::IMPROVED_STAB_BI_CG) {
    KSPSetType(ksp, KSPIBCGS);
  } else {
    print_error("*** Error: Detected unknown PETSc KSP type.\n");
    exit_mpi();
  }

  if(lin_input.pc == LinearSolverData::PETSC_PC_DEFAULT) {
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
    else if(lin_input.pc == LinearSolverData::BLOCK_JACOBI)
      PCSetType(pc, PCBJACOBI);
    else if(lin_input.pc == LinearSolverData::MG)
      PCSetType(pc, PCMG);
    else { 
      print_error("*** Error: Detected unknown PETSc KSP preconditioner type.\n");
      exit_mpi();
    }
  }

  if(strcmp(lin_input.options_file, "")) {
    FILE *file = fopen(lin_input.options_file, "r");
    if(file == NULL) {
      print_error("*** Error: Cannot open PETSc options file %s.\n", lin_input.options_file);
      exit_mpi();
    }
    fclose(file); // just make sure the file can be opened...
    PetscOptionsInsert(NULL, NULL, NULL, lin_input.options_file);
  }

  KSPSetFromOptions(ksp); //overrides any options specified above

  rnorm_history.resize(1000); //1000 entries should be more than enough
  KSPSetResidualHistory(ksp, rnorm_history.data(), rnorm_history.size(), PETSC_TRUE); //reset for each Solve

  // Set up log file
  equation_name = equation_name_;
  log_filename = lin_input.logfile;
  if(!log_filename.empty()) {
    FILE *file = fopen(log_filename.c_str(), "w");
    if(file == NULL) {
      print_error("*** Error: Unable to open file %s for printing the log of linear system solver.\n");
      exit_mpi();
    }
    if(equation_name.empty())
      print(file, "## Solution of linear systems - Log file.\n");
    else
      print(file, "## Solution of linear systems (%s) - Log file.\n", equation_name.c_str());
  
    fclose(file);
    MPI_Barrier(comm);
  }
  write_log_to_screen = lin_input.write_log_to_screen == LinearSolverData::YES;

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
LinearSystemSolver::GetSolverType(string *ksp_type, string *pc_type)
{
  if(ksp_type) {
    KSPType ksp_type_;
    KSPGetType(ksp, &ksp_type_);
    *ksp_type = (string)ksp_type_;
  }

  if(pc_type) {
    PC pc;
    KSPGetPC(ksp, &pc);
    PCType pc_type_;
    PCGetType(pc, &pc_type_);
    *pc_type = (string)pc_type_;
  }
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
LinearSystemSolver::ComputeResidual(SpaceVariable3D &b, SpaceVariable3D &x,
                                    SpaceVariable3D &res)
{
  ApplyLinearOperator(x, res); //res = Ax 
  res.AXPlusBY(-1.0, 1.0, b); //res = -1.0*res + b = b - Ax
}

//-----------------------------------------------------

void
LinearSystemSolver::UsePreviousPreconditioner(bool reuse_or_not)
{
  KSPSetReusePreconditioner(ksp, reuse_or_not ? PETSC_TRUE : PETSC_FALSE);
}

//-----------------------------------------------------

bool
LinearSystemSolver::Solve(SpaceVariable3D &b, SpaceVariable3D &x,
                          LinearSolverConvergenceReason *reason, int *numIts, std::vector<double> *rnorm)
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
  x.SyncLocalToGlobal(); //update the "localVec" of x to match xx
  // ---------------------------------------------------


  KSPConvergedReason ksp_code; //positive if convergence; negative if diverged
  KSPGetConvergedReason(ksp, &ksp_code);
  bool success = ksp_code>0;

  if(reason) { //user requested convergence/divergence reason
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


  int nIts = 0;
  KSPGetIterationNumber(ksp, &nIts);
  if(numIts) //user requested output of number of iterations
    *numIts = nIts;

  // log
  if(rnorm || !log_filename.empty() || write_log_to_screen) {//need residual norm history
    int nEntries(0);
    KSPGetResidualHistory(ksp, NULL, &nEntries);

    if(rnorm) {
      rnorm->resize(nEntries, 0);
      for(int i=0; i<nEntries; i++) //copy data instead of passing rnorm_history (safer)
        (*rnorm)[i] = rnorm_history[i];
    }

    if(write_log_to_screen) {
      if(success) {
        if(equation_name.empty())
          print("  o Linear solver converged (It. %d).\n", nIts);
        else
          print("  o Linear solver for %s converged (It. %d).\n", equation_name.c_str(), nIts);
      } else {
        if(equation_name.empty())
          print_warning("  o Linear solver failed to converged (It. %d).\n", nIts);
        else
          print_warning("  o Linear solver for %s failed to converged (It. %d).\n",
                        equation_name.c_str(), nIts);
      }
      for(int i=0; i<nEntries; i++)
        print("    > It. %d: residual = %e.\n", i+1, rnorm_history[i]);
    }

    if(!log_filename.empty()) {
      FILE *file = fopen(log_filename.c_str(), "a");
      if(!file) {
        print_error("*** Error: Unable to open file %s for printing the log of linear system solver.\n");
        exit_mpi();
      }
      if(success) 
        print(file, "  o Linear solver converged (It. %d).\n", nIts);
      else
        print(file, "  o Linear solver failed to converged (It. %d).\n", nIts);

      for(int i=0; i<nEntries; i++)
        print(file, "    > It. %d: residual = %e.\n", i+1, rnorm_history[i]);

      fclose(file);
      MPI_Barrier(comm);
    }
  }


  return success;
}

//-----------------------------------------------------
//-----------------------------------------------------







