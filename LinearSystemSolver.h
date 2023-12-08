/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _LINEAR_SYSTEM_SOLVER_H_
#define _LINEAR_SYSTEM_SOLVER_H_
#include<petscksp.h>
#include<LinearOperator.h>

/*********************************************************************
 * class LinearSystemSolver is responsible for solving large-scale linear
 * systems Ax = b, where x and b are SpaceVariable3D. The actual work is
 * done by PETSc. In some sense, LinearSystemSolver is just
 * a wrapper over PETSc/KSP. 
 *********************************************************************
*/

class LinearSystemSolver : public LinearOperator {

  KSP ksp;
  std::vector<double> rnorm_history; //!< stores the history of residual norm for the "Solve"

public:

  enum ConvergenceReason {NONE = 0, CONVERGED_REL_TOL = 1, CONVERGED_ABS_TOL = 2, CONVERGED_OTHER = 3,
                          DIVERGED_ITS = 4, DIVERGED_DTOL = 5, DIVERGED_OTHER = 6};


  LinearSystemSolver(MPI_Comm &comm_, DM &dm_, LinearSolverData &lin_input);
  ~LinearSystemSolver(); 
  void Destroy();

  void SetLinearOperator(std::vector<RowEntries>& row_entries);

  bool Solve(SpaceVariable3D &b, SpaceVariable3D &x, //!< x: both input (initial guess) & output (solution)
             ConvergenceReason *reason = NULL, int *numIts = NULL, std::vector<double> *rnorm = NULL) ;

  void GetTolerances(double *rtol, double *abstol, double *dtol, int *maxits); //!< set NULL to params not needed

private:

  void SetTolerances(LinearSolverData &lin_input);

};

#endif
