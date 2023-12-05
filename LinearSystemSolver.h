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

public:

  LinearSystemSolver(MPI_Comm &comm_, DM &dm_, LinearSolverData &lin_input);
  ~LinearSystemSolver(); 
  void Destroy();

  void SetLinearOperator(std::vector<RowEntries>& row_entries);

  int Solve(SpaceVariable3D &b, SpaceVariable3D &x); //!< x: both input (initial guess) & output (solution)

  void GetTolerances(double *rtol, double *abstol, double *dtol, int *maxits); //!< set NULL to params not needed

private:

  void SetTolerances(LinearSolverData &lin_input);

};

#endif
