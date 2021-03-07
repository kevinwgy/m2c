#ifndef _EXACT_RIEMANN_SOLVER_BASE_H_
#define _EXACT_RIEMANN_SOLVER_BASE_H_

#include <FluxFcnBase.h>
/*****************************************************************************************
 * Base class for solving one-dimensional, single- or two-material Riemann problems
 *****************************************************************************************/

class ExactRiemannSolverBase {

protected:
  FluxFcnBase&         ff; //!< this is the base that does not implement any numerical flux function
  vector<VarFcnBase*>& vf; //!< although fluxFcn has vf, it's convenient to keep another copy here

  int maxIts_main, maxIts_shock;
  int numSteps_rarefaction;
  double tol_main, tol_shock, tol_rarefaction;

public:
  ExactRiemannSolverBase(FluxFcnBase &ff_, std::vector<VarFcnBase*> &vf_, IoData &iod);
  virtual ~FluxFcnBase() {}

  virtual void ComputeRiemannSolution(int dir/*0~x,1~y,2~z*/, double *Vm, int idm, /*"left" state*/, 
                                      double *Vp, int idp, /*"right" state*/, 
                                      double *V /*solution at xi = 0 (i.e. x=0) */);

protected: //internal functions

  virtual void ComputeRhoUStarLeft(double rhol, double ul, double pl, double ps, int idl/*inputs*/, 
                   double rhols0, double rhols1/*initial guesses for Hugo. eq.*/,
                   double &rhols, double uls/*outputs*/,
                   bool *trans_rare = NULL, double *Vrare_x0 = NULL/*filled only if found a trans. rarefaction*/); 

  virtual void ComputeRhoUStarRight(double rhor, double ur, double pr, double ps, int idr/*inputs*/,
                   double rhors0, double rhors1/*initial guesses for Hugo. eq.*/,
                   double &rhors, double urs/*outputs*/,
                   bool *trans_rare = NULL, double *Vrare_x0 = NULL/*filled only if found a trans. rarefaction*/);

};

#endif
