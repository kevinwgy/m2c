#ifndef _EXACT_RIEMANN_SOLVER_BASE_H_
#define _EXACT_RIEMANN_SOLVER_BASE_H_

#include <VarFcnBase.h>
#include <vector>
/*****************************************************************************************
 * Base class for solving one-dimensional, single- or two-material Riemann problems
 *****************************************************************************************/
#define PRINT_RIEMANN_SOLUTION 0

class ExactRiemannSolverBase {

protected:
  vector<VarFcnBase*>& vf; 

  int maxIts_main, maxIts_shock;
  int numSteps_rarefaction;
  double tol_main, tol_shock, tol_rarefaction;

public:
  ExactRiemannSolverBase(std::vector<VarFcnBase*> &vf_, ExactRiemannSolverData &iod_riemann);
  virtual ~ExactRiemannSolverBase() {}

  virtual void ComputeRiemannSolution(int dir/*0~x,1~y,2~z*/, double *Vm, int idm /*"left" state*/, 
                                      double *Vp, int idp /*"right" state*/, 
                                      double *V, int &id /*solution at xi = 0 (i.e. x=0) */,
                                      double *Vsm /*left 'star' solution*/,
                                      double *Vsp /*right 'star' solution*/);

#if PRINT_RIEMANN_SOLUTION == 1
  vector<vector<double> > sol1d;
#endif

protected: //internal functions

//! Nested class / functor: Hugoniot eq. (across a shock wave) as a function of rho_K* (K = l,r)
  struct HugoniotEquation {

    HugoniotEquation(VarFcnBase* vf_, double rho_, double p_, double ps_)
                    : vf(vf_), rho(rho_), p(p_), ps(ps_) { 
      e = vf->GetInternalEnergyPerUnitMass(rho_,p_);
      pavg = 0.5*(p_ + ps_);
      one_over_rho = 1.0/rho_;
    }

    inline double operator() (double rhos) {
      es = vf->GetInternalEnergyPerUnitMass(rhos, ps);
      return es - e  + pavg*(1.0/rhos - one_over_rho);
    }

    private:
    VarFcnBase* vf;
    double rho, p, e, ps, es, pavg, one_over_rho;
  };

  virtual void ComputeRhoUStar(int wavenumber /*1 or 3*/,
                   double rho, double u, double p, double ps, int id/*inputs*/,
                   double rhos0, double rhos1/*initial guesses for Hugo. eq.*/,
                   double &rhos, double &us/*outputs*/,
                   bool *trans_rare = NULL, double *Vrare_x0 = NULL/*filled only if found tran rf*/,
                   bool rough_estimate = false/*set to true for the initial guesses*/);

  virtual void Rarefaction_OneStepRK4(int wavenumber/*1 or 3*/, int id,
                   double rho_0, double u_0, double p_0 /*start state*/,
                   double drho /*step size*/,
                   double &rho, double &u, double &p, double &xi /*output*/);
};

#endif
