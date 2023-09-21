/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _EXACT_RIEMANN_SOLVER_BASE_H_
#define _EXACT_RIEMANN_SOLVER_BASE_H_

#include <VarFcnBase.h>
#include <vector>

/*****************************************************************************************
 * Base class for solving one-dimensional, single- or two-material Riemann problems
 *****************************************************************************************/
//#define PRINT_RIEMANN_SOLUTION 0

class ExactRiemannSolverBase {

protected:
  vector<VarFcnBase*>& vf; 
  ExactRiemannSolverData& iod_riemann;  

  int maxIts_main, maxIts_bracket, maxIts_shock;
  int numSteps_rarefaction;
  double tol_main;
  double tol_shock;
  double tol_rarefaction; // non-dimensional, for rarefaction end points
  double min_pressure, failure_threshold, pressure_at_failure;
  std::vector<std::vector<double> > integrationPath1; // first index: 1-pressure, 2-density, 3-velocity
  std::vector<std::vector<double> > integrationPath3;

  bool surface_tension; // an indicator of whether consider surface tension

public:

  ExactRiemannSolverBase(std::vector<VarFcnBase*> &vf_, ExactRiemannSolverData &iod_riemann_);

  virtual ~ExactRiemannSolverBase() {}

  virtual double GetSurfaceTensionCoefficient();

  virtual int ComputeRiemannSolution(double *dir/*unit normal*/, double *Vm, int idm /*"left" state*/, 
                                     double *Vp, int idp /*"right" state*/,
                                     double *Vs, int &id /*solution at xi = 0 (i.e. x=0) */,
                                     double *Vsm /*left 'star' solution*/,
                                     double *Vsp /*right 'star' solution*/,
                                     double curvature = 0.0);

  virtual void PrintStarRelations(double rhol, double ul, double pl, int idl,
                          double rhor, double ur, double pr, int idr,
                          double pmin, double pmax, double dp);

  virtual int ComputeOneSidedRiemannSolution(double *dir/*unit normal towards interface/wall*/, 
                                             double *Vm, int idm /*left state*/,
                                             double *Ustar, /*interface/wall velocity (3D)*/
                                             double *Vs, int &id, /*solution at xi = 0 (i.e. x=0), id = -1 if invalid*/
                                             double *Vsm /*left 'star' solution*/);

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

  virtual bool FindInitialInterval(double rhol, double ul, double pl, double el, double cl, int idl,
           double rhor, double ur, double pr, double er, double cr, int idr, /*inputs*/
           double &p0, double &rhol0, double &rhor0, double &ul0, double &ur0,
           double &p1, double &rhol1, double &rhor1, double &ul1, double &ur1/*outputs*/);

  virtual bool FindInitialFeasiblePoints(double rhol, double ul, double pl, double el, double cl, int idl,
           double rhor, double ur, double pr, double er, double cr, int idr, /*inputs*/
           double &p0, double &rhol0, double &rhor0, double &ul0, double &ur0,
           double &p1, double &rhol1, double &rhor1, double &ul1, double &ur1/*outputs*/);

  virtual int FindInitialFeasiblePointsByAcousticTheory(double rhol, double ul, double pl, double el, double cl, int idl,
           double rhor, double ur, double pr, double er, double cr, int idr, /*inputs*/
           double &p0, double &rhol0, double &rhor0, double &ul0, double &ur0,
           double &p1, double &rhol1, double &rhor1, double &ul1, double &ur1/*outputs*/);

  virtual bool ComputeRhoUStar(int wavenumber /*1 or 3*/,
		   std::vector<std::vector<double>>& integrationPath /*3 by n, first index: 1-pressure, 2-density, 3-velocity*/,
                   double rho, double u, double p, double ps, int id/*inputs*/,
                   double rhos0, double rhos1/*initial guesses for Hugo. eq.*/,
                   double &rhos, double &us/*outputs*/,
                   bool *trans_rare = NULL, double *Vrare_x0 = NULL/*filled only if found tran rf*/);

  virtual bool Rarefaction_OneStepRK4(int wavenumber/*1 or 3*/, int id,
                            double rho_0, double u_0, double p_0 /*start state*/, 
                            double dp /*step size*/,
                            double &rho, double &u, double &p, double &xi /*output*/,
                            double & uErr, double & rhoErr /*output: absolute error in us*/);

  virtual void FinalizeSolution(double *dir, double *Vm, double *Vp,
           double rhol, double ul, double pl, int idl,
           double rhor, double ur, double pr, int idr,
           double rhol2, double rhor2, double u2, double p2,
           bool trans_rare, double Vrare_x0[3], /*inputs*/
           double *Vs, int &id, double *Vsm, double *Vsp /*outputs*/);


  //! For one-sided Riemann problem
  
  //! In the case of a shock, need two initial guesses of pressure
  bool FindInitialIntervalOneSided(double rhol, double ul, double pl, double el, double cl, int idl, double ustar,
           double &p0, double &rhol0, double &ul0, double &p1, double &rhol1, double &ul1/*outputs*/);

  bool FindInitialFeasiblePointsOneSided(double rhol, double ul, double pl, double el, double cl, int idl, double ustar,
           double &p0, double &rhol0, double &ul0, double &p1, double &rhol1, double &ul1/*outputs*/);

  int FindInitialFeasiblePointsOneSidedByAcousticTheory(double rhol, double ul,
           double pl, double el, double cl, int idl, double ustar,
           double &p0, double &rhol0, double &ul0, double &p1, double &rhol1, double &ul1/*outputs*/);

  //! Integrate the isentropic relations to the wall velocity us.
  bool ComputeOneSidedRarefaction(double rho, double u, double p, double e,
                                  double c, int id, double us/*inputs*/,
                                  double &rhos, double &ps/*outputs*/,
                                  bool *trans_rare, double *Vrare_x0/*filled only if found tran rf*/);

  void FinalizeOneSidedSolution(double *dir, double *Vm, double rhol, double ul, double pl, int idl,
                                double rhol2, double u2/*ustar*/, double p2,
                                bool trans_rare, double Vrare_x0[3], /*inputs*/
                                double *Vs, int &id, double *Vsm /*outputs*/);


};


/*****************************************************************************************
 * A derived class that contains an older version of the Riemann solver that does not apply
 * adapative Runge-Kutta for integrating the isentropic relations in the case of rarefaction.
 * Slower, but *may* be more robust. Used as a fail-safe procedure.
 *****************************************************************************************/

class ExactRiemannSolverNonAdaptive: public ExactRiemannSolverBase {

public:
  ExactRiemannSolverNonAdaptive(std::vector<VarFcnBase*> &vf_, ExactRiemannSolverData &iod_riemann_) : ExactRiemannSolverBase(vf_, iod_riemann_) {};

  int ComputeRiemannSolution(double *dir/*unit normal*/, double *Vm, int idm /*"left" state*/, 
                             double *Vp, int idp /*"right" state*/, 
                             double *Vs, int &id /*solution at xi = 0 (i.e. x=0) */,
                             double *Vsm /*left 'star' solution*/,
                             double *Vsp /*right 'star' solution*/,
                             double curvature = 0.0);

protected:
  bool ComputeRhoUStar(int wavenumber /*1 or 3*/,
		   std::vector<std::vector<double>>& integrationPath /*3 by n, first index: 1-pressure, 2-density, 3-velocity*/,
                   double rho, double u, double p, double ps, int id/*inputs*/,
                   double rhos0, double rhos1/*initial guesses for Hugo. eq.*/,
                   double &rhos, double &us/*outputs*/,
                   bool *trans_rare = NULL, double *Vrare_x0 = NULL/*filled only if found tran rf*/);

  bool Rarefaction_OneStepRK4(int wavenumber/*1 or 3*/, int id,
                            double rho_0, double u_0, double p_0 /*start state*/, 
                            double dp /*step size*/,
                            double &rho, double &u, double &p, double &xi /*output*/);

 
};

#endif
