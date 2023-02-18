/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _NON_IDEAL_SAHA_EQUATION_SOLVER_H_
#define _NON_IDEAL_SAHA_EQUATION_SOLVER_H_

#include<SahaEquationSolver.h>

//----------------------------------------------------------------
// Class NonIdealSahaEquationSolver is responsible for solving the 
// non-ideal Saha equation for one material (w/ a fixed id)
// Input: v and id 
// Output: Zav, ne, nh, alphas.
// (A dummy solver is defined for materials not undergoing ionization)
//----------------------------------------------------------------

class NonIdealSahaEquationSolver : public SahaEquationSolver {

  //! constants
  double eps0; //!< vacuum permittivity;
  double pi;
  double factor_deltaI, factor_LambB;

  //! variables for *temporary* use
  std::vector<std::vector<double> > f; //!< f_{r,j}, j=0,...,elem.size-1, r=0,...,rmax(j), (f[0][j]=0, not used)
  std::vector<std::vector<double> > alpha; //!< alpha_{r,j}, same dimensions as f

public:

  NonIdealSahaEquationSolver(MaterialIonizationModel& iod_ion_mat_, IoData& iod_, VarFcnBase* vf_, MPI_Comm* comm);

  ~NonIdealSahaEquationSolver();

  void Solve(double* v, double& zav, double& nh, double& ne, std::map<int, std::vector<double> >& alpha_rj,
             double* lambD = NULL);

protected:

  // computes the depression of ionization energy, for a given Debye length lambD and temperature T
  double ComputeDeltaI(int r, int j, double T, double nh, double zav, double one_over_lambD);

  // returns Zej (for a given j), fills "f" if zav!=0. Also fills "alpha" if compute_alpha == true
  double ComputeStateForElement(int j, double T, double nh, double zav, 
                                double one_over_lambD, bool compute_alpha);


  //! nested class / functor: nonlinear equation for lambD (Debye length). This is the master equation
  //! note that the independent variable is actually 1/lambD, not lambD (which can be inf)
  class LambDEquation { 
    double factor_lambD;
    double T, nh;
    NonIdealSahaEquationSolver& saha;
    double *zav_ptr; //stores the value obtained from last call to operator()
  public:
    LambDEquation(NonIdealSahaEquationSolver& saha_, double T_, double nh_, double *zav_ptr_);
    ~LambDEquation() {}
    double operator() (double one_over_lambD);
    double ComputeRHS(double one_over_lambD, double zav);

  private:
    //! nested class / functor: nonlinear equation for Zav, given lambD
    class ZavEquation {
      double T, nh, one_over_lambD;
      NonIdealSahaEquationSolver& saha;
    public:
      ZavEquation(NonIdealSahaEquationSolver &saha_, double T_, double nh_, double one_over_lambD_);
      ~ZavEquation() {}
      double operator() (double zav) {return zav - ComputeRHS(zav);}
    private:
      double ComputeRHS(double zav); //!< compute the right-hand-side of the Zav equation
    };
  };


};



#endif
