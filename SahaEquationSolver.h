/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _SAHA_EQUATION_SOLVER_H_
#define _SAHA_EQUATION_SOLVER_H_

#include<VarFcnBase.h>
#include<AtomicIonizationData.h>

//----------------------------------------------------------------
// Class SahaEquationSolver is responsible for solving the ideal
// Saha equation for one material (w/ a fixed id)
// Input: v and id 
// Output: Zav, ne, nh, alphas.
// (A dummy solver is defined for materials not undergoing ionization)
//----------------------------------------------------------------

class SahaEquationSolver {

protected: //should be accessible for derived classes

  // constants
  double h; //!< planck_constant; 
  double e; //!< electron_charge;
  double me; //!< electron_mass;
  double kb; //!< boltzmann_constant;

  int max_mean_atomic_number;

  double molar_mass; //the material is possibly a mixture

  double Tmin; //!< min temperature specified by user for ionization cal.

  // IoData
  MaterialIonizationModel* iod_ion_mat;

  std::vector<AtomicIonizationData> elem; //!< chemical elements / species

  VarFcnBase* vf;

public:

  SahaEquationSolver(IoData& iod, VarFcnBase* vf_); //!< creates a dummy solver

  SahaEquationSolver(MaterialIonizationModel& iod_ion_mat_, IoData& iod_, VarFcnBase* vf_, MPI_Comm* comm);

  virtual ~SahaEquationSolver();

  virtual void Solve(double* v, double& zav, double& nh, double& ne, std::map<int, std::vector<double> >& alpha_rj,
                     double* lambD = NULL);

  int GetNumberOfElements() {return elem.size();}

protected:

  //! nested class / functor: nonlinear equation for Zav
  class ZavEquation {
    std::vector<std::vector<double> > fprod; // f_{r,j}(T,...)*f_{r-1,j}(T,...)*...*f_{0,j}(T,...)
    std::vector<AtomicIonizationData>& elem;
  public:
    ZavEquation(double kb, double T, double nh, double me, double h, std::vector<AtomicIonizationData>& elem_);
    ~ZavEquation() {}
    double operator() (double zav) {return zav - ComputeRHS(zav);}
    double GetFProd(int r, int j) {assert(j<(int)fprod.size() && r<(int)fprod[j].size()); return fprod[j][r];}
    double GetZej(double zav, int j) {assert(j<(int)fprod.size()); return ComputeRHS_ElementJ(zav, j);}
  private:
    double ComputeRHS(double zav); //!< compute the right-hand-side of the Zav equation
    double ComputeRHS_ElementJ(double zav, int j);
  };

};



#endif
