#ifndef _SAHA_EQUATION_SOLVER_H_
#define _SAHA_EQUATION_SOLVER_H_

#include<VarFcnBase.h>

//----------------------------------------------------------------
// Class SahaEquationSolver is responsible for solving the ideal
// or non-ideal Saha equation for one material (w/ a fixed id)
// Input: v and id 
// Output: Zav, ne, nh, alphas.
// (A dummy solver is defined for materials not undergoing ionization)
//----------------------------------------------------------------

struct ElementData {

  double molar_fraction;
  int atomic_number;
  std::vector<double> I; //!< ionization energy
  std::vector<std::vector<double> > g; //!< angular momentum
  std::vector<std::vector<double> > E; //!< excitation energy
  int rmax; //max charge state

  ElementData() : molar_fraction(-1.0), atomic_number(-1) {}
  ~ElementData() {}

};


class SahaEquationSolver {

  // constants
  double h; //planck_constant; 
  double e; //electron_charge;
  double me; //electron_mass;
  double kb; //boltzmann_constant;

  int max_atomic_number;

  // IoData
  MaterialIonizationModel* iod_ion_mat;

  std::vector<ElementData> elem; //chemical elements / species

  VarFcnBase* vf;

public:

  SahaEquationSolver(IoData& iod, VarFcnBase* vf_); //creates a dummy solver

  SahaEquationSolver(MaterialIonizationModel& iod_ion_mat_, IoData& iod_, VarFcnBase* vf_);

  ~SahaEquationSolver();

  void Solve(double* v, double& zav, double& nh, double& ne, std::map<int, std::vector<double> >& alpha_rj);

  int GetNumberOfElements() {return elem.size();}

protected:

  void GetDataInFile(std::fstream& file, vector<double> &X, int MaxCount, bool non_negative);

  //! nested class / functor: nonlinear equation for Zav
  class ZavEquation {
    double kbT; //kb*T
    double nh; //p/(kb*T)
    double fcore; //(2*pi*me*kb*T/(h*h))^(3/2)
    std::vector<ElementData>& elem;
  public:
    ZavEquation(double kb, double T, double p, double me, double h, std::vector<ElementData>& elem_);
    ~ZavEquation() {}
    double operator() (double zav) {return zav - ComputeRHS(zav);}
  private:
    double ComputeRHS(double zav); //!< compute the right-hand-side of the Zav equation
  };

};



#endif
