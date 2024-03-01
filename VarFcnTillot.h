/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _VAR_FCN_TILLOT_H_
#define _VAR_FCN_TILLOT_H_

#include <VarFcnBase.h>
#include <polynomial_equations.h>
#include <ordinary_differential_equations.h>
#include <vector>
#include <algorithm> //std::lower_bound
#include <boost/math/tools/roots.hpp>

/********************************************************************************
 * This class is the VarFcn class for the Tillotson equation of state (EOS)
 * Only elementary functions are declared and/or defined here.
 * All arguments must be pertinent to only a single grid node or a single
 * state.
 *
 * The Tillotson EOS presented in Brundage, 2013 is implemented here. For details,
 * see KW's notes.
 *
 *   p      : pressure
 *   e      : internal energy per unit mass.
 *   rho    : density
 *   eta    : rho/rho0
 *   mu     : eta - 1 
 *   omega  : rho0/rho - 1
 *   chi    : 1/(e/(e0*eta*eta)+1)
 *
 *   rho0   : density in the ambient state
 *   e0     : internal energy (per unit mass) in the ambient state
 *   a, b   : non-D model parameters. (a+b/4: Gruneisen param in ambient state; a+b: ... at low T (e=0))
 *   A, B   : model parameters. dim: [force]/[length]^2. (A: bulk modulus at p = 0 and e = 0)
 *   alpha  : non-D model parameters (in the formula for "hot expanded states")
 *   beta   : non-D model parameters (in the formula for "hot expanded states")
 *   rhoIV  : "incipient vaporization" density (constant model param.)
 *   eIV    : "incipient vaporization" internal energy (constant model param.)
 *   eCV    : "complete vaporization" internal energy (constant model param.) (eCV-eIV: latent heat per unit mass)
 *
 *   cv     : specific heat at constant volume
 *   cp     : specific heat at constant pressure 
 *   h0     : specific enthalpy corresponding to T0 (only used if T is calculated using cp)
 *   T0     : reference temperature at rho0 and e0 (w/ linear cv law), rho0 and h0 (w/ linear cp law), or
 *            rho0 and e=ecold (w/ T(rho,e) law) 
 *   temperature_depends_on_density: whether T depends on both rho and e, or just e.
 *
 *   Note: We apply the convention that one should first check Cases 1,2,3 before calling Case 4 ("p1|2") functions.
 *         Some Case 4 functions may crash (assertion failure) if it is not applied in the correct domain.
 *
 *   References: KW's notes.
 *   
 ********************************************************************************/

class VarFcnTillot : public VarFcnBase {

private:

  double rho0, e0, a, b, A, B, alpha, beta;
  double rhoIV, eIV, eCV;
  double cv, T0;
  bool temperature_depends_on_density;

  bool use_cp;
  double cp, h0;

  double elat; //!< eCV-eIV;
  double invcv; //!< 1/cv;
  double invcp; //!< 1/cp;

  double tol; //!< convergence tolerance, non-D.

  //! Stores the rho-e trajectory of the interatomic/intermolecular potential ("cold internal energy") for 
  //! temperature-related calculations. Only constructed & used if temperature_depends_on_density == true.
  std::vector<double> ecold_plus, rho_plus; //!< for rho>=rho0
  std::vector<double> ecold_minus, rho_minus; //!< for rho<=rho0
  int Nmax = (int)1e7; //!< max size of ecold_plus & ecold_minus, hard-coded for the moment

public:

  VarFcnTillot(MaterialModelData &data);
  ~VarFcnTillot() {}

  inline double GetPressure(double rho, double e) {return (this->*GetPressureCase[GetCaseWithRhoE(rho,e)])(rho,e);}

  double GetInternalEnergyPerUnitMass(double rho, double p);

  double GetDensity(double p, double e);

  inline double GetDpdrho(double rho, double e) {return (this->*GetDpdrhoCase[GetCaseWithRhoE(rho,e)])(rho,e);}

  inline double GetBigGamma(double rho, double e) {return (this->*GetGammaCase[GetCaseWithRhoE(rho,e)])(rho,e);}

  double GetTemperature(double rho, double e);

  inline double GetReferenceTemperature() {return T0;} //!< reference temperature (ambient state)

  inline double GetReferenceInternalEnergyPerUnitMass() {return e0;} //!< ambient state

  double GetInternalEnergyPerUnitMassFromTemperature(double rho, double T);

  double GetInternalEnergyPerUnitMassFromEnthalpy(double rho, double h);

private:

  void SetupColdEnergyTrajectory();
  void ExtendColdEnergyTrajectoryUpwards(double rhomax);
  void ExtendColdEnergyTrajectoryDownwards(double rhomin);

  //! Get ecold from trajectory. Extend the trajectory if needed.
  double GetColdEnergy(double rho);

  inline double GetChiWithEta(double eta, double e) {return 1.0/(e/(e0*eta*eta)+1.0);}
  inline double GetChiWithOmega(double omega, double e) {return 1.0/(e/e0*(omega+1.0)*(omega+1)+1.0);}

  //! Determine the case id, given rho and e. Returns 0, 1, 2, or 3 for Case 1, 2, 3, and 1|2
  int GetCaseWithRhoE(double rho, double e) {
    if(rho<=0.0) {
      fprintf(stdout,"\033[0;31m*** Error: VarFcnTillot::GetCaseWithRhoE detected non-positive rho (%e).\033[0m\n",
              rho);
      exit(-1);
    }
    if(e<0.0 && verbose>=1) {
      fprintf(stdout,"\033[0;35mWarning: In VarFcnTillot::GetCaseWithRhoE detected negative e (%e).\033[0m\n", e);
    }
    if(rho>=rho0) return 0; // Case 1 
    if(e>=eCV)    return 1; // Case 2
    if(rho<rhoIV) return 2; // Case 3
    if(e<=eIV)    return 0; // Case 1
    return 3; // Case 1|2
  }
   
  
  /********************************
   *            Case 1
   *******************************/
  double GetPressure1(double rho, double e) {
    double eta = rho/rho0;
    double mu  = eta - 1.0;
    return (a + b*GetChiWithEta(eta,e))*rho*e + (A + B*mu)*mu;
  }

  double GetGamma1(double rho, double e) {
    double chi = GetChiWithEta(rho/rho0, e);
    return a + b*chi*chi;
  }

  double GetDpdrho1(double rho, double e) {
    double eta = rho/rho0;
    double mu  = eta - 1.0;
    double chi = GetChiWithEta(eta, e);
    return a*e + b*e*chi*chi*(1.0 + 3.0*e/(e0*eta*eta)) + (A + 2.0*B*mu)/rho0;
  }
         
  double GetInternalEnergyPerUnitMass1(double rho, double p);

  /********************************
   *            Case 2
   *******************************/
  double GetPressure2(double rho, double e) {
    double mu    = rho/rho0 - 1.0;
    double omega = rho0/rho - 1.0;
    double rho_e = rho*e;
    return a*rho_e + (b*rho_e*GetChiWithOmega(omega,e) + A*mu*exp(-beta*omega))*exp(-alpha*omega*omega);
  }
 
  double GetGamma2(double rho, double e) {
    double omega = rho0/rho - 1.0;
    double chi   = GetChiWithOmega(omega, e);
    return a + b*chi*chi*exp(-alpha*omega*omega);  
  }

  double GetDpdrho2(double rho, double e) {
    double omega = rho0/rho - 1.0;
    double chi   = GetChiWithOmega(omega, e);
    double omom  = omega*(1.0+omega);
    return a*e + A/rho0*(1.0 - omom*(beta+2.0*alpha*omega))*exp(-omega*(beta+alpha*omega))
               + b*e*chi*chi*(1.0 + 2.0*alpha*omom + e/e0*(1.0+omega)*(1.0+omega)*(3.0+2.0*alpha*omom))*exp(-alpha*omega*omega);
  }

  double GetInternalEnergyPerUnitMass2(double rho, double p);



  /********************************
   *            Case 3 
   *******************************/
  double GetPressure3(double rho, double e) {
    double eta = rho/rho0;
    return (a + b*GetChiWithEta(eta,e))*rho*e + A*(eta - 1.0);
  }

  double GetGamma3(double rho, double e) { //same as GetGamma1
    double chi = GetChiWithEta(rho/rho0, e);
    return a + b*chi*chi;
  }

  double GetDpdrho3(double rho, double e) {
    double eta = rho/rho0;
    double chi = GetChiWithEta(eta, e);
    return a*e + b*e*chi*chi*(1.0 + 3.0*e/(e0*eta*eta)) + A/rho0;
  }

  double GetInternalEnergyPerUnitMass3(double rho, double p);





  /********************************
   *           Case 1|2 
   *******************************/
  inline double GetPressure12(double rho, double e) {
    return ((eCV-e)*GetPressure1(rho,e) + (e-eIV)*GetPressure2(rho,e))/elat;
  }

  inline double GetGamma12(double rho, double e) {
    return (  (GetPressure2(rho,e)-GetPressure1(rho,e))/rho 
            + (eCV-e)*GetGamma1(rho,e) + (e-eIV)*GetGamma2(rho,e) )/elat;
  }

  inline double GetDpdrho12(double rho, double e) {
    return ((eCV-e)*GetDpdrho1(rho,e) + (e-eIV)*GetDpdrho2(rho,e))/elat;
  }

  double GetInternalEnergyPerUnitMass12(double rho, double p);



  /********************************
   * Function Arrays (All Cases)      
   *******************************/
  typedef double (VarFcnTillot::*DoubleFunction) (double rho, double e);
  DoubleFunction GetPressureCase[4] = {&VarFcnTillot::GetPressure1, &VarFcnTillot::GetPressure2,
                                       &VarFcnTillot::GetPressure3, &VarFcnTillot::GetPressure12};
  DoubleFunction GetGammaCase[4]    = {&VarFcnTillot::GetGamma1, &VarFcnTillot::GetGamma2,
                                       &VarFcnTillot::GetGamma3, &VarFcnTillot::GetGamma12};
  DoubleFunction GetDpdrhoCase[4]   = {&VarFcnTillot::GetDpdrho1, &VarFcnTillot::GetDpdrho2,
                                       &VarFcnTillot::GetDpdrho3, &VarFcnTillot::GetDpdrho12};





  /*************************************
   * Structs for Numerical Root-Finding
   ************************************/
  
  //! Solve f(e) = p - p1|2(rho,e) = 0 for e, given rho and p. (Divided by p to non-dimensionalize.)
  struct SpecificEnergyEquationP12 {
    SpecificEnergyEquationP12(double rho_, double p_, VarFcnTillot *vf_) 
        : vf(vf_), rho(rho_), p(p_), pfabs(fabs(p_)) {}
    inline double operator() (double e) {
      return pfabs==0.0 ?  -(vf->GetPressure12(rho, e)) : (p - vf->GetPressure12(rho, e))/pfabs;}

  private:
    double rho, p, pfabs;
    VarFcnTillot *vf;
  };   


  //! Solve f(e) = e + p(rho,e)/rho - (h0 + cp*(T-T0)) = 0 for e, given rho and T.
  //! (Divided by h0+cp(T-T0) to nondimensionalize.)
  struct EnthalpyTemperatureEquation {
    EnthalpyTemperatureEquation(double rho_, double T_, VarFcnTillot *vf_)
        : vf(vf_), rho(rho_), T(T_), enthalpy(vf_->h0 + vf_->cp*(T_ - vf_->T0)) {}
    double enthalpy; //!< public member
    inline double operator() (double e) {
      return enthalpy==0.0 ? e + vf->GetPressure(rho,e)/rho
                           : (e + vf->GetPressure(rho,e)/rho - enthalpy)/fabs(enthalpy);}
    
  private:
    double rho, T;
    VarFcnTillot *vf;
  };


  //! Solve f(e) = e + p(rho,e)/rho - h = 0 for e, given rho and h
  //! (Divide by h to nondimensionalize.)
  struct EnthalpyEquation {
    EnthalpyEquation(double rho_, double h_, VarFcnTillot *vf_)
        : vf(vf_), rho(rho_), h(h_), hfabs(fabs(h_)) {}
    inline double operator() (double e) {
      return hfabs==0.0 ? e + vf->GetPressure(rho,e)/rho : (e + vf->GetPressure(rho,e)/rho - h)/hfabs;}

  private:
    double rho, h, hfabs;
    VarFcnTillot *vf;
  };


  //! Solve f(rho) = p(rho,e) - p = 0 for rho, given p and e
  //! (Divided by p to nondimensionalize)
  struct DensityEquation {
    DensityEquation(double p_, double e_, VarFcnTillot *vf_) : p(p_), pfabs(fabs(p_)), e(e_), vf(vf_) {}
    inline double operator() (double rho) {
      return pfabs==0.0 ? vf->GetPressure(rho,e) : (vf->GetPressure(rho,e)-p)/pfabs;}

    private:
      double p, pfabs, e;
      VarFcnTillot *vf;
  };

};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

VarFcnTillot::VarFcnTillot(MaterialModelData &data) : VarFcnBase(data)
{

  if(data.eos != MaterialModelData::TILLOTSON){
    fprintf(stdout, "\033[0;31m*** Error: MaterialModelData is not of type Tillotson.\033[0m\n");
    exit(-1);
  }

  type = TILLOTSON;

  
  rho0   = data.tillotModel.rho0;
  e0     = data.tillotModel.e0;
  a      = data.tillotModel.a;
  b      = data.tillotModel.b;
  A      = data.tillotModel.A;
  B      = data.tillotModel.B;
  alpha  = data.tillotModel.alpha;
  beta   = data.tillotModel.beta;

  rhoIV  = data.tillotModel.rhoIV;
  eIV    = data.tillotModel.eIV;
  eCV    = data.tillotModel.eCV;

  cv     = data.tillotModel.cv;
  T0     = data.tillotModel.T0;
  temperature_depends_on_density = data.tillotModel.temperature_depends_on_density 
                                     == TillotsonModelData::YES;

  cp     = data.tillotModel.cp;
  h0     = data.tillotModel.h0;

  tol = 1.0e-8; //!< TODO: hard-coded for the moment, non-dimensional.

  if(temperature_depends_on_density) {
    use_cp = false; //we use e-T relation to determine the thermal vibration energy
    SetupColdEnergyTrajectory(); //!< fill ecold_plus, rho_plus, ecold_minus, rho_minus
  } 
  else {
    use_cp = (cp>0 && cv<=0.0) ? true : false;
  }

  invcv = cv==0.0 ? 0.0 : 1.0/cv;
  invcp = cp==0.0 ? 0.0 : 1.0/cp;

  elat = eCV-eIV;

  if(elat<=0.0) {
    fprintf(stdout, "\033[0;31m*** Error: VarFcnTillot detected eCV (%e) <= eIV (%e).\033[0m\n", eCV, eIV);
    exit(-1);
  }
  if(rhoIV>=rho0) {
    fprintf(stdout, "\033[0;31m*** Error: VarFcnTillot detected rhoIV (%e) >= rho0 (%e).\033[0m\n", rhoIV, rho0);
    exit(-1);
  }
  if(e0<=0.0) {
    fprintf(stdout, "\033[0;31m*** Error: VarFcnTillot detected non-positive e0 (%e).\033[0m\n", e0);
    exit(-1);
  }
  if(rhoIV<=0.0 || rho0<=0.0) {
    fprintf(stdout, "\033[0;31m*** Error: VarFcnTillot detected non-positive rhoIV (%e) or rho0 (%e).\033[0m\n", rhoIV, rho0);
    exit(-1);
  }

}

//------------------------------------------------------------------------------

void
VarFcnTillot::SetupColdEnergyTrajectory()
{
  assert(temperature_depends_on_density); //otherwise, no need to do this

  // Build the trajectory above rho0. 
  rho_plus.push_back(rho0);
  ecold_plus.push_back(0.0); //ecold(rho0) = 0.0; 
  ExtendColdEnergyTrajectoryUpwards(2.0*rho0); 

  // Build the trajectory below rho0.
  rho_minus.push_back(rho0); 
  ecold_minus.push_back(0.0); //ecold(rho0) = 0.0
  ExtendColdEnergyTrajectoryDownwards(0.5*rho0);
}

//------------------------------------------------------------------------------

void
VarFcnTillot::ExtendColdEnergyTrajectoryUpwards(double rhomax)
{
  if(rhomax<=rho_plus.back())
    return;

  // sets up the ODE
  auto fun = [&](double *ec, double rho, double *dec_drho) {
    *dec_drho = GetPressure(rho, *ec)/(rho*rho); return;};

  // determine the initial step size
  double drho0;
  int mysize = rho_plus.size();
  if(mysize>2) {
    drho0 = rho_plus[mysize-2] - rho_plus[mysize-3]; //rho_plus[mysize-1]-rho_plus[mysize-2] may not
                                                     //be the desired step size, as the last element
                                                     //of rho_plus is pre-specified/fixed.
  } else
    drho0 = 1e-3*rho0;

  // integration
  double ecold_max;
  std::vector<double> rho_new, ecold_new;
  int err = MathTools::runge_kutta_45(fun, 1, rho_plus[mysize-1], &ecold_plus[mysize-1], drho0, rhomax,
                                      &ecold_max, [](double*){return false;}, 
                                      tol/100.0, Nmax-mysize, &rho_new, &ecold_new); //smaller tol, because we will interpolate
  assert(!err);

  int new_size = rho_new.size();
  assert(new_size == (int)ecold_new.size());

  rho_plus.reserve(mysize + new_size - 1);
  ecold_plus.reserve(mysize + new_size - 1);
  for(int i=1; i<new_size; i++) {
    rho_plus.push_back(rho_new[i]);
    ecold_plus.push_back(ecold_new[i]);
  }
}

//------------------------------------------------------------------------------

void
VarFcnTillot::ExtendColdEnergyTrajectoryDownwards(double rhomin)
{
  if(rhomin>=rho_minus.back())
    return;

  assert(rhomin>0.0);

  // sets up the ODE
  auto fun = [&](double *ec, double rho, double *dec_drho) {
    *dec_drho = GetPressure(rho, *ec)/(rho*rho); return;};

  // determine the initial step size (<0 in this case)
  double drho0; 
  int mysize = rho_minus.size();
  if(mysize>2) {
    drho0 = rho_minus[mysize-2] - rho_minus[mysize-3]; //rho_minus[mysize-1]-rho_minus[mysize-2] may not
                                                       //be the desired step size, as the last element
                                                       //of rho_minus is pre-specified/fixed.
  } else
    drho0 = -1e-3*rho0;

  // integration
  double ecold_min;
  std::vector<double> rho_new, ecold_new;
  int err = MathTools::runge_kutta_45(fun, 1, rho_minus[mysize-1], &ecold_minus[mysize-1], drho0, rhomin,
                                      &ecold_min, [](double*){return false;}, 
                                      tol/100.0, Nmax-mysize, &rho_new, &ecold_new); //smaller tol, because we will interpolate
  assert(!err);

  int new_size = rho_new.size();
  assert(new_size == (int)ecold_new.size());

  rho_minus.reserve(mysize + new_size - 1);
  ecold_minus.reserve(mysize + new_size - 1);
  for(int i=1; i<new_size; i++) {
    rho_minus.push_back(rho_new[i]);
    ecold_minus.push_back(ecold_new[i]);
  }
}

//------------------------------------------------------------------------------

double
VarFcnTillot::GetColdEnergy(double rho)
{
  int ind, index;
  if(rho>=rho0) { // get from ecold_plus
    ind = std::lower_bound(rho_plus.begin(), rho_plus.end(), rho) - rho_plus.begin();
    if(ind == (int)rho_plus.size()) {
      assert(rho_plus.back()<rho);
      ExtendColdEnergyTrajectoryUpwards(rho);
      ind = rho_plus.size()-1;
    }
    if(rho_plus[ind]==rho) {
      index = ind;
      return ecold_plus[index];
    } else { //linear interpolation
      index = ind-1; 
      return ((rho_plus[ind]-rho)*ecold_plus[index] + (rho-rho_plus[index])*ecold_plus[ind])/
             (rho_plus[ind]-rho_plus[index]);
    }  
  }
  else { // get from ecold_minus
    ind = std::lower_bound(rho_minus.begin(), rho_minus.end(), rho,
                               [](double a, double b){return a>b;}) - rho_minus.begin();
    if(ind == (int)rho_minus.size()) {
      assert(rho_minus.back()>rho);
      ExtendColdEnergyTrajectoryDownwards(rho);
      ind = rho_minus.size()-1;
    }
    if(rho_minus[ind]==rho) {
      index = ind;
      return ecold_minus[index];
    } else { //linear interpolation
      index = ind-1; 
      return ((rho_minus[ind]-rho)*ecold_minus[index] + (rho-rho_minus[index])*ecold_minus[ind])/
             (rho_minus[ind]-rho_minus[index]);
    }  
  }

  return 0.0; //shouldn't get here.
}

//------------------------------------------------------------------------------

double
VarFcnTillot::GetInternalEnergyPerUnitMass1(double rho, double p)
{
  assert(rho>=rhoIV);

  double eta = rho/rho0;
  double mu  = eta - 1.0;
  double inv_e0_eta_eta = 1.0/(e0*eta*eta);

  double c0 = -p + A*mu + B*mu*mu;
  double c1 = c0*inv_e0_eta_eta + (a+b)*rho;
  double c2 = a*rho*inv_e0_eta_eta;

  double e_root_1, e_root_2;
  int num_real_roots = MathTools::quadratic_equation_solver(c2, c1, c0, e_root_1, e_root_2);
   
  if(num_real_roots==2) { // e_root_1 > e_root_2 (see quadratic_equation_solver)
    if(e_root_1<0.0 && verbose==1) {
      fprintf(stdout, "\033[0;35mWarning: VarFcnTillot(P1) got two negative e's (%e,%e) for rho = %e, p = %e. "
              "Taking the first one.\033[0m\n", e_root_1, e_root_2, rho, p);
      return e_root_1;
    }
    else if(e_root_2>0.0 && verbose==1) {
      fprintf(stdout, "\033[0;35mWarning: VarFcnTillot(P1) got two positive e's (%e,%e) for rho = %e, p = %e. "
              "Taking the first one.\033[0m\n", e_root_1, e_root_2, rho, p);
      return e_root_1;
    }
    return e_root_1;
  }
  else if(num_real_roots==1) {
    if(e_root_1<0.0 && verbose==1) {
      fprintf(stdout, "\033[0;35mWarning: VarFcnTillot(P1) got a negative e (%e) for rho = %e, p = %e. \033[0m\n",
              e_root_1, rho, p);
      return e_root_1;
    }
    return e_root_1;
  }
  else { // num_real_roots = 0
    fprintf(stdout, "\033[0;35mWarning: VarFcnTillot(P1) failed to find an e for rho = %e, p = %e. "
            "Solution (real number) does not exist. Returning e = %e (1.0e-6*e0).\033[0m\n", 
            rho, p, 1.0e-6*e0);
    return 1.0e-6*e0;
  }

  return 0; //will not reach here
}

//------------------------------------------------------------------------------

double
VarFcnTillot::GetInternalEnergyPerUnitMass2(double rho, double p)
{

  assert(rho<rho0);

  double eta   = rho/rho0;
  double omega = rho0/rho - 1.0;
  double mu    = eta - 1.0;
  double inv_e0_eta_eta = 1.0/(e0*eta*eta);

  double c0 = -p + A*mu*exp(-beta*omega - alpha*omega*omega);
  double c1 = c0*inv_e0_eta_eta + (a + b*exp(-alpha*omega*omega))*rho;
  double c2 = a*rho*inv_e0_eta_eta;

  double e_root_1, e_root_2;
  int num_real_roots = MathTools::quadratic_equation_solver(c2, c1, c0, e_root_1, e_root_2);
   
  if(num_real_roots==2) { // e_root_1 > e_root_2 (see quadratic_equation_solver)
    if(e_root_1<0.0 && verbose>=1) {
      fprintf(stdout, "\033[0;35mWarning: VarFcnTillot(P2) got two negative e's (%e,%e) for rho = %e, p = %e. "
              "Taking the first one.\033[0m\n", e_root_1, e_root_2, rho, p);
      return e_root_1;
    }
    else if(e_root_2>0.0 && verbose>=1) {
      fprintf(stdout, "\033[0;35mWarning: VarFcnTillot(P2) got two positive e's (%e,%e) for rho = %e, p = %e. "
              "Taking the first one.\033[0m\n", e_root_1, e_root_2, rho, p);
      return e_root_1;
    }
    return e_root_1;
  }
  else if(num_real_roots==1) {
    if(e_root_1<0.0 && verbose>=1) {
      fprintf(stdout, "\033[0;35mWarning: VarFcnTillot(P2) got a negative e (%e) for rho = %e, p = %e. \033[0m\n",
              e_root_1, rho, p);
      return e_root_1;
    }
    return e_root_1;
  }
  else { // num_real_roots = 0
    fprintf(stdout, "\033[0;35mWarning: VarFcnTillot(P2) failed to find an e for rho = %e, p = %e. "
            "Solution (real number) does not exist. Returning e = %e (1.0e-6*e0).\033[0m\n",
            rho, p, 1.0e-6*e0);
    return 1.0e-6*e0;
  }

  return 0; //will not reach here
}

//------------------------------------------------------------------------------

double
VarFcnTillot::GetInternalEnergyPerUnitMass3(double rho, double p)
{

  assert(rho<rhoIV);

  //Almost the same as Case 1, except that B is set to 0.
  
  double eta = rho/rho0;
  double mu  = eta - 1.0;
  double inv_e0_eta_eta = 1.0/(e0*eta*eta);

  double c0 = -p + A*mu; //This is the only difference from Case 1
  double c1 = c0*inv_e0_eta_eta + (a+b)*rho;
  double c2 = a*rho*inv_e0_eta_eta;

  double e_root_1, e_root_2;
  int num_real_roots = MathTools::quadratic_equation_solver(c2, c1, c0, e_root_1, e_root_2);
   
  if(num_real_roots==2) { // e_root_1 > e_root_2 (see quadratic_equation_solver)
    if(e_root_1<0.0 && verbose>=1) {
      fprintf(stdout, "\033[0;35mWarning: VarFcnTillot(P3) got two negative e's (%e,%e) for rho = %e, p = %e. "
              "Taking the first one.\033[0m\n", e_root_1, e_root_2, rho, p);
      return e_root_1;
    }
    else if(e_root_2>0.0 && verbose>=1) {
      fprintf(stdout, "\033[0;35mWarning: VarFcnTillot(P3) got two positive e's (%e,%e) for rho = %e, p = %e. "
              "Taking the first one.\033[0m\n", e_root_1, e_root_2, rho, p);
      return e_root_1;
    }
    return e_root_1;
  }
  else if(num_real_roots==1) {
    if(e_root_1<0.0 && verbose>=1) {
      fprintf(stdout, "\033[0;35mWarning: VarFcnTillot(P3) got a negative e (%e) for rho = %e, p = %e. \033[0m\n",
              e_root_1, rho, p);
      return e_root_1;
    }
    return e_root_1;
  }
  else { // num_real_roots = 0
    fprintf(stdout, "\033[0;35mWarning: VarFcnTillot(P3) failed to find an e for rho = %e, p = %e. "
            "Solution (real number) does not exist. Returning e = %e (1.0e-6*e0).\033[0m\n",
            rho, p, 1.0e-6*e0);
    return 1.0e-6*e0;
  }

  return 0; //will not reach here
}

//------------------------------------------------------------------------------

double
VarFcnTillot::GetInternalEnergyPerUnitMass12(double rho, double p)
{
  // Note: Assumes that the solution e is between eIV and eCV!

  assert(rho<rho0 && rho>=rhoIV);

  SpecificEnergyEquationP12 equation(rho, p, this);

  double e_low  = eIV;
  double e_high = eCV; // naturally, this should be a bracketing interval.
  double f_low  = equation(e_low);
  double f_high = equation(e_high);
  if(f_low*f_high>0) {
    fprintf(stdout,"\033[0;31m*** Error: VarFcnTillot::GetInternalEnergyPerUnitMass12 called w. incorrect inputs."
                   " rho = %e, p = %e. f(%e) = %e, f(%e) = %e.\033[0m\n",
            rho, p, e_low, f_low, e_high, f_high);
    exit(-1);
  }

  boost::uintmax_t maxit = 500; //!< "maxit" is both an input and an output!
  double tolerance = tol*std::min(e_low, elat); //elat = e_high - e_low
  std::pair<double,double> sol = boost::math::tools::toms748_solve(equation, e_low, e_high, f_low, f_high,
                                 [=](double rr0, double rr1){return fabs(rr1-rr0)<tolerance;},
                                 maxit);

  return 0.5*(sol.first + sol.second);

}

//------------------------------------------------------------------------------

double
VarFcnTillot::GetInternalEnergyPerUnitMass(double rho, double p)
{
  if(rho<=0.0) {
    fprintf(stdout,"\033[0;31m*** Error: VarFcnTillot::GetInternalEnergyPerUnitMass detected negative rho (%e).\033[0m\n",
            rho);
    exit(-1);
  }

  double e;

  if(rho>=rho0) {// Case 1
    e = GetInternalEnergyPerUnitMass1(rho,p);
  } 
  else if(rho<rhoIV) { // Case 2 or 3
    e = GetInternalEnergyPerUnitMass3(rho,p);
    if(e>=eCV) {// Case 2
      e = GetInternalEnergyPerUnitMass2(rho,p);
      if(e<eCV) {
        fprintf(stdout,"\033[0;31m*** Error: VarFcnTillot::GetInternalEnergyPerUnitMass failed for "
                       "rho = %e, p = %e.\033[0m\n", rho, p);
        exit(-1);
      }
    }
  }
  else { // rho is between rhoIV and rho0 ==> Case 1, 2, or 1|2
    e = GetInternalEnergyPerUnitMass1(rho,p);
    if(e>eIV) { 
      e = GetInternalEnergyPerUnitMass2(rho,p);
      if(e<eCV) 
        e = GetInternalEnergyPerUnitMass12(rho,p);
    }
  }

  if(verbose>=1 && e<0.0) {
    fprintf(stdout,"\033[0;35mWarning: VarFcnTillot::GetInternalEnergyPerUnitMass got negative "
                   "e (%e) for rho = %e, p = %e.\033[0m\n", e, rho, p);
  }

  return e;

}

//------------------------------------------------------------------------------

double
VarFcnTillot::GetDensity(double p, double e)
{

  DensityEquation equation(p, e, this);

  // find bracketing interval
  double rho_low = 1.0e-2*rho0;  //density must be positive!
  double rho_high = 1.0e2*rho0;
  double f_low(1.0), f_high(1.0);
  bool found_low(false), found_high(false);
  for(int iter=0; iter<20; iter++) {
    f_low = equation(rho_low);
    if(f_low<=0.0) {
      found_low = true;
      break;
    }
    rho_high = rho_low; //rho_low is actually an upperbound
    f_high = f_low;
    found_high = true; 

    rho_low /= 10.0;
  }
  if(!found_low) {
    rho_low *= 10.0;
    if(verbose>=1)
      fprintf(stdout,"\033[0;35mWarning: In VarFcnTillot::GetDensity, setting rho to %e,"
                     " for p = %e, e = %e (fun = %e).\033[0m\n", rho_low, p, e, f_low);
    return rho_low;
  }

  if(!found_high) {
    for(int iter=0; iter<20; iter++) {
      f_high = equation(rho_high);
      if(f_high>=0.0) {
        found_high = true;
        break;
      }
      rho_low = rho_high;
      f_low = f_high;
      found_low = true; //not really needed.
      rho_high *= 10.0;
    }
    if(!found_high) {
      rho_high /= 10.0;
      if(verbose>=1)
        fprintf(stdout,"\033[0;35mWarning: In VarFcnTillot::GetDensity, setting rho to %e,"
                       " for p = %e, e = %e (fun = %e).\033[0m\n", rho_high, p, e, f_high);
      return rho_high;
    }
  }

  //fprintf(stdout,"rho_low = %e, rho_high = %e.\n", rho_low, rho_high);

  boost::uintmax_t maxit = 500; //!< "maxit" is both an input and an output!
  double tolerance = tol*std::min(fabs(rho_low), fabs(rho_high-rho_low));
  std::pair<double,double> sol = boost::math::tools::toms748_solve(equation, rho_low, rho_high, f_low, f_high,
                                 [=](double rr0, double rr1){return fabs(rr1-rr0)<tolerance;},
                                 maxit);

  return 0.5*(sol.first + sol.second);

}


//------------------------------------------------------------------------------

double
VarFcnTillot::GetTemperature(double rho, double e)
{

  if(!temperature_depends_on_density) {
    if(use_cp) {
      double p = GetPressure(rho, e);
      return T0 + invcp*(e + p/rho - h0);
    } else
      return T0 + invcv*(e-e0);
  }

  // Implementing the method described in Brundage 2013
  return T0 + invcv*(e - GetColdEnergy(rho));

}

//------------------------------------------------------------------------------

double
VarFcnTillot::GetInternalEnergyPerUnitMassFromTemperature(double rho, double T)
{

  if(temperature_depends_on_density)
    return GetColdEnergy(rho) + cv*(T-T0); //Done

  if(!use_cp)  
    return e0 + cv*(T-T0); //Done

  // use_cp == true

  EnthalpyTemperatureEquation equation(rho, T, this);

  // find bracketing interval
  double e_low = equation.enthalpy>0 ? 1.0e-2*equation.enthalpy : 1.0e-2*e0; //try positive values first
  double e_high = equation.enthalpy>0 ? equation.enthalpy : 10.0*e0;
  double f_low(1.0), f_high(1.0);
  bool found_low(false), found_high(false);
  for(int iter=0; iter<15; iter++) {
    f_low = equation(e_low);
    if(f_low<=0.0) {
      found_low = true;
      break;
    }
    // if f_low>0, it can be the upperbound
    e_high = e_low;
    f_high = f_low;
    found_high = true; 

    e_low /= 10.0;
  }
  if(!found_low) { //try negative values...
    e_low = -1e-5*e0;
    for(int iter=0; iter<15; iter++) {
      f_low = equation(e_low);
      if(f_low<=0.0) {
        found_low = true;
        break;
      }
      e_low *= 10.0;
    }
    if(verbose>=1 && !found_low) {
      e_low /= 10.0;
      fprintf(stdout,"\033[0;35mWarning: In VarFcnTillot::GetInternalEnergyPerUnitMassFromTemperature, setting e to %e,"
                     " for rho = %e, T = %e (fun = %e).\033[0m\n", e_low, rho, T, f_low);
      return e_low;
    }
  }

  if(!found_high) { //e_high is always set to a positive number...
    for(int iter=0; iter<20; iter++) {
      f_high = equation(e_high);
      if(f_high>=0.0) {
        found_high = true;
        break;
      }
      e_low = e_high;
      f_low = f_high;
      found_low = true; //not really needed
      e_high *= 10.0;
    }
    if(!found_high) {
      e_high /= 10.0;
      if(verbose>=1)
        fprintf(stdout,"\033[0;35mWarning: In VarFcnTillot::GetInternalEnergyPerUnitMassFromTemperature, setting e to %e,"
                       " for rho = %e, T = %e (fun = %e).\033[0m\n", e_high, rho, T, f_high);
      return e_high;
    }
  }

  //fprintf(stdout,"e_low = %e, e_high = %e.\n", e_low, e_high);

  boost::uintmax_t maxit = 500; //!< "maxit" is both an input and an output!
  double tolerance = tol*std::min(fabs(e_low), fabs(e_high-e_low));
  std::pair<double,double> sol = boost::math::tools::toms748_solve(equation, e_low, e_high, f_low, f_high,
                                 [=](double rr0, double rr1){return fabs(rr1-rr0)<tolerance;},
                                 maxit);

  return 0.5*(sol.first + sol.second);

}

//------------------------------------------------------------------------------

double
VarFcnTillot::GetInternalEnergyPerUnitMassFromEnthalpy(double rho, double h)
{

  EnthalpyEquation equation(rho, h, this);

  // find bracketing interval
  double e_low = h>0 ? 1.0e-2*h: 1.0e-2*e0; //first try positive values.
  double e_high = h>0 ? h : 10.0*e0;
  double f_low(1.0), f_high(1.0);
  bool found_low(false), found_high(false);
  for(int iter=0; iter<15; iter++) {
    f_low = equation(e_low);
    if(f_low<=0.0) {
      found_low = true;
      break;
    }
    e_high = e_low;
    f_high = f_low;
    found_high = true; 

    e_low /= 10.0;
  }
  if(!found_low) {
    e_low = -1e-5*e0; //try negative values.
    for(int iter=0; iter<15; iter++) {
      f_low = equation(e_low);
      if(f_low<=0.0) {
        found_low = true;
        break;
      }
      e_low *= 10.0;
    }
    if(!found_low && verbose>=1) {
      e_low /= 10.0;
      if(verbose>=1)
        fprintf(stdout,"\033[0;35mWarning: In VarFcnTillot::GetInternalEnergyPerUnitMassFromEnthalpy, setting e to %e,"
                       " for rho = %e, h = %e (fun = %e).\033[0m\n", e_low, rho, h, f_low);
      return e_low;
    }
  }

  if(!found_high) {
    for(int iter=0; iter<20; iter++) {
      f_high = equation(e_high);
      if(f_high>=0.0) {
        found_high = true;
        break;
      }
      e_low = e_high;
      f_low = f_high;
      found_low = true; //not really needed
      e_high *= 10.0;
    }
    if(!found_high) {
      e_high /= 10.0;
      if(verbose>=1)
        fprintf(stdout,"\033[0;35mWarning: In VarFcnTillot::GetInternalEnergyPerUnitMassFromEnthalpy, setting e to %e,"
                       " for rho = %e, h = %e (fun = %e).\033[0m\n", e_high, rho, h, f_high);
      return e_high;
    }
  }

  //fprintf(stdout,"e_low = %e, e_high = %e.\n", e_low, e_high);

  boost::uintmax_t maxit = 500; //!< "maxit" is both an input and an output!
  double tolerance = tol*std::min(fabs(e_low), fabs(e_high-e_low));
  std::pair<double,double> sol = boost::math::tools::toms748_solve(equation, e_low, e_high, f_low, f_high,
                                 [=](double rr0, double rr1){return fabs(rr1-rr0)<tolerance;},
                                 maxit);

  return 0.5*(sol.first + sol.second);

}

//------------------------------------------------------------------------------


//------------------------------------------------------------------------------





#endif
