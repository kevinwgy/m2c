/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _VAR_FCN_TILLOT_H_
#define _VAR_FCN_TILLOT_H_

#include <VarFcnBase.h>
#include <polynomial_equations.h>
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
 *   T0     : temperature at rho0 and e0
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

  inline double GetChiWithEta(double eta, double e) {return 1.0/(e/(e0*eta*eta)+1.0);}
  inline double GetChiWithOmega(double omega, double e) {return 1.0/(e/e0*(omega+1.0)*(omega+1)+1.0);}

  //! Determine the case id, given rho and e. Returns 0, 1, 2, or 3 for Case 1, 2, 3, and 1|2
  int GetCaseWithRhoE(double rho, double e) {
    if(rho<=0.0 || e<0.0) {
      fprintf(stdout,"\033[0;31m*** Error: VarFcnTillot::GetCaseWithRhoE detected negative rho (%e) or e (%e).\n\033[0m",
              rho, e);
      exit(-1);
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
    return a*rho_e + (b*rho_e*GetChiWithOmega(omega,e) + A*mu*exp(-beta*mu))*exp(-alpha*omega*omega);
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
    double mu  = eta - 1.0;
    double chi = GetChiWithEta(eta, e);
    return a*e + b*e*chi*chi*(1.0 + 3.0*e/(e0*eta*eta)) + A/rho0;
  }

  double GetInternalEnergyPerUnitMass3(double rho, double p);





  /********************************
   *           Case 1|2 
   *******************************/
  inline double GetPressure12(double rho, double e) {
    return ((eCV-e)*GetPressure1(rho,e) + (e-eIV)*GetPressure2(rho,e))/(eCV-eIV);
  }

  inline double GetGamma12(double rho, double e) {
    return (  (GetPressure2(rho,e)-GetPressure1(rho,e))/rho 
            + (eCV-e)*GetGamma1(rho,e) + (e-eIV)*GetGamma2(rho,e) )/(eCV-eIV);
  }

  inline double GetDpdrho12(double rho, double e) {
    return ((eCV-e)*GetDpdrho1(rho,e) + (e-eIV)*GetDpdrho2(rho,e))/(eCV-eIV);
  }

  double GetInternalEnergyPerUnitMass12(double rho, double p);

  //! (p - p1|2(rho,e))/p = 0 for e, given rho and p. (Divided by p to non-dimensionalize.)
  struct SpecificEnergyEquationP12 {
    SpecificEnergyEquationP12(double rho_, double p_, VarFcnTillot *vf_) 
        : vf(vf_), rho(rho_), p(p_) {}
    inline double operator() (double e) {
      return p==0 ? -(vf->GetPressure12(rho, e)) : 1.0 - (vf->GetPressure12(rho, e))/p;}

  private:
    double rho, p;
    VarFcnTillot *vf;
  }    


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

};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

VarFcnTillot::VarFcnTillot(MaterialModelData &data) : VarFcnBase(data)
{

  if(data.eos != MaterialModelData::TILLOTSON){
    fprintf(stdout, "\033[0;31m*** Error: MaterialModelData is not of type Tillotson.\n\033[0m");
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

  // Brundage 2013 did not clearly say how to integrate e_c(rho) in the "expansion" region...
  // Based on what he said, we would be evaluating p1, p2, p3, and p1|2 with negative e.
  // I am not gonna implement it until I understand it... (KW)
  if(temperature_depends_on_density) {
    fprintf(stdout, "\033[0;31m*** Error: VarFcnTillot does not support "
                    "\"TemperatureDependsOnDensity\" at the moment.\n\033[0m");
    exit(-1); 
  }

  cp     = data.tillotModel.cp;
  h0     = data.tillotModel.h0;

  if(temperature_depends_on_density) {
    use_cp = false;
  } else {
    use_cp = (cp>0 && cv<=0.0) ? true : false;
  }

  invcv = cv==0.0 ? 0.0 : 1.0/cv;
  invcp = cp==0.0 ? 0.0 : 1.0/cp;

  elat = eCV-eIV;

  tol = 1.0e-8; //!< TODO: hard-coded for the moment, non-dimensional.

  if(elat<=0.0) {
    fprintf(stdout, "\033[0;31m*** Error: VarFcnTillot detected eCV (%e) <= eIV (%e).\n\033[0m", eCV, eIV);
    exit(-1);
  }
  if(rhoIV>=rho0) {
    fprintf(stdout, "\033[0;31m*** Error: VarFcnTillot detected rhoIV (%e) >= rho0 (%e).\n\033[0m", rhoIV, rho0);
    exit(-1);
  }
  if(rhoIV<=0.0 || rho0<=0.0) {
    fprintf(stdout, "\033[0;31m*** Error: VarFcnTillot detected non-positive rhoIV (%e) or rho0 (%e).\n\033[0m", rhoIV, rho0);
    exit(-1);
  }

}

//------------------------------------------------------------------------------

double
VarFcnTillot::GetInternalEnergyPerUnitMass1(double rho, double p)
{
  assert(rho>=rhoIV);

  double eta = rho/rho0;
  double mu  = eta - 1.0;
  double inv_e0_eta_eta = 1.0/(e0*eta*eta);

  c0 = -p + A*mu + B*mu*mu
  c1 = c0*inv_e0_eta_eta + (a+b)*rho;
  c2 = a*rho*inv_e0_eta_eta;

  double e_root_1, e_root_2;
  int num_real_roots = MathTools::quadratic_equation_solver(c2, c1, c0, e_root_1, e_root_2);
   
  if(num_real_roots==2) { // e_root_1 > e_root_2 (see quadratic_equation_solver)
    if(e_root_1<0.0) {
      fprintf(stdout, "\033[0;31m*** Error: VarFcnTillot(P1) got two negative e's (%e,%e) for rho = %e, p = %e.\n\033[0m",
              e_root_1, e_root_2, rho, p);
      exit(-1);
    }
    else if(e_root_2>0.0) {
      fprintf(stdout, "\033[0;31m*** Error: VarFcnTillot(P1) got two positive e's (%e,%e) for rho = %e, p = %e.\n\033[0m",
              e_root_1, e_root_2, rho, p);
      exit(-1);
    }
    return e_root_1;
  }
  else if(num_real_roots==1) {
    if(e_root_1<0.0) {
      fprintf(stdout, "\033[0;31m*** Error: VarFcnTillot(P1) got a negative e (%e) for rho = %e, p = %e. \n\033[0m",
              e_root_1, rho, p);
      exit(-1);
    }
    return e_root_1;
  }
  else { // num_real_roots = 0
    fprintf(stdout, "\033[0;31m*** Error: VarFcnTillot(P1) failed to find an e for rho = %e, p = %e. "
            "Solution does not exist.\n\033[0m", rho, p);
    exit(-1);
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

  c0 = -p + A*mu*exp(-beta*omega - alpha*omega*omega);
  c1 = c0*inv_e0_eta_eta + (a + b*exp(-alpha*omega*omega))*rho;
  c2 = a*rho*inv_e0_eta_eta;

  double e_root_1, e_root_2;
  int num_real_roots = MathTools::quadratic_equation_solver(c2, c1, c0, e_root_1, e_root_2);
   
  if(num_real_roots==2) { // e_root_1 > e_root_2 (see quadratic_equation_solver)
    if(e_root_1<0.0) {
      fprintf(stdout, "\033[0;31m*** Error: VarFcnTillot(P2) got two negative e's (%e,%e) for rho = %e, p = %e. \n\033[0m",
              e_root_1, e_root_2, rho, p);
      exit(-1);
    }
    else if(e_root_2>0.0) {
      fprintf(stdout, "\033[0;31m*** Error: VarFcnTillot(P2) got two positive e's (%e,%e) for rho = %e, p = %e. \n\033[0m",
              e_root_1, e_root_2, rho, p);
      exit(-1);
    }
    return e_root_1;
  }
  else if(num_real_roots==1) {
    if(e_root_1<0.0) {
      fprintf(stdout, "\033[0;31m*** Error: VarFcnTillot(P2) got a negative e (%e) for rho = %e, p = %e. \n\033[0m",
              e_root_1, rho, p);
      exit(-1);
    }
    return e_root_1;
  }
  else { // num_real_roots = 0
    fprintf(stdout, "\033[0;31m*** Error: VarFcnTillot(P2) failed to find an e for rho = %e, p = %e. "
            "Solution does not exist.\n\033[0m", rho, p);
    exit(-1);
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

  c0 = -p + A*mu; //This is the only difference from Case 1
  c1 = c0*inv_e0_eta_eta + (a+b)*rho;
  c2 = a*rho*inv_e0_eta_eta;

  double e_root_1, e_root_2;
  int num_real_roots = MathTools::quadratic_equation_solver(c2, c1, c0, e_root_1, e_root_2);
   
  if(num_real_roots==2) { // e_root_1 > e_root_2 (see quadratic_equation_solver)
    if(e_root_1<0.0) {
      fprintf(stdout, "\033[0;31m*** Error: VarFcnTillot(P3) got two negative e's (%e,%e) for rho = %e, p = %e. \n\033[0m",
              e_root_1, e_root_2, rho, p);
      exit(-1);
    }
    else if(e_root_2>0.0) {
      fprintf(stdout, "\033[0;31m*** Error: VarFcnTillot(P3) got two positive e's (%e,%e) for rho = %e, p = %e. \n\033[0m",
              e_root_1, e_root_2, rho, p);
      exit(-1);
    }
    return e_root_1;
  }
  else if(num_real_roots==1) {
    if(e_root_1<0.0) {
      fprintf(stdout, "\033[0;31m*** Error: VarFcnTillot(P3) got a negative e (%e) for rho = %e, p = %e. \n\033[0m",
              e_root_1, rho, p);
      exit(-1);
    }
    return e_root_1;
  }
  else { // num_real_roots = 0
    fprintf(stdout, "\033[0;31m*** Error: VarFcnTillot(P3) failed to find an e for rho = %e, p = %e. "
            "Solution does not exist.\n\033[0m", rho, p);
    exit(-1);
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
  double f_high = equation(f_high);
  if(f_low*f_high>0) {
    fprintf(stdout,"\033[0;31m*** Error: VarFcnTillot::GetInternalEnergyPerUnitMass12 called w. incorrect inputs."
                   " rho = %e, p = %e. f(%e) = %e, f(%e) = %e.\n\033[0m",
            rho, p, e_low, f_low, e_high, f_high);
    exit(-1);
  }

  boost::uintmax_t maxit = 500; //!< "maxit" is both an input and an output!
  double tolerance = tol*e_low;
  std::pair<double,double> sol = toms748_solve(equation, e_low, e_high, f_low, f_high,
                                            [=](double rr0, double rr1){return fabs(rr1-rr0)<tolerance;},
                                            maxit);

  return 0.5*(sol.first + sol.second);

}

//------------------------------------------------------------------------------

double
VarFcnTillot::GetInternalEnergyPerUnitMass(double rho, double p)
{
  if(rho<=0.0) {
    fprintf(stdout,"\033[0;31m*** Error: VarFcnTillot::GetInternalEnergyPerUnitMass detected negative rho (%e).\n\033[0m",
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
                       "rho = %e, p = %e.\n\033[0m", rho, p);
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

  if(e<0) {
    fprintf(stdout,"\033[0;31m*** Error: VarFcnTillot::GetInternalEnergyPerUnitMass got negative "
                       "e (%e) for rho = %e, p = %e.\n\033[0m", e, rho, p);
    exit(-1);
  }

  return e;

}

//------------------------------------------------------------------------------

double
VarFcnTillot::GetDensity(double p, double e)
{
  //TODO: This function is not really needed at the moment. Therefore, it is not implemented.

  fprintf(stdout,"\033[0;31m*** Error: VarFcnTillot::GetDensity has not been implemented.\n\033[0m");
  exit(-1);

  return 0.0;
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

  // Brundage 2013 did not clearly say how to integrate e_c(rho) in the "expansion" region...
  // Based on what he said, we would be evaluating p1, p2, p3, and p1|2 with negative e.
  // I am not gonna implement it until I understand it... (KW)
  return 0.0;

}

//------------------------------------------------------------------------------

double
VarFcnTillot::GetInternalEnergyPerUnitMassFromTemperature(double rho, double T)
{

  assert(!temperature_depends_on_density);

  if(use_cp) {
    double e_plus_p_over_rho = h0 + cp*(T-T0);

  } else
    return e0 + cv*(T-T0);
}

//------------------------------------------------------------------------------





#endif
