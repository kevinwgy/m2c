/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _VAR_FCN_ANEOS_EX1_H_
#define _VAR_FCN_ANEOS_EX1_H_

#include<VarFcnANEOSBase.h>
#include<polylogarithm_function.h>
#include<tuple>
#include<boost/math/tools/roots.hpp>
#include<boost/math/interpolators/cubic_b_spline.hpp>  //spline interpolation

extern double avogadro_number;

/****************************************************************************
 * This class is an example of ANEOS, given in Section 2.2 of J.J. Sanchez (Sandia),
 * CMAME 2021.
 * e(rho,T) = e_c(rho) + e_l(rho,T) + \Delta e.
 * F(rho,T) = e_c(rho) + F_l(rho,T) + \Delta e. (specific Helmholtz free energy)
 * e_c is modeled by the Birch-Murnaghan equation.
 * e_l is modeled using the Debye model. An analytical formula is used to
 * calculate/approximate the Debye function, which would otherwise require
 * numerical integration.
 ***************************************************************************/

class VarFcnANEOSEx1 : public VarFcnANEOSBase {

private:
  
  //! Physical parameters
  double r0; //!< reference density at 0 Kelvin
  double b0; //!< reference bulk modulus
  double b0prime; //!< derivative of b0 w.r.t. pressure
  double delta_e; //!< energy shift constant (often set to 0)
  double R_over_w; //!< gas constant (avogadro*boltzmann) divided by molar mass
  double T0; //!< reference temperature (ambient state)
  double e0; //!< reference internal energy (ambient state) --- only for output 
  double Gamma0; //!< reference Gruneisen parameter 
  double rho0; //!< reference density (ambient state)

  double pi4_over_15; //!< pi^4/15

  //! Calculated values (latest few)
  std::vector<std::tuple<double,double,double> > rho_e_T;
  std::vector<std::tuple<double,double,double,double> > rho_e_p_T;

  //! cubic splines for interpolating the polylogarithm functions involved in Debye function
  double splines_expmx_min, splines_expmx_max; //!< min and max of exp(-x) of the sample points
  std::vector<double> Li2, Li3, Li4; //!< sample values 
  boost::math::cubic_b_spline<double> *spline_Li2, *spline_Li3, *spline_Li4;
  
  //! Numerical parameters
  double tol_Debye; //!< non-dimensional error tolerance (hard-coded for the moment)

public:

  VarFcnANEOSEx1(MaterialModelData &data);
  ~VarFcnANEOSEx1();

  double GetPressure(double rho, double e);
  double GetInternalEnergyPerUnitMass(double rho, double p);
  double GetDensity(double p, double e);
  //double GetDpdrho(double rho, double e);   //!< Using finite difference, defined in base class
  //double GetBigGamma(double rho, double e); //!< Using finite difference, defined in base class
  double GetTemperature(double rho, double e);
  inline double GetReferenceTemperature() {return T0;} //!< reference temperature (ambient state)
  inline double GetReferenceInternalEnergyPerUnitMass() {return e0;} //!< ambient state 
  double GetInternalEnergyPerUnitMassFromTemperature(double rho, double T);
  double GetInternalEnergyPerUnitMassFromEnthalpy(double rho, double h);

private:

  //! calculate Debye temperature Theta(rho)
  inline double ComputeDebyeTemperature(double rho) {return T0*pow(rho/rho0, Gamma0);}

  //! calculate Theta'(rho)
  inline double ComputeDebyeTemperatureDerivative(double rho) {
    return T0*Gamma0/rho0*pow(rho/rho0, Gamma0-1.0);}

  //! calculate e_{cold} using the Birch-Murnaghan EOS.
  inline double ComputeColdSpecificEnergy(double rho) {
    double x = pow(rho/r0, 2.0/3.0) - 1.0;
    return 1.125*b0/r0*x*x*(0.5*x*(b0prime-4.0) + 1.0);} //Birch-Murnaghan eq.

  //! calculate d(e_{cold})/d(\rho)
  inline double ComputeColdSpecificEnergyDerivative(double rho) {
    double x = pow(rho/r0, 2.0/3.0) - 1.0;
    return 1.5*(0.75*(b0prime-4.0)*x*x + x)*b0/(r0*r0)*pow(rho/r0, -1.0/3.0);}
  
  //! calculate e_l
  inline double ComputeThermalSpecificEnergy(double rho, double T) {
    double Theta = ComputeDebyeTemperature(rho);   assert(T>0);
    return R_over_w*(1.125*Theta + 3.0*T*EvaluateDebyeFunction(Theta/T));}

  //! calculate F_l(rho,T)
  double ComputeThermalSpecificHelmholtz(double rho, double T) {
    double Theta = ComputeDebyeTemperature(rho);   assert(T>0);
    double Theta_over_T = Theta/T;
    return R_over_w*(1.125*Theta + 3.0*T*log(1.0-exp(-Theta_over_T)) - T*EvaluateDebyeFunction(Theta_over_T));}

  //! calculates dF_l(rho,T)/drho
  double ComputeThermalSpecificHelmholtzDerivativeRho(double rho, double T) {
    double ThetaPrime = ComputeDebyeTemperatureDerivative(rho);
    double Theta_over_T = ComputeDebyeTemperature(rho)/T;
    assert(Theta_over_T>0);
    return R_over_w*ThetaPrime*(1.125 + 3.0/Theta_over_T*EvaluateDebyeFunction(Theta_over_T));}

  //! calculate dF_l(rho,T)/dT
  inline double ComputeThermalSpecificHelmholtzDerivativeT(double rho, double T) {
    double Theta_over_T = ComputeDebyeTemperature(rho)/T;   assert(Theta_over_T>0.0);
    return R_over_w*(3.0*log(1.0 - exp(-Theta_over_T)) - 4.0*EvaluateDebyeFunction(Theta_over_T));}

  //! evaluate Debye function D(x)
  inline double EvaluateDebyeFunction(double x) {
    return (spline_Li4 && spline_Li3 && spline_Li2) ? EvaluateDebyeFunctionByInterpolation(x)
                                                    : EvaluateDebyeFunctionOnTheFly(x);}

  //! evaluate D'(x) (Note: Can also be done by differentiating the splines (i.e. spline.prime(x)))
  inline double EvaluateDebyeFunctionDerivative(double x) {
    double expmx = exp(-x);
    assert(expmx != 1.0);
    return -3.0/x*EvaluateDebyeFunction(x) + 3.0*expmx/(1.0-expmx);}

  //! evaluate Debye function D(x) on the fly
  double EvaluateDebyeFunctionOnTheFly(double x);
  
  //! evaluate Debye function D(x) by spline interpolation
  double EvaluateDebyeFunctionByInterpolation(double x);

  //! build spline interpolation for Debye function D(x)
  void InitializeInterpolationForDebyeFunction(double expmx_min, double expmx_max, int sample_size);

  //! Update rho_e_T
  void Update_rho_e_T(double rho, double e, double T) {
    for(int i=0; i<(int)rho_e_T.size()-1; i++)
      rho_e_T[i] = rho_e_T[i+1];
    rho_e_T.back() = std::make_tuple(rho,e,T);
  }

  //! Update rho_e_p_T
  void Update_rho_e_p_T(double rho, double e, double p, double T) {
    for(int i=0; i<(int)rho_e_p_T.size()-1; i++)
      rho_e_p_T[i] = rho_e_p_T[i+1];
    rho_e_p_T.back() = std::make_tuple(rho,e,p,T);
  } 


protected:

  //---------------------------------------------------------------------
  //! Solve p(rho,T)/(rho*rho) - e_cold'(rho) - dF_l(rho,T)/drho = 0 for T
  struct ThermalHelmholtzDerivativeRhoEquation {

    ThermalHelmholtzDerivativeRhoEquation(double rho_, double dFl_drho_, VarFcnANEOSEx1 *vf_) 
        : vf(vf_), rho(rho_), dFl_drho(dFl_drho_) {
      Theta = vf->ComputeDebyeTemperature(rho);    
      ThetaPrime = vf->ComputeDebyeTemperatureDerivative(rho);
    }

    inline double operator() (double T) { //!< similar to ComputeThermalSpecificHelmholtzDerivativeRho
      assert(T>0);
      double Theta_over_T = Theta/T;
      return dFl_drho - (vf->R_over_w)*ThetaPrime*(1.125 + 3.0/Theta_over_T*vf->EvaluateDebyeFunction(Theta_over_T)); 
    }

  private:
    double rho, dFl_drho, Theta, ThetaPrime;
    VarFcnANEOSEx1 *vf;

  };

  //---------------------------------------------------------------------
  //! Solve el - e_l(rho,T) = 0, given rho and el.
  struct ThermalSpecificEnergyEquation { 

    ThermalSpecificEnergyEquation(double rho_, double el_, VarFcnANEOSEx1 *vf_) : vf(vf_), rho(rho_), el(el_) {
      Theta = vf->ComputeDebyeTemperature(rho);
      Theta_times_9eights = 1.125*Theta;
    }

    inline double operator() (double T) { //!< similar to ComputeThermalSpecificEnergy
      return el - vf->R_over_w*(Theta_times_9eights + 3.0*T*vf->EvaluateDebyeFunction(Theta/T));
    }

  private:
    double rho, el, Theta, Theta_times_9eights;
    VarFcnANEOSEx1 *vf;

  };


  //---------------------------------------------------------------------
  //! Solve (h - p(rho,T) - rho/e(rho,T) = 0) for T    
  struct SpecificEnthalpyEquation {

    SpecificEnthalpyEquation(double rho_, double h_, VarFcnANEOSEx1 *vf_) : vf(vf_), rho(rho_), h(h_) {
      Theta = vf->ComputeDebyeTemperature(rho);
      ThetaPrime = vf->ComputeDebyeTemperatureDerivative(rho);
      assert(Theta != 0.0);
      ThetaPrime_over_Theta = ThetaPrime/Theta;
      e_cold = vf->ComputeColdSpecificEnergy(rho);
      e_cold_prime = vf->ComputeColdSpecificEnergyDerivative(rho);
    }

    double operator() (double T) {
      assert(T>0.0);
      double Theta_over_T = Theta/T;
      double el = (vf->R_over_w)*(1.125*Theta + 3.0*T*(vf->EvaluateDebyeFunction(Theta_over_T)));
      double p = rho*rho*(e_cold_prime + ThetaPrime_over_Theta*el);
      double e = e_cold + el + vf->delta_e;
      assert(e!=0.0);
      return h - p - rho/e;
    }

  private:
    double rho, h, Theta, ThetaPrime, ThetaPrime_over_Theta, e_cold, e_cold_prime;
    VarFcnANEOSEx1 *vf;

  };


};


//---------------------------------------------------------------------
//---------------------------------------------------------------------
//! Constructor
VarFcnANEOSEx1::VarFcnANEOSEx1(MaterialModelData &data) : VarFcnANEOSBase(data)
{

  if(data.eos != MaterialModelData::ANEOS_BIRCH_MURNAGHAN_DEBYE) {
    fprintf(stdout, "*** Error: MaterialModelData is not of type ANEOS_BIRCH_MURNAGHAN_DEBYE.\n");
    exit(-1);
  }

  type = ANEOS_BIRCH_MURNAGHAN_DEBYE;

  r0 = data.abmdModel.zeroKelvinDensity;
  b0 = data.abmdModel.b0;
  b0prime = data.abmdModel.b0prime;
  delta_e = data.abmdModel.delta_e;
  R_over_w = avogadro_number*data.abmdModel.boltzmann_constant/data.abmdModel.molar_mass;
  T0 = data.abmdModel.T0; 
  e0 = data.abmdModel.e0; 
  Gamma0 = data.abmdModel.Gamma0;
  rho0 = data.abmdModel.rho0;

  pi4_over_15 = pow(2.0*acos(0.0), 4)/15.0;

  tol_Debye = 1.0e-8;

  if(data.abmdModel.debye_evaluation == ANEOSBirchMurnaghanDebyeModelData::CUBIC_SPLINE_INTERPOLATION)
    InitializeInterpolationForDebyeFunction(1.0e-20, 0.999999, 15000);
  else {
    spline_Li2 = NULL;
    spline_Li3 = NULL;
    spline_Li4 = NULL;
    splines_expmx_min = DBL_MAX;
    splines_expmx_max = DBL_MAX;
  }

  // store the latest 4 solutions (can be changed)
  rho_e_T.resize(4);
  rho_e_p_T.resize(4);

}

//---------------------------------------------------------------------
//! Destructor
VarFcnANEOSEx1::~VarFcnANEOSEx1()
{
  if(spline_Li2) delete spline_Li2;
  if(spline_Li3) delete spline_Li3;
  if(spline_Li4) delete spline_Li4;
}

//---------------------------------------------------------------------
//! Compute pressure: p = rho*rho*(\partial F(rho,T))/(\partial rho)
double 
VarFcnANEOSEx1::GetPressure(double rho, double e) 
{
  double T = GetTemperature(rho, e);
  return rho*rho*(ComputeColdSpecificEnergyDerivative(rho) + 
                  ComputeThermalSpecificHelmholtzDerivativeRho(rho,T));
}

//---------------------------------------------------------------------
//! Compute e from rho and p
double
VarFcnANEOSEx1::GetInternalEnergyPerUnitMass(double rho, double p)
{

  // Check storage
  for(auto&& mytuple : rho_e_p_T)
    if(std::get<0>(mytuple) == rho && std::get<2>(mytuple) == p)
      return std::get<1>(mytuple);

  // -------------------------------------------------------------------------------
  // solve a nonlinear equation (p(rho,T)/(rho*rho) - e_cold'(rho) - dF_l(rho,T)/drho = 0) to find T
  // -------------------------------------------------------------------------------
  assert(rho>0.0);
  double e_cold_prime = ComputeColdSpecificEnergyDerivative(rho);
  double dFl_drho = p/(rho*rho) - e_cold_prime;
  //fprintf(stdout,"e_cold_prime = %e, dFl_drho = %e.\n", e_cold_prime, dFl_drho);
  ThermalHelmholtzDerivativeRhoEquation equation(rho, dFl_drho, this);

  // find bracketing interval
  double T_low(500.0), T_high(5000.0);
  double f_low  = equation(T_low);
  double f_high = equation(T_high);
  int iter = 0;
  while(f_low*f_high>0.0) {
    //fprintf(stdout,"f(%e) = %e, f(%e) = %e.\n", T_low, f_low, T_high, f_high);
    if(++iter>10)
      break;
    T_low  /= 2.0;
    T_high *= 4.0;
    f_low  = equation(T_low);
    f_high = equation(T_high);
  }
  if(iter>10) {
    fprintf(stdout,"\033[0;31m*** Error: In VarFcnANEOSEx1, unable to locate temperature between %e and %e."
                   " (rho = %e, p = %e).\033[0m\n", T_low, T_high, rho, p);
    exit(-1);
  }

  boost::uintmax_t maxit = 500; //!< "maxit" is both an input and an output!
  double tol = tol_Debye*T_low; 
  std::pair<double,double> sol = boost::math::tools::toms748_solve(equation, T_low, T_high, f_low, f_high,
                                 [=](double rr0, double rr1){return fabs(rr1-rr0)<tol;},
                                 maxit); 

  double T = 0.5*(sol.first + sol.second);

  // Calculates e from rho and T
  double e = ComputeColdSpecificEnergy(rho) + ComputeThermalSpecificEnergy(rho,T) + delta_e;

  Update_rho_e_p_T(rho, e, p, T);

  return e;
}

//---------------------------------------------------------------------
//! Calculates rho from p and e
double
VarFcnANEOSEx1::GetDensity(double p, double e) 
{

  // Check storage
  for(auto&& mytuple : rho_e_p_T)
    if(std::get<1>(mytuple) == e && std::get<2>(mytuple) == p)
      return std::get<0>(mytuple);
  
  //TODO: This function is not really needed at the moment. Therefore, it is not implemented

  fprintf(stdout,"\033[0;31m*** Error: Implementation of VarFcnANEOSEx1::GetDensity is incomplete.\n\033[0m");
  exit(-1); 
  return 0.0;
}

//---------------------------------------------------------------------
//! Temperature law, defined naturally by the EOS. (ANEOS is a complete EOS.)
double
VarFcnANEOSEx1::GetTemperature(double rho, double e)
{
  
  // Check storage
  for(auto&& mytuple : rho_e_T)
    if(std::get<0>(mytuple) == rho && std::get<1>(mytuple) == e)
      return std::get<2>(mytuple);
  for(auto&& mytuple : rho_e_p_T)
    if(std::get<0>(mytuple) == rho && std::get<1>(mytuple) == e)
      return std::get<3>(mytuple);

  // -------------------------------
  // solve a nonlinear equation (e - e_cold(rho) - delta_e - el(rho,T) = 0) to find T
  // -------------------------------
  double el = e - delta_e - ComputeColdSpecificEnergy(rho);
  ThermalSpecificEnergyEquation equation(rho, el, this);
  // find bracketing interval
  double T_low = 100.0; 
  double f_low(1.0);
  bool found = false;
  for(int iter=0; iter<5; iter++) {
    f_low = equation(T_low);
    if(f_low>=0.0) {
      found = true;
      break;
    }
    T_low /= 4.0; 
  }
  if(!found) {
    T_low *= 4.0;
    if(verbose>=1)
      fprintf(stdout,"\033[0;35mWarning: In VarFcnANEOSEx1::GetTemperature, setting temperature to %e,"
                     " for rho = %e, e = %e (fun = %e).\033[0m\n", T_low, rho, e, f_low);
    return T_low;
  } 

  double T_high = 2000.0; 
  double f_high(1.0);
  found = false;
  for(int iter=0; iter<8; iter++) {
    f_high = equation(T_high);
    if(f_high<=0.0) {
      found = true;
      break;
    }
    T_high *= 4.0; 
  }
  if(!found) {
    T_high /= 4.0;
    if(verbose>=1)
      fprintf(stdout,"\033[0;35mWarning: In VarFcnANEOSEx1::GetTemperature, setting temperature to %e,"
                     " for rho = %e, e = %e (fun = %e).\033[0m\n", T_high, rho, e, f_high);
    return T_high;
  } 

  boost::uintmax_t maxit = 500; //!< "maxit" is both an input and an output!
  double tol = tol_Debye*T_low; 
  std::pair<double,double> sol = boost::math::tools::toms748_solve(equation, T_low, T_high, f_low, f_high,
                                 [=](double rr0, double rr1){return fabs(rr1-rr0)<tol;},
                                 maxit); 

  double T = 0.5*(sol.first + sol.second);
  Update_rho_e_T(rho, e, T);

  return T;
}

//---------------------------------------------------------------------
//! Computes e from rho and T
double
VarFcnANEOSEx1::GetInternalEnergyPerUnitMassFromTemperature(double rho, double T)
{

  // Check storage
  for(auto&& mytuple : rho_e_T)
    if(std::get<0>(mytuple) == rho && std::get<2>(mytuple) == T)
      return std::get<1>(mytuple);
  for(auto&& mytuple : rho_e_p_T)
    if(std::get<0>(mytuple) == rho && std::get<3>(mytuple) == T)
      return std::get<1>(mytuple);

  // Computation
  double e = ComputeColdSpecificEnergy(rho) + ComputeThermalSpecificEnergy(rho,T) + delta_e;

  Update_rho_e_T(rho, e, T);

  return e;
}

//---------------------------------------------------------------------
//! Computes e from rho and h
double
VarFcnANEOSEx1::GetInternalEnergyPerUnitMassFromEnthalpy(double rho, double h)
{

  // ---------------------------------------------------------------------
  // Solve a nonlinear equation (h - p(rho,T) - rho/e(rho,T) = 0) for T    
  // ---------------------------------------------------------------------
  SpecificEnthalpyEquation equation(rho, h, this);

  // find bracketing interval
  double T_low(100.0), T_high(2000.0);
  double f_low  = equation(T_low);
  double f_high = equation(T_high);
  int iter = 0;
  while(f_low*f_high>0.0) {
    if(++iter>10)
      break;
    T_low  /= 2.0;
    T_high *= 4.0;
    f_low  = equation(T_low);
    f_high = equation(T_high);
  }
  if(iter>10) {
    fprintf(stdout,"\033[0;31m*** Error: In VarFcnANEOSEx1, unable to locate temperature between %e and %e."
                   " (rho = %e, h = %e).\033[0m\n", T_low, T_high, rho, h);
    exit(-1);
  }

  boost::uintmax_t maxit = 500; //!< "maxit" is both an input and an output!
  double tol = tol_Debye*T_low; 
  std::pair<double,double> sol = boost::math::tools::toms748_solve(equation, T_low, T_high, f_low, f_high,
                                 [=](double rr0, double rr1){return fabs(rr1-rr0)<tol;},
                                 maxit); 

  double T = 0.5*(sol.first + sol.second);

  return  ComputeColdSpecificEnergyDerivative(rho) + ComputeThermalSpecificEnergy(rho,T) + delta_e;
}

//---------------------------------------------------------------------
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------

//---------------------------------------------------------------------
//! evaluate Debye function D(x) on the fly
double 
VarFcnANEOSEx1::EvaluateDebyeFunctionOnTheFly(double x) 
{
  assert(x>0.0); 
  double expmx = exp(-x);
  return 3.0/(x*x*x)*(-6.0*MathTools::polylogarithm_function(4, expmx, 100, tol_Debye)
                      - 6.0*x*MathTools::polylogarithm_function(3, expmx, 100, tol_Debye)
                      - 3.0*x*x*MathTools::polylogarithm_function(2, expmx, 100, tol_Debye)
                      + x*x*x*log(1.0 - expmx) + pi4_over_15);

}

//---------------------------------------------------------------------
//! evaluate Debye function D(x) by spline interpolation
double
VarFcnANEOSEx1::EvaluateDebyeFunctionByInterpolation(double x) 
{
  assert(x>0.0); 
  double expmx = exp(-x);
  if(expmx<splines_expmx_min || expmx>splines_expmx_max)
    return EvaluateDebyeFunctionOnTheFly(x);
  return 3.0/(x*x*x)*(-6.0*(*spline_Li4)(expmx) - 6.0*x*(*spline_Li3)(expmx) - 3.0*x*x*(*spline_Li2)(expmx)
                      + x*x*x*log(1.0 - expmx) + pi4_over_15);
}

//---------------------------------------------------------------------
//! build spline interpolation for Debye function D(x)
void 
VarFcnANEOSEx1::InitializeInterpolationForDebyeFunction(double expmx_min, double expmx_max, 
                                                        int sample_size)
{
  splines_expmx_min = expmx_min;
  splines_expmx_max = expmx_max;
  double delta = (splines_expmx_max - splines_expmx_min)/(sample_size-1.0);
  Li4.resize(sample_size);
  Li3.resize(sample_size);
  Li2.resize(sample_size);
  double expmx;
  for(int i=0; i<sample_size; i++) {
    expmx = expmx_min + i*delta; 
    Li4[i] = MathTools::polylogarithm_function(4, expmx, 100, tol_Debye);
    Li3[i] = MathTools::polylogarithm_function(3, expmx, 100, tol_Debye);
    Li2[i] = MathTools::polylogarithm_function(2, expmx, 100, tol_Debye);
  }
  spline_Li4 = new boost::math::cubic_b_spline<double>(Li4.begin(), Li4.end(), splines_expmx_min, delta);
  spline_Li3 = new boost::math::cubic_b_spline<double>(Li3.begin(), Li3.end(), splines_expmx_min, delta);
  spline_Li2 = new boost::math::cubic_b_spline<double>(Li2.begin(), Li2.end(), splines_expmx_min, delta);
}

//---------------------------------------------------------------------

#endif
