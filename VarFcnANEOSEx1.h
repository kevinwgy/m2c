#ifndef _VAR_FCN_ANEOS_EX1_H_
#define _VAR_FCN_ANEOS_EX1_H_

#include<VarFcnANEOSBase.h>
#include<polylogarithm_function.h>

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
  
  double tol_Debye; //!< error tolerance (hard-coded for the moment)

  double r0; //!< reference density at 0 Kelvin
  double b0; //!< reference bulk modulus
  double b0prime; //!< derivative of b0 w.r.t. pressure
  double b0_times_9eighths; //!< b0*9/8
  double b0prime_minus_4; //!< b0prime - 4
  double delta_e; //!< energy shift constant (often set to 0)

  std::vector<tuple<double,double,double> > rho_e_T;

public:
  Don't forget to destroy the pointers!

  //----- EOS-Specific Functions -----//

  // Compute pressure: p = rho*rho*(\partial F(rho,T))/(\partial rho)
  inline double GetPressure(double rho, e) const {
    double T = GetTemperature(rho, e);

    I AM HERE



  }

  //! temperature law, defined naturally by the EOS. (ANEOS is a complete EOS.)
  double GetTemperature(double rho, double e) const{
    
    double T(0.0);
    if(GetTemperatureFromHistory(rho,e,T))
      return T;

    double el = e - delta_e - ComputeZeroKelvinSpecificEnergy(rho);

    // -------------------------------
    // solve a nonlinear equation to find T
    // -------------------------------
    ThermalVibrationEnergyEquation equation(rho, el, this);
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
        fprintf(stderr,"\033[0;35mWarning: In VarFcnANEOSEx1::GetTemperature, setting temperature to %e,"
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
        fprintf(stderr,"\033[0;35mWarning: In VarFcnANEOSEx1::GetTemperature, setting temperature to %e,"
                       " for rho = %e, e = %e (fun = %e).\033[0m\n", T_high, rho, e, f_high);
      return T_high;
    } 

    boost::uintmax_t maxit = 500; //!< "maxit" is both an input and an output!
    double tol = 10.0*tol_Debye; 
    std::pair<double,double> sol = toms748_solve(equation, T_low, T_high, f_low, f_high,
                                              [=](double r0, double r1){return fabs(r1-r0)<tol;},
                                              maxit); 

    T = 0.5*(sol.first + sol.second);
    UpdateTemperatureHistory(rho, e, T);

    return T;
  }


  //! temperature law, defined separately for each EOS
  virtual double GetReferenceTemperature() const{
    fprintf(stderr,"\033[0;31m*** Error:  GetReferenceTemperature Function not defined\n\033[0m");
    exit(-1); return 0.0;}

  //! temperature law, defined separately for each EOS
  virtual double GetInternalEnergyPerUnitMassFromTemperature(double rho, double T) const{
    fprintf(stderr,"\033[0;31m*** Error:  GetInternalEnergyPerUnitMassFromTemperature Function not defined\n\033[0m");
    exit(-1); return 0.0;}

  //! calculate e from rho and h
  virtual double GetInternalEnergyPerUnitMassFromEnthalpy(double rho, double h) const{
    fprintf(stderr,"\033[0;31m*** Error:  GetInternalEnergyPerUnitMassFromEnthalpy Function not defined\n\033[0m");
    exit(-1); return 0.0;}

  //checks that the Euler equations are still hyperbolic
  virtual bool CheckState(double rho, double p, bool silence = false) const{
    if(!std::isfinite(rho) || !std::isfinite(p)) {
      if(!silence)
        fprintf(stderr, "*** Error: CheckState failed. rho = %e, p = %e.\n\033[0m", rho, p);
      return true;
    }
    if(rho <= 0.0) {
      if(!silence && verbose>1)
        fprintf(stdout, "Warning: Negative density or violation of hyperbolicity. rho = %e, p = %e.\n", rho, p);
      return true;
    }
    double e = GetInternalEnergyPerUnitMass(rho,p);
    double c2 = GetDpdrho(rho, e) + p/rho*GetBigGamma(rho, e);
    if(c2<=0){
      if(!silence && verbose>1)
        fprintf(stdout, "Warning: Negative density or violation of hyperbolicity. rho = %e, p = %e.\n", rho, p);
      return true;
    }
    return false;
  }

  //checks that the Euler equations are still hyperbolic
  virtual bool CheckState(double *V, bool silence = false) const{
    if(!std::isfinite(V[0]) || !std::isfinite(V[1]) || !std::isfinite(V[2]) || !std::isfinite(V[3]) || !std::isfinite(V[4])){
      if(!silence)
        fprintf(stderr, "\033[0;31m*** Error: CheckState failed. V = %e %e %e %e %e\n\033[0m", V[0], V[1], V[2], V[3], V[4]);
      return true;
    }
    return CheckState(V[0], V[4]); 
  }
 
  //check for phase transitions
  virtual bool CheckPhaseTransition(int id/*id of the other phase*/) const{
    return false; //by default, phase transition is not allowed/considered
  }

  //----- Transformation Operators -----//
  virtual void ConservativeToPrimitive(double *U, double *V); 
  virtual void PrimitiveToConservative(double *V, double *U);

  //----- General Functions -----//
  inline int GetType() const{ return type; }

  virtual double ComputeSoundSpeed(double rho, double e);
  virtual double ComputeSoundSpeedSquare(double rho, double e); //!< this one does not crash on negative c^2
  virtual double ComputeMachNumber(double *V);
  virtual double ComputeEnthalpyPerUnitMass(double rho, double p); //!< h = e + p/rho
  virtual double ComputeTotalEnthalpyPerUnitMass(double *V); //!< H = 1/rho*(E + p)

  // Clipping
  virtual bool ClipDensityAndPressure(double *V, double *U = 0);


protected:

  //solve e_l(rho,T) = ... Debye model ... for T, given rho and el.
  struct ThermalVibrationEnergyEquation { 

    ThermalVibrationEnergyEquation(double rho_, double el_, VarFcnANEOSEx1 *vf_) : vf(vf_), rho(rho_), el(el_) {
      Theta = vf->T0*pow(rho/(vf->rho0), vf->Gamma0); 
      Theta_times_9eights = 1.125*Theta;
    }

    inline double operator() (double T) {
      return el - vf->R_over_w*(Theta_times_9eights + 3.0*T*vf->EvaluateDebyeFunction(Theta/T));
    }

  private:
    double rho, el, Theta, Theta_times_9eights;
    VarFcnANEOSEx1 *vf;

  };


private:

  //! calculates e_{cold} using the Birch-Murnaghan EOS.
  inline double ComputeZeroKelvinSpecificEnergy(double rho) { 
    double x = pow(rho/r0, 2.0/3.0) - 1.0;
    return b0_times_9eighths*x*x*(0.5*x*b0prime_minus_4 + 1.0); //Birch-Murnaghan eq.
  }

  //! evaluates Debye function D(x)
  inline double EvaluateDebyeFunction(double x) {
    return (spline_Li4 && spline_Li3 && spline_Li2) ? EvaluateDebyeFunctionByInterpolation(x)
                                                    : EvaluateDebyeFunctionOnTheFly(x);
  }

  //! evaluates Debye function D(x) on the fly
  inline double EvaluateDebyeFunctionOnTheFly(double x) {
    assert(x>0.0); 
    double expmx = exp(-x);
    return 3.0/(x*x*x)*(-6.0*MathTools::polylogarithm_function(4, expmx, 100, tol_Debye)
                        - 6.0*x*MathTools::polylogarithm_function(3, expmx, 100, tol_Debye)
                        - 3.0*x*x*MathTools::polylogarithm_function(2, expmx, 100, tol_Debye)
                        + x*x*x*log(1.0 - expmx) + pi4_over_15);

  }

  //! evaluates Debye function D(x) by spline interpolation
  inline double EvaluateDebyeFunctionByInterpolation(double x) {
    assert(x>0.0); 
    double expmx = exp(-x);
    if(expmx<splines_expmx_min || expmx>splines_expmx_max)
      return EvaluateDebyeFunctionOnTheFly(x);
    return 3.0/(x*x*x)*(-6.0*(*spline_Li4)(expmx) - 6*x*(*spline_Li3)(expmx) - 3*x*x*(*spline_Li2)(expm;
  }
 
  //! build spline interpolation for Debye function D(x)
  void InitializeInterpolationForDebyeFunction(double expmx_min, double expmx_max, int sample_size) {
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
    spline_Li4 = new boost::math::cubic_b_spline<double>(Li4.begin(), Li4.end(), splines_expmx_min, delta)
    spline_Li3 = new boost::math::cubic_b_spline<double>(Li3.begin(), Li3.end(), splines_expmx_min, delta)
    spline_Li2 = new boost::math::cubic_b_spline<double>(Li2.begin(), Li2.end(), splines_expmx_min, delta)
  }

  //! Get temperature from storage
  bool GetTemperatureFromHistory(double rho, double e, double &T) {
    for(auto&& mytuple : rho_e_T)
      if(std::get<0>(mytuple) == rho && std::get<1>(mytuple) == e) {
        T = std::get<2>(mytuple);
        return true;
      }
    return false;
  }

  //! Update rho_e_T
  void UpdateTemperatureHistory(double rho, double e, double T) {
    for(int i=0; i<rho_e_T.size()-1; i++)
      rho_e_T[i] = rho_e_T[i+1];
    rho_e_T.back() = std::make_tuple(rho,e,T);
  }

};


#endif
