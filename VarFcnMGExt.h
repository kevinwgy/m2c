/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _VAR_FCN_MG_EXT_H
#define _VAR_FCN_MG_EXT_H

#include <VarFcnBase.h>
#include <ordinary_differential_equations.h>
//#include <boost/math/interpolators/cubic_b_spline.hpp>  //spline interpolation
#include <fstream>
#include <cassert>

extern int verbose;

/********************************************************************************
 * This class is the VarFcn class for the *extended* Mie-Gruneisen equation of state (EOS)
 * that applies different equations for compression and tension, see A. Robinson (SNL Report, 2019).
 * It also supports the temperature law presented there, which is consistent with the first law of
 * thermodynamics.
 *
 * Only elementary functions are declared and/or defined here.
 * All arguments must be pertinent to only a single grid node or a single
 * state.
 *
 * Pressure equation: 
 *      / rho0*c0*c0*eta*(1-Gamma0/2*eta)/(1-s*eta)^2     + rho0*Gamma0*(e-e0) ... if eta>=0,
 *  p = | rho0*c0*c0*eta*(1-Gamma0/2*eta)                 + rho0*Gamma0*(e-e0) ... if eta_min<=eta<0,
 *      \ rho0*c0*c0*eta_min*(1-Gamma0/2*(2*eta-eta_min)) + rho0*Gamma0*(e-e0) ... if eta<eta_min.
 *
 * Temperature equation:
 *  - Option 0: T = (e - e_R(eta))/cv + T_R(eta), where
 *
 *               / c0*c0*eta*eta/(2(1-s*eta)^2) + e0                               ... if eta>=0,
 *    e_R(eta) = | c0*c0*eta*eta/2              + e0                               ... if eta_min<=eta<0,
 *               \ c0*c0*eta_min*eta_min/2      + e0 + c0*c0*eta_min*(eta-eta_mim) ... if eta<eta_min.
 *
 *               / exp(Gamm0*eta)*(T0 + c0*c0*s/cv*Int(0,eta)[exp(Gamm0*z)*z*z/(1-s*z)^3] ... if eta>=0,
 *    T_R(eta) = | T0*exp(Gamma0*eta) ... if eta_min<=eta<0,
 *               \ T0*exp(Gamma0*eta) ... if eta<eta_min.
 *
 *  - Option 1 (Simplified): T = (e - e0)/cv + T0;
 *
 *  - Option 2 (Simplified): T = (e + p/rho - h0)/cp + T0;
 *
 *  - Activation Mechanism
 *    o Tlaw = OriginalCv   ==> Option 0,
 *    o Tlaw = SimplifiedCv ==> Option 1,
 *    o Tlaw = SimplifiedCp ==> Option 2.
 *
 *
 * Variables and Parameters: 
 *
 *   e      : internal energy per unit mass.
 *   eta    : 1 - rho0/rho (volumetric strain)
 *   rho    : density
 *   p      : pressure
 *
 *   rho0   : ref. density
 *   c0     : bulk speed of sound
 *   Gamma0 : Gruneisen coefficient at ref. state
 *   s      : slope of shock Hugoniot
 *   e0     : internal energy at ref. state 
 *   T0     : temperature (K) at ref. state
 *   eta_min: a negative value (tension) that triggers special treatment (not "min eta"!)
 *   cv     : specific heat capacity at constant volume
 *   cp     : specific heat capacity at constant pressure
 *
 * References: KW's notes,  Allen Robinson (SNL)'s technical report (2019)
 ********************************************************************************/

class VarFcnMGExt : public VarFcnBase {

private:

  double rho0;
  double c0;
  double Gamma0;
  double s;
  double e0;

  double eta_min; //!< a negative value that triggers special treatment (not "min eta"!)
  double p_min; //!< p_min = rho0*c0*c0*eta_min; (NOTE: This is NOT the "pmin" in the base class.)

  double eta_max; //!< 1/s - tol; (too much compression when reaching it. May get division-by-zero.)

  double rho0_c0_c0;    //!< rho0*c0*c0
  double Gamma0_over_2; //!< Gamma0/2
  double Gamma0_rho0;   //!< Gamma0*rho0
  double er_min;        //!< 0.5*c0*c0*eta_min*eta_min + e0

  enum TemperatureLaw {ORIGINAL_CV = 0, SIMPLIFIED_CV = 1, SIMPLIFIED_CP = 2} Tlaw;
  double T0; //!< ref. temperature

  double cv; //!< specific heat at constant volume
  double invcv;

  double cp;   //!< specific heat at constant pressure
  double invcp;
  double h0;   //!< ref. enthalpy corresponding to T0

  // -----------------------------
  //! Calculate Tr(eta) (eta>0) by uniform sampling & spline integration
  double delta_eta;
  std::vector<double> Trs;
  //boost::math::cubic_b_spline<double>* Tr_spline;  //!< Tried. Not better (& slower) than linear interp
  // -----------------------------

public:
  VarFcnMGExt(MaterialModelData &data);
  ~VarFcnMGExt() {
    //if(Tr_spline) delete Tr_spline;
  }

  //! ----- EOS-Specific Functions -----
  double GetPressure(double rho, double e) {
    double eta = 1.0 - rho0/rho;
    if(eta>=0.0) 
      return rho0_c0_c0*eta*(1.0 - Gamma0_over_2*eta)/((1.0-s*eta)*(1.0-s*eta)) + Gamma0_rho0*(e-e0);
    else if(eta>=eta_min)
      return rho0_c0_c0*eta*(1.0 - Gamma0_over_2*eta) + Gamma0_rho0*(e-e0);
    return rho0_c0_c0*eta_min*(1.0 - Gamma0_over_2*(2.0*eta-eta_min)) + Gamma0_rho0*(e-e0);
  }  

  inline double GetInternalEnergyPerUnitMass(double rho, double p) {
    double eta = 1.0 - rho0/rho;
    return (p - GetPr(eta))/Gamma0_rho0 + GetEr(eta); 
  }

  double GetDensity(double p, double e);

  inline double GetDpdrho(double rho, [[maybe_unused]] double e){
    double eta = 1.0 - rho0/rho;
    return (GetDPrDeta(eta) - Gamma0_rho0*GetDErDeta(eta))*GetDetaDrho(rho);
  }

  inline double GetBigGamma(double rho, [[maybe_unused]] double e) {return Gamma0_rho0/rho;}

  double GetTemperature(double rho, double e);

  inline double GetReferenceTemperature() {return T0;}

  inline double GetReferenceInternalEnergyPerUnitMass() {return e0;}

  double GetInternalEnergyPerUnitMassFromTemperature(double rho, double T);

  inline double GetInternalEnergyPerUnitMassFromEnthalpy(double rho, double h);

private:

  void SetupTrInterpolation();

  inline double GetPr(double eta) {
    if(eta>=0.0) {
      if(eta>=eta_max) {
        if(verbose>=1)
          fprintf(stdout,"\033[0;35mWarning: (Ext-MG EOS) Exceeding the compression limit "
                  "(eta=%e, eta_max=%e).\n",eta, eta_max);
        eta = eta_max;
      }
      double f = 1.0/(1.0-s*eta);
      return rho0_c0_c0*eta*f*f;
    }
    else if(eta>=eta_min) 
      return rho0_c0_c0*eta;

    return p_min;
  }

  inline double GetEr(double eta) {
    if(eta>=0.0) {
      if(eta>=eta_max) {
        if(verbose>=1)
          fprintf(stdout,"\033[0;35mWarning: (Ext-MG EOS) Exceeding the compression limit "
                  "(eta=%e, eta_max=%e).\n",eta, eta_max);
        eta = eta_max;
      }
      double f = 1.0/(1.0-s*eta);
      return 0.5*c0*c0*eta*eta*f*f + e0;
    }
    else if(eta>=eta_min)
      return 0.5*c0*c0*eta*eta + e0;

    return er_min + p_min/rho0*(eta-eta_min);
  }

  inline double GetTr(double eta) {
    if(eta>=0.0) {
      assert(eta<1.0);
      if(eta>=eta_max) {
        if(verbose>=1)
          fprintf(stdout,"\033[0;35mWarning: (Ext-MG EOS) Exceeding the compression limit "
                  "(eta=%e, eta_max=%e).\n",eta, eta_max);
        eta = eta_max;
      }

      double division = eta/delta_eta; 
      int N = int(division); //Note: cast-to-int is different from "floor" for negative numbers. Not relevant
                             //      here though. cast-to-int is much faster than "floor"!
      assert(N+1<(int)Trs.size());
      double remainder = division - (double)N;
      return Trs[N]*(1.0-remainder)+Trs[N+1]*remainder;

      //return (*Tr_spline)(eta);
    }
    return T0*exp(Gamma0*eta);
  }

  inline double GetDetaDrho(double rho) {return rho0/(rho*rho);}

  inline double GetDPrDeta(double eta) {
    if(eta>=0.0) {
      if(eta>=eta_max) {
        if(verbose>=1)
          fprintf(stdout,"\033[0;35mWarning: (Ext-MG EOS) Exceeding the compression limit "
                  "(eta=%e, eta_max=%e).\n",eta, eta_max);
        eta = eta_max;
      }
      double f = 1.0/(1.0-s*eta);
      return rho0_c0_c0*(1+s*eta)*f*f*f;
    }
    else if(eta>=eta_min)
      return rho0_c0_c0;

    return 0.0;
  }

  inline double GetDErDeta(double eta) {
    if(eta>=0.0) {
      if(eta>=eta_max) {
        if(verbose>=1)
          fprintf(stdout,"\033[0;35mWarning: (Ext-MG EOS) Exceeding the compression limit "
                  "(eta=%e, eta_max=%e).\n",eta, eta_max);
        eta = eta_max;
      }
      double f = 1.0/(1.0-s*eta);
      return c0*c0*eta*f*f*f;
    }
    else if(eta>=eta_min)
      return c0*c0*eta;

    return p_min/rho0;
  }

};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

VarFcnMGExt::VarFcnMGExt(MaterialModelData &data) : VarFcnBase(data)
{

  if(data.eos != MaterialModelData::EXTENDED_MIE_GRUNEISEN){
    fprintf(stdout, "*** Error: MaterialModelData is not of type Extended Mie-Gruneisen.\n");
    exit(-1);
  }

  type   = EXTENDED_MIE_GRUNEISEN;

  rho0   = data.mgextModel.rho0;
  c0     = data.mgextModel.c0;
  Gamma0 = data.mgextModel.Gamma0;
  s      = data.mgextModel.s;
  e0     = data.mgextModel.e0;

  eta_min= data.mgextModel.eta_min;
  p_min  = rho0*c0*c0*eta_min;

  eta_max = 1.0/s - 1e-6;

  cv = data.mgextModel.cv;
  invcv = cv==0.0 ? 0.0 : 1.0/cv;
  T0 = data.mgextModel.T0;

  cp = data.mgextModel.cp;
  invcp = cp==0.0 ? 0.0 : 1.0/cp;
  h0 = data.mgextModel.h0;

  Tlaw = ORIGINAL_CV;
  if(data.mgextModel.Tlaw == ExtendedMieGruneisenModelData::SIMPLIFIED_CV)
    Tlaw = SIMPLIFIED_CV;
  else if(data.mgextModel.Tlaw == ExtendedMieGruneisenModelData::SIMPLIFIED_CP)
    Tlaw = SIMPLIFIED_CP;


  // define some frequently used constants
  rho0_c0_c0 = rho0*c0*c0;
  Gamma0_over_2 = 0.5*Gamma0;
  Gamma0_rho0 = Gamma0*rho0;

  er_min = 0.5*c0*c0*eta_min*eta_min + e0;


  if(rho0<=0.0 || c0<=0.0 || Gamma0<=0.0 || s<=0.0) {
    fprintf(stdout, "*** Error: VarFcnMGExt detected non-positive rho0 (%e), c0 (%e), "
                    "Gamma0 (%e), or s (%e).\n", rho0, c0, Gamma0, s);
    exit(-1);
  }

  if(eta_min>=0.0) {
    fprintf(stdout, "*** Error: VarFcnMGExt detected positive or zero eta_min (%e).\n",
                    eta_min);
    exit(-1);
  }


  // build the splines for interpolation
  SetupTrInterpolation();
   
}

//------------------------------------------------------------------------------

void VarFcnMGExt::SetupTrInterpolation()
{
  //Tr_spline = NULL;
  Trs.clear();

  double eta0 = 0.0;
  [[maybe_unused]] double sol0 = 0.0; 
  double eta1 = std::min(eta_max, 1.0);
  [[maybe_unused]] double sol1;
  int N = 500001;
  delta_eta = (eta1 - eta0)/(N-1.0);  //small enough?

  // sets up the integration
  auto fun = [&]([[maybe_unused]] double *Int, double z, double *Integrand) {
    double c = 1.0/(1.0 - s*z);
    *Integrand = exp(Gamma0*z)*z*z*c*c*c; return;};

  int err = MathTools::runge_kutta_4(fun, 1, eta0, &sol0, delta_eta, eta1, &sol1,
                                     [](double*){return false;}, &Trs);
  assert(!err);

  if((int)Trs.size()==N+1) //the last one must be created due to round-off. we just drop it.
    Trs.pop_back();
  assert((int)Trs.size()==N);

  double coeff = c0*c0*s*invcv;
  double eta;
  for(int i=0; i<N; i++) {
    eta =  eta0 + i*delta_eta;
    Trs[i] = exp(Gamma0*eta)*(T0 + coeff*Trs[i]);
  }

  //Tr_spline = new boost::math::cubic_b_spline<double>(Trs.begin(), Trs.end(), eta0, delta_eta);
}

//------------------------------------------------------------------------------

double VarFcnMGExt::GetDensity(double p, double e)
{

  // ---------------------------------------------------------------------------
  // (dp/drho)|(e=const) is always positive. Therefore, we first evaluate p(rho0,e). 
  // If it is lower than p, then, rho>rho0 (i.e., eta>0). Otherwise, we check the
  // other two formulas.
  // -----------------------------
  double ptrial = GetPressure(rho0,e);
  if(ptrial==p) //so lucky...
    return rho0;

  if(ptrial<p) {
    // -----------------------------
    // Check for eta>=0 
    // -----------------------------
    //solving a quadratic equation a*eta^2 + b*eta + c = 0 for eta ==> rho = rho0/(1-eta)
    double c = p - Gamma0_rho0*(e - e0);
    double a = c*s*s + Gamma0_over_2*rho0_c0_c0;
    double b = -(2.0*c*s + rho0_c0_c0);

    if(a==0) //linear equation: eta = -c/b, rho = rho0/(1-eta)
      return rho0/(1 + c/b);

    double b2m4ac = b*b - 4.0*a*c;
    if(b2m4ac<0) {
      fprintf(stdout, "*** Error: The Extended M-G EOS is invalid for the given p(%e) and e(%e) "
              "--- unable to solve it for rho.\n", p, e);
      exit(-1);
    }

    b2m4ac = sqrt(b2m4ac);
    double eta1 = (-b + b2m4ac)/(2.0*a);
    double eta2 = (-b - b2m4ac)/(2.0*a);
    if(eta1>eta2)
      std::swap(eta1,eta2); //eta2 should be the bigger one

    if(eta1>=1.0 || eta2<0) { //both are negative, or both are positive
      fprintf(stdout, "*** Error: Cannot determine the solution (rho) of the Extended M-G EOS "
              "(eta1 = %e, eta2 = %e). \n", eta1, eta2);
      exit(-1);
    }

    if(eta2>=1.0) { //eta1 should be the valid root.
      if(eta1<0.0) {
        fprintf(stdout, "*** Error: Cannot find the solution (rho) of the Extended M-G EOS "
                "(eta1 = %e, eta2 = %e). \n", eta1, eta2);
        exit(-1);
      }   
      return rho0/(1.0 - eta1);
    }

    // now, eta2 must be between 0 and 1 (i.e., valid)
    
    if(eta1>=0.0 && eta1 != eta2) { //both eta1 and eta2 can be valid...
      fprintf(stdout, "*** Error: Cannot find the solution (rho) of the Extended M-G EOS "
              "(eta1 = %e, eta2 = %e). \n", eta1, eta2);
      exit(-1);
    }

    return rho0/(1.0 - eta2);
  }

  //If we get here, it means ptrial>p ==> eta<0

  // --------------------------------------------------------
  // Check for eta<eta_min (because it is a linear equation)
  // --------------------------------------------------------
  double c = (p - Gamma0_rho0*(e-e0))/rho0_c0_c0;
  double eta = (1.0-c/eta_min)/Gamma0 + 0.5*eta_min;
  if(eta<=eta_min)
    return rho0/(1.0-eta);


  //If we get here, it means eta_min < eta < 0

  //solving a quadratic equation a*eta^2 + b*eta + c = 0 for eta ==> rho = rho0/(1-eta)
  double b2m4ac = 1.0 - 2.0*Gamma0*c;
  if(b2m4ac<0) {
    fprintf(stdout, "*** Error: The Extended M-G EOS is invalid for the given p(%e) and e(%e) "
            "--- unable to solve it for rho (eta_min<eta<0).\n", p, e);
    exit(-1);
  }
  b2m4ac = sqrt(b2m4ac);
  eta = (1.0 - sqrt(b2m4ac))/Gamma0; //eta has to be the one with (-). The other root is positive

  if(eta>0 || eta<eta_min) {
    fprintf(stdout, "*** Error: Cannot determine the solution (rho) of the Extended M-G EOS "
            "(eta = %e). \n", eta);
    exit(-1);
  }

  return rho0/(1.0 - eta);
 
}

//------------------------------------------------------------------------------

double VarFcnMGExt::GetTemperature(double rho, double e)
{
  if(Tlaw == SIMPLIFIED_CV)
    return T0 + invcv*(e-e0);

  if(Tlaw == ORIGINAL_CV) {
    double eta = 1.0 - rho0/rho;
    return GetTr(eta) + invcv*(e - GetEr(eta));
  }

  // SimplifiedCp
  double p = GetPressure(rho, e);
  return T0 + invcp*(e + p/rho - h0);
}

//------------------------------------------------------------------------------

double VarFcnMGExt::GetInternalEnergyPerUnitMassFromTemperature(double rho, double T) 
{
  if(Tlaw == SIMPLIFIED_CV)
    return e0 + cv*(T-T0);

  if(Tlaw == ORIGINAL_CV) {
    double eta = 1.0 - rho0/rho;
    return GetEr(eta) + cv*(T - GetTr(eta));
  }

  // SimplifiedCp
  double eta = 1.0 - rho0/rho;
  return (h0 + cp*(T-T0) + (Gamma0_rho0*GetEr(eta) - GetPr(eta))/rho)/(1.0 + Gamma0_rho0/rho);
}

//------------------------------------------------------------------------------

inline
double VarFcnMGExt::GetInternalEnergyPerUnitMassFromEnthalpy(double rho, double h) 
{
  double eta = 1.0 - rho0/rho;
  return (h + (Gamma0_rho0*GetEr(eta) - GetPr(eta))/rho)/(1.0 + Gamma0_rho0/rho);
}

//------------------------------------------------------------------------------

#endif
