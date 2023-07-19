/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _VAR_FCN_MG_EXT_H
#define _VAR_FCN_MG_EXT_H

#include <VarFcnBase.h>
#include <fstream>
#include <cassert>

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
  double p_min; //!< p_min = rho0*c0*c0*eta_min;

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

public:
  VarFcnMGExt(MaterialModelData &data);
  ~VarFcnMGExt() {}

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

  inline double GetTemperature(double rho, double e);

  inline double GetReferenceTemperature() {return T0;}

  inline double GetReferenceInternalEnergyPerUnitMass() {return e0;}

  inline double GetInternalEnergyPerUnitMassFromTemperature(double rho, double T);

  inline double GetInternalEnergyPerUnitMassFromEnthalpy(double rho, double h);

protected:

  inline double GetPr(double eta) {
    if(eta>=0.0) {
      double f = 1.0/(1.0-s*eta);  f = f*f;  assert(isfinite(f));
      return rho0_c0_c0*eta*f;
    }
    else if(eta>=eta_min) 
      return rh0_c0_c0*eta;

    return p_min;
  }

  inline double GetEr(double eta) {
    if(eta>=0.0) {
      double f = 1.0/(1.0-s*eta);  f = f*f;  assert(isfinite(f));
      return 0.5*c0*c0*eta*eta*f + e0;
    }
    else if(eta>=eta_min)
      return 0.5*c0*c0*eta*eta + e0;

    return er_min + pmin/rho0*(eta-eta_min);
  }

  inline double GetTr(double eta) {
    if(eta>=0.0) {
      assert(eta<1.0);
      return (*Tr_spline)(eta);
    }
    return T0*exp(Gamma0*eta);
  }

  inline double GetDetaDrho(double rho) {return rho0/(rho*rho);}

  inline double GetDPrDeta(double eta) {
    if(eta>=0.0) {
      double f = 1.0/(1.0-s*eta);  f = f*f*f;  assert(isfinite(f));
      return rho0_c0_c0*(1+s*eta)*f;
    }
    else if(eta>=eta_min)
      return rh0_c0_c0;

    return 0.0;
  }

  inline double GetDErDeta(double eta) {
    if(eta>=0.0) {
      double f = 1.0/(1.0-s*eta);  f = f*f*f;  assert(isfinite(f));
      return c0*c0*eta*f;
    }
    else if(eta>=eta_min)
      return c0*c0*eta;

    return p_min/rho0;
  }

};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

VarFcnMGExt::VarFcnMGExt(MaterialModelData &data) : VarFcnBase(data) {

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
}

//------------------------------------------------------------------------------

double VarFcnMGExt::GetDensity(double p, double e) {

  // -----------------------------
  // Step 1: Check for eta>=0 
  // -----------------------------
  //solving a quadratic equation a*eta^2 + b*eta + c = 0 for eta ==> rho = rho0/(1-eta)
  double c = p - Gamma0_rho0*(e - e0);
  double a = c*s*s + Gamma0_over_2*rho0_c0_c0;
  double b = -(2.0*c*s + rho0_c0_c0);

  rho = 0.5*rho0; //give it a tensile rho!
  if(a==0) {//linear equation: eta = -c/b, rho = rho0/(1-eta)
    rho = rho0/(1 + c/b);
  } else {
    double b2m4ac = b*b - 4.0*a*c;
    if(b2m4ac<0) 
      goto FAILED_STEP1;

    b2m4ac = sqrt(b2m4ac);
    double rho1 = rho0/(1.0 - (-b + b2m4ac)/(2.0*a));
    double rho2 = rho0/(1.0 - (-b - b2m4ac)/(2.0*a));

    if(rho1>rho2)
      std::swap(rho1,rho2); //rho2 should be the bigger one

    if(rho1>0) { //both are are positive
      fprintf(stdout, "*** Error: Cannot determine the solution (rho) of the M-G EOS "
              "(rho1 = %e, rho2 = %e). \n", rho1, rho2);
      exit(-1);
    }

    rho = rho2; 
  }

  if(rho>=rho0) //good!
    return rho;

FAILED_STEP1:

  // -----------------------------
  //Step 2: Check for eta>=eta_min
  // -----------------------------

  I AM HERE
}

//------------------------------------------------------------------------------

inline 
double VarFcnMGExt::GetTemperature(double rho, double e) {

  if(use_cp) {
    double p = GetPressure(rho, e);
    return T0 + invcp*(e + p/rho - h0);
  } else
    return T0 + invcv*(e-e0);
}

//------------------------------------------------------------------------------

inline
double VarFcnMGExt::GetInternalEnergyPerUnitMassFromTemperature(double rho, double T) 
{
  if(use_cp) {
    double eta = 1.0 - rho0/rho;
    double first_term = rho0_c0_c0*eta*(1.0 - Gamma0_over_2*eta)/((1.0-s*eta)*(1.0-s*eta));
    return (h0 + cp*(T-T0) + (Gamma0_rho0*e0 - first_term)/rho)/(Gamma0_rho0/rho + 1.0);
  } else
    return e0 + cv*(T-T0);
}

//------------------------------------------------------------------------------

inline
double VarFcnMGExt::GetInternalEnergyPerUnitMassFromEnthalpy(double rho, double h) 
{
  double eta = 1.0 - rho0/rho;
  double first_term = rho0_c0_c0*eta*(1.0 - Gamma0_over_2*eta)/((1.0-s*eta)*(1.0-s*eta));
  return (h + (Gamma0_rho0*e0 - first_term)/rho)/(Gamma0_rho0/rho + 1.0);
}

//------------------------------------------------------------------------------

#endif
