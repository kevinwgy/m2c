/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _VAR_FCN_MG_H
#define _VAR_FCN_MG_H

#include <VarFcnBase.h>
#include <fstream>

/********************************************************************************
 * This class is the VarFcn class for the Mie-Gruneisen equation of state (EOS). It
 * applies the same equation --- see below --- for compression and tension. It only
 * supports simplified linear temperature laws (cv or cp specified).
 *
 * NOTE: In general, the user should try to apply the extended version, implemented in
 *       VarFcnMGExt.h
 *
 * Only elementary functions are declared and/or defined here.
 * All arguments must be pertinent to only a single grid node or a single
 * state.
 *
 * EOS: Pressure = rho0*c0*c0*eta*(1-Gamma0/2*eta)/(1-s*eta)^2 + rho0*Gamma0*(e-e0)
 *
 * where
 *
 *   e      : internal energy per unit mass.
 *   eta    : 1 - rho0/rho
 *   rho    : density
 *
 *   rho0   : ref. density
 *   c0     : bulk speed of sound
 *   Gamma0 : Gruneisen coefficient at ref. state
 *   s      : slope of shock Hugoniot
 *   e0     : internal energy at ref. state (including temperature T0)
 *
 * References: Shafquat Islam's report (01/2021), Allen Robinson's technical report (2019)
 *
 *   Note: default temperature law is de = cv*dT or dh = cp*dT. 
 ********************************************************************************/

class VarFcnMG : public VarFcnBase {

private:
  double rho0;
  double c0;
  double Gamma0;
  double s;
  double e0;

  double rho0_c0_c0;    //!< rho0*c0*c0
  double Gamma0_over_2; //!< Gamma0/2
  double Gamma0_rho0;   //!< Gamma0*rho0

  double cv; //!< specific heat at constant volume
  double invcv;
  double T0; //!< ref. temperature

  bool use_cp; //!< whether cp (instead of cv) is used in the temperature law
  double cp;   //!< specific heat at constant pressure
  double invcp;
  double h0;   //!< ref. enthalpy corresponding to T0

public:
  VarFcnMG(MaterialModelData &data);
  ~VarFcnMG() {}

  //! ----- EOS-Specific Functions -----
  inline double GetPressure(double rho, double e) {
    double eta = 1.0 - rho0/rho;
    return rho0_c0_c0*eta*(1.0 - Gamma0_over_2*eta)/((1.0-s*eta)*(1.0-s*eta)) + Gamma0_rho0*(e-e0);}

  inline double GetInternalEnergyPerUnitMass(double rho, double p) {
    double eta = 1.0 - rho0/rho;
    return (p - rho0_c0_c0*eta*(1.0 - Gamma0_over_2*eta)/((1.0-s*eta)*(1.0-s*eta)))/Gamma0_rho0 + e0;
  }

  double GetDensity(double p, double e);

  inline double GetDpdrho(double rho, [[maybe_unused]] double e){
    double eta = 1.0 - rho0/rho;
    double S = 1.0 - s*eta;
    return rho0_c0_c0*(1.0 + (s - Gamma0)*eta)/(S*S*S)*rho0/(rho*rho);
  }

  inline double GetBigGamma(double rho, [[maybe_unused]] double e) {return Gamma0_rho0/rho;}

  double GetTemperature(double rho, double e);

  inline double GetReferenceTemperature() {return T0;}

  inline double GetReferenceInternalEnergyPerUnitMass() {return e0;}

  double GetInternalEnergyPerUnitMassFromTemperature(double rho, double T);

  double GetInternalEnergyPerUnitMassFromEnthalpy(double rho, double h);

};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

VarFcnMG::VarFcnMG(MaterialModelData &data) : VarFcnBase(data) {

  if(data.eos != MaterialModelData::MIE_GRUNEISEN){
    fprintf(stdout, "*** Error: MaterialModelData is not of type Mie-Gruneisen.\n");
    exit(-1);
  }

  type   = MIE_GRUNEISEN;

  rho0   = data.mgModel.rho0;
  c0     = data.mgModel.c0;
  Gamma0 = data.mgModel.Gamma0;
  s      = data.mgModel.s;
  e0     = data.mgModel.e0;

  cv = data.mgModel.cv;
  invcv = cv==0.0 ? 0.0 : 1.0/cv;
  T0 = data.mgModel.T0;

  cp = data.mgModel.cp;
  invcp = cp==0.0 ? 0.0 : 1.0/cp;
  h0 = data.mgModel.h0;

  use_cp = (cp>0 && cv<=0.0) ? true : false;

  rho0_c0_c0 = rho0*c0*c0;
  Gamma0_over_2 = 0.5*Gamma0;
  Gamma0_rho0 = Gamma0*rho0;

  if(rho0<=0.0 || c0<=0.0 || Gamma0<=0.0 || s<=0.0) {
    fprintf(stdout, "*** Error: VarFcnMG detected non-positive rho0 (%e), c0 (%e), "
                    "Gamma0 (%e), or s (%e).\n", rho0, c0, Gamma0, s);
    exit(-1);
  }
}

//------------------------------------------------------------------------------

double VarFcnMG::GetDensity(double p, double e) {

  //solving a quadratic equation a*eta^2 + b*eta + c = 0 for eta ==> rho = rho0/(1-eta)
  double c = p - Gamma0_rho0*(e - e0);
  double a = c*s*s + Gamma0_over_2*rho0_c0_c0;
  double b = -(2.0*c*s + rho0_c0_c0);

  if(a==0) {//linear equation: eta = -c/b, rho = rho0/(1-eta)

    return rho0/(1 + c/b);

  } else {

    double b2m4ac = b*b - 4.0*a*c;
    if(b2m4ac<0) {
      fprintf(stdout, "*** Error: The M-G EOS is invalid for the given p(%e) and e(%e) --- unable to solve it for rho.\n",
              p, e);
      exit(-1);
    }

    b2m4ac = sqrt(b2m4ac);
    double rho1 = rho0/(1.0 - (-b + b2m4ac)/(2.0*a));
    double rho2 = rho0/(1.0 - (-b - b2m4ac)/(2.0*a));

    if(rho1>rho2)
      std::swap(rho1,rho2); //rho2 should be the bigger one

    if(rho1>0 || rho2<0) { //both are negative, or both are positive
      fprintf(stdout, "*** Error: Cannot determine the solution (rho) of the M-G EOS (rho1 = %e, rho2 = %e). \n", 
              rho1, rho2);
      exit(-1);
    }
    
    return rho2; 

  }
}

//------------------------------------------------------------------------------

double VarFcnMG::GetTemperature(double rho, double e) {

  if(use_cp) {
    double p = GetPressure(rho, e);
    return T0 + invcp*(e + p/rho - h0);
  } else
    return T0 + invcv*(e-e0);
}

//------------------------------------------------------------------------------

double VarFcnMG::GetInternalEnergyPerUnitMassFromTemperature(double rho, double T) 
{
  if(use_cp) {
    double eta = 1.0 - rho0/rho;
    double first_term = rho0_c0_c0*eta*(1.0 - Gamma0_over_2*eta)/((1.0-s*eta)*(1.0-s*eta));
    return (h0 + cp*(T-T0) + (Gamma0_rho0*e0 - first_term)/rho)/(Gamma0_rho0/rho + 1.0);
  } else
    return e0 + cv*(T-T0);
}

//------------------------------------------------------------------------------

double VarFcnMG::GetInternalEnergyPerUnitMassFromEnthalpy(double rho, double h) 
{
  double eta = 1.0 - rho0/rho;
  double first_term = rho0_c0_c0*eta*(1.0 - Gamma0_over_2*eta)/((1.0-s*eta)*(1.0-s*eta));
  return (h + (Gamma0_rho0*e0 - first_term)/rho)/(Gamma0_rho0/rho + 1.0);
}

//------------------------------------------------------------------------------

#endif
