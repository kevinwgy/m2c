/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _VAR_FCN_VanDerWaals_H
#define _VAR_FCN_VanDerWaals_H

#include <VarFcnBase.h>
#include <fstream>

/********************************************************************************
 * This class is the VarFcn class for the van der Waals EOS.
 * Only elementary functions are declared and/or defined here.
 * All arguments must be pertinent to only a single grid node or a single
 * state.
 *
 * Note: This implementation assumes a constant cv. For vdW, this means cp is not a 
 *       constant. The specific heat ratio, gamma, is also not a constant. (By first law of
 *       thermodynamics.)
 *
 * EOS: Pressure = rho*Rs/(1-b*rho)*(T0 + (e+a*(rho-rho0))/cv) - a*rho*rho
 *      Temperature = T0 + (e + a*(rho-rho0))/cv
 *      (where Rs is the specific gas constant)
 *
 * Parameters:
 *      - User inputs (specified in simulation input file) are based on the common
 *        expression of the EOS: p = RT/(V-bb) - aa/V^2, where R is the molar gas constant,
 *        V is molar volume.
 *        o R: molar gass constant [energy]/([Kelvin].[mol])
 *        o M: material's molar mass [mass]/[mol]
 *        o aa: molecular attraction [energy].[length]^3/[mol]^2
 *        o bb: min. molar volume [length]^3/[mol]
 *        o (rho0, T0): reference density [mass]/[length]^3 and temperature [Kelvin]
 *          (reference pressure p0 is determined by rho0 and T0)
 *        o cv: specific heat capacity at constant volume [energy]/([mass].[Kelvin])
 *
 *      - IMPORTANT: The functions in this file do not directly operate on R, aa, and bb.
 *                   They are converts from per mol to per mass in the constructor.
 *      - Inside the code:
 *        o Rs = R/M
 *        o a = aa/M
 *        o b = bb/M
 ********************************************************************************/
class VarFcnVanDerWaals : public VarFcnBase {

private:

  double Rs;
  double a, b;
  double rho0; //!< ref. density 
  double T0; //!< ref. temperature
  double cv; //!< specific heat at constant volume
  double invcv;

public:
  VarFcnVanDerWaals(MaterialModelData &data);
  ~VarFcnVanDerWaals() {}

  //! ----- EOS-Specific Functions -----
  inline double GetPressure(double rho, double e) {
    return rho*Rs/(1.0-b*rho)*(T0 + (e + a*(rho-rho0))/cv) - a*rho*rho;}

  inline double GetInternalEnergyPerUnitMass(double rho, double p) {
    return cv*(1.0-b*rho)*(p+a*rho*rho)/(rho*Rs) - cv*T0 - a*(rho-rho0);

  inline double GetDensity(double p, double e) {return (p+gam*Pstiff)/(gam1*e);}

  double GetDpdrho(double rho, double e);

  inline double GetBigGamma(double rho, [[maybe_unused]] double e) {return Rs/(cv*(1.0-b*rho));}

  inline double GetTemperature(double rho, double e) {return T0 + (e+a*(rho-rho0))/cv;}

  inline double GetReferenceTemperature() {return T0;}
  inline double GetReferenceInternalEnergyPerUnitMass() {return 0.0;}

  inline double GetInternalEnergyPerUnitMassFromTemperature(double rho, double T) {
    return cv*(T-T0) - a*(rho-rho0);
  }
  
  double GetInternalEnergyPerUnitMassFromEnthalpy(double rho, double h);

  //! Verify hyperbolicity (i.e. c^2 > 0) , rho>0 and 1/rho>b
  inline bool CheckState(double rho, double p, bool silence = false) {
    VarFcnBase::CheckState(rho, p, silence);
    if(1.0/rho<=b) {
      if(!silence)
        fprintf(stdout, "*** Error: CheckState failed (vdW). rho = %e, p = %e, b = %e.\n", 
                rho, p, b);
      return true;
    }
    return false;
  }

};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

VarFcnVanDerWaals::VarFcnVanDerWaals(MaterialModelData &data) : VarFcnBase(data) {

  if(data.eos != MaterialModelData::VAN_DER_WAALS){
    fprintf(stdout, "*** Error: MaterialModelData is not of type VAN_DER_WAALS.\n");
    exit(-1);
  }

  type = VAN_DER_WAALS;

  double molar_mass = data.vdwModel.M;
  Rs = data.vdwModel.R/molar_mass;
  a  = data.vdwModel.aa/molar_mass;
  b  = data.vdwModel.bb/molar_mass;
  
  T0 = data.vdwModel.T0;
  rho0 = data.vdwModel.rho0;

  cv = data.vdwModel.cv;
  invcv = cv==0.0 ? 0.0 : 1.0/cv;
}

//------------------------------------------------------------------------------

double
VarFcnVanDerWaals::GetDpdrho(double rho, double e)
{
  double beta = 1.0 - b*rho;
  double T = GetTemperature(rho, e);
  return Rs*T/(beta*beta) + a*Rs*rho/(cv*beta) - 2*a*rho;
}

//------------------------------------------------------------------------------

double
VarFcnVanDerWaals::GetInternalEnergyPerUnitMassFromEnthalpy(double rho, double h)
{
  double T = (h + cv*T0 + 2.0*a*rho - a*rho0)/(cv + Rs/(1.0-b*rho));
  return cv*(T-T0) - a*(rho-rho0);
}

//------------------------------------------------------------------------------


#endif
