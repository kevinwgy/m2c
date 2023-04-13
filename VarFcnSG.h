/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _VAR_FCN_SG_H
#define _VAR_FCN_SG_H

#include <VarFcnBase.h>
#include <fstream>

/********************************************************************************
 * This class is the VarFcn class for the Stiffened Gas equation of state (EOS)
 * Only elementary functions are declared and/or defined here.
 * All arguments must be pertinent to only a single grid node or a single
 * state.
 *
 * EOS: Pressure = (gam - 1)*Density*e - gam*Pc
 * where
 *   e  : internal energy per unit mass.
 *   Pc : pressure constant. (named Pstiff below)
 *
 * Note: Three different temperature laws have been implemented. 
 *       Method 1 assumes de = cv*dT, and T is a function of only e. 
 *       Method 2 assumes dh = cp*dT, and T is a function of only h.
 *       Method 3 assumes de = cv*dT, but T is a function of rho and e. (See K.W. notes)
 *
 *   For a perfect gas, Method 1 leads to dh = cp*dT, where cp = gamma*cv is the
 *   specific heat at constant pressure. For a general stiffened gas,
 *   dh =/= cp*dT! See KW's notes. (One could have assumed dh = cp*dT, but then
 *   de =/= cv*dT.)
 *
 *   If      cv>0  && rho0>0 ==> Method 3
 *   Else if cv<=0 && cp>0   ==> Method 2
 *   Else                    ==> Method 1
 ********************************************************************************/
class VarFcnSG : public VarFcnBase {

private:
  double gam;
  double Pstiff;

  double invgam;  //!< 1/gamma
  double gam1;    //!< gamma-1
  double invgam1; //!< 1/(gamma-1)

  double cv; //!< specific heat at constant volume
  double invcv;
  double T0; //!< ref. temperature
  double e0; //!< ref. internal energy corresponding to T0
  double rho0; //!< ref. density (for T = T(rho,e))
  bool use_cv_advanced;

  bool use_cp; //!< whether cp (instead of cv) is used in the temperature law
  double cp;   //!< specific heat at constant pressure
  double invcp;
  double h0;   //!< ref. enthalpy corresponding to T0

public:
  VarFcnSG(MaterialModelData &data);
  ~VarFcnSG() {}

  //! ----- EOS-Specific Functions -----
  inline double GetPressure(double rho, double e) {return gam1*rho*e - gam*Pstiff;}
  inline double GetInternalEnergyPerUnitMass(double rho, double p) {return (p+gam*Pstiff)/(gam1*rho);}
  inline double GetDensity(double p, double e) {return (p+gam*Pstiff)/(gam1*e);}
  inline double GetDpdrho([[maybe_unused]] double rho, double e) {return gam1*e;}
  inline double GetBigGamma([[maybe_unused]] double rho, [[maybe_unused]] double e) {return gam1;}

  inline double GetTemperature(double rho, double e) {
    if(use_cv_advanced) { //Method 3
      return invcv*(e + Pstiff/rho) + pow(rho/rho0, gam1)*(T0 - invcv*(e0 + Pstiff/rho0));
    } else if(use_cp) { //Method 2
      double p = GetPressure(rho, e);
      return T0 + invcp*(e + p/rho - h0);
    } else //Method 1
      return T0 + invcv*(e-e0);
  }

  inline double GetReferenceTemperature() {return T0;}
  inline double GetReferenceInternalEnergyPerUnitMass() {return e0;}

  inline double GetInternalEnergyPerUnitMassFromTemperature(double rho, double T) {
    if(use_cv_advanced) {
      return cv*T - Pstiff/rho - pow(rho/rho0, gam1)*(cv*T0 - (e0 + Pstiff/rho0));
    } else if(use_cp) 
      return invgam*(h0 + cp*(T-T0)) + Pstiff/rho;
    else
      return e0 + cv*(T-T0);
  }
  
  inline double GetInternalEnergyPerUnitMassFromEnthalpy(double rho, double h) {return invgam*h+Pstiff/rho;}


  //! Verify hyperbolicity (i.e. c^2 > 0): Report error if rho < 0 or p + Pstiff < 0 (Not p + gamma*Pstiff). 
  inline bool CheckState(double rho, double p, bool silence = false) {
    if(!std::isfinite(rho) || !std::isfinite(p)) {
      if(!silence)
        fprintf(stdout, "*** Error: CheckState failed. rho = %e, p = %e.\n", rho, p);
      return true;
    }
    if(rho <= 0.0 || p+Pstiff <= 0.0){
      if(!silence && verbose>1)
        fprintf(stdout, "Warning: Negative density or violation of hyperbolicity. rho = %e, p = %e.\n", rho, p);
      return true;
    }
    return false;
  }

};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

VarFcnSG::VarFcnSG(MaterialModelData &data) : VarFcnBase(data) {

  if(data.eos != MaterialModelData::STIFFENED_GAS){
    fprintf(stdout, "*** Error: MaterialModelData is not of type STIFFENED_GAS.\n");
    exit(-1);
  }

  type = STIFFENED_GAS;

  gam = data.sgModel.specificHeatRatio;
  invgam = (gam==0.0) ? 0.0 : 1.0/gam;
  gam1 = gam -1.0;
  invgam1 = 1.0/gam1;
  Pstiff = data.sgModel.pressureConstant;

  cv = data.sgModel.cv;
  invcv = cv==0.0 ? 0.0 : 1.0/cv;
  T0 = data.sgModel.T0;
  e0 = data.sgModel.e0;

  cp = data.sgModel.cp;
  invcp = cp==0.0 ? 0.0 : 1.0/cp;
  h0 = data.sgModel.h0;

  rho0 = data.sgModel.rho0;

  use_cv_advanced = (cv>0 && rho0>0) ? true : false;
  use_cp = (cp>0 && cv<=0.0) ? true : false;
    
}

//------------------------------------------------------------------------------

#endif
