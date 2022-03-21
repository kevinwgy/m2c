#ifndef _VAR_FCN_SG_H
#define _VAR_FCN_SG_H

#include <VarFcnBase.h>
#include <fstream>

/********************************************************************************
 * This class is the VarFcn class for the Stiffened Gas EOS in Euler
 * Equations. Only elementary functions are declared and/or defined here.
 * All arguments must be pertinent to only a single grid node or a single
 * state.
 *
 * EOS: Pressure = (gam - 1)*Density*e - gam*Pc
 * where
 *   e  : internal energy per unit mass.
 *   Pc : pressure constant. (named Pstiff below)
 *
 *   The default temperature law is de = cv*dT, where cv is assumed to be a constant
 *   For a perfect gas, this leads to dh = cp*dT, where cp = gamma*cv is the
 *   specific heat at constant pressure. For a general stiffened gas,
 *   dh =/= cp*dT! See KW's notes. (One could have assumed dh = cp*dT, but then
 *   de =/= cv*dT.)
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

  bool use_cp; //!< whether cp (instead of cv) is used in the temperature law
  double cp;   //!< specific heat at constant pressure
  double invcp;
  double h0;   //!< ref. enthalpy corresponding to T0

public:
  VarFcnSG(MaterialModelData &data);
  ~VarFcnSG() {}

  //! ----- EOS-Specific Functions -----
  inline double GetPressure(double rho, double e) const {return gam1*rho*e - gam*Pstiff;}
  inline double GetInternalEnergyPerUnitMass(double rho, double p) const {return (p+gam*Pstiff)/(gam1*rho);}
  inline double GetDensity(double p, double e) const {return (p+gam*Pstiff)/(gam1*e);}
  inline double GetDpdrho(double rho, double e) const{return gam1*e;}
  inline double GetBigGamma(double rho, double e) const {return gam1;}

  inline double GetTemperature(double rho, double e) const {
    if(use_cp) {
      double p = GetPressure(rho, e);
      return T0 + invcp*(e + p/rho - h0);
    } else
      return T0 + invcv*(e-e0);
  }

  inline double GetReferenceTemperature() const {return T0;}

  inline double GetInternalEnergyPerUnitMassFromTemperature(double rho, double T) const {
    if(use_cp) 
      return invgam*(h0 + cp*(T-T0)) + Pstiff/rho;
    else
      return e0 + cv*(T-T0);
  }
  
  inline double GetInternalEnergyPerUnitMassFromEnthalpy(double rho, double h) const {return invgam*h+Pstiff/rho;}


  //! Verify hyperbolicity (i.e. c^2 > 0): Report error if rho < 0 or p + Pstiff < 0 (Not p + gamma*Pstiff). 
  inline bool CheckState(double rho, double p, bool silence = false) const{
    if(!isfinite(rho) || !isfinite(p)) {
      if(!silence)
        fprintf(stderr, "*** Error: CheckState failed. rho = %e, p = %e.\n", rho, p);
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
inline
VarFcnSG::VarFcnSG(MaterialModelData &data) : VarFcnBase(data) {

  if(data.eos != MaterialModelData::STIFFENED_GAS){
    fprintf(stderr, "*** Error: MaterialModelData is not of type GAS\n");
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

  use_cp = (cp>0 && cv<=0.0) ? true : false;
    
}

//------------------------------------------------------------------------------

#endif
