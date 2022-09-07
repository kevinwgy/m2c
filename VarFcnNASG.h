#ifndef _VAR_FCN_NASG_H
#define _VAR_FCN_NASG_H

#include <VarFcnBase.h>
#include <fstream>

/********************************************************************************
 * This class is the VarFcn class for the Noble-Abel Stiffened Gas EOS in Euler
 * Equations. Only elementary functions are declared and/or defined here.
 * All arguments must be pertinent to only a single grid node or a single
 * state.
 *
 * EOS: Pressure = (gam - 1)*(e - q)/(1/rho - b) - gam*Pc
 * where
 *   e  : internal energy per unit mass.
 *   Pc : pressure constant. 
 *   gam: specific heat ratio
 *   b  : specific volume constant
 *   q  : specific internal energy constant
 *   q' : specific entropy constant
 *
 *   cv : specific heat at constant volume
 *   cp = gam*cv
 *
 *   Temperature law: T = (e - q - Pc*(1/rho - b))/cv;
 *                    Equivalently, T = (h - q - b*p)/cp
 *   In general, de =/= cv*dT, and dh =/= cp*dT.
 ********************************************************************************/
class VarFcnNASG : public VarFcnBase {

private:
  double gam;
  double pc;
  double b; 
  double q; 
  double qprime;
  double cv; //!< specific heat at constant volume

  double invgam;  //!< 1/gamma
  double gam1;    //!< gamma-1
  double invgam1; //!< 1/(gamma-1)
  double invcv; //!< 1/cv
  double cp; //!< gam*cv
  double invcp; //!< 1/cp


public:
  VarFcnNASG(MaterialModelData &data);
  ~VarFcnNASG() {}

  //! ----- EOS-Specific Functions -----
  inline double GetPressure(double rho, double e) {return gam1*(e-q)/(1.0/rho - b) - gam*pc;}
  inline double GetInternalEnergyPerUnitMass(double rho, double p) {return (p+gam*pc)*(1.0/rho-b)/gam1 + q;}
  inline double GetDensity(double p, double e) {return 1.0/(gam1*(e-q)/(p+gam*pc) + b);}
  inline double GetDpdrho(double rho, double e) {double V = 1.0/rho; return gam1*V*V*(e-q)/((V-b)*(V-b));}
  inline double GetBigGamma(double rho, double e) {double V = 1.0/rho; return gam1*V/(V-b);}

  inline double GetTemperature(double rho, double e) {
    I AM HERE
  }

  inline double GetReferenceTemperature() {return T0;}
  inline double GetReferenceInternalEnergyPerUnitMass() {return e0;}

  inline double GetInternalEnergyPerUnitMassFromTemperature(double rho, double T) {
    if(use_cv_advanced) {
      return cv*T - pc/rho - pow(rho/rho0, gam1)*(cv*T0 - (e0 + pc/rho0));
    } else if(use_cp) 
      return invgam*(h0 + cp*(T-T0)) + pc/rho;
    else
      return e0 + cv*(T-T0);
  }
  
  inline double GetInternalEnergyPerUnitMassFromEnthalpy(double rho, double h) {return invgam*h+pc/rho;}


  //! Verify hyperbolicity (i.e. c^2 > 0): Report error if rho < 0 or p + pc < 0 (Not p + gamma*pc). 
  inline bool CheckState(double rho, double p, bool silence = false) {
    if(!std::isfinite(rho) || !std::isfinite(p)) {
      if(!silence)
        fprintf(stderr, "*** Error: CheckState failed. rho = %e, p = %e.\n", rho, p);
      return true;
    }
    if(rho <= 0.0 || p+pc <= 0.0){
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
VarFcnNASG::VarFcnNASG(MaterialModelData &data) : VarFcnBase(data) {

  if(data.eos != MaterialModelData::NOBLE_ABEL_STIFFENED_GAS){
    fprintf(stderr, "*** Error: MaterialModelData is not of type NOBLE_ABEL_STIFFENED_GAS.\n");
    exit(-1);
  }

  type = NOBLE_ABEL_STIFFENED_GAS;

  gam    = data.nasgModel.specificHeatRatio;
  pc     = data.nasgModel.pressureConstant;
  q      = data.nasgModel.energyConstant;
  qprime = data.nasgModel.entropyConstant;
  b      = data.nasgModel.volumeConstant;

  cv     = data.nasgModel.cv;

  invgam = (gam==0.0) ? 0.0 : 1.0/gam;
  gam1 = gam -1.0;
  invgam1 = 1.0/gam1;
  invcv = cv==0.0 ? 0.0 : 1.0/cv;
  cp = gam*cv;
  invcp = cp==0.0 ? 0.0 : 1.0/cp;

}

//------------------------------------------------------------------------------

#endif
