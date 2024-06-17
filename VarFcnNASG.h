/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _VAR_FCN_NASG_H
#define _VAR_FCN_NASG_H

#include <VarFcnBase.h>
#include <fstream>

/********************************************************************************
 * This class is the VarFcn class for the Noble-Abel Stiffened Gas equation of
 * state (EOS). Only elementary functions are declared and/or defined here.
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
  double cv; //!< specific heat at constant volume
  double bigC; //!< Constant of integration in temperature law. Maybe set to 0; (New: Dimension is [volume]/[mass])

  double gam1;    //!< gamma-1
  double invgam1; //!< 1/(gamma-1)
  double gam_pc; //!< gam*pc
  double invcv; //!< 1/cv


public:
  VarFcnNASG(MaterialModelData &data);
  ~VarFcnNASG() {}

  //! ----- EOS-Specific Functions -----
  inline double GetPressure(double rho, double e) {return gam1*(e-q)/(1.0/rho - b) - gam_pc;}
  inline double GetInternalEnergyPerUnitMass(double rho, double p) {return invgam1*(p+gam_pc)*(1.0/rho-b) + q;}
  inline double GetDensity(double p, double e) {return 1.0/(gam1*(e-q)/(p+gam_pc) + b);}
  inline double GetDpdrho(double rho, double e) {double V = 1.0/rho; return gam1*V*V*(e-q)/((V-b)*(V-b));}
  inline double GetBigGamma(double rho, [[maybe_unused]] double e) {return gam1/(1.0 - b*rho);}

  inline double GetTemperature(double rho, double e) {double V = 1.0/rho; return invcv*(V-b)*((e - q)/(V-b) - pc - pc*pow((bigC/(V-b)),gam)/gam1);}

  inline double GetReferenceTemperature() {return 0.0;}
  inline double GetReferenceInternalEnergyPerUnitMass() {return 0.0;}

  inline double GetInternalEnergyPerUnitMassFromTemperature(double rho, double T) {double V = 1.0/rho; return (cv*T/(V-b) + pc + pc*pow((bigC/(V-b)),gam)/gam1) * (V-b) + q;}
  
  inline double GetInternalEnergyPerUnitMassFromEnthalpy(double rho, double h) {
    double V = 1.0/rho;  return ((h+gam_pc*V)*(V-b)+V*gam1*q)/(gam*V-b);}


  //! Verify hyperbolicity (i.e. c^2 > 0): Report error if rho < 0 or p + pc < 0 (Not p + gamma*pc). 
  inline bool CheckState(double rho, double p, bool silence = false) {
    if(!std::isfinite(rho) || !std::isfinite(p)) {
      if(!silence)
        fprintf(stdout, "*** Error: CheckState failed. rho = %e, p = %e.\n", rho, p);
      return true;
    }
    if(rho <= 0.0 || p+pc <= 0.0){ //if p+pc<=0, c^2<=0
      if(!silence && verbose>1)
        fprintf(stdout, "Warning: Negative density or violation of hyperbolicity. rho = %e, p = %e.\n", rho, p);
      return true;
    }
    return false;
  }

};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

VarFcnNASG::VarFcnNASG(MaterialModelData &data) : VarFcnBase(data) {

  if(data.eos != MaterialModelData::NOBLE_ABEL_STIFFENED_GAS){
    fprintf(stdout, "*** Error: MaterialModelData is not of type NOBLE_ABEL_STIFFENED_GAS.\n");
    exit(-1);
  }

  type = NOBLE_ABEL_STIFFENED_GAS;

  gam    = data.nasgModel.specificHeatRatio;
  pc     = data.nasgModel.pressureConstant;
  q      = data.nasgModel.energyConstant;
  b      = data.nasgModel.volumeConstant;

  cv     = data.nasgModel.cv;
  bigC   = data.nasgModel.bigC;
  gam1 = gam -1.0;
  invgam1 = 1.0/gam1;
  gam_pc = gam*pc;
  invcv = cv==0.0 ? 0.0 : 1.0/cv;

}

//------------------------------------------------------------------------------

#endif
