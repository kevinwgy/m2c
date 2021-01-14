#ifndef _VAR_FCN_SGEULER_H
#define _VAR_FCN_SGEULER_H

#include <VarFcnBase.h>
#include <fstream>
#include <Utils.h>

/********************************************************************************
 * This class is the VarFcn class for the Stiffened Gas EOS in Euler
 * Equations. Only elementary functions are declared and/or defined here.
 * All arguments must be pertinent to only a single grid node or a single
 * state.
 *
 * lay-out of the base class is:
 *  - 1 -  Transformation Operators
 *  - 2 -  General Functions
 *  - 3 -  Equations of State Parameters
 *  - 4 -  EOS related functions
 *
 *--------------------------------------------------------------------------
 *
 * EOS: Pressure = (gam - 1)*Density*e - gam*Pc
 * where
 *   e  : internal energy per unit mass.
 *   Pc : pressure constant. (named Pstiff below)
 *
 *   Note that the complete EOS is given by the above and h = cp * T
 *   where cp is constant. For a perfect gas, this leads to epsilon = cv * T
 *   but for a stiffened gas, epsilon = cv * T does not hold.
 *   For a stiffened gas, choosing epsilon = cv * T would lead to a non-constant cp...
 *
 ********************************************************************************/
class VarFcnSGEuler : public VarFcnBase {

private:
  double gam;
  double Pstiff;

  double gam1; //!< gamma-1
  double invgam1; //!< 1/(gamma-1)

public:
  VarFcnSGEuler(FluidModelData &data);
  ~VarFcnSGEuler() {}

  //! ----- EOS-Specific Functions -----
  inline double GetPressure(double rho, double e) const {return gam1*rho*e - gam*Pstiff;}
  inline double GetInternalEnergyPerUnitMass(double rho, double p) const {return (p+gam*Pstiff)/(gam1*rho);}
  inline double GetDensity(double p, double e) const {return (p+gam*Pstiff)/(gam1*e);}
  inline double GetDpdrho(double rho, double e) const{return gam1*e;}
  inline double GetBigGamma(double rho, double e) const {return gam1;}

  //! Verify hyperbolicity (i.e. c^2 > 0): Report error if rho < 0 or p + Pstiff < 0 (Not p + gamma*Pstiff). 
  inline bool CheckState(double *V) const{
    if(V[0] <= 0.0 || V[4]+Pstiff <= 0.0){
      fprintf(stdout, "Warning:  found negative density (%e) or negative pressure (p = %e, Pstiff = %e).\n", V[0], V[4], Pstiff);
      return true;
    }
    return false;
  }

};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
inline
VarFcnSGEuler::VarFcnSGEuler(FluidModelData &data) : VarFcnBase(data) {

  if(data.fluid != FluidModelData::STIFFENED_GAS){
    fprintf(stderr, "*** Error: FluidModelData is not of type GAS\n");
    exit(1);
  }

  if(data.gasModel.type == GasModelData::STIFFENED)
    type = STIFFENEDGAS;
  else
    fprintf(stdout, "*** Error: VarFcnSGEuler::type is undefined since data.gasModel.type = %d\n", data.gasModel.type);

  gam = data.gasModel.specificHeatRatio;
  gam1 = gam -1.0;
  invgam1 = 1.0/gam1;
  Pstiff = data.gasModel.pressureConstant;

}

//------------------------------------------------------------------------------

#endif
