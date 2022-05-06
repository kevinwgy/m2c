#ifndef _VAR_FCN_DUMMY_H
#define _VAR_FCN_DUMMY_H

#include <VarFcnBase.h>
#include <fstream>

/********************************************************************************
 * This class is the VarFcn class for an 'inactive' material, often assigned to
 * nodes that are inside a solid body or 'covered' by an embedded surface with
 * finite, but small, thickness.
 * This dummy material should only be instantiated internally, not by the user.
 ********************************************************************************/
class VarFcnDummy : public VarFcnBase {

private:
  double p0, rho0, e0, T0;
public:
  VarFcnDummy(StateVariable& sv) : VarFcnBase(sv), 
                                   p0(sv.pressure), rho0(sv.density), T0(sv.temperature),
                                   e0(sv.internal_energy_per_mass) {}
  ~VarFcnDummy() {}

  //! ----- EOS-Specific Functions -----
  inline double GetPressure(double rho, double e) const {return p0;}
  inline double GetInternalEnergyPerUnitMass(double rho, double p) const {return e0;}
  inline double GetDensity(double p, double e) const {return rho0;}
  inline double GetDpdrho(double rho, double e) const{return 0.0;}
  inline double GetBigGamma(double rho, double e) const {return 0.0;}
  inline double GetTemperature(double rho, double e) const {return T0;}
  inline double GetReferenceTemperature() const {return T0;}
  inline double GetInternalEnergyPerUnitMassFromTemperature(double rho, double T) const {return e0;}
  inline double GetInternalEnergyPerUnitMassFromEnthalpy(double rho, double h) const {return e0;}
  inline bool   CheckState(double rho, double p, bool silence = false) const {return false;}

};

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

#endif
