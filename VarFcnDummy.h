#ifndef _VAR_FCN_DUMMY_H
#define _VAR_FCN_DUMMY_H

#include <VarFcnBase.h>
#include <fstream>

/********************************************************************************
 * This class is the VarFcn class for an 'inactive' material, often assigned to
 * nodes that are inside an embedded solid body or 'covered' by an embedded surface 
 * with finite, but small, thickness.
 * This dummy material should only be instantiated internally, not by the user.
 ********************************************************************************/
class VarFcnDummy : public VarFcnBase {

private:
  double p0, rho0, e0, T0;
public:
  VarFcnDummy(StateVariable& sv) : VarFcnBase(sv), 
                                   p0(sv.pressure), rho0(sv.density), T0(sv.temperature),
                                   e0(sv.internal_energy_per_mass) {type = DUMMY;}
  ~VarFcnDummy() {}

  //! ----- EOS-Specific Functions -----
  inline double GetPressure(double rho, double e) const {return p0;}
  inline double GetInternalEnergyPerUnitMass(double rho, double p) const {return e0;}
  inline double GetDensity(double p, double e) const {return rho0;}
  inline double GetDpdrho(double rho, double e) const{return 0.0;}
  inline double GetBigGamma(double rho, double e) const {return 0.0;}
  inline double GetTemperature(double rho, double e) const {return T0;}
  inline double GetReferenceTemperature() const {return T0;}
  inline double GetReferenceInternalEnergyPerUnitMass() const {return e0;}
  inline double GetInternalEnergyPerUnitMassFromTemperature(double rho, double T) const {return e0;}
  inline double GetInternalEnergyPerUnitMassFromEnthalpy(double rho, double h) const {return e0;}
  inline bool   CheckState(double rho, double p, bool silence = false) const {return false;}
  inline bool   CheckState(double *V, bool silence = false) const {return false;}
  inline bool   CheckPhaseTransition(int id/*id of the other phase*/) const {return false;}

  //! Overwrite the calculations done in the base class
  inline void   ConservativeToPrimitive(double *U, double *V) {
    V[0] = rho0; V[1] = V[2] = V[3] = 0.0; V[4] = p0;}
  inline void   PrimitiveToConservative(double *V, double *U) {
    U[0] = rho0; U[1] = U[2] = U[3] = 0.0; U[4] = rho0*e0;}
  inline double ComputeSoundSpeed(double rho, double e) {return DBL_MIN;}
  inline double ComputeSoundSpeedSquare(double rho, double e) {return DBL_MIN;}
  inline double ComputeMachNumber(double *V) {return 0.0;}
  inline double ComputeEnthalpyPerUnitMass(double rho, double p) {return e0;}
  inline double ComputeTotalEnthalpyPerUnitMass(double *V) {return e0;}
  inline bool   ClipDensityAndPressure(double *V, double *U) {
    V[0] = rho0; V[1] = V[2] = V[3] = 0.0; V[4] = p0;
    if(U) {U[0] = rho0; U[1] = U[2] = U[3] = 0.0; U[4] = rho0*e0;}
    return false;
  }

};

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

#endif
