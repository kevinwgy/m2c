/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

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
  VarFcnDummy(StateVariable& sv) : VarFcnBase(), 
                                   p0(sv.pressure), rho0(sv.density), T0(sv.temperature),
                                   e0(sv.internal_energy_per_mass) {type = DUMMY;}
  ~VarFcnDummy() {}

  //! ----- EOS-Specific Functions -----
  inline double GetPressure([[maybe_unused]] double rho, [[maybe_unused]] double e) {return p0;}
  inline double GetInternalEnergyPerUnitMass([[maybe_unused]] double rho, [[maybe_unused]] double p) {return e0;}
  inline double GetDensity([[maybe_unused]] double p, [[maybe_unused]] double e) {return rho0;}
  inline double GetDpdrho([[maybe_unused]] double rho, [[maybe_unused]] double e) {return 0.0;}
  inline double GetBigGamma([[maybe_unused]] double rho, [[maybe_unused]] double e) {return 0.0;}
  inline double GetTemperature([[maybe_unused]] double rho, [[maybe_unused]] double e) {return T0;}
  inline double GetReferenceTemperature() {return T0;}
  inline double GetReferenceInternalEnergyPerUnitMass() {return e0;}
  inline double GetInternalEnergyPerUnitMassFromTemperature([[maybe_unused]] double rho, [[maybe_unused]] double T) {return e0;}
  inline double GetInternalEnergyPerUnitMassFromEnthalpy([[maybe_unused]] double rho, [[maybe_unused]] double h) {return e0;}
  inline bool   CheckState([[maybe_unused]] double rho, [[maybe_unused]] double p, [[maybe_unused]] bool silence = false) {return false;}
  inline bool   CheckState([[maybe_unused]] double *V, [[maybe_unused]] bool silence = false) {return false;}
  inline bool   CheckPhaseTransition([[maybe_unused]] int id/*id of the other phase*/) {return false;}

  //! Overwrite the calculations done in the base class
  inline void   ConservativeToPrimitive([[maybe_unused]] double *U, double *V) {
    V[0] = rho0; V[1] = V[2] = V[3] = 0.0; V[4] = p0;}
  inline void   PrimitiveToConservative([[maybe_unused]] double *V, double *U) {
    U[0] = rho0; U[1] = U[2] = U[3] = 0.0; U[4] = rho0*e0;}
  inline double ComputeSoundSpeed([[maybe_unused]] double rho, [[maybe_unused]] double e) {return DBL_MIN;}
  inline double ComputeSoundSpeedSquare([[maybe_unused]] double rho, [[maybe_unused]] double e) {return DBL_MIN;}
  inline double ComputeMachNumber([[maybe_unused]] double *V) {return 0.0;}
  inline double ComputeEnthalpyPerUnitMass([[maybe_unused]] double rho, [[maybe_unused]] double p) {return e0;}
  inline double ComputeTotalEnthalpyPerUnitMass([[maybe_unused]] double *V) {return e0;}
  inline bool   ClipDensityAndPressure(double *V, double *U) {
    V[0] = rho0; V[1] = V[2] = V[3] = 0.0; V[4] = p0;
    if(U) {U[0] = rho0; U[1] = U[2] = U[3] = 0.0; U[4] = rho0*e0;}
    return false;
  }

};

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

#endif
