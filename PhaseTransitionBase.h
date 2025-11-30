/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _PHASE_TRANSITION_BASE_H_
#define _PHASE_TRANSITION_BASE_H_

#include <VarFcnBase.h>
#include <IoData.h>
#include <cassert>

/****************************************************************************
 * Classes defined in this file are responsible for determining whether 
 * phase/material transition should occur for a given state. 
 * One instantiation should be created for each pair of phases/materials 
 * (fromID, toID). References to their VarFcns are stored in this class.
 ***************************************************************************/

extern int verbose;

class PhaseTransitionBase {

protected:

  MaterialTransitionData::KineticsModel kinetics_model;

  int fromID, toID;

  VarFcnBase &vf1;
  VarFcnBase &vf2;

  double pmin, pmax, Tmin, Tmax;
  double latent_heat;

public:

  PhaseTransitionBase(MaterialTransitionData &data, VarFcnBase &vf1_, VarFcnBase &vf2_) 
     : fromID(data.from_id), toID(data.to_id),
       vf1(vf1_), vf2(vf2_),
       Tmin(data.temperature_lowerbound), Tmax(data.temperature_upperbound),
       pmin(data.pressure_lowerbound), pmax(data.pressure_upperbound),
       latent_heat(data.latent_heat),
       kinetics_model(data.kinetics)
  { }

  virtual ~PhaseTransitionBase() {}

  inline int FromID() {return fromID;}
  inline int ToID() {return toID;}
  inline double GetLatentHeat() {return latent_heat;}

  //-----------------------------------------------------------------------------
  // Get the amount of latent heat (lambda) that should be deposited over time dt
  //-----------------------------------------------------------------------------
  virtual double GetDeltaLambda(double lambda, [[maybe_unused]] double dt) {
    return std::max(0.0, lambda);
  }

  //-----------------------------------------------------------------------------
  // Deposits delta lambda after phase transition; update pressure using *vf2*
  //-----------------------------------------------------------------------------
  virtual void DepositDeltaLambdaAfterTransition(double *v, double &lambda, double dt,
                                                 double *delta_lam = NULL) {
    double e    = vf2.GetInternalEnergyPerUnitMass(v[0], v[4]);
    double dlam = GetDeltaLambda(lambda, dt);
    e      += dlam;
    lambda -= dlam;
    v[4] = vf2.GetPressure(v[0], e);

    if(delta_lam)
      *delta_lam = dlam;
  }
     
  //-----------------------------------------------------------------------------
  // Handles energy exchange between e and lambda, and checks for a phase
  // transition. If a transition is detected, computes the amount of lambda
  // converted to e (delta_lam) based on the kinetics model and timestep dt.
  // When no transition occurs, delta_lam (if provided) is set to zero.
  // Returns true if a phase transition is detected; false otherwise.
  //-----------------------------------------------------------------------------
  virtual bool Transition(double *v, double &lambda, double dt, double *delta_lam = NULL);

};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

inline bool
PhaseTransitionBase::Transition(double *v, double &lambda, double dt, double *delta_lam)
{
  // check temperature
  double e = vf1.GetInternalEnergyPerUnitMass(v[0], v[4]);
  double T = vf1.GetTemperature(v[0], e);

  if(delta_lam)
    *delta_lam = 0.0; //reset

  if(T<=Tmax) { //no phase transition, but if lambda (latent heat reservoir) is non-zero, should pour it 
                //back to raise temperature
    if(lambda>0) {
      double e_vap = vf1.GetInternalEnergyPerUnitMassFromTemperature(v[0], Tmax); 
      double de = e_vap - e;

      double dlam = std::min(de, lambda);

      lambda -= dlam;
      e      += dlam;

      v[4] = vf1.GetPressure(v[0], e);
    }

    return false;

  } else { // excessive heat should go to lambda. Then, check if latent heat is reached

    double e_vap = vf1.GetInternalEnergyPerUnitMassFromTemperature(v[0], Tmax);
    double de = e - e_vap;

    lambda += de;
    e      -= de; //= e_vap

    if(lambda >= latent_heat) {//DETECTED PHASE TRANSITION
      double dlam = GetDeltaLambda(lambda, dt);
      e      += dlam;
      lambda -= dlam;
      v[4] = vf2.GetPressure(v[0], e);
      if(delta_lam)
        *delta_lam = dlam;
      return true;
    } 
    else {
      v[4] = vf1.GetPressure(v[0], e);
      return false;
    }
  }

  // TODO: may also check pressure in future
}


#endif
