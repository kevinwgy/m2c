#ifndef _PHASE_TRANSITION_H_
#define _PHASE_TRANSITION_H_

#include <VarFcnBase.h>
#include <IoData.h>

/****************************************************************************
 * Classes defined in this file are responsible for determining whether 
 * phase/material transition should occur for a given state. It is NOT 
 * responsible for changing/updating the state variable when phase/material 
 * transition occurs. It is also not responsible for updating the level set 
 * function (phi). 
 * One instantiation should be created for each pair of phases/materials 
 * (fromID, toID). References to their VarFcns are stored in this class.
 ***************************************************************************/

extern int verbose;

class PhaseTransitionBase {

public:

  enum Type {BASE = 0} type;

  int fromID;
  int toID;
  VarFcnBase &vf1;
  VarFcnBase &vf2;

  double pmin, pmax, Tmin, Tmax;
  double latent_heat;
  double emax, emin; //internal energy that would trigger phase transition

  PhaseTransitionBase(MaterialTransitionData &data, VarFcnBase &vf1_, VarFcnBase &vf2_) 
     : fromID(data.from_id), toID(data.to_id),
       vf1(vf1_), vf2(vf2_),
       Tmin(data.temperature_lowerbound), Tmax(data.temperature_upperbound),
       pmin(data.pressure_lowerbound), pmax(data.pressure_upperbound),
       latent_heat(data.latent_heat)
  {
    type = BASE; //no other choices at the moment

    emax = vf1.GetInternalEnergyPerUnitMassFromTemperature(0.0, Tmax) + latent_heat;
    emin = vf1.GetInternalEnergyPerUnitMassFromTemperature(0.0, Tmin) - latent_heat;
  }

  virtual ~PhaseTransitionBase() {}

  virtual bool Transition(double *v) {

    // check temperature
    double e = vf1.GetInternalEnergyPerUnitMass(v[0], v[4]);
    if(e<=emin || e>emax)
      return true;

    // check pressure
    // TODO
    
    return false;
  }

};


#endif
