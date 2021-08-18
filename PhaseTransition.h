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

  PhaseTransitionBase(MaterialTransitionData &data, VarFcnBase &vf1_, VarFcnBase &vf2_) 
     : fromID(data.from_id), toID(data.to_id),
       vf1(vf1_), vf2(vf2_),
       Tmin(data.temperature_lowerbound), Tmax(data.temperature_upperbound),
       pmin(data.pressure_lowerbound), pmax(data.pressure_upperbound)
  {
    type = BASE; //no other choices at the moment
  }

  virtual ~PhaseTransitionBase() {}

  virtual bool Transition(double *v) {
    // check pressure
    if(v[0]<=pmin || v[0]>pmax)
      return true;
    // check temperature
    double e = vf1.GetInternalEnergyPerUnitMass(v[0], v[4]);
    double T = vf1.GetTemperature(v[0], e);
    if(T<=Tmin || T>Tmax)
      return true;

    return false;
  }

};


#endif
