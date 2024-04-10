/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _PHASE_TRANSITION_H_
#define _PHASE_TRANSITION_H_

#include <VarFcnBase.h>
#include <IoData.h>
#include <cassert>

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
  //double emax, emin; //internal energy that would trigger phase transition

  PhaseTransitionBase(MaterialTransitionData &data, VarFcnBase &vf1_, VarFcnBase &vf2_) 
     : fromID(data.from_id), toID(data.to_id),
       vf1(vf1_), vf2(vf2_),
       Tmin(data.temperature_lowerbound), Tmax(data.temperature_upperbound),
       pmin(data.pressure_lowerbound), pmax(data.pressure_upperbound),
       latent_heat(data.latent_heat)
  {
    type = BASE; //no other choices at the moment

//    emax = vf1.GetInternalEnergyPerUnitMassFromTemperature(0.0, Tmax) + latent_heat;
//    emin = vf1.GetInternalEnergyPerUnitMassFromTemperature(0.0, Tmin) - latent_heat;
  }

  virtual ~PhaseTransitionBase() {}

  inline int FromID() {return fromID;}
  inline int ToID() {return toID;}


  virtual bool Transition(double *v, double &lambda) {

    // check temperature
    double e = vf1.GetInternalEnergyPerUnitMass(v[0], v[4]);
    double T = vf1.GetTemperature(v[0], e);

    if(T<=Tmax) { //no phase transition, but if lambda (latent heat reservoir) is non-zero, should pour it 
                  //back to raise temperature
      if(lambda>0) {
        double e_vap = vf1.GetInternalEnergyPerUnitMassFromTemperature(v[0], Tmax); 
        //assert(e_vap >= e); //this may be violated due to round-off error
        double de = e_vap - e;

        double dlam = std::min(de, lambda);

        lambda -= dlam;
        e      += dlam;

        v[4] = vf1.GetPressure(v[0], e);
      }

      return false;

    } else { // excessive heat should go to lambda. Then, check if latent heat is reached

      double e_vap = vf1.GetInternalEnergyPerUnitMassFromTemperature(v[0], Tmax);
      //assert(e_vap < e);  //this may be violated due to round-off error

      double de = e - e_vap;

      lambda += de;
      e      -= de; //= e_vap

      v[4] = vf1.GetPressure(v[0], e);

      if(lambda >= latent_heat) {//DETECTED PHASE TRANSITION

        //double h = e + v[4]/v[0]; //enthalpy before phase transition (using vf1)
        //h += lambda; //enthalpy after phase transition
        e += lambda; //Correction. Check Xuning JFM
        lambda = 0.0;
        //e = vf2.GetInternalEnergyPerUnitMassFromEnthalpy(v[0], h); //rho is always fixed
        v[4] = vf2.GetPressure(v[0], e);
        return true;

      } else
        return false;

    }

    // TODO: may also check pressure in future
    
  }

};


#endif
