/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _PHASE_TRANSITION_CONST_RATE_H_
#define _PHASE_TRANSITION_CONST_RATE_H_

#include<PhaseTransitionBase.h>

/****************************************************************************
 * A derived class for a constant phase transition rate (simplest kinetics model)
 ***************************************************************************/

class PhaseTransitionConstRate : public PhaseTransitionBase {

  double energy_transfer_rate; //!< [latent_heat]/[time]

public:

  PhaseTransitionConstRate(MaterialTransitionData &data, VarFcnBase &vf1_, VarFcnBase &vf2_) 
     : PhaseTransitionBase(data, vf1_, vf2_)
  {
    assert(kinetics_model == MaterialTransitionData::CONSTANT);
    if(data.constant_model.duration<=0.0) {
      fprintf(stdout, "*** Error: Duration of phase transition (kinetics model) should be positive (%e).\n",
              data.constant_model.duration);
      exit(-1);
    }
    energy_transfer_rate = latent_heat/data.constant_model.duration;
  }

  ~PhaseTransitionConstRate() {}

  double GetDeltaLambda(double lambda, double dt) {
    if(lambda<=0.0)
      return 0.0;
    return std::min(lambda, energy_transfer_rate*dt);
  }

};


#endif
