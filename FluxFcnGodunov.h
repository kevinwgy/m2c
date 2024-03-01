/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _FLUX_FCN_GODUNOV_H_
#define _FLUX_FCN_GODUNOV_H_

#include <FluxFcnBase.h>
#include <ExactRiemannSolverBase.h>

/****************************************************************************************
 * The Godunov flux, based on solving the exact Riemann problem
 ***************************************************************************************/

class FluxFcnGodunov : public FluxFcnBase {

public:

  FluxFcnGodunov(std::vector<VarFcnBase*> &varFcn, IoData &iod) 
      : FluxFcnBase(varFcn), riemann(varFcn,iod.exact_riemann) { } 
    
  ~FluxFcnGodunov() {}

  inline void ComputeNumericalFluxAtCellInterface(int dir /*0~x, 1~y, 2~z*/, double *Vm/*minus*/, 
                                                  double *Vp/*plus*/, int id, double *flux/*F,G,or H*/);

private:

  ExactRiemannSolverBase riemann;

};

//----------------------------------------------------------------------------------------

inline
void FluxFcnGodunov::ComputeNumericalFluxAtCellInterface(int dir, double *Vm, double *Vp, int id, double *flux)
{

  int nDOF = 5;
  double Vmid[nDOF], Vsm[nDOF], Vsp[nDOF];
  int midid;

  Vec3D normal(0.0,0.0,0.0);
  normal[dir] = 1.0;

  riemann.ComputeRiemannSolution(normal, Vm, id, Vp, id, Vmid, midid, Vsm, Vsp);
 
  if(dir==0) 
    EvaluateFluxFunction_F(Vmid, id, flux);
  else if(dir == 1)
    EvaluateFluxFunction_G(Vmid, id, flux);
  else if(dir == 2)
    EvaluateFluxFunction_H(Vmid, id, flux);
  else {
    fprintf(stdout, "*** Error: Dir. (%d) not recognized.\n", dir);
    exit(-1);
  }

}

//----------------------------------------------------------------------------------------


#endif
