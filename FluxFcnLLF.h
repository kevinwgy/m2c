/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _FLUX_FCN_LLF_H_
#define _FLUX_FCN_LLF_H_

#include <FluxFcnBase.h>
#include <math.h>
using std::fabs;
using std::max;
/****************************************************************************************
 * The local Lax-Friedrichs Flux (a.k.a. Rusanov's method), for general EOS. Ref: Leveque, Chap 12.5
 ***************************************************************************************/

class FluxFcnLLF : public FluxFcnBase {

public:

  FluxFcnLLF(std::vector<VarFcnBase*> &varFcn, [[maybe_unused]] IoData &iod) : FluxFcnBase(varFcn) { } 
    
  ~FluxFcnLLF() {}

  inline void ComputeNumericalFluxAtCellInterface(int dir /*0~x, 1~y, 2~z*/, double *Vm/*minus*/, 
                                                  double *Vp/*plus*/, int id, double *flux/*F,G,or H*/);

  inline void ComputeNumericalFluxAtMaterialInterface(int dir/*0~x,1~y,2~z*/, double *Vminus/*left*/, int idm,
                                                       double *Vplus/*right*/, int idp, double *F);

private:

  inline void ComputeMaxEigenvalue(int dir /*0~x, 1~y, 2~z*/, double *Vm, double *Vp, int id, double &a);

  inline void ComputeMaxEigenvalueAtMaterialInterface(int dir /*0~x, 1~y, 2~z*/, double *Vm, int idm,
                                                      double *Vp, int idp, double &a);

};

//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
// Internal function, Used only within this file
//----------------------------------------------------------------------------------------

inline 
void FluxFcnLLF::ComputeMaxEigenvalue(int dir /*0~x, 1~y, 2~z*/, double *Vm, double *Vp, int id, double &a)
{
  
  double velocity_m, velocity_p;
  if     (dir==0) {velocity_m = Vm[1];  velocity_p = Vp[1];} //x --> u;
  else if(dir==1) {velocity_m = Vm[2];  velocity_p = Vp[2];} //y --> v;
  else            {velocity_m = Vm[3];  velocity_p = Vp[3];} //z --> w;

  double e_m = vf[id]->GetInternalEnergyPerUnitMass(Vm[0], Vm[4]);
  double e_p = vf[id]->GetInternalEnergyPerUnitMass(Vp[0], Vp[4]);
  double c_m = vf[id]->ComputeSoundSpeedSquare(Vm[0], e_m);

  if(c_m<0) {
    fprintf(stdout,"*** Error: c^2 (square of sound speed) = %e in LLF flux function. Vm = %e, %e, %e, %e, %e, ID = %d.\n",
            c_m, Vm[0], Vm[1], Vm[2], Vm[3], Vm[4], id);
    exit(-1);
  } else
    c_m = sqrt(c_m);

  double c_p = vf[id]->ComputeSoundSpeedSquare(Vp[0], e_p);

  if(c_p<0) {
    fprintf(stdout,"*** Error: c^2 (square of sound speed) = %e in LLF flux function. Vp = %e, %e, %e, %e, %e, ID = %d.\n",
            c_p, Vp[0], Vp[1], Vp[2], Vp[3], Vp[4], id);
    exit(-1);
  } else
    c_p = sqrt(c_p);

  a = max( fabs(velocity_m) + c_m, 
           fabs(velocity_p) + c_p );
}

//----------------------------------------------------------------------------------------

inline
void FluxFcnLLF::ComputeNumericalFluxAtCellInterface(int dir, double *Vm, double *Vp, int id, double *flux)
{
  // Compute lambda, alpha, and R
  int nDOF = 5;
  double a;
  ComputeMaxEigenvalue(dir/*0~x,1~y,2~z*/, Vm, Vp, id, a);

  //LLF flux
  double fm[5], fp[5];
  if(dir==0) {
    EvaluateFluxFunction_F(Vm, id, fm);
    EvaluateFluxFunction_F(Vp, id, fp);
  } else if (dir==1) {
    EvaluateFluxFunction_G(Vm, id, fm);
    EvaluateFluxFunction_G(Vp, id, fp);
  } else if (dir==2) {
    EvaluateFluxFunction_H(Vm, id, fm);
    EvaluateFluxFunction_H(Vp, id, fp);
  } else {
    fprintf(stdout, "*** Error: Dir. (%d) not recognized.\n", dir);
    exit(-1);
  }

  double Um[5], Up[5];
  vf[id]->PrimitiveToConservative(Vm, Um);
  vf[id]->PrimitiveToConservative(Vp, Up);

  for(int i=0; i<nDOF; i++)
    flux[i] = 0.5*( (fm[i]+fp[i]) - a*(Up[i]-Um[i]) );
}


//----------------------------------------------------------------------------------------
// Internal function, Used only within this file
//----------------------------------------------------------------------------------------

inline
void FluxFcnLLF::ComputeMaxEigenvalueAtMaterialInterface(int dir, double *Vm, int idm,
                                                         double *Vp, int idp, double &a)
{
  
  double velocity_m, velocity_p;
  if     (dir==0) {velocity_m = Vm[1];  velocity_p = Vp[1];} //x --> u;
  else if(dir==1) {velocity_m = Vm[2];  velocity_p = Vp[2];} //y --> v;
  else            {velocity_m = Vm[3];  velocity_p = Vp[3];} //z --> w;

  double e_m = vf[idm]->GetInternalEnergyPerUnitMass(Vm[0], Vm[4]);
  double e_p = vf[idp]->GetInternalEnergyPerUnitMass(Vp[0], Vp[4]);
  double c_m = vf[idm]->ComputeSoundSpeedSquare(Vm[0], e_m);

  if(c_m<0) {
    fprintf(stdout,"*** Error: c^2 (square of sound speed) = %e in LLF flux function. Vm = %e, %e, %e, %e, %e, ID = %d.\n",
            c_m, Vm[0], Vm[1], Vm[2], Vm[3], Vm[4], idm);
    exit(-1);
  } else
    c_m = sqrt(c_m);

  double c_p = vf[idp]->ComputeSoundSpeedSquare(Vp[0], e_p);

  if(c_p<0) {
    fprintf(stdout,"*** Error: c^2 (square of sound speed) = %e in LLF flux function. Vp = %e, %e, %e, %e, %e, ID = %d.\n",
            c_p, Vp[0], Vp[1], Vp[2], Vp[3], Vp[4], idp);
    exit(-1);
  } else
    c_p = sqrt(c_p);

  a = max( fabs(velocity_m) + c_m, 
           fabs(velocity_p) + c_p );
}

//----------------------------------------------------------------------------------------

inline
void FluxFcnLLF::ComputeNumericalFluxAtMaterialInterface(int dir, double *Vm, int idm,
                                                         double *Vp, int idp, double *flux)
{
  // Compute lambda, alpha, and R
  int nDOF = 5;
  double a;
  ComputeMaxEigenvalueAtMaterialInterface(dir/*0~x,1~y,2~z*/, Vm, idm, Vp, idp, a);

  //LLF flux
  double fm[5], fp[5];
  if(dir==0) {
    EvaluateFluxFunction_F(Vm, idm, fm);
    EvaluateFluxFunction_F(Vp, idp, fp);
  } else if (dir==1) {
    EvaluateFluxFunction_G(Vm, idm, fm);
    EvaluateFluxFunction_G(Vp, idp, fp);
  } else if (dir==2) {
    EvaluateFluxFunction_H(Vm, idm, fm);
    EvaluateFluxFunction_H(Vp, idp, fp);
  } else {
    fprintf(stdout, "*** Error: Dir. (%d) not recognized.\n", dir);
    exit(-1);
  }

  double Um[5], Up[5];
  vf[idm]->PrimitiveToConservative(Vm, Um);
  vf[idp]->PrimitiveToConservative(Vp, Up);

  for(int i=0; i<nDOF; i++)
    flux[i] = 0.5*( (fm[i]+fp[i]) - a*(Up[i]-Um[i]) );
}
#endif
