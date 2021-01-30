#ifndef _FLUX_FCN_LLF_H_
#define _FLUX_FCN_LLF_H_

#include <FluxFcnBase.h>
#include <Vector5D.h>
#include <math.h>
using std::fabs;
using std::max;
/****************************************************************************************
 * The local Lax-Friedrichs Flux (a.k.a. Rusanov's method), for general EOS. Ref: Leveque, Chap 12.5
 ***************************************************************************************/

class FluxFcnLLF : public FluxFcnBase {

public:

  FluxFcnLLF(VarFcnBase *varFcn, IoData &iod) : FluxFcnBase(varFcn) { } 
    
  ~FluxFcnLLF() {}

  inline void ComputeNumericalFluxAtCellInterface(int dir /*0~x, 1~y, 2~z*/, double *Vm/*minus*/, 
                                                  double *Vp/*plus*/, double *flux/*F,G,or H*/);

private:
  inline void ComputeMaxEigenvalue(int dir /*0~x, 1~y, 2~z*/, double *Vm, double *Vp, double &a);

};

//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
// Internal function, Used only within this file
//----------------------------------------------------------------------------------------

inline 
void FluxFcnLLF::ComputeMaxEigenvalue(int dir /*0~x, 1~y, 2~z*/, double *Vm, double *Vp, double &a)
{
  
  double velocity_m, velocity_p;
  if     (dir==0) {velocity_m = Vm[1];  velocity_p = Vp[1];} //x --> u;
  else if(dir==1) {velocity_m = Vm[2];  velocity_p = Vp[2];} //y --> v;
  else            {velocity_m = Vm[3];  velocity_p = Vp[3];} //z --> w;

  double e_m = vf->GetInternalEnergyPerUnitMass(Vm[0], Vm[4]);
  double e_p = vf->GetInternalEnergyPerUnitMass(Vp[0], Vp[4]);
  double c_m = vf->ComputeSoundSpeed(Vm[0], e_m);
  double c_p = vf->ComputeSoundSpeed(Vp[0], e_p);

  a = max( fabs(velocity_m) + c_m, 
           fabs(velocity_p) + c_p );
}

//----------------------------------------------------------------------------------------

inline
void FluxFcnLLF::ComputeNumericalFluxAtCellInterface(int dir, double *Vm, double *Vp, double *flux)
{
  // Compute lambda, alpha, and R
  int nDOF = 5;
  double a;
  ComputeMaxEigenvalue(dir/*0~x,1~y,2~z*/, Vm, Vp, a);

  //LLF flux
  double fm[5], fp[5];
  EvaluateFluxFunction_F(Vm,fm);
  EvaluateFluxFunction_F(Vp,fp);

  double Um[5], Up[5];
  vf->PrimitiveToConservative(Vm, Um);
  vf->PrimitiveToConservative(Vp, Up);

  for(int i=0; i<nDOF; i++)
    flux[i] = 0.5*( (fm[i]+fp[i]) - a*(Up[i]-Um[i]) );
}

//----------------------------------------------------------------------------------------


#endif
