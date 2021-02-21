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

  FluxFcnLLF(std::vector<VarFcnBase*> &varFcn, IoData &iod) : FluxFcnBase(varFcn) { } 
    
  ~FluxFcnLLF() {}

  inline void ComputeNumericalFluxAtCellInterface(int dir /*0~x, 1~y, 2~z*/, double *Vm/*minus*/, 
                                                  double *Vp/*plus*/, int id, double *flux/*F,G,or H*/);

private:
  inline void ComputeMaxEigenvalue(int dir /*0~x, 1~y, 2~z*/, double *Vm, double *Vp, int id, double &a);

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
  double c_m = vf[id]->ComputeSoundSpeed(Vm[0], e_m);
  double c_p = vf[id]->ComputeSoundSpeed(Vp[0], e_p);

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
    print_error("Error: Dir. (%d) not recognized.\n", dir);
    exit_mpi();
  }

  double Um[5], Up[5];
  vf[id]->PrimitiveToConservative(Vm, Um);
  vf[id]->PrimitiveToConservative(Vp, Up);

  for(int i=0; i<nDOF; i++)
    flux[i] = 0.5*( (fm[i]+fp[i]) - a*(Up[i]-Um[i]) );
}

//----------------------------------------------------------------------------------------


#endif
