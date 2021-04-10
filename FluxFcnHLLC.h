#ifndef _FLUX_FCN_HLLC_H_
#define _FLUX_FCN_HLLC_H_

#include <FluxFcnBase.h>
#include <math.h>
using std::fabs;
/****************************************************************************************
 * Harten-Lax-van Leer-Contact (HLLC) flux for general EOS. Ref: Toro's book, Chapter 10
 ***************************************************************************************/

class FluxFcnHLLC : public FluxFcnBase {

private:

  double eps; //a small tolerance in Hu et al.'s paper to avoid divid-by-zero. (non-dimensional)

public:

  FluxFcnHLLC(std::vector<VarFcnBase*> &varFcn, IoData &iod) : FluxFcnBase(varFcn) 
    {eps = 1e-10;} //!< Hard-coded for the moment. If needed, can be made a user input (IoData)
    
  ~FluxFcnHLLC() {}

  inline void ComputeNumericalFluxAtCellInterface(int dir /*0~x, 1~y, 2~z*/, double *Vm/*minus*/, 
                                                  double *Vp/*plus*/, int id, double *flux/*F,G,or H*/);

private:
  //! Internal functions
  inline double Average(double w1, double f1, double w2, double f2);

  inline void ComputeMinMaxWaveSpeedsByRoeAverage(int dir /*0~x, 1~y, 2~z*/, double *Vm, double *Vp, int id,
                                                  double &Sm, double &Sp);

  inline void ComputeFstar(int dir /*0~x, 1~y, 2~z*/, double *V, double S, double Sstar, int id, 
                           double *Fstar);
};

//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
// Internal functions, Used only within this file
//----------------------------------------------------------------------------------------

inline
double FluxFcnHLLC::Average(double w1, double f1, double w2, double f2)
{
  return (w1*f1 + w2*f2)/(w1 + w2); 
}

//----------------------------------------------------------------------------------------

inline 
void FluxFcnHLLC::ComputeMinMaxWaveSpeedsByRoeAverage(int dir /*0~x, 1~y, 2~z*/, double *Vm, double *Vp, int id,
                                                       double &Sm, double &Sp)
{
  //Based on Hu et al., 2009 and Toro's book Chapter 10. Also see Kevin's notes
  
  // 1. Compute the intermediate state variables (with "hat")
  double rho_hat, u_hat, v_hat, w_hat, H_hat, c_hat;
  
  rho_hat = sqrt(Vm[0]*Vp[0]);

  double sqrt_rhom = sqrt(Vm[0]);
  double sqrt_rhop = sqrt(Vp[0]);

  double Hm = vf[id]->ComputeTotalEnthalpyPerUnitMass(Vm);
  double Hp = vf[id]->ComputeTotalEnthalpyPerUnitMass(Vp);
  H_hat = Average(sqrt_rhom, Hm, sqrt_rhop, Hp);

  // now calculate c_hat step by step
  double drho = Vp[0] - Vm[0];
  double du   = Vp[1] - Vm[1];
  double dv   = Vp[2] - Vm[2];
  double dw   = Vp[3] - Vm[3];
  double dp   = Vp[4] - Vm[4];

  double diffvelo;
  if     (dir==0) diffvelo = du; //x --> du;
  else if(dir==1) diffvelo = dv; //y --> dv;
  else            diffvelo = dw; //z --> dw;

  double diff = diffvelo/(sqrt_rhom+sqrt_rhop); 
  double p_over_rho_hat = Average(sqrt_rhom, Vm[4]/Vm[0], sqrt_rhop, Vp[4]/Vp[0])
                        + 0.5*diff*diff;

  double em   = vf[id]->GetInternalEnergyPerUnitMass(Vm[0],Vm[4]);
  double ep   = vf[id]->GetInternalEnergyPerUnitMass(Vp[0],Vp[4]);
  double de   = ep - em;
  double e_hat = Average(sqrt_rhom, em, sqrt_rhop, ep);

  double w_rho = drho*drho/(rho_hat*rho_hat);
  double w_e   = de*de/(e_hat*e_hat);
  double denominator =  w_rho + w_e + eps;

  double dpdrho_roeavg = Average(sqrt_rhom, vf[id]->GetDpdrho(Vm[0],em), 
                                 sqrt_rhop, vf[id]->GetDpdrho(Vp[0],ep));
  double Gamma_roeavg = Average(sqrt_rhom, vf[id]->GetBigGamma(Vm[0],em), 
                                sqrt_rhop, vf[id]->GetBigGamma(Vp[0],ep));

  double dpdrho_hat_numerator_term1 = (w_e + eps)*dpdrho_roeavg;
  double dpdrho_hat_numerator_term2 = (dp - Gamma_roeavg*rho_hat*de)*drho/(rho_hat*rho_hat); //avoid dividing by drho (possibly 0)
  double dpdrho_hat = (dpdrho_hat_numerator_term1 + dpdrho_hat_numerator_term2)/denominator;

  double Gamma_hat_numerator_term1 = (w_rho + eps)*Gamma_roeavg; 
  double Gamma_hat_numerator_term2 = (dp - dpdrho_roeavg*drho)*de/(rho_hat*e_hat*e_hat); //avoid dividing by de (possibly 0)
  double Gamma_hat = (Gamma_hat_numerator_term1 + Gamma_hat_numerator_term2)/denominator;

  double c_hat_square = dpdrho_hat + Gamma_hat*p_over_rho_hat;
  if(c_hat_square <= 0) {
    fprintf(stderr,"Warning: The artificial state in the generalized Roe flux function loses hyperbolicity (c_hat_square = %e). Setting c_hat = 0.\n", c_hat_square);
    c_hat = eps;
    c_hat_square = c_hat*c_hat;
  } else
    c_hat = sqrt(c_hat_square);
    
  double cm = vf[id]->ComputeSoundSpeedSquare(Vm[0], em);

  if(cm<0) {
    fprintf(stderr,"*** Error: c^2 (square of sound speed) = %e in HLLC flux function. Vm = %e, %e, %e, %e, %e, ID = %d.\n",
            cm, Vm[0], Vm[1], Vm[2], Vm[3], Vm[4], id);
    exit_mpi();
  } else
    cm = sqrt(cm);

  double cp = vf[id]->ComputeSoundSpeedSquare(Vp[0], ep);

  if(cp<0) {
    fprintf(stderr,"*** Error: c^2 (square of sound speed) = %e in HLLC flux function. Vp = %e, %e, %e, %e, %e, ID = %d.\n",
            cp, Vp[0], Vp[1], Vp[2], Vp[3], Vp[4], id);
    exit_mpi();
  } else
    cp = sqrt(cp);


  // 2. Calculate Sm and Sp, the minimum and maximum wave speeds.
  switch (dir) {
    case 0: //x-dir, i.e. for F
      u_hat = Average(sqrt_rhom, Vm[1], sqrt_rhop, Vp[1]);
      Sm = std::min(Vm[1] - cm, u_hat - c_hat);
      Sp = std::max(Vp[1] + cp, u_hat + c_hat); 
      break;
    case 1: //y-dir, i.e. for G
      v_hat = Average(sqrt_rhom, Vm[2], sqrt_rhop, Vp[2]);
      Sm = std::min(Vm[2] - cm, v_hat - c_hat);
      Sp = std::max(Vp[2] + cp, v_hat + c_hat); 
      break;
    case 2: //z-dir, i.e. for H
      w_hat = Average(sqrt_rhom, Vm[3], sqrt_rhop, Vp[3]);
      Sm = std::min(Vm[3] - cm, w_hat - c_hat);
      Sp = std::max(Vp[3] + cp, w_hat + c_hat); 
      break;
    default:
      fprintf(stderr,"*** Error: Incorrect use of function FluxFcnHLLC::ComputeMinMaxWaveSpeedsByRoeAverage.\n");
      exit_mpi();
  }
}

//----------------------------------------------------------------------------------------

inline 
void FluxFcnHLLC::ComputeFstar(int dir /*0~x, 1~y, 2~z*/, double *V, double S, double Sstar, int id,
                               double *Fstar)
{
  int nDOF = 5;
  double velo = V[dir+1];
  double rhostar = V[0]*(S - velo)/(S - Sstar);

  double U[nDOF];
  vf[id]->PrimitiveToConservative(V, U);

  // compute Ustar
  double Ustar[nDOF];
  Ustar[0] = 1.0;
  Ustar[1] = (dir==0) ? Sstar : V[1];
  Ustar[2] = (dir==1) ? Sstar : V[2];
  Ustar[3] = (dir==2) ? Sstar : V[3];
  Ustar[4] = U[4]/V[0] + (Sstar - velo)*(Sstar + V[4]/(V[0]*(S-velo)));

  for(int i=0; i<nDOF; i++)
    Ustar[i] *= rhostar;

  // compute Fstar
  switch (dir) {
    case 0 :
      EvaluateFluxFunction_F(V,id,Fstar); break;
    case 1 :
      EvaluateFluxFunction_G(V,id,Fstar); break;
    case 2 :
      EvaluateFluxFunction_H(V,id,Fstar); break;
  }

  for(int i=0; i<nDOF; i++)
    Fstar[i] += S*(Ustar[i] - U[i]);
}

//----------------------------------------------------------------------------------------

inline
void FluxFcnHLLC::ComputeNumericalFluxAtCellInterface(int dir, double *Vm, double *Vp, int id,
                                                      double *flux)
{
  // estimates of min and max wave speeds
  double Sm = 0.0, Sp = 0.0;
  ComputeMinMaxWaveSpeedsByRoeAverage(dir, Vm, Vp, id, Sm, Sp);

  // calculate the HLLC flux function
  if (0.0 <= Sm) {
    switch (dir) {
      case 0:
        EvaluateFluxFunction_F(Vm, id, flux); break;
      case 1:
        EvaluateFluxFunction_G(Vm, id, flux); break;
      case 2:
        EvaluateFluxFunction_H(Vm, id, flux); break;
    }
  } else if (0.0 >= Sp) {
    switch (dir) {
      case 0:
        EvaluateFluxFunction_F(Vp, id, flux); break;
      case 1:
        EvaluateFluxFunction_G(Vp, id, flux); break;
      case 2:
        EvaluateFluxFunction_H(Vp, id, flux); break;
    }
  } else {
    // now we need to find the "star" states
    double velo_m = Vm[dir+1];
    double velo_p = Vp[dir+1];
    double Sstar  = (Vp[4] - Vm[4] + Vm[0]*velo_m*(Sm - velo_m) - Vp[0]*velo_p*(Sp - velo_p))
                  / (Vm[0]*(Sm - velo_m) - Vp[0]*(Sp - velo_p));

    if (Sstar >= 0.0)
      ComputeFstar(dir, Vm, Sm, Sstar, id, flux);
    else if (Sstar <= 0.0)
      ComputeFstar(dir, Vp, Sp, Sstar, id, flux);
    else {
      fprintf(stderr,"*** Error: Logic error in FluxFcnHLLC::ComputeNumericalFluxAtCellInterface.\n");
      exit_mpi();
    }
  }
}

//----------------------------------------------------------------------------------------


#endif
