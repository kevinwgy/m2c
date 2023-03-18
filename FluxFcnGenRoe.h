/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _FLUX_FCN_GEN_ROE_H_
#define _FLUX_FCN_GEN_ROE_H_

#include <FluxFcnBase.h>
#include <Vector5D.h>
#include <math.h>
using std::fabs;

/****************************************************************************************
 * Roe's Flux (Roe-Pike Flux) for general EOS. Ref: Hu, Aadams, Iaccarino, JCP 2009
 ***************************************************************************************/

class FluxFcnGenRoe : public FluxFcnBase {

private:

  double del; //a small tolerance in Harten's entropy fix (non-dimensional, multiplied to the max char. speed)
  double eps; //a small tolerance in Hu et al.'s paper to avoid divid-by-zero. (non-dimensional)

public:

  FluxFcnGenRoe(std::vector<VarFcnBase*> &varFcn, IoData &iod) : FluxFcnBase(varFcn) {
    del = iod.schemes.ns.delta; eps = 1e-10;} //!< Hard-coded for the moment. If needed, can be made a user input (IoData)
    
  ~FluxFcnGenRoe() {}

  inline void ComputeNumericalFluxAtCellInterface(int dir /*0~x, 1~y, 2~z*/, double *Vm/*minus*/, 
                                                  double *Vp/*plus*/, int id, double *flux/*F,G,or H*/);

private:
  //! Internal functions
  inline double Average(double w1, double f1, double w2, double f2);

  inline void ComputeLambdaAlphaR(int dir /*0~x, 1~y, 2~z*/, double *Vm, double *Vp, int id,
                                  double &lam1, double &lam2, double &lam3, double &lam4, double &lam5,
                                  double &a1, double &a2, double &a3, double &a4, double &a5,
                                  double *r1, double *r2, double *r3, double *r4, double *r5);

  //! Harten's entropy fix. See LeVeque's book, Section 15.3.5.
  inline void ApplyEntropyFix(double &lam, double tol); //lam will be modified
                      
};

//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
// Internal functions, Used only within this file
//----------------------------------------------------------------------------------------

inline
double FluxFcnGenRoe::Average(double w1, double f1, double w2, double f2)
{
  return (w1*f1 + w2*f2)/(w1 + w2); 
}

//----------------------------------------------------------------------------------------

inline 
void FluxFcnGenRoe::ComputeLambdaAlphaR(int dir /*0~x, 1~y, 2~z*/, double *Vm, double *Vp, int id,
                        double &lam1, double &lam2, double &lam3, double &lam4, double &lam5,
                        double &a1, double &a2, double &a3, double &a4, double &a5,
                        double *r1, double *r2, double *r3, double *r4, double *r5)
{
  //Based on Hu et al., 2009. Also see Kevin's notes
  
  // 1. Compute the intermediate state variables (with "hat")
  double rho_hat, u_hat, v_hat, w_hat, H_hat, c_hat;
  
  rho_hat = sqrt(Vm[0]*Vp[0]);

  double sqrt_rhom = sqrt(Vm[0]);
  double sqrt_rhop = sqrt(Vp[0]);

  u_hat = Average(sqrt_rhom, Vm[1], sqrt_rhop, Vp[1]);
  v_hat = Average(sqrt_rhom, Vm[2], sqrt_rhop, Vp[2]);
  w_hat = Average(sqrt_rhom, Vm[3], sqrt_rhop, Vp[3]);

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
    fprintf(stdout,"Warning: The artificial state in the generalized Roe flux function loses hyperbolicity (c_hat_square = %e). Setting c_hat = %e.\n", c_hat_square, eps);
    fprintf(stdout,"Vm = %e %e %e %e %e, Vp = %e %e %e %e %e, dir = %d, id = %d \n",
            Vm[0], Vm[1], Vm[2], Vm[3], Vm[4], Vp[0], Vp[1], Vp[2], Vp[3], Vp[4], dir, id);
    c_hat = eps;
    c_hat_square = c_hat*c_hat;
  } else
    c_hat = sqrt(c_hat_square);
    

  // 2. Now, we calculate lambda, alpha, and R
  //    (lambda and R follow take the *form* of the exact eigenvalues & eigenvectors of the Jacobian matrix.)
  switch (dir) {
    case 0: //x-dir, i.e. for F
      EvaluateEigensOfJacobian_F(u_hat, v_hat, w_hat, e_hat, c_hat, H_hat, Gamma_hat, dpdrho_hat, 
                                 lam1, lam2, lam3, lam4, lam5, r1, r2, r3, r4, r5);
      a1 = (dp - rho_hat*c_hat*du)/(2*c_hat_square);
      a2 = -dp/c_hat_square + drho; 
      a3 = rho_hat*dv;
      a4 = rho_hat*dw;
      a5 = (dp + rho_hat*c_hat*du)/(2*c_hat_square);
      break;
    case 1: //y-dir, i.e. for G
      EvaluateEigensOfJacobian_G(u_hat, v_hat, w_hat, e_hat, c_hat, H_hat, Gamma_hat, dpdrho_hat, 
                                 lam1, lam2, lam3, lam4, lam5, r1, r2, r3, r4, r5);
      a1 = (dp - rho_hat*c_hat*dv)/(2*c_hat_square);
      a2 = rho_hat*du;
      a3 = -dp/c_hat_square + drho;
      a4 = rho_hat*dw;
      a5 = (dp + rho_hat*c_hat*dv)/(2*c_hat_square);
      break;
    case 2: //z-dir, i.e. for H
      EvaluateEigensOfJacobian_H(u_hat, v_hat, w_hat, e_hat, c_hat, H_hat, Gamma_hat, dpdrho_hat, 
                                 lam1, lam2, lam3, lam4, lam5, r1, r2, r3, r4, r5);
      a1 = (dp - rho_hat*c_hat*dw)/(2*c_hat_square);
      a2 = rho_hat*du;
      a3 = rho_hat*dv;
      a4 = -dp/c_hat_square + drho;
      a5 = (dp + rho_hat*c_hat*dw)/(2*c_hat_square);
      break;
    default:
      fprintf(stdout,"*** Error: Incorrect use of function FluxFcnGenRoe::ComputeLambdaAlphaR.\n");
      exit_mpi();
  }
}

//----------------------------------------------------------------------------------------

inline
void FluxFcnGenRoe::ComputeNumericalFluxAtCellInterface(int dir, double *Vm, double *Vp, int id, double *flux)
{
  // Compute lambda, alpha, and R
  int nDOF = 5;
  double lam[nDOF], a[nDOF];
  Vec5D r[nDOF];
  ComputeLambdaAlphaR(dir/*0~x,1~y,2~z*/, Vm, Vp, id,
                      lam[0], lam[1], lam[2], lam[3], lam[4], a[0], a[1], a[2], a[3], a[4], 
                      r[0], r[1], r[2], r[3], r[4]);

  // Entropy fix
  double tol = lam[0];
  for(int i=1; i<5; i++)
    if(fabs(lam[i])>tol)
      tol = fabs(lam[i]);
  tol *= del;
  for(int i=0; i<5; i++)
    ApplyEntropyFix(lam[i], tol);

  //Roe's flux
  double fm[5], fp[5];
  if(dir==0) {
    EvaluateFluxFunction_F(Vm, id, fm);
    EvaluateFluxFunction_F(Vp, id, fp);
  } else if(dir==1) {
    EvaluateFluxFunction_G(Vm, id, fm);
    EvaluateFluxFunction_G(Vp, id, fp);
  } else { //dir = 2
    EvaluateFluxFunction_H(Vm, id, fm);
    EvaluateFluxFunction_H(Vp, id, fp);
  }

  for(int i=0; i<5; i++) {
    flux[i] = 0.5*(fm[i]+fp[i]);
    for(int p=0; p<nDOF; p++)
      flux[i] -= 0.5*fabs(lam[p])*a[p]*r[p][i];
  }
}

//----------------------------------------------------------------------------------------
// Harten's entropy fix. See LeVeque's book, Section 15.3.5.
inline 
void FluxFcnGenRoe::ApplyEntropyFix(double &lam, double tol)
{
  if(fabs(lam)<tol)
    lam = (lam*lam + tol*tol)/(2*tol);
}

//----------------------------------------------------------------------------------------


#endif
