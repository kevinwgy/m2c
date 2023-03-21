/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _FLUX_FCN_BASE_H_
#define _FLUX_FCN_BASE_H_

#include <VarFcnBase.h>
#include <vector>

/*****************************************************************************************
 * Base class for calculating fluxes within the interiors of all the material subdomains, 
 * Fluxes across material interfaces are handled ELSEWHERE.
 *****************************************************************************************/

class FluxFcnBase {

protected:
  std::vector<VarFcnBase*> &vf;

public:
  FluxFcnBase(std::vector<VarFcnBase*> &varFcn) : vf(varFcn) { }
  virtual ~FluxFcnBase() { }

  inline void EvaluateFluxFunction_F(double *V, int id, double *F);
  inline void EvaluateFluxFunction_G(double *V, int id, double *G);
  inline void EvaluateFluxFunction_H(double *V, int id, double *H);

  /** Evaluate eigenvalues and eigenvectors of the Jacobian matrices, using V and the EOS. */
  inline void EvaluateEigensOfJacobian_F(double *V, int id, //input state
                                         double &lam1, double &lam2, double &lam3, double &lam4, double &lam5,
                                         double *r1 = 0, double *r2 = 0, double *r3 = 0, double *r4 = 0, double *r5 = 0);
  inline void EvaluateEigensOfJacobian_G(double *V, int id, //input state
                                         double &lam1, double &lam2, double &lam3, double &lam4, double &lam5,
                                         double *r1 = 0, double *r2 = 0, double *r3 = 0, double *r4 = 0, double *r5 = 0);
  inline void EvaluateEigensOfJacobian_H(double *V, int id, //input state
                                         double &lam1, double &lam2, double &lam3, double &lam4, double &lam5,
                                         double *r1 = 0, double *r2 = 0, double *r3 = 0, double *r4 = 0, double *r5 = 0);
  
  /** Evaluate eigenvalues and eigenvectors of the Jacobian matrices, with all the variables specified explicitly.
    * The input variables can, but do not have to, satisfy the EOS. This is useful for evaluating some numerical flux functions which
    * simplify borrows the form of the Jacobian matrix */
  inline void EvaluateEigensOfJacobian_F(double u, double v, double w, double e, double c, double H, double Gamma, double dpdrho, 
                                         double &lam1, double &lam2, double &lam3, double &lam4, double &lam5,
                                         double *r1 = 0, double *r2 = 0, double *r3 = 0, double *r4 = 0, double *r5 = 0);
  inline void EvaluateEigensOfJacobian_G(double u, double v, double w, double e, double c, double H, double Gamma, double dpdrho, 
                                         double &lam1, double &lam2, double &lam3, double &lam4, double &lam5,
                                         double *r1 = 0, double *r2 = 0, double *r3 = 0, double *r4 = 0, double *r5 = 0);
  inline void EvaluateEigensOfJacobian_H(double u, double v, double w, double e, double c, double H, double Gamma, double dpdrho, 
                                         double &lam1, double &lam2, double &lam3, double &lam4, double &lam5,
                                         double *r1 = 0, double *r2 = 0, double *r3 = 0, double *r4 = 0, double *r5 = 0);

  inline void EvaluateMaxEigenvalues(double *V, int id, double &lam_f_max, double &lam_g_max, double &lam_h_max);

  /** Functions related to characteristic variables (see KW's notes)*/
  inline void PrimitiveToPrimitiveCharacteristic(int dir/*0~x/F,1~y/G,2~z/H*/, double *V0, double *dV, int id, double *dW);
  inline void PrimitiveCharacteristicToPrimitive(int dir/*0~x/F,1~y/G,2~z/H*/, double *V0, double *dW, int id, double *dV);
  inline void ConservativeToConservativeCharacteristic(int dir/*0~x/F,1~y/G,2~z/H*/, double *V0, double *dU, int id, double *dW);
  inline void ConservativeCharacteristicToConservative(int dir/*0~x/F,1~y/G,2~z/H*/, double *V0, double *dW, int id, double *dU);


  /** The following function(s) depend on the numerical method for flux calculation.
    * Should be defined in derived classes. */
  virtual void ComputeNumericalFluxAtCellInterface([[maybe_unused]] int dir/*0~x,1~y,2~z*/, [[maybe_unused]] double *Vminus/*left*/, 
                                                   [[maybe_unused]] double *Vplus/*right*/, [[maybe_unused]] int id, 
                                                   [[maybe_unused]] double *F) {
    print_error("*** Error: ComputeNumericalFluxAtCellInterface function not defined.\n"); 
    exit_mpi();}

  //! Computes flux across material interface. Currently, only available in the case of Local Lax-Friedrichs (LLF)
  virtual void ComputeNumericalFluxAtMaterialInterface([[maybe_unused]] int dir/*0~x,1~y,2~z*/, [[maybe_unused]] double *Vminus/*left*/,
                                                       [[maybe_unused]] int idm, [[maybe_unused]] double *Vplus/*right*/,
                                                       [[maybe_unused]] int idp, [[maybe_unused]] double *F) {
    print_error("*** Error: ComputeNumericalFluxAtMaterialInterface function not defined.\n"); 
    exit_mpi();}

};

//----------------------------------------------------------------------------------------

inline 
void FluxFcnBase::EvaluateFluxFunction_F(double *V, int id, double *F)
{
  double H = vf[id]->ComputeTotalEnthalpyPerUnitMass(V); //H = 1/rho*(E+p) = e + 1/2||u||^2 + p/rho
  double rhou = V[0]*V[1];

  F[0] = rhou; //rho*u 
  F[1] = rhou*V[1] + V[4]; //rho*u*u + p
  F[2] = rhou*V[2]; //rho*u*v
  F[3] = rhou*V[3]; //rho*u*w
  F[4] = rhou*H; //rho*H*u 
}

//------------------------------------------------------------------------------

inline 
void FluxFcnBase::EvaluateFluxFunction_G(double *V, int id, double *G)
{
  double H = vf[id]->ComputeTotalEnthalpyPerUnitMass(V); 
  double rhov = V[0]*V[2];

  G[0] = rhov; //rho*v
  G[1] = rhov*V[1]; //rho*v*u
  G[2] = rhov*V[2] + V[4]; //rho*v*v+p
  G[3] = rhov*V[3]; //rho*v*w
  G[4] = rhov*H; //rho*H*v
}

//------------------------------------------------------------------------------

inline
void FluxFcnBase::EvaluateFluxFunction_H(double *V, int id, double *HH)
{
  double H = vf[id]->ComputeTotalEnthalpyPerUnitMass(V); 
  double rhow = V[0]*V[3];

  HH[0] = rhow; //rho*w
  HH[1] = rhow*V[1]; //rho*w*u
  HH[2] = rhow*V[2]; //rho*w*v
  HH[3] = rhow*V[3] + V[4]; //rho*w*w + p
  HH[4] = rhow*H; //rho*H*w
}

//------------------------------------------------------------------------------

inline 
void FluxFcnBase::EvaluateEigensOfJacobian_F(double *V, int id,
                                double &lam1, double &lam2, double &lam3, double &lam4, double &lam5,
                                double *r1, double *r2, double *r3, double *r4, double *r5)
{
  double u = V[1];
  double v = V[2];
  double w = V[3];
  double e = vf[id]->GetInternalEnergyPerUnitMass(V[0],V[4]);
  double c = vf[id]->ComputeSoundSpeedSquare(V[0], e);

  if(c<0) {
    fprintf(stdout,"*** Error: c^2 (square of sound speed) = %e in EvaluateEigensOfJacobian_F. V = %e, %e, %e, %e, %e, ID = %d.\n",
            c, V[0], V[1], V[2], V[3], V[4], id);            
    exit_mpi();
  } else
    c = sqrt(c);

  double H = vf[id]->ComputeTotalEnthalpyPerUnitMass(V); 
  double Gamma = vf[id]->GetBigGamma(V[0], e);
  double dpdrho = vf[id]->GetDpdrho(V[0], e);

  EvaluateEigensOfJacobian_F(u,v,w,e,c,H,Gamma,dpdrho,lam1,lam2,lam3,lam4,lam5,r1,r2,r3,r4,r5);
}

//------------------------------------------------------------------------------

inline 
void FluxFcnBase::EvaluateEigensOfJacobian_F(double u, double v, double w, double e, double c, double H, double Gamma, double dpdrho,
                      double &lam1, double &lam2, double &lam3, double &lam4, double &lam5,
                      double *r1, double *r2, double *r3, double *r4, double *r5)
{
  double vel2 = u*u + v*v + w*w;

  //eigenvalues
  lam1 = u - c;
  lam2 = lam3 = lam4 = u;
  lam5 = u + c;

  if(!r1 || !r2 || !r3 || !r4 || !r5) //no need to output eigenvectors
    return;

  //eigenvectores
  r1[0] = 1;        r2[0] = 1;                        r3[0] = 0;  r4[0] = 0;  r5[0] = 1;
  r1[1] = u - c;    r2[1] = u;                        r3[1] = 0;  r4[1] = 0;  r5[1] = u + c;
  r1[2] = v;        r2[2] = v;                        r3[2] = 1;  r4[2] = 0;  r5[2] = v;
  r1[3] = w;        r2[3] = w;                        r3[3] = 0;  r4[3] = 1;  r5[3] = w;
  r1[4] = H - u*c;  r2[4] = e+0.5*vel2-dpdrho/Gamma;  r3[4] = v;  r4[4] = w;  r5[4] = H + u*c;

}

//------------------------------------------------------------------------------

inline
void FluxFcnBase::EvaluateEigensOfJacobian_G(double *V, int id,
                      double &lam1, double &lam2, double &lam3, double &lam4, double &lam5,
                      double *r1, double *r2, double *r3, double *r4, double *r5)
{
  double u = V[1];
  double v = V[2];
  double w = V[3];
  double e = vf[id]->GetInternalEnergyPerUnitMass(V[0],V[4]);
  double c = vf[id]->ComputeSoundSpeedSquare(V[0], e);

  if(c<0) {
    fprintf(stdout,"*** Error: c^2 (square of sound speed) = %e in EvaluateEigensOfJacobian_G. V = %e, %e, %e, %e, %e, ID = %d.\n",
            c, V[0], V[1], V[2], V[3], V[4], id);            
    exit_mpi();
  } else
    c = sqrt(c);

  double H = vf[id]->ComputeTotalEnthalpyPerUnitMass(V); 
  double Gamma = vf[id]->GetBigGamma(V[0], e);
  double dpdrho = vf[id]->GetDpdrho(V[0], e);

  EvaluateEigensOfJacobian_G(u,v,w,e,c,H,Gamma,dpdrho,lam1,lam2,lam3,lam4,lam5,r1,r2,r3,r4,r5);
}

//------------------------------------------------------------------------------

inline
void FluxFcnBase::EvaluateEigensOfJacobian_G(double u, double v, double w, double e, double c, double H, double Gamma, double dpdrho,
                      double &lam1, double &lam2, double &lam3, double &lam4, double &lam5,
                      double *r1, double *r2, double *r3, double *r4, double *r5)
{
  double vel2 = u*u + v*v + w*w;

  //eigenvalues
  lam1 = v - c;
  lam2 = lam3 = lam4 = v;
  lam5 = v + c;

  if(!r1 || !r2 || !r3 || !r4 || !r5) //no need to output eigenvectors
    return;

  //eigenvectores
  r1[0] = 1;        r2[0] = 0;  r3[0] = 1;                        r4[0] = 0;  r5[0] = 1;
  r1[1] = u;        r2[1] = 1;  r3[1] = u;                        r4[1] = 0;  r5[1] = u;
  r1[2] = v - c;    r2[2] = 0;  r3[2] = v;                        r4[2] = 0;  r5[2] = v + c;
  r1[3] = w;        r2[3] = 0;  r3[3] = w;                        r4[3] = 1;  r5[3] = w;
  r1[4] = H - v*c;  r2[4] = u;  r3[4] = e+0.5*vel2-dpdrho/Gamma;  r4[4] = w;  r5[4] = H + v*c;

}

//------------------------------------------------------------------------------

inline
void FluxFcnBase::EvaluateEigensOfJacobian_H(double *V, int id,
                      double &lam1, double &lam2, double &lam3, double &lam4, double &lam5,
                      double *r1, double *r2, double *r3, double *r4, double *r5)
{
  double u = V[1];
  double v = V[2];
  double w = V[3];
  double e = vf[id]->GetInternalEnergyPerUnitMass(V[0],V[4]);
  double c = vf[id]->ComputeSoundSpeedSquare(V[0], e);

  if(c<0) {
    fprintf(stdout,"*** Error: c^2 (square of sound speed) = %e in EvaluateEigensOfJacobian_H. V = %e, %e, %e, %e, %e, ID = %d.\n",
            c, V[0], V[1], V[2], V[3], V[4], id);            
    exit_mpi();
  } else
    c = sqrt(c);

  double H = vf[id]->ComputeTotalEnthalpyPerUnitMass(V); 
  double Gamma = vf[id]->GetBigGamma(V[0], e);
  double dpdrho = vf[id]->GetDpdrho(V[0], e);

  EvaluateEigensOfJacobian_H(u,v,w,e,c,H,Gamma,dpdrho,lam1,lam2,lam3,lam4,lam5,r1,r2,r3,r4,r5);
}

//------------------------------------------------------------------------------

inline
void FluxFcnBase::EvaluateEigensOfJacobian_H(double u, double v, double w, double e, double c, double H, double Gamma, double dpdrho,
                      double &lam1, double &lam2, double &lam3, double &lam4, double &lam5,
                      double *r1, double *r2, double *r3, double *r4, double *r5)
{
  double vel2 = u*u + v*v + w*w;

  //eigenvalues
  lam1 = w - c;
  lam2 = lam3 = lam4 = w;
  lam5 = w + c;

  if(!r1 || !r2 || !r3 || !r4 || !r5) //no need to output eigenvectors
    return;

  //eigenvectores
  r1[0] = 1;        r2[0] = 0;  r3[0] = 0;  r4[0] = 1;                        r5[0] = 1;
  r1[1] = u;        r2[1] = 1;  r3[1] = 0;  r4[1] = u;                        r5[1] = u;
  r1[2] = v;        r2[2] = 0;  r3[2] = 1;  r4[2] = v;                        r5[2] = v;
  r1[3] = w - c;    r2[3] = 0;  r3[3] = 0;  r4[3] = w;                        r5[3] = w + c;
  r1[4] = H - w*c;  r2[4] = u;  r3[4] = v;  r4[4] = e+0.5*vel2-dpdrho/Gamma;  r5[4] = H + w*c;

}

//------------------------------------------------------------------------------

inline 
void FluxFcnBase::EvaluateMaxEigenvalues(double *V, int id, double &lam_f_max, double &lam_g_max, double &lam_h_max)
{
  double e = vf[id]->GetInternalEnergyPerUnitMass(V[0],V[4]);
  double c = vf[id]->ComputeSoundSpeedSquare(V[0], e);

  if(c<0) {
    fprintf(stdout,"*** Error: c^2 (square of sound speed) = %e in EvaluateMaxEigenvalues. V = %e, %e, %e, %e, %e, ID = %d.\n",
            c, V[0], V[1], V[2], V[3], V[4], id);            
    exit_mpi();
  } else
    c = sqrt(c);


  lam_f_max = std::fabs(V[1]) + c;
  lam_g_max = std::fabs(V[2]) + c;
  lam_h_max = std::fabs(V[3]) + c;
}

//------------------------------------------------------------------------------
// Given a ref state (V0) and a difference in primitive state variables (dV), compute
// the difference in characterstic varaibles (dW). See KW's notes.
inline
void FluxFcnBase::PrimitiveToPrimitiveCharacteristic(int dir/*0~x/F,1~y/G,2~z/H*/, double *V0, double *dV, 
                                                     int id, double *dW)
{
  //Step 1. Multiply dV by Q0 (dU/dV evaluated at V0)
  double e0 = vf[id]->GetInternalEnergyPerUnitMass(V0[0],V0[4]);
  double c0 = vf[id]->ComputeSoundSpeedSquare(V0[0], e0);
  if(c0<0) {
    fprintf(stdout,"*** Error: c0^2 (square of sound speed) = %e in PrimitiveToPrimitiveCharacteristic. "
                   "V0 = %e, %e, %e, %e, %e, ID = %d.\n", c0, V0[0], V0[1], V0[2], V0[3], V0[4], id);            
    exit_mpi();
  } else
    c0 = sqrt(c0);

  double dpdrho0 = vf[id]->GetDpdrho(V0[0], e0);
  double Gamma0 = vf[id]->GetBigGamma(V0[0], e0);
  double kin0 = 0.5*(V0[1]*V0[1]+V0[2]*V0[2]+V0[3]*V0[3]);
  double Q0_51 = e0 + kin0 - dpdrho0/Gamma0;
  double Tmp[5];
  Tmp[0] = dV[0];
  Tmp[1] = V0[1]*dV[0] + V0[0]*dV[1];
  Tmp[2] = V0[2]*dV[0] + V0[0]*dV[2];
  Tmp[3] = V0[3]*dV[0] + V0[0]*dV[3];
  Tmp[4] = Q0_51*dV[0] + V0[0]*V0[1]*dV[1] + V0[0]*V0[2]*dV[2] + V0[0]*V0[3]*dV[3] + 1.0/Gamma0*dV[4];
  
  //Step 2. Multiply by R0^{-1}
  double un0 = V0[dir+1]; //normal vel.
  double n[3] = {0,0,0}; n[dir] = 1.0; //normal direction 
  double R0inv[5][5];
  double alpha0 = V0[0]*Gamma0/(V0[4]*Gamma0 + V0[0]*dpdrho0);
  double beta0 = (e0 - kin0 + V0[4]/V0[0])*alpha0; 

  R0inv[0][0] = 0.5*(1.0 - beta0 + un0/c0);
  for(int i=0; i<3; i++)
    R0inv[0][1+i] = -0.5*(alpha0*V0[1+i] + n[i]/c0);
  R0inv[0][4] = alpha0/2.0;
  for(int i=0; i<3; i++)
    R0inv[1+i][0] = -V0[1+i] + (beta0+un0)*n[i];
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      R0inv[1+i][1+j] = (i==j) + alpha0*n[i]*V0[1+j] - n[i]*n[j];
  for(int i=0; i<3; i++)
    R0inv[1+i][4] = -alpha0*n[i];
  R0inv[4][0] = 0.5*(1.0 - beta0 - un0/c0);
  for(int i=0; i<3; i++)
    R0inv[4][1+i] = -0.5*(alpha0*V0[1+i] - n[i]/c0);
  R0inv[4][4] = alpha0/2.0;  

  for(int i=0; i<5; i++) {
    dW[i] = 0.0;
    for(int j=0; j<5; j++)
      dW[i] += R0inv[i][j]*Tmp[j];
  } 

}

//------------------------------------------------------------------------------
// The reverse operation of the previous function (primitive -> primitive characteristic)
inline
void FluxFcnBase::PrimitiveCharacteristicToPrimitive(int dir/*0~x/F,1~y/G,2~z/H*/, double *V0, double *dW, 
                                                     int id, double *dV)
{
  //Step 1. Multiply dW by R0
  double un0 = V0[dir+1]; //normal vel.
  double n[3] = {0,0,0}; n[dir] = 1.0; //normal direction 
  double e0 = vf[id]->GetInternalEnergyPerUnitMass(V0[0],V0[4]);
  double c0 = vf[id]->ComputeSoundSpeedSquare(V0[0], e0);
  if(c0<0) {
    fprintf(stdout,"*** Error: c0^2 (square of sound speed) = %e in PrimitiveCharacteristicToPrimitive. "
                   "V0 = %e, %e, %e, %e, %e, ID = %d.\n", c0, V0[0], V0[1], V0[2], V0[3], V0[4], id);            
    exit_mpi();
  } else
    c0 = sqrt(c0);

  double dpdrho0 = vf[id]->GetDpdrho(V0[0], e0);
  double Gamma0 = vf[id]->GetBigGamma(V0[0], e0);
  double kin0 = 0.5*(V0[1]*V0[1]+V0[2]*V0[2]+V0[3]*V0[3]);
  double H0 = e0 + kin0 + V0[4]/V0[0];

  double R0[5][5];
  R0[0][0] = 1.0;
  for(int i=0; i<3; i++)
    R0[0][1+i] = n[i];
  R0[0][4] = 1.0;
  for(int i=0; i<3; i++)
    R0[1+i][0] = V0[1+i] - c0*n[i];
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      R0[1+i][1+j] = (i==j) + V0[1+i]*n[j] - n[i]*n[j];
  for(int i=0; i<3; i++)
    R0[1+i][4] = V0[1+i] + c0*n[i];
  R0[4][0] = H0 - un0*c0;
  double param = e0 + kin0 - dpdrho0/Gamma0 - un0;
  for(int i=0; i<3; i++)
    R0[4][1+i] = V0[1+i] + param*n[i];
  R0[4][4] = H0 + un0*c0;
  
  double Tmp[5];
  for(int i=0; i<5; i++) {
    Tmp[i] = 0.0;
    for(int j=0; j<5; j++)
      Tmp[i] += R0[i][j]*dW[j];
  } 

  //Step 2. Multiply by Qbar^{-1}
  dV[0] = Tmp[0];
  dV[1] = -V0[1]/V0[0]*Tmp[0] + Tmp[1]/V0[0];
  dV[2] = -V0[2]/V0[0]*Tmp[0] + Tmp[2]/V0[0];
  dV[3] = -V0[3]/V0[0]*Tmp[0] + Tmp[3]/V0[0];
  dV[4] = ((-e0 + kin0)*Gamma0 + dpdrho0)*Tmp[0] 
        - Gamma0*(V0[1]*Tmp[1]+V0[2]*Tmp[2]+V0[3]*Tmp[3] - Tmp[4]);

}

//------------------------------------------------------------------------------

inline
void FluxFcnBase::ConservativeToConservativeCharacteristic(int dir/*0~x/F,1~y/G,2~z/H*/, double *V0, double *dU,
                                                           int id, double *dW)
{
  //Multiply dU by R0^{-1}
  double e0 = vf[id]->GetInternalEnergyPerUnitMass(V0[0],V0[4]);
  double c0 = vf[id]->ComputeSoundSpeedSquare(V0[0], e0);
  if(c0<0) {
    fprintf(stdout,"*** Error: c0^2 (square of sound speed) = %e in ConservativeToConservativeCharacteristic. "
                   "V0 = %e, %e, %e, %e, %e, ID = %d.\n", c0, V0[0], V0[1], V0[2], V0[3], V0[4], id);            
    exit_mpi();
  } else
    c0 = sqrt(c0);

  double dpdrho0 = vf[id]->GetDpdrho(V0[0], e0);
  double Gamma0 = vf[id]->GetBigGamma(V0[0], e0);
  double kin0 = 0.5*(V0[1]*V0[1]+V0[2]*V0[2]+V0[3]*V0[3]);
  double un0 = V0[dir+1]; //normal vel.
  double n[3] = {0,0,0}; n[dir] = 1.0; //normal direction 
  double alpha0 = V0[0]*Gamma0/(V0[4]*Gamma0 + V0[0]*dpdrho0);
  double beta0 = (e0 - kin0 + V0[4]/V0[0])*alpha0; 

  double R0inv[5][5];
  R0inv[0][0] = 0.5*(1.0 - beta0 + un0/c0);
  for(int i=0; i<3; i++)
    R0inv[0][1+i] = -0.5*(alpha0*V0[1+i] + n[i]/c0);
  R0inv[0][4] = alpha0/2.0;
  for(int i=0; i<3; i++)
    R0inv[1+i][0] = -V0[1+i] + (beta0+un0)*n[i];
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      R0inv[1+i][1+j] = (i==j) + alpha0*n[i]*V0[1+j] - n[i]*n[j];
  for(int i=0; i<3; i++)
    R0inv[1+i][4] = -alpha0*n[i];
  R0inv[4][0] = 0.5*(1.0 - beta0 - un0/c0);
  for(int i=0; i<3; i++)
    R0inv[4][1+i] = -0.5*(alpha0*V0[1+i] - n[i]/c0);
  R0inv[4][4] = alpha0/2.0;  

  for(int i=0; i<5; i++) {
    dW[i] = 0.0;
    for(int j=0; j<5; j++)
      dW[i] += R0inv[i][j]*dU[j];
  } 
}

//------------------------------------------------------------------------------

inline
void FluxFcnBase::ConservativeCharacteristicToConservative(int dir/*0~x/F,1~y/G,2~z/H*/, double *V0, double *dW,
                                                           int id, double *dU)
{
  //Multiply dW by R0
  double un0 = V0[dir+1]; //normal vel.
  double n[3] = {0,0,0}; n[dir] = 1.0; //normal direction 
  double e0 = vf[id]->GetInternalEnergyPerUnitMass(V0[0],V0[4]);
  double c0 = vf[id]->ComputeSoundSpeedSquare(V0[0], e0);
  if(c0<0) {
    fprintf(stdout,"*** Error: c0^2 (square of sound speed) = %e in ConservativeCharacteristicToConservative. "
                   "V0 = %e, %e, %e, %e, %e, ID = %d.\n", c0, V0[0], V0[1], V0[2], V0[3], V0[4], id);            
    exit_mpi();
  } else
    c0 = sqrt(c0);

  double dpdrho0 = vf[id]->GetDpdrho(V0[0], e0);
  double Gamma0 = vf[id]->GetBigGamma(V0[0], e0);
  double kin0 = 0.5*(V0[1]*V0[1]+V0[2]*V0[2]+V0[3]*V0[3]);
  double H0 = e0 + kin0 + V0[4]/V0[0];

  double R0[5][5];
  R0[0][0] = 1.0;
  for(int i=0; i<3; i++)
    R0[0][1+i] = n[i];
  R0[0][4] = 1.0;
  for(int i=0; i<3; i++)
    R0[1+i][0] = V0[1+i] - c0*n[i];
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      R0[1+i][1+j] = (i==j) + V0[1+i]*n[j] - n[i]*n[j];
  for(int i=0; i<3; i++)
    R0[1+i][4] = V0[1+i] + c0*n[i];
  R0[4][0] = H0 - un0*c0;
  double param = e0 + kin0 - dpdrho0/Gamma0 - un0;
  for(int i=0; i<3; i++)
    R0[4][1+i] = V0[1+i] + param*n[i];
  R0[4][4] = H0 + un0*c0;
  
  for(int i=0; i<5; i++) {
    dU[i] = 0.0;
    for(int j=0; j<5; j++)
      dU[i] += R0[i][j]*dW[j];
  } 

}

//------------------------------------------------------------------------------


#endif
