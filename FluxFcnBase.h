#ifndef _FLUX_FCN_BASE_H_
#define _FLUX_FCN_BASE_H_

#include <VarFcnBase.h>

//----------------------------------------------------------------------------------------

class FluxFcnBase {

protected:
  VarFcnBase *vf;

public:
  FluxFcnBase(IoData VarFcnBase *varFcn);
  virtual ~FluxFcnBase() { vf = 0; }

  inline void EvaluateFluxFunction_F(double *V, double *F);
  inline void EvaluateFluxFunction_G(double *V, double *G);
  inline void EvaluateFluxFunction_H(double *V, double *H);

  /** Evaluate eigenvalues and eigenvectors of the Jacobian matrices, using V and the EOS. */
  inline void EvaluateEigensOfJacobian_F(double *V, //input state
                                         double &lam1, double &lam2, double &lam3, double &lam4, double &lam5,
                                         double *r1 = 0, double *r2 = 0, double *r3 = 0, double *r4 = 0, double *r5 = 0);
  inline void EvaluateEigensOfJacobian_G(double *V, //input state
                                         double &lam1, double &lam2, double &lam3, double &lam4, double &lam5,
                                         double *r1 = 0, double *r2 = 0, double *r3 = 0, double *r4 = 0, double *r5 = 0);
  inline void EvaluateEigensOfJacobian_H(double *V, //input state
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

  inline void EvaluateMaxEigenvalues(double *V, double &lam_f_max, &lam_g_max, &lam_h_max);

  /** The following function(s) depend on the numerical method for flux calculation.
    * Should be defined in derived classes. */
  virtual void ComputeNumericalFluxAtCellInterface(int dir/*0~x,1~y,2~z*/,double *Vminus/*left*/, double *Vplus/*right*/, double *F) {
    print_error("*** Error: ComputeNumericalFluxAtCellInterface function not defined.\n"); 
    exit_mpi();}

};

//----------------------------------------------------------------------------------------

inline
FluxFcnBase::FluxFcnBase(VarFcnBase *varFcn) : vf(varFcn)
{
  
}

//------------------------------------------------------------------------------

inline 
void FluxFcnBase::EvaluateFluxFunction_F(double *V, double *F)
{
  double H = vf->ComputeTotalEnthalpyPerUnitMass(V); //H = 1/rho*(E+p) = e + 1/2||u||^2 + p/rho
  double rhou = V[0]*V[1];

  F[0] = rhou; //rho*u 
  F[1] = rhou*V[1] + V[4]; //rho*u*u + p
  F[2] = rhou*V[2]; //rho*u*v
  F[3] = rhou*V[3]; //rho*u*w
  F[4] = rhou*H; //rho*H*u 
}

//------------------------------------------------------------------------------

inline 
void FluxFcnBase::EvaluateFluxFunction_G(double *V, double *G)
{
  double H = vf->ComputeTotalEnthalpyPerUnitMass(V); 
  double rhov = V[0]*V[2];

  G[0] = rhov; //rho*v
  G[1] = rhov*V[1]; //rho*v*u
  G[2] = rhov*V[2] + V[4]; //rho*v*v+p
  G[3] = rhov*V[3]; //rho*v*w
  G[4] = rhov*H; //rho*H*v
}

//------------------------------------------------------------------------------

inline
void FluxFcnBase::EvaluateFluxFunction_H(double *V, double *H)
{
  double H = vf->ComputeTotalEnthalpyPerUnitMass(V); 
  double rhow = V[0]*V[3];

  H[0] = rhow; //rho*w
  H[1] = rhow*V[1]; //rho*w*u
  H[2] = rhow*V[2]; //rho*w*v
  H[3] = rhow*V[3] + V[4]; //rho*w*w + p
  H[4] = rhow*H; //rho*H*w
}

//------------------------------------------------------------------------------

inline 
void FluxFcnBase::EvaluateEigensOfJacobian_F(double *V, 
                                double &lam1, double &lam2, double &lam3, double &lam4, double &lam5,
                                double *r1, double *r2, double *r3, double *r4, double *r5)
{
  double u = V[1];
  double v = V[2];
  double w = V[3];
  double e = vf->GetInternalEnergyPerUnitMass(V[0],V[4]);
  double c = vf->ComputeSoundSpeed(V[0], e);
  double H = vf->ComputeTotalEnthalpyPerUnitMass(V); 
  double Gamma = vf->GetBigGamma(V[0], e);
  double dpdrho = vf->GetDpdrho(V[0], e);

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
void FluxFcnBase::EvaluateEigensOfJacobian_G(double *V, 
                      double &lam1, double &lam2, double &lam3, double &lam4, double &lam5,
                      double *r1, double *r2, double *r3, double *r4, double *r5)
{
  double u = V[1];
  double v = V[2];
  double w = V[3];
  double e = vf->GetInternalEnergyPerUnitMass(V[0],V[4]);
  double c = vf->ComputeSoundSpeed(V[0], e);
  double H = vf->ComputeTotalEnthalpyPerUnitMass(V); 
  double Gamma = vf->GetBigGamma(V[0], e);
  double dpdrho = vf->GetDpdrho(V[0], e);

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
void FluxFcnBase::EvaluateEigensOfJacobian_H(double *V, 
                      double &lam1, double &lam2, double &lam3, double &lam4, double &lam5,
                      double *r1, double *r2, double *r3, double *r4, double *r5)
{
  double u = V[1];
  double v = V[2];
  double w = V[3];
  double e = vf->GetInternalEnergyPerUnitMass(V[0],V[4]);
  double c = vf->ComputeSoundSpeed(V[0], e);
  double H = vf->ComputeTotalEnthalpyPerUnitMass(V); 
  double Gamma = vf->GetBigGamma(V[0], e);
  double dpdrho = vf->GetDpdrho(V[0], e);

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
  r1[0] = 1;        r2[0] = 0;  r3[0] = 0;  r3[0] = 1;                        r5[0] = 1;
  r1[1] = u;        r2[1] = 1;  r3[1] = 0;  r3[1] = u;                        r5[1] = u;
  r1[2] = v;        r2[2] = 0;  r3[2] = 1;  r3[2] = v;                        r5[2] = v;
  r1[3] = w - c;    r2[3] = 0;  r3[3] = 0;  r3[3] = w;                        r5[3] = w + c;
  r1[4] = H - w*c;  r2[4] = u;  r3[4] = v;  r3[4] = e+0.5*vel2-dpdrho/Gamma;  r5[4] = H + w*c;

}

//------------------------------------------------------------------------------

inline 
void FluxFcnBase::EvaluateMaxEigenvalues(double *V, double &lam_f_max, &lam_g_max, &lam_h_max)
{
  double e = vf->GetInternalEnergyPerUnitMass(V[0],V[4]);
  double c = vf->ComputeSoundSpeed(V[0], e);

  lam_f_max = std::fabs(V[1]) + c;
  lam_g_max = std::fabs(V[2]) + c;
  lam_h_max = std::fabs(V[3]) + c;
}

//------------------------------------------------------------------------------

#endif
