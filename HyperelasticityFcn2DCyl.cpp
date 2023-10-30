/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<HyperelasticityFcn2DCyl.h>
#include<linear_algebra.h>
#include<cassert>

//----------------------------------------------------------------------

void
HyperelasticityFcnBase2DCyl::EvaluateHyperelasticFluxFunction_F(double* flux, double* F, double* V,
                                                                bool deviator_only)
{
  if(type == NONE) {
    for(int i=0; i<5; i++)
      flux[i] = 0.0;
    return;
  }

  double sigma[3], sigma_phiphi;
  GetCauchyStressTensor(F, V, sigma, sigma_phiphi);

  if(deviator_only) {
    double p = -1.0/3.0*(sigma[0] + sigma[2] + sigma_phiphi); //hydrostatic pressure
    sigma[0] += p;
    sigma[2] += p;
    sigma_phiphi += p;
  }

  flux[0] = 0.0;
  flux[1] = sigma[0]; //xx
  flux[2] = sigma[1]; //xy
  flux[3] = 0.0;
  flux[4] = V[1]*sigma[0] + V[2]*sigma[1];
}

//----------------------------------------------------------------------

void
HyperelasticityFcnBase2DCyl::EvaluateHyperelasticFluxFunction_G(double* flux, double* F, double* V,
                                                                bool deviator_only)
{
  if(type == NONE) {
    for(int i=0; i<5; i++)
      flux[i] = 0.0;
    return;
  }

  double sigma[3], sigma_phiphi;
  GetCauchyStressTensor(F, V, sigma, sigma_phiphi);

  if(deviator_only) {
    double p = -1.0/3.0*(sigma[0] + sigma[2] + sigma_phiphi); //hydrostatic pressure
    sigma[0] += p;
    sigma[2] += p;
    sigma_phiphi += p;
  }

  flux[0] = 0.0;
  flux[1] = sigma[1]; //xy
  flux[2] = sigma[2]; //yy
  flux[3] = 0.0;
  flux[4] = V[1]*sigma[1] + V[2]*sigma[2];
}

//----------------------------------------------------------------------

void
HyperelasticityFcnBase2DCyl::ConvertPK2ToCauchy(double* P, double *F, double J, double *sigma)
{
  //Note: This function takes P2D, F2D, J (not J2D!), and outputs sigma_2D (2x2 matrix, dim=3 by symm)
  //Note: (Danger!) Do NOT use F2x2, M2x2, N2x2 --- they could be "active" elsewhere
  double MM2x2[4], NN2x2[4];
  assert(J>0.0);
  MathTools::LinearAlgebra::CalculateCTimesMatrixA2x2(1.0/J, F, MM2x2); //MM2x2 = 1/J*F
  MathTools::LinearAlgebra::CalculateMatrixMatrixProduct2x2(MM2x2, P, NN2x2); //NN2x2 = MM2x2*P 
  MathTools::LinearAlgebra::CalculateABTranspose2x2(NN2x2, F, MM2x2); //MM2x2 = NN2x2*F'

  sigma[0] = MM2x2[0];
  sigma[1] = MM2x2[1];
  sigma[2] = MM2x2[3];
}

//----------------------------------------------------------------------

HyperelasticityFcnSaintVenantKirchhoff2DCyl::
HyperelasticityFcnSaintVenantKirchhoff2DCyl(HyperelasticityModelData &hyper, VarFcnBase &vf_)
    : HyperelasticityFcnBase2DCyl(vf_)
{
  type = SAINTVENANT_KIRCHHOFF;

  double EE = hyper.youngs_modulus;
  double nu = hyper.poissons_ratio;

  lambda = (EE*nu)/((1.0+nu)*(1.0-2.0*nu)); //first Lame constant
  mu     = EE/(2.0*(1.0+nu)); //second Lame constant (Shear modulus)

  assert(lambda>0.0);
  assert(mu>0.0);
}

//----------------------------------------------------------------------

void
HyperelasticityFcnSaintVenantKirchhoff2DCyl::
GetCauchyStressTensor(double *F, [[maybe_unused]] double *V, double *sigma, double &sigma_phiphi)
{
  // Note: F = [dz/dZ  0  dz/dR;  0  r/R  0;  dr/dZ  0  dr/dR]; //"x = z", "y = r"; also, column-first
  F2x2[0] = F[0];  F2x2[2] = F[6];
  F2x2[1] = F[2];  F2x2[3] = F[8];
  double r_over_R = F[4], r_over_R_2 = r_over_R*r_over_R;

  MathTools::LinearAlgebra::CalculateATransposeA2x2(F2x2,M2x2); //M2x2 = C2D = F2D'F2D

  double J = r_over_R*sqrt(MathTools::LinearAlgebra::CalculateDeterminant2x2(M2x2)); //J = r/R*sqrt(|C2D|)
  assert(J>0.0);

  M2x2[0] -= 1.0;
  M2x2[3] -= 1.0;
  MathTools::LinearAlgebra::CalculateCTimesMatrixA2x2(0.5, M2x2, M2x2); //M2x2=E2D=1/2(C2D-I), Green strain
  double lambda_trace = lambda*(MathTools::LinearAlgebra::CalculateMatrixTrace2x2(M2x2) +
                                0.5*(r_over_R_2 - 1.0));

  //calculates the second Piola-Kirchhoff stress tensor: P2D = lambda*tr(E)*I2D + 2*mu*E2D
  MathTools::LinearAlgebra::CalculateCTimesMatrixA2x2(2.0*mu, M2x2, M2x2);
  M2x2[0] += lambda_trace;
  M2x2[3] += lambda_trace;
  
  //convert to sigma2D (dim:3)
  ConvertPK2ToCauchy(M2x2, F2x2, J, sigma);

  sigma_phiphi = 1.0/J*r_over_R_2*(lambda_trace + mu*(r_over_R_2 - 1.0));

}

//----------------------------------------------------------------------

HyperelasticityFcnModifiedSaintVenantKirchhoff2DCyl::
HyperelasticityFcnModifiedSaintVenantKirchhoff2DCyl(HyperelasticityModelData &hyper, VarFcnBase &vf_)
    : HyperelasticityFcnBase2DCyl(vf_)
{
  type = MODIFIED_SAINTVENANT_KIRCHHOFF;

  double EE = hyper.youngs_modulus;
  double nu = hyper.poissons_ratio;

  kappa = EE/(3.0*(1.0-2.0*nu)); //bulk modulus
  mu    = EE/(2.0*(1.0+nu)); //second Lame constant (Shear modulus)

  assert(kappa>0.0);
  assert(mu>0.0);
}

//----------------------------------------------------------------------

void
HyperelasticityFcnModifiedSaintVenantKirchhoff2DCyl::
GetCauchyStressTensor(double *F, [[maybe_unused]] double *V, double *sigma, double &sigma_phiphi)
{
  // Note: F = [dz/dZ  0  dz/dR;  0  r/R  0;  dr/dZ  0  dr/dR]; //"x = z", "y = r"; also, column-first
  F2x2[0] = F[0];  F2x2[2] = F[6];
  F2x2[1] = F[2];  F2x2[3] = F[8];
  double r_over_R = F[4], r_over_R_2 = r_over_R*r_over_R;

  MathTools::LinearAlgebra::CalculateATransposeA2x2(F2x2,N2x2); //N = C2D = F2D'F2D
  for(int i=0; i<4; i++)
    M2x2[i] = N2x2[i];  //M = N = C2D
  M2x2[0] -= 1.0;
  M2x2[3] -= 1.0;     //M = C2D - I2D
  MathTools::LinearAlgebra::CalculateCTimesMatrixA2x2(0.5, M2x2, M2x2); //M = E2D 

  //calculates Cinv
  double Cinv[4], J;
  MathTools::LinearAlgebra::CalculateMatrixInverseAndDeterminant2x2(N2x2, Cinv, &J);
  J = r_over_R*sqrt(J);
  assert(J>0.0);

  //calculates the second Piola-Kirchhoff stress tensor
  MathTools::LinearAlgebra::CalculateMatrixC1APlusC2B2x2(kappa*log(J), Cinv, 2.0*mu, M2x2, N2x2);
  
  //convert to sigma_2D (dim:3)
  ConvertPK2ToCauchy(N2x2, F2x2, J, sigma);

  //get sigma_phiphi
  sigma_phiphi = 1.0/J*(kappa*log(J) + mu*r_over_R_2*(r_over_R_2 - 1));
}

//----------------------------------------------------------------------

HyperelasticityFcnNeoHookean2DCyl::
HyperelasticityFcnNeoHookean2DCyl(HyperelasticityModelData &hyper, VarFcnBase &vf_)
    : HyperelasticityFcnBase2DCyl(vf_)
{
  type = NEO_HOOKEAN;

  double EE = hyper.youngs_modulus;
  double nu = hyper.poissons_ratio;

  kappa = EE/(3.0*(1.0-2.0*nu)); //bulk modulus
  mu    = EE/(2.0*(1.0+nu)); //second Lame constant (Shear modulus)

  assert(kappa>0.0);
  assert(mu>0.0);
}

//----------------------------------------------------------------------

void
HyperelasticityFcnNeoHookean2DCyl::
GetCauchyStressTensor(double *F, [[maybe_unused]] double *V, double *sigma, double &sigma_phiphi)
{
  // Note: F = [dz/dZ  0  dz/dR;  0  r/R  0;  dr/dZ  0  dr/dR]; //"x = z", "y = r"; also, column-first
  F2x2[0] = F[0];  F2x2[2] = F[6];
  F2x2[1] = F[2];  F2x2[3] = F[8];
  double r_over_R = F[4], r_over_R_2 = r_over_R*r_over_R;

  MathTools::LinearAlgebra::CalculateAATranspose2x2(F2x2,N2x2); //N = B2D = F2D*F2D': left Cauchy-Green

  // calculate I1, I3 (principal invariants)
  double I1 = MathTools::LinearAlgebra::CalculateFirstPrincipalInvariant2x2(N2x2) + r_over_R_2;
  double I3 = r_over_R_2*MathTools::LinearAlgebra::CalculateSecondPrincipalInvariant2x2(N2x2);
  assert(I1>0.0);
  assert(I3>0.0);

  double J  = sqrt(I3);
  double Jf = pow(J, -2.0/3.0);
  I1 *= Jf;

  double factor1 = mu/J*Jf;
  double factor2 = kappa*(J-1.0) - 1.0/(3.0*J)*(mu*I1);

  MathTools::LinearAlgebra::CalculateCTimesMatrixA2x2(factor1, N2x2, N2x2);

  N2x2[0] += factor2;
  N2x2[3] += factor2;

  sigma[0] = N2x2[0];
  sigma[1] = N2x2[1];
  sigma[2] = N2x2[3];

  sigma_phiphi = factor1*r_over_R_2 + factor2;
}

//----------------------------------------------------------------------

HyperelasticityFcnMooneyRivlin2DCyl::
HyperelasticityFcnMooneyRivlin2DCyl(HyperelasticityModelData &hyper, VarFcnBase &vf_)
    : HyperelasticityFcnBase2DCyl(vf_)
{
  type = MOONEY_RIVLIN;

  double EE = hyper.youngs_modulus;
  double nu = hyper.poissons_ratio;

  C01 = hyper.C01;
  kappa = EE/(3.0*(1.0-2.0*nu)); //bulk modulus
  double mu = EE/(2.0*(1.0+nu)); //second Lame constant (Shear modulus)
  C10 = 0.5*mu - C01;

  assert(kappa>0.0);
  assert(mu>0.0);
}

//----------------------------------------------------------------------

void
HyperelasticityFcnMooneyRivlin2DCyl::
GetCauchyStressTensor(double *F, [[maybe_unused]] double *V, double *sigma, double &sigma_phiphi)
{
  // Note: F = [dz/dZ  0  dz/dR;  0  r/R  0;  dr/dZ  0  dr/dR]; //"x = z", "y = r"; also, column-first
  F2x2[0] = F[0];  F2x2[2] = F[6];
  F2x2[1] = F[2];  F2x2[3] = F[8];
  double r_over_R = F[4], r_over_R_2 = r_over_R*r_over_R;

  MathTools::LinearAlgebra::CalculateAATranspose2x2(F2x2,N2x2); //N = B2D = F2D*F2D': left Cauchy-Green
  MathTools::LinearAlgebra::CalculateMatrixMatrixProduct2x2(N2x2,N2x2,M2x2); //M = B2D*B2D

  // calculate I1, I2, I3 (principal invariants)
  double N2x2_trc = MathTools::LinearAlgebra::CalculateMatrixTrace2x2(N2x2);
  double N2x2_det = MathTools::LinearAlgebra::CalculateDeterminant2x2(N2x2);
  double I1 = N2x2_trc + r_over_R_2;
  double I2 = N2x2_det + r_over_R_2*N2x2_trc;
  double I3 = r_over_R_2*N2x2_det;
  assert(I1>0.0);
  assert(I2>0.0);
  assert(I3>0.0);

  double J  = sqrt(I3);
  double Jf = pow(J, -2.0/3.0);
  I1 *= Jf;
  I2 *= Jf*Jf;

  double factor1 = 2.0/J*Jf*(C10+I1*C01);
  double factor2 = kappa*(J-1.0) - 2.0/(3.0*J)*(C10*I1 + 2.0*C01*I2);
  double factor3 = -2.0/J*Jf*Jf*C01;

  MathTools::LinearAlgebra::CalculateMatrixC1APlusC2B2x2(factor1, N2x2, factor3, M2x2, N2x2);

  N2x2[0] += factor2;
  N2x2[3] += factor2;

  sigma[0] = N2x2[0];
  sigma[1] = N2x2[1];
  sigma[2] = N2x2[3];

  sigma_phiphi = factor1*r_over_R_2 + factor3*r_over_R_2*r_over_R_2 + factor2;
}


//----------------------------------------------------------------------





