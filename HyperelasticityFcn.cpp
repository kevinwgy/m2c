/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 *j**********************************************************************/

#include<HyperelasticityFcn.h>
#include<linear_algebra.h>
#include<cassert>

//----------------------------------------------------------------------

void
HyperelasticityFcnBase::EvaluateHyperelasticFluxFunction_F(double* flux, double* F, double* V,
                                                           bool deviator_only)
{
  if(type == NONE) {
    for(int i=0; i<5; i++)
      flux[i] = 0.0;
    return;
  }

  double sigma[6];
  GetCauchyStressTensor(F, V, sigma);

  if(deviator_only) {
    double p = -1.0/3.0*(sigma[0] + sigma[3] + sigma[5]); //hydrostatic pressure
    sigma[0] += p;
    sigma[3] += p;
    sigma[5] += p;
  }

  flux[0] = 0.0;
  flux[1] = sigma[0]; //xx
  flux[2] = sigma[1]; //xy
  flux[3] = sigma[2]; //xz
  flux[4] = V[1]*sigma[0] + V[2]*sigma[1] + V[3]*sigma[2];
}

//----------------------------------------------------------------------

void
HyperelasticityFcnBase::EvaluateHyperelasticFluxFunction_G(double* flux, double* F, double* V,
                                                           bool deviator_only)
{
  if(type == NONE) {
    for(int i=0; i<5; i++)
      flux[i] = 0.0;
    return;
  }

  double sigma[6];
  GetCauchyStressTensor(F, V, sigma);

  if(deviator_only) {
    double p = -1.0/3.0*(sigma[0] + sigma[3] + sigma[5]); //hydrostatic pressure
    sigma[0] += p;
    sigma[3] += p;
    sigma[5] += p;
  }

  flux[0] = 0.0;
  flux[1] = sigma[1]; //xy
  flux[2] = sigma[3]; //yy
  flux[3] = sigma[4]; //yz
  flux[4] = V[1]*sigma[1] + V[2]*sigma[3] + V[3]*sigma[4];
}

//----------------------------------------------------------------------

void
HyperelasticityFcnBase::EvaluateHyperelasticFluxFunction_H(double* flux, double* F, double* V,
                                                           bool deviator_only)
{
  if(type == NONE) {
    for(int i=0; i<5; i++)
      flux[i] = 0.0;
    return;
  }

  double sigma[6];
  GetCauchyStressTensor(F, V, sigma);

  if(deviator_only) {
    double p = -1.0/3.0*(sigma[0] + sigma[3] + sigma[5]); //hydrostatic pressure
    sigma[0] += p;
    sigma[3] += p;
    sigma[5] += p;
  }

  flux[0] = 0.0;
  flux[1] = sigma[2]; //xz
  flux[2] = sigma[4]; //yz
  flux[3] = sigma[5]; //zz
  flux[4] = V[1]*sigma[2] + V[2]*sigma[4] + V[3]*sigma[5];
}


//----------------------------------------------------------------------

void
HyperelasticityFcnBase::ConvertPK2ToCauchy(double* P, double *F, double J, double *sigma)
{
  //Note: (Danger!) Should not use M3x3 and N3x3 in this function. They could be active elsewhere!
  
  assert(J>0.0);
  MathTools::LinearAlgebra::CalculateCTimesMatrixA3x3(1.0/J, F, MM3x3); //M = 1/J*F
  MathTools::LinearAlgebra::CalculateMatrixMatrixProduct3x3(MM3x3, P, NN3x3); //N = M*P 
  MathTools::LinearAlgebra::CalculateABTranspose3x3(NN3x3, F, MM3x3); //M = N*F'

  sigma[0] = MM3x3[0];
  sigma[1] = MM3x3[1];
  sigma[2] = MM3x3[2];
  sigma[3] = MM3x3[4];
  sigma[4] = MM3x3[5];
  sigma[5] = MM3x3[8];
}

//----------------------------------------------------------------------

HyperelasticityFcnSaintVenantKirchhoff::
HyperelasticityFcnSaintVenantKirchhoff(HyperelasticityModelData &hyper, VarFcnBase &vf_)
    : HyperelasticityFcnBase(vf_)
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
HyperelasticityFcnSaintVenantKirchhoff::GetCauchyStressTensor(double *F, [[maybe_unused]] double *V, double *sigma)
{

  MathTools::LinearAlgebra::CalculateATransposeA3x3(F,M3x3); //M(C) = F'F: right Cauchy-Green def. tensor
  M3x3[0] -= 1.0;
  M3x3[4] -= 1.0;
  M3x3[8] -= 1.0;     //M = M - I
  MathTools::LinearAlgebra::CalculateCTimesMatrixA3x3(0.5, M3x3, M3x3); //M(E) = 1/2(C-I): Green strain
  double lambda_trace = lambda*MathTools::LinearAlgebra::CalculateMatrixTrace3x3(M3x3);

  //calculates the second Piola-Kirchhoff stress tensor: P = lambda*tr(E)*I + 2*mu*E
  MathTools::LinearAlgebra::CalculateCTimesMatrixA3x3(2.0*mu, M3x3, M3x3);
  M3x3[0] += lambda_trace;
  M3x3[4] += lambda_trace;
  M3x3[8] += lambda_trace;
  
  //convert to sigma (dim:6)
  double J = MathTools::LinearAlgebra::CalculateDeterminant3x3(F);
  assert(J>0.0);
  ConvertPK2ToCauchy(M3x3, F, J, sigma);

}

//----------------------------------------------------------------------

HyperelasticityFcnModifiedSaintVenantKirchhoff::
HyperelasticityFcnModifiedSaintVenantKirchhoff(HyperelasticityModelData &hyper, VarFcnBase &vf_)
    : HyperelasticityFcnBase(vf_)
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
HyperelasticityFcnModifiedSaintVenantKirchhoff::GetCauchyStressTensor(double *F, 
                                                                      [[maybe_unused]] double *V, double *sigma)
{

  MathTools::LinearAlgebra::CalculateATransposeA3x3(F,N3x3); //N(C) = F'F: right Cauchy-Green def. tensor
  for(int i=0; i<9; i++)
    M3x3[i] = N3x3[i];  //M = N
  M3x3[0] -= 1.0;
  M3x3[4] -= 1.0;
  M3x3[8] -= 1.0;     //M = M - I
  MathTools::LinearAlgebra::CalculateCTimesMatrixA3x3(0.5, M3x3, M3x3); //M(E) = 1/2(C-I): Green strain

  //calculates Cinv
  double Cinv[9], J;
  MathTools::LinearAlgebra::CalculateMatrixInverseAndDeterminant3x3(N3x3, Cinv, &J);
  assert(J>0);
  J = sqrt(J);

  //calculates the second Piola-Kirchhoff stress tensor: P = kappa*ln(J)*inv(C) + 2*mu*E
  MathTools::LinearAlgebra::CalculateMatrixC1APlusC2B3x3(kappa*log(J), Cinv, 2.0*mu, M3x3, N3x3);
  
  //convert to sigma (dim:6)
  ConvertPK2ToCauchy(N3x3, F, J, sigma);

}

//----------------------------------------------------------------------

HyperelasticityFcnNeoHookean::
HyperelasticityFcnNeoHookean(HyperelasticityModelData &hyper, VarFcnBase &vf_)
    : HyperelasticityFcnBase(vf_)
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
HyperelasticityFcnNeoHookean::GetCauchyStressTensor(double *F, [[maybe_unused]] double *V, double *sigma)
{

  MathTools::LinearAlgebra::CalculateAATranspose3x3(F,N3x3); //N(B) = FF': left Cauchy-Green def. tensor

  // calculate I1, I2, I3 (principal invariants)
  double I1 = MathTools::LinearAlgebra::CalculateFirstPrincipalInvariant3x3(N3x3);
  double I3 = MathTools::LinearAlgebra::CalculateThirdPrincipalInvariant3x3(N3x3);
  assert(I1>0.0);
  assert(I3>0.0);

  double J  = sqrt(I3);
  double Jf = pow(J, -2.0/3.0);
  I1 *= Jf;

  double factor1 = mu/J*Jf;
  double factor2 = kappa*(J-1.0) - 1.0/(3.0*J)*(mu*I1);

  MathTools::LinearAlgebra::CalculateCTimesMatrixA3x3(factor1, N3x3, N3x3);

  N3x3[0] += factor2;
  N3x3[4] += factor2;
  N3x3[8] += factor2;

  sigma[0] = N3x3[0];
  sigma[1] = N3x3[1];
  sigma[2] = N3x3[2];
  sigma[3] = N3x3[4];
  sigma[4] = N3x3[5];
  sigma[5] = N3x3[8];

}

//----------------------------------------------------------------------

HyperelasticityFcnMooneyRivlin::
HyperelasticityFcnMooneyRivlin(HyperelasticityModelData &hyper, VarFcnBase &vf_)
    : HyperelasticityFcnBase(vf_)
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
HyperelasticityFcnMooneyRivlin::GetCauchyStressTensor(double *F, [[maybe_unused]] double *V, double *sigma)
{

  MathTools::LinearAlgebra::CalculateAATranspose3x3(F,N3x3); //N(B) = FF': left Cauchy-Green def. tensor
  MathTools::LinearAlgebra::CalculateMatrixMatrixProduct3x3(N3x3,N3x3,M3x3); //M = B*B

  // calculate I1, I2, I3 (principal invariants)
  double I1 = MathTools::LinearAlgebra::CalculateFirstPrincipalInvariant3x3(N3x3);
  double I2 = MathTools::LinearAlgebra::CalculateSecondPrincipalInvariant3x3(N3x3);
  double I3 = MathTools::LinearAlgebra::CalculateThirdPrincipalInvariant3x3(N3x3);
  assert(I1>0.0);
  assert(I2>0.0);
  assert(I3>0.0);

  double J  = sqrt(I3);
  double Jf = pow(J, -2.0/3.0);
  I1 *= Jf;
  I2 *= Jf*Jf;

  double factor1 = 2.0/J*Jf*(C10+I1*C01);
  double factor3 = -2.0/J*Jf*Jf*C01;
  double factor2 = kappa*(J-1.0) - 2.0/(3.0*J)*(C10*I1 + 2.0*C01*I2);

  MathTools::LinearAlgebra::CalculateMatrixC1APlusC2B3x3(factor1, N3x3, factor3, M3x3, N3x3);

  N3x3[0] += factor2;
  N3x3[4] += factor2;
  N3x3[8] += factor2;

  sigma[0] = N3x3[0];
  sigma[1] = N3x3[1];
  sigma[2] = N3x3[2];
  sigma[3] = N3x3[4];
  sigma[4] = N3x3[5];
  sigma[5] = N3x3[8];

}


//----------------------------------------------------------------------





