/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _HYPERELASTICITY_FCN_H_
#define _HYPERELASTICITY_FCN_H_

#include<VarFcnBase.h>

/****************************************************
 * Class HyperelasticityFcnBase and the derived classes are
 * responsible for calculating the hyperelastic stress tensor.
 * Each material has a separate 'HyperelasticityFcn'
 *
 * Note 1: The input 'F' is assumed to be the *elastic* deformation
 *         gradient. This may not be the same as the complete
 *         deformation gradient, when there is plasticity 
 *         (or other things). 
 * Note 2: When combining EOS and hyperelasticity models, a
 *         method is to compute pressure from EOS, and *only*
 *         the deviatoric stress from hyperelasticity. In this
 *         class, the "GetCauchyStress" function computes the
 *         complete stress tensor. But the "EvaluateFlux" function
 *         allows the user to make a decision between the full
 *         tensor and the deviatoric part.
 * Note 3: Matrices follow the 'column-major' / 'column-first'
 *         convention. For example, A[1] is the A(2,1) entry.
 * Note 4: The Cauchy stress tensor is symmetric, so only 6
 *         entries are stored: sigma[0] = sigma_xx,
 *         sigma[1] = sigma_xy, sigma[2] = sigma_xz,
 *         sigma[3] = sigma_yy, sigma[4] = sigma_yz,
 *         sigma[5] = sigma_zz
 ****************************************************/

//---------------------------------------------------------------------------------
//
class HyperelasticityFcnBase {

protected:

  VarFcnBase &vf;

  double M3x3[9], N3x3[9]; //!< for temporary use
  double MM3x3[9], NN3x3[9]; //!< for ConvertPK2ToCauchy only

public:

  enum Type {NONE = 0, SAINTVENANT_KIRCHHOFF = 1, MODIFIED_SAINTVENANT_KIRCHHOFF = 2,
             NEO_HOOKEAN = 3, MOONEY_RIVLIN = 4} type;

  HyperelasticityFcnBase(VarFcnBase &vf_) : vf(vf_), type(NONE) {}
  virtual ~HyperelasticityFcnBase() {}

  virtual void GetCauchyStressTensor([[maybe_unused]] double *F, [[maybe_unused]] double *V, double *sigma) { //!< V: state var.
    for(int i=0; i<6; i++)
      sigma[i] = 0.0;
  }

  //! compute the flux function
  void EvaluateHyperelasticFluxFunction_F(double* flux/*output*/, double* F, double* V/*state var.*/,
                                          bool deviator_only = true);
  void EvaluateHyperelasticFluxFunction_G(double* flux/*output*/, double* F, double* V/*state var.*/,
                                          bool deviator_only = true);
  void EvaluateHyperelasticFluxFunction_H(double* flux/*output*/, double* F, double* V/*state var.*/,
                                          bool deviator_only = true);

protected:

  void ConvertPK2ToCauchy(double* P, double *F, double J, double *sigma); //!< sigma[0-5] (because of symmetry)

};

//---------------------------------------------------------------------------------

class HyperelasticityFcnSaintVenantKirchhoff : public HyperelasticityFcnBase {

  double lambda, mu; //first and second Lame constants

public:

  HyperelasticityFcnSaintVenantKirchhoff(HyperelasticityModelData &hyper, VarFcnBase &vf_);
  ~HyperelasticityFcnSaintVenantKirchhoff() {}
  
  void GetCauchyStressTensor(double *F, double *V, double *sigma);

};

//---------------------------------------------------------------------------------

class HyperelasticityFcnModifiedSaintVenantKirchhoff : public HyperelasticityFcnBase {

  double kappa, mu; //bulk modulus and the second Lame constant (i.e. shear modulus)

public:

  HyperelasticityFcnModifiedSaintVenantKirchhoff(HyperelasticityModelData &hyper, VarFcnBase &vf_);
  ~HyperelasticityFcnModifiedSaintVenantKirchhoff() {}
  
  void GetCauchyStressTensor(double *F, double *V, double *sigma);

};

//---------------------------------------------------------------------------------

class HyperelasticityFcnNeoHookean: public HyperelasticityFcnBase {

  double kappa, mu; //bulk modulus and the second Lame constant (i.e. shear modulus)

public:

  HyperelasticityFcnNeoHookean(HyperelasticityModelData &hyper, VarFcnBase &vf_);
  ~HyperelasticityFcnNeoHookean() {}
  
  void GetCauchyStressTensor(double *F, double *V, double *sigma);

};

//---------------------------------------------------------------------------------

class HyperelasticityFcnMooneyRivlin: public HyperelasticityFcnBase {

  double kappa, C01, C10; //kappa: bulk modulus, C01+C10 = mu/2

public:

  HyperelasticityFcnMooneyRivlin(HyperelasticityModelData &hyper, VarFcnBase &vf_);
  ~HyperelasticityFcnMooneyRivlin() {}
  
  void GetCauchyStressTensor(double *F, double *V, double *sigma);

};

//---------------------------------------------------------------------------------


#endif
