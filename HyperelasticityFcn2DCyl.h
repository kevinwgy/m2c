/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _HYPERELASTICITY_FCN_2DCYL_H_
#define _HYPERELASTICITY_FCN_2DCYL_H_

#include<HyperelasticityFcn.h>

/****************************************************
 * Class HyperelasticityFcnBase2DCyl and the derived classes
 * are children of HyperelasticityFcnBase. These classes are
 * designed for problems with cylindrical symmetry solved on
 * a 2D mesh. The "x" and "y" coordinates of the 2D mesh
 * correspond to the "z" and "r" cylindrical coordinates,
 * respectively. This class contains functions that
 * compute the fluxes associated with "\sigma_{2D}$ 
 * --- see KW's notes. This class also contains functions that
 * compute the source terms associated with the hoop stress. 
 *
 * Note 1: The input matrix F (deformation gradient) is constructed
 *       as follows:
 *       F[0] = dz/dZ;   F[3] = 0.0;   F[6] = dz/dR;
 *       F[1] = 0.0;     F[4] = r/R;   F[7] = 0.0;
 *       F[2] = dr/dZ;   F[5] = 0.0;   F[8] = dr/dR;
 *
 * Note 2: When combining EOS and hyperelasticity models, a
 *         method is to compute pressure from EOS, and *only*
 *         the deviatoric stress from hyperelasticity. In this
 *         class, the "GetCauchyStress" function computes the
 *         complete stress tensor. But the "EvaluateFlux" function
 *         allows the user to make a decision between the full
 *         tensor and the deviatoric part.
 * Note 3: Matrices follow the 'column-major' / 'column-first'
 *         convention. For example, A[1] is the A(2,1) entry.
 * Note 4: The Cauchy stress tensor (\sigma_{2D}) is symmetric, so
 *         only 3 entries are stored: sigma[0] = sigma_zz,
 *         sigma[1] = sigma_zr, sigma[2] = sigma_rr
 ****************************************************/

//---------------------------------------------------------------------------------
//
class HyperelasticityFcnBase2DCyl : public HyperelasticityFcnBase {

protected:

  double F2x2[4], M2x2[4], N2x2[4]; //!< for temporary use (avoid overriding!)
  double MM2x2[4], NN2x2[4]; //!< for ConvertPK2ToCauchy only.

public:

  HyperelasticityFcnBase2DCyl(VarFcnBase &vf_) : HyperelasticityFcnBase(vf_) {}
  virtual ~HyperelasticityFcnBase2DCyl() {}

  virtual void GetCauchyStressTensor([[maybe_unused]] double *F, [[maybe_unused]] double *V, 
                                     double *sigma, double &sigma_phiphi) {
    for(int i=0; i<3; i++)
      sigma[i] = 0.0;

    sigma_phiphi = 0.0; //hoop stress
  }

  //! compute the flux function
  void EvaluateHyperelasticFluxFunction_F(double* flux/*output*/, double* F, double* V/*state var.*/,
                                          bool deviator_only = true);
  void EvaluateHyperelasticFluxFunction_G(double* flux/*output*/, double* F, double* V/*state var.*/,
                                          bool deviator_only = true);

protected:

  //! PK2->Cauchy. Input: P2D, F2D, J (not J2D). Output: sigma_2D[0,1,2] (2x2 matrix, w/ symmetry)
  void ConvertPK2ToCauchy(double* P, double *F, double J, double *sigma); //!< sigma[0-2] (by symmetry)

};

//---------------------------------------------------------------------------------

class HyperelasticityFcnSaintVenantKirchhoff2DCyl : public HyperelasticityFcnBase2DCyl {

  double lambda, mu; //first and second Lame constants

public:

  HyperelasticityFcnSaintVenantKirchhoff2DCyl(HyperelasticityModelData &hyper, VarFcnBase &vf_);
  ~HyperelasticityFcnSaintVenantKirchhoff2DCyl() {}
  
  void GetCauchyStressTensor(double *F, double *V, double *sigma, double &sigma_phiphi);

};

//---------------------------------------------------------------------------------

class HyperelasticityFcnModifiedSaintVenantKirchhoff2DCyl : public HyperelasticityFcnBase2DCyl {

  double kappa, mu; //bulk modulus and the second Lame constant (i.e. shear modulus)

public:

  HyperelasticityFcnModifiedSaintVenantKirchhoff2DCyl(HyperelasticityModelData &hyper, VarFcnBase &vf_);
  ~HyperelasticityFcnModifiedSaintVenantKirchhoff2DCyl() {}
  
  void GetCauchyStressTensor(double *F, double *V, double *sigma, double &sigma_phiphi);

};

//---------------------------------------------------------------------------------

class HyperelasticityFcnNeoHookean2DCyl: public HyperelasticityFcnBase2DCyl {

  double kappa, mu; //bulk modulus and the second Lame constant (i.e. shear modulus)

public:

  HyperelasticityFcnNeoHookean2DCyl(HyperelasticityModelData &hyper, VarFcnBase &vf_);
  ~HyperelasticityFcnNeoHookean2DCyl() {}
  
  void GetCauchyStressTensor(double *F, double *V, double *sigma, double &sigma_phiphi);

};

//---------------------------------------------------------------------------------

class HyperelasticityFcnMooneyRivlin2DCyl: public HyperelasticityFcnBase2DCyl {

  double kappa, C01, C10; //kappa: bulk modulus, C01+C10 = mu/2

public:

  HyperelasticityFcnMooneyRivlin2DCyl(HyperelasticityModelData &hyper, VarFcnBase &vf_);
  ~HyperelasticityFcnMooneyRivlin2DCyl() {}
  
  void GetCauchyStressTensor(double *F, double *V, double *sigma, double &sigma_phiphi);

};

//---------------------------------------------------------------------------------


#endif
