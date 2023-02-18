/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _HEAT_DIFFUSION_OPERATOR_H_
#define _HEAT_DIFFUSION_OPERATOR_H_

#include <IoData.h>
#include <SpaceVariable.h>
#include <GradientCalculatorBase.h>
#include <VarFcnBase.h>
#include <Interpolator.h>
#include <HeatDiffusionFcn.h>
#include <memory>

class EmbeddedBoundaryDataSet;

/***************************************************
 * class HeatDiffusionOperator calculates the diffusive 
 * heat fluxes (on the LEFT hand side of the Navier-Stokes
 * equations) at cell interfaces, 
 * and adds it to the total flux function
 **************************************************/

class HeatDiffusionOperator
{
  MeshData& iod_mesh;

  EquationsData& iod_eqs;
  
  //! Variable function
  vector<VarFcnBase*>& varFcn; //!< each material has a varFcn
   
  //! Heat diffusion function (one for each material)
  vector<HeatDiffusionFcnBase*> heatdiffFcn;

  //! Mesh info
  SpaceVariable3D &coordinates;
  SpaceVariable3D &delta_xyz;
  SpaceVariable3D &volume;

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain

  //! interpolator
  InterpolatorBase &interpolator;

  //! gradient calculator
  GradientCalculatorBase &grad;

  //! internal variables
  SpaceVariable3D T;//!< temperature
  SpaceVariable3D dTdx_i_minus_half;  //!< dT/dx at i +/- 1/2   
  SpaceVariable3D dTdy_j_minus_half;  //!< dT/dy at j +/- 1/2 
  SpaceVariable3D dTdz_k_minus_half;  //!< dT/dz at k +/- 1/2 
  SpaceVariable3D dTdr; //Used for cylindrical or Spherical coordinates


public:

  HeatDiffusionOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, MeshData &iod_mesh_, EquationsData &iod_eqs_,
                        vector<VarFcnBase*> &varFcn_,
                        SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_, SpaceVariable3D &volume_,
                        InterpolatorBase &interpolator_, GradientCalculatorBase &grad_);

  ~HeatDiffusionOperator();
  //! destroy internal variables
  void Destroy();

  //! calculate diffusion fluxes (on the left-hand-side of the N-S equations; add them to R) 
  void AddDiffusionFluxes(SpaceVariable3D &V, SpaceVariable3D &ID, 
                          vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS,
                          SpaceVariable3D &R);

  //! Add symmetry Term induced by heat diffusion if needed
  void AddSymmetryDiffusionTerms(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &R);

private:
  
  void AddSphericalSymmetryDiffusionTerms(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &R);

  void AddCylindricalSymmetryDiffusionTerms(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &R);

};

#endif
