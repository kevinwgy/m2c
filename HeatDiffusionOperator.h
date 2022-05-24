#ifndef _HEAT_DIFFUSION_OPERATOR_H_
#define _HEAT_DIFFUSION_OPERATOR_H_

#include <GradientCalculatorBase.h>
#include <VarFcnBase.h>
#include <Interpolator.h>
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
  EquationsData& iod_eqs;
  
  //! Variable function
  vector<VarFcnBase*>& varFcn; //!< each material has a varFcn

  //! Mesh info
  SpaceVariable3D &coordinates;
  SpaceVariable3D &delta_xyz;

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


public:

  HeatDiffusionOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, EquationsData &iod_eqs_,
                        vector<VarFcnBase*> &varFcn_,
                        SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_,
                        InterpolatorBase &interpolator_, GradientCalculatorBase &grad_);

  ~HeatDiffusionOperator();
  //! destroy internal variables
  void Destroy();

  //! calculate diffusion fluxes (on the left-hand-side of the N-S equations; add them to R) 
  void AddDiffusionFluxes(SpaceVariable3D &V, SpaceVariable3D &ID, 
                          vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS,
                          SpaceVariable3D &R);

};

#endif
