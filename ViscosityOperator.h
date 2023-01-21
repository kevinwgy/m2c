#ifndef _VISCOSITY_OPERATOR_H_
#define _VISCOSITY_OPERATOR_H_

#include <GradientCalculatorBase.h>
#include <Interpolator.h>
#include <ViscoFcn.h>
#include <memory>

class EmbeddedBoundaryDataSet;

/***************************************************
 * class ViscosityOperator calculates the viscous
 * fluxes (on the LEFT hand side of the Navier-Stokes
 * equations) at cell interfaces, 
 * and adds it to the total flux function
 **************************************************/

class ViscosityOperator
{
  EquationsData& iod_eqs;
  
  //! Variable function
  vector<VarFcnBase*>& varFcn; //!< each material has a varFcn

  //! Viscosity function (one for each material)
  vector<ViscoFcnBase*> visFcn;

  //! Mesh info
  SpaceVariable3D &coordinates;
  SpaceVariable3D &delta_xyz;

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain

  //! interpolator
  InterpolatorBase &interpolator;

  //! gradient calculator
  GradientCalculatorBase &grad;

  //! internal variables -- storing data at cell interfaces
  SpaceVariable3D V_i_minus_half;//!< velocity (dim = 3) 
  SpaceVariable3D V_j_minus_half;//!< velocity (dim = 3) 
  SpaceVariable3D V_k_minus_half;//!< velocity (dim = 3) 
  SpaceVariable3D dVdx_i_minus_half;  //!< du/dx, dv/dx, dw/dx at i +/- 1/2 (dim = 3)  
  SpaceVariable3D dVdx_j_minus_half;  //!< du/dx, dv/dx, dw/dx at j +/- 1/2 (dim = 3)
  SpaceVariable3D dVdx_k_minus_half;  //!< du/dx, dv/dx, dw/dx at k +/- 1/2 (dim = 3)
  SpaceVariable3D dVdy_i_minus_half;  //!< du/dy, dv/dy, dw/dy at i +/- 1/2 (dim = 3)  
  SpaceVariable3D dVdy_j_minus_half;  //!< du/dy, dv/dy, dw/dy at j +/- 1/2 (dim = 3)
  SpaceVariable3D dVdy_k_minus_half;  //!< du/dy, dv/dy, dw/dy at k +/- 1/2 (dim = 3)
  SpaceVariable3D dVdz_i_minus_half;  //!< du/dz, dv/dz, dw/dz at i +/- 1/2 (dim = 3)  
  SpaceVariable3D dVdz_j_minus_half;  //!< du/dz, dv/dz, dw/dz at j +/- 1/2 (dim = 3)
  SpaceVariable3D dVdz_k_minus_half;  //!< du/dz, dv/dz, dw/dz at k +/- 1/2 (dim = 3)

  //! internal variables for cylindrical symmetry (x~axial, y~radial)
  bool cylindrical_symmetry;
  SpaceVariable3D *dudx, *dvdx, *dudy, *dvdy; //!< at nodes
  SpaceVariable3D *Lam, *Mu, *dLamdx, *dLamdy;

public:

  ViscosityOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, EquationsData &iod_eqs_,
                    vector<VarFcnBase*> &varFcn_,
                    SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_,
                    InterpolatorBase &interpolator_, GradientCalculatorBase &grad_);

  ~ViscosityOperator();

  //! calculate diffusion fluxes (on the left-hand-side of the N-S equations; add them to R) 
  void AddDiffusionFluxes(SpaceVariable3D &V, SpaceVariable3D &ID, 
                          vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS,
                          SpaceVariable3D &R);

  //! destroy internal variables
  void Destroy();

private:

  void PopulateGhostNodes(SpaceVariable3D &V, SpaceVariable3D &ID,
                          vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS);

};

#endif
