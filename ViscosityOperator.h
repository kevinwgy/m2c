/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _VISCOSITY_OPERATOR_H_
#define _VISCOSITY_OPERATOR_H_

#include <GradientCalculatorBase.h>
#include <GhostFluidOperator.h>
#include <GlobalMeshInfo.h>
#include <Interpolator.h>
#include <ViscoFcn.h>
#include <memory>

class EmbeddedBoundaryDataSet;
struct Vec5D;

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
  vector<ViscoFcnBase*> visFcn; //!< undefined for INACTIVE_MATERIAL_ID

  //! Mesh info
  SpaceVariable3D &coordinates;
  SpaceVariable3D &delta_xyz;
  SpaceVariable3D &volume;

  GlobalMeshInfo& global_mesh;

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain
  int NX, NY, NZ;

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


  //! Ghost fluid method
  GhostFluidOperator *gfo;
  SpaceVariable3D *Velog; //!< a copy of velocity with ghost nodes populated
  
  // ---------------------------------------------------------------------
  //! internal variables for cylindrical symmetry (x~axial, y~radial)
  bool cylindrical_symmetry;
  SpaceVariable3D *DDXm, *DDXp, *DDYm, *DDYp; //!< biased 3rd order FD at nodes
  SpaceVariable3D *Lam, *Mu;
  SpaceVariable3D *scalarG2; //!< with 2 ghost layers

  //! Class for calculating spatial gradient (FDM)
  GradientCalculatorBase *grad_minus; //left-biased FD
  GradientCalculatorBase *grad_plus;  //right-biased FD
  // ---------------------------------------------------------------------


public:

  ViscosityOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_,
                    vector<VarFcnBase*> &varFcn_, GlobalMeshInfo &global_mesh_,
                    SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_,
                    SpaceVariable3D &volume_,
                    InterpolatorBase &interpolator_, GradientCalculatorBase &grad_,
                    bool with_embedded_boundary);

  ~ViscosityOperator();

  inline vector<ViscoFcnBase*>& GetViscoFcns() {return visFcn;}

  //! calculate diffusion fluxes (on the left-hand-side of the N-S equations; add them to R) 
  void AddDiffusionFluxes(SpaceVariable3D &V, SpaceVariable3D &ID, 
                          vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS,
                          SpaceVariable3D &R);

  //! destroy internal variables
  void Destroy();

private:

  //! enforce cylindrical symmetry (on the left-hand-side of the N-S equations; added to res)
  void AddCylindricalSymmetryTerms(Vec5D*** v5, double*** id, Vec3D*** dxyz, 
                                   vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS,
                                   Vec5D*** res);

  double CalculateLocalDiv2D(Vec5D*** v, double*** id, Vec3D*** coords, int i, int j, int k);


};

#endif
