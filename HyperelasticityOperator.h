/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _VISCOELASTICITY_OPERATOR_H_
#define _VISCOELASTICITY_OPERATOR_H_
#include<ReferenceMapOperator.h>
#include<HyperelasticityFcn.h>
#include<Interpolator.h>

class EmbeddedBoundaryDataSet;

/*********************************************************
 * class HyperelasticityOperator is responsible for
 * computing the elastic and hyperelastic stresses for
 * solid and solid-like materials.
 * Note: This class treats 3D domains and 2D + cylindrical
 *       symmetry separately. This is slightly different
 *       from how convection and viscosity fluxes are calculated.
 *********************************************************/

class HyperelasticityOperator
{
  MPI_Comm& comm;
  IoData &iod;

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain

  bool cylindrical_symmetry;

  //! global mesh
  GlobalMeshInfo& global_mesh;

  //! Variable function
  vector<VarFcnBase*>& varFcn; //!< each material has a varFcn
 
  //! Solver of the reference map
  ReferenceMapOperator refmap;

  //! Hyperelasticity models (one for each material)
  vector<HyperelasticityFcnBase*> hyperFcn;

  //! Gradient calculator
  GradientCalculatorBase &grad;

  //! Interpolator
  InterpolatorBase &interpolator;

  //! Deformation gradient (dim = 9, i.e. 3x3 matrix)
  SpaceVariable3D F;

  //! Volumetric deformation (sqrt(|F^T*F|), rho0/rho)
  SpaceVariable3D J;

  //! Internal variables (for temporary use), dim = 3
  SpaceVariable3D Var1, Var2, Var3;


  //! internal variables -- storing data at cell interfaces
  SpaceVariable3D V_i_minus_half;//!< velocity (dim = 3)
  SpaceVariable3D V_j_minus_half;//!< velocity (dim = 3)
  SpaceVariable3D V_k_minus_half;//!< velocity (dim = 3)
  SpaceVariable3D dXidx_i_minus_half;  //!< dxi_x/dx, dxi_y/dx, dxi_z/dx at i +/- 1/2 (dim = 3)
  SpaceVariable3D dXidx_j_minus_half;  //!< dxi_x/dx, dxi_y/dx, dxi_z/dx at j +/- 1/2 (dim = 3)
  SpaceVariable3D dXidx_k_minus_half;  //!< dxi_x/dx, dxi_y/dx, dxi_z/dx at k +/- 1/2 (dim = 3)
  SpaceVariable3D dXidy_i_minus_half;  //!< dxi_x/dy, dxi_y/dy, dxi_z/dy at i +/- 1/2 (dim = 3)
  SpaceVariable3D dXidy_j_minus_half;  //!< dxi_x/dy, dxi_y/dy, dxi_z/dy at j +/- 1/2 (dim = 3)
  SpaceVariable3D dXidy_k_minus_half;  //!< dxi_x/dy, dxi_y/dy, dxi_z/dy at k +/- 1/2 (dim = 3)
  SpaceVariable3D dXidz_i_minus_half;  //!< dxi_x/dz, dxi_y/dz, dxi_z/dz at i +/- 1/2 (dim = 3)
  SpaceVariable3D dXidz_j_minus_half;  //!< dxi_x/dz, dxi_y/dz, dxi_z/dz at j +/- 1/2 (dim = 3)
  SpaceVariable3D dXidz_k_minus_half;  //!< dxi_x/dz, dxi_y/dz, dxi_z/dz at k +/- 1/2 (dim = 3)


public:

  HyperelasticityOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_,
                          vector<VarFcnBase*> &varFcn_,
                          SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_,
                          GlobalMeshInfo &global_mesh_,
                          InterpolatorBase &interpolator_, GradientCalculatorBase &grad_,
                          std::vector<GhostPoint> &ghost_nodes_inner_,
                          std::vector<GhostPoint> &ghost_nodes_outer_);
  ~HyperelasticityOperator();

  void Destroy();

  void InitializeReferenceMap(SpaceVariable3D &Xi);

  void ApplyBoundaryConditionsToReferenceMap(SpaceVariable3D &Xi);

  void ComputeReferenceMapResidual(SpaceVariable3D &V, SpaceVariable3D &Xi, SpaceVariable3D &R);

  void ComputeDeformationGradientAtNodes(SpaceVariable3D &Xi); //!< Only within the physical domain

  void AddHyperelasticityFluxes(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &Xi,
                                vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS,
                                SpaceVariable3D &R);

private:

  void ComputeDeformGradAtNodes3D(SpaceVariable3D &Xi); //!< Only within the physical domain

  void ComputeDeformGradAtNodes2DCylindrical(SpaceVariable3D &Xi); //!< Only within the physical domain

  void AddFluxes3D(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &Xi,
                   vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS,
                   SpaceVariable3D &R);

  void AddFluxes2DCylindrical(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &Xi,
                              vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS,
                              SpaceVariable3D &R);

};

#endif
