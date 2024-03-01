/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _VISCOELASTICITY_OPERATOR_H_
#define _VISCOELASTICITY_OPERATOR_H_
#include<ReferenceMapOperator.h>
#include<HyperelasticityFcn2DCyl.h>
#include<Interpolator.h>

class EmbeddedBoundaryDataSet;
class Vec5D;

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

  SpaceVariable3D &coordinates; //!< mesh coords.
  SpaceVariable3D &volume; //!< volume of control volumes
  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain

  bool cylindrical_symmetry;

  vector<bool> deviator_only; //!< For each material, only apply the stress deviator?

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
                          SpaceVariable3D &volume_, GlobalMeshInfo &global_mesh_,
                          InterpolatorBase &interpolator_, GradientCalculatorBase &grad_,
                          std::vector<GhostPoint> &ghost_nodes_inner_,
                          std::vector<GhostPoint> &ghost_nodes_outer_);
  ~HyperelasticityOperator();

  void Destroy();

  void InitializeReferenceMap(SpaceVariable3D &Xi);

  void ApplyBoundaryConditionsToReferenceMap(SpaceVariable3D &Xi);

  void ComputeReferenceMapResidual(SpaceVariable3D &V, SpaceVariable3D &Xi, SpaceVariable3D &R, double time);

  void ComputeDeformationGradientAtNodes(SpaceVariable3D &Xi); //!< Only within the physical domain

  void ComputePrincipalStresses(SpaceVariable3D &Xi, SpaceVariable3D &V, SpaceVariable3D &ID,
                                SpaceVariable3D &PS);

  //! This is a function customized for "probe output". Be careful if you use it elsewhere!
  void ComputePrincipalStressesAtProbes(SpaceVariable3D &Xi, SpaceVariable3D &ID,
           std::vector<Int3> &ijk, std::vector<std::pair<int, std::array<bool,8> > > &ijk_valid,
           std::vector<Vec3D> &trilinear_coords, Vec5D*** v, std::vector<Vec3D> &sol);
           
  //! "f" is the 3x3 (i.e. double[9]) deformation gradient
  Vec3D ComputePrincipalStressesAtPoint(double *f, Vec5D &v, int id);

  void AddHyperelasticityFluxes(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &Xi,
                                vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS,
                                SpaceVariable3D &R);

  //! Test cases
  void PrescribeVelocityForTesting(SpaceVariable3D &V, double time);

private:

  void ComputeDeformGradAtNodes3D(SpaceVariable3D &Xi); //!< Only within the physical domain

  void ComputeDeformGradAtNodes2DCylindrical(SpaceVariable3D &Xi); //!< Only within the physical domain

  void ComputePrincipalStresses3D(SpaceVariable3D &Xi, SpaceVariable3D &V, SpaceVariable3D &ID,
                                  SpaceVariable3D &PS);

  void ComputePrincipalStresses2DCylindrical(SpaceVariable3D &Xi, SpaceVariable3D &V,
                                             SpaceVariable3D &ID, SpaceVariable3D &PS);

  Vec3D ComputePrincipalStressesAtPoint3D(double *f, Vec5D &v, int id);

  Vec3D ComputePrincipalStressesAtPoint2DCylindrical(double *f, Vec5D &v, int id);

  void AddFluxes3D(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &Xi,
                   vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS,
                   SpaceVariable3D &R);

  void AddFluxes2DCylindrical(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &Xi,
                              vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS,
                              SpaceVariable3D &R);

  void AddCylindricalSourceTerms(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &Xi,
                                 vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS,
                                 SpaceVariable3D &R);

};

#endif
