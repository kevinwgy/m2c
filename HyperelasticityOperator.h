#ifndef _VISCOELASTICITY_OPERATOR_H_
#define _VISCOELASTICITY_OPERATOR_H_
#include<ReferenceMapOperator.h>
#include<Interpolator.h>

/*********************************************************
 * class HyperelasticityOperator is responsible for
 * computing the elastic and hyperelastic stresses for
 * solid and solid-like materials.
 *********************************************************/

class HyperelasticityOperator
{
  MPI_Comm& comm;
  IoData &iod;

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain

  //! global mesh
  GlobalMeshInfo& global_mesh;

  //! Solver of the reference map
  ReferenceMapOperator refmap;

  //! Gradient calculator
  GradientCalculatorBase &grad;

  //! Interpolator
  InterpolatorBase &interpolator;

  //! Deformation gradient (dim = 9, i.e. 3x3 matrix)
  SpaceVariable3D F;

  //! Determinant of deformation gradient
  SpaceVariable3D DetF;

  //! Internal variables (for temporary use), dim = 3
  SpaceVariable3D Var1, Var2, Var3;

public:

  HyperelasticityOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_, SpaceVariable3D &coordinates_,
                          SpaceVariable3D &delta_xyz_, GlobalMeshInfo &global_mesh_,
                          InterpolatorBase &interpolator_, GradientCalculatorBase &grad_,
                          std::vector<GhostPoint> &ghost_nodes_inner_,
                          std::vector<GhostPoint> &ghost_nodes_outer_);
  ~HyperelasticityOperator();

  void Destroy();

  void InitializeReferenceMap(SpaceVariable3D &Xi);

  void ApplyBoundaryConditionsToReferenceMap(SpaceVariable3D &Xi);

  void ComputeReferenceMapResidual(SpaceVariable3D &V, SpaceVariable3D &Xi, SpaceVariable3D &R);


private :

  void ComputeDeformationGradient(SpaceVariable3D &Xi); //!< Only within the physical domain

};

#endif
