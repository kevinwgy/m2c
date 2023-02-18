/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _REFERENCE_MAP_OPERATOR_H_
#define _REFERENCE_MAP_OPERATOR_H_

#include<IoData.h>
#include<GradientCalculatorBase.h>
#include<GhostPoint.h>
#include<GlobalMeshInfo.h>

/*********************************************************
 * class ReferenceMapOperator is responsible for
 * computing the function that maps current configuration
 * to the original configuration.
 * Similar to LevelSetOperator, except that the reference
 * map is a vector.
 *********************************************************/

class ReferenceMapOperator
{
  MPI_Comm& comm;
  IoData &iod;
  
  //! Mesh info
  SpaceVariable3D& coordinates;
  GlobalMeshInfo& global_mesh;

  //! ghost nodes
  vector<GhostPoint> ghost_nodes_inner; //!< ghost nodes inside the physical domain (shared with other subd)
  vector<GhostPoint> ghost_nodes_outer; //!< ghost nodes outside the physical domain
  vector<int>        ghost_nodes_outer_tag; //!< -1: corner (not used) 0: slip wall/symm, 1: stick, 2: open

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain

  //! Class for calculating spatial gradient (FDM)
  GradientCalculatorBase *grad_minus; //left-biased FD
  GradientCalculatorBase *grad_plus;  //right-biased FD

  //! Internal variables
  SpaceVariable3D Xil, Xir, Xib, Xit, Xik, Xif;
  SpaceVariable3D vectorG2; //2 ghost layers

public:

  ReferenceMapOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_, SpaceVariable3D &coordinates_,
                       SpaceVariable3D &delta_xyz_, GlobalMeshInfo& global_mesh_,
                       std::vector<GhostPoint> &ghost_nodes_inner_,
                       std::vector<GhostPoint> &ghost_nodes_outer_);  
  ~ReferenceMapOperator();

  void Destroy();

  void SetInitialCondition(SpaceVariable3D &Xi);

  void ApplyBoundaryConditions(SpaceVariable3D &Xi);

  void ComputeResidual(SpaceVariable3D &V, SpaceVariable3D &Xi, SpaceVariable3D &R);

private:

  void TagExternalGhostNodes(); //fill 'ghost_nodes_outer_tag'

};

#endif
