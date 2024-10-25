/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _GHOST_FLUID_OPERATOR_H_
#define _GHOST_FLUID_OPERATOR_H_

#include<SpaceVariable.h>
#include<NeighborCommunicator.h>
#include<memory> //unique_ptr

class EmbeddedBoundaryDataSet;
class EmbeddedBoundaryFormula;

/**************************************************************
* Class GhostFluidOperator is responsible for populating ghost
* nodes in ``inactive'' regions next to an embedded boundary.
* It is a "process-on-order" type of class that requires other
* classes to provide the interface tracking information and the
* state variables to be operated on.
**************************************************************/

class GhostFluidOperator {

  MPI_Comm &comm;
  int mpi_rank, mpi_size;

  NeighborCommunicator neicomm;

  GlobalMeshInfo &global_mesh;

  SpaceVariable3D Tag;
  std::vector<Int3> ghosts_ijk;

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain
  int NX, NY, NZ;

  //! Ghost nodes --- not limited to inactive nodes
  std::vector<std::pair<Int3, EmbeddedBoundaryFormula> > ghostNodes1; //!< ghosts whose image is inside the subdomain

  std::vector<std::vector<Int3> > ghostNodes2; //!< ghosts whose image is inside another subdomain
  std::vector<int> ghostNodes2_sender; //!< owner (subD id.) of each node in ghostNodes2

  std::vector<std::vector<std::pair<Int3, EmbeddedBoundaryFormula> > > friendsGhostNodes;
  std::vector<int> friendsGhostNodes_receiver;

  //! Values at ghost nodes 
  std::vector<double> scalar_gn1; //!< ghostNodes1
  std::vector<std::vector<double> > scalar_gn2; //!< ghostNodes2

public:

  GhostFluidOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, GlobalMeshInfo &global_mesh_);

  ~GhostFluidOperator();

  void Destroy();

  //! call neicomm to do the job (layer: graph distance, dist: physical distance)
  int SetupCustomNeighborCommunicator(int layer, double dist = 0.0);

  int PopulateInactiveNodesForVisco(SpaceVariable3D &V, SpaceVariable3D &ID,
                                    std::vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS,
                                    SpaceVariable3D &Vgf);

  int PopulateInactiveNodesForInco(SpaceVariable3D &V3, SpaceVariable3D &ID,
                                   std::vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS);

private:

  //! Find ghost nodes. Update Tag and ghosts_ijk.
  int TagInactiveNodes(double*** id, int nLayers); //!< returns the total number of ghosts

  //! Find projection points (not unique! just find one.). Also finds the image by mirroring
  void ProjectGhostsToInterface(std::vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS, 
                                std::vector<Vec3D> &ghosts_xyz, std::vector<Vec3D> &ghosts_xp,
                                std::vector<Vec3D> &ghosts_vp, std::vector<Int3> &images_ijk,
                                std::vector<Vec3D> &images_xi, std::vector<Vec3D> &images_xyz,
                                std::vector<int> &ghosts_tag);

  //! Find a *valid* image for each ghost node. "ghosts_tag" is completed. Some images are changed.
  void FindImagesForGhosts(double*** id, std::vector<Vec3D> &ghosts_xyz,
                           std::vector<Int3> &images_ijk, std::vector<Vec3D> &images_xi,
                           std::vector<Vec3D> &images_xyz, std::vector<int> &ghosts_tag);

};


#endif
