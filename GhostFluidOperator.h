#ifndef _GHOST_FLUID_OPERATOR_H_
#define _GHOST_FLUID_OPERATOR_H_

#include<SpaceVariable.h>
#include<NeighborCommunicator.h>
#include<GlobalMeshInfo.h>
#include<vector>
#include<memory> //unique_ptr

class EmbeddedBoundaryDataSet;

/**************************************************************
* Class GhostFluidOperator is responsible for populating ghost
* nodes in ``inactive'' regions next to an embedded boundary.
* It is a "process-on-order" type of class that requires other
* classes to provide the interface tracking information and the
* state variables to be operated on.
**************************************************************/

class GhostFluidOperator {

  MPI_Comm &comm;
  int rank, size;

  NeighborCommunicator *neicomm_ptr;

  GlobalMeshInfo &global_mesh;

  SpaceVariable3D Tag;

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain
  int NX, NY, NZ;


public:

  GhostFluidOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, GlobalMeshInfo &global_mesh_);

  ~GhostFluidOperator();

  void Destroy();

  int PopulateGhostNodesForViscosityOperator(SpaceVariable3D &V, SpaceVariable3D &ID,
                                             std::vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS,
                                             SpaceVariable3D &Vgf);


};


#endif
