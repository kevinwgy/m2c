#ifndef _GHOST_FLUID_OPERATOR_H_
#define _GHOST_FLUID_OPERATOR_H_

#include<NeighborCommunicator.h>
#include<GlobalMeshInfo.h>
#include<GhostPoint.h>
#include<EmbeddedBoundaryFormula.h>
#include<vector>

class SpaceVariable3D;
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

public:

  GhostFluidOperator(MPI_Comm &comm_, GlobalMeshInfo &global_mesh_);

  ~GhostFluidOperator();

  void Destroy();

  int PopulateGhostNodesForViscosityOperator(SpaceVariable3D &V, SpaceVariable3D &ID,
                                             std::vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS,
                                             SpaceVariable3D &Vgf);


};


#endif
