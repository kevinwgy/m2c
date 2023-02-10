#ifndef _GHOST_FLUID_OPERATOR_H_
#define _GHOST_FLUID_OPERATOR_H_

#include<NeighborCommunicator.h>
#include<GhostPoint.h>
#include<EmbeddedBoundaryFormula.h>

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
  NeighborCommunicator neighbor_comm;

public:

  GhostFluidOperator(xxx)










};


#endif
