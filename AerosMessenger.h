#ifndef _AEROS_MESSENGER_H_
#define _AEROS_MESSENGER_H_

#include<mpi.h>
#include<IoData.h>
#include<Vector3D.h>

// -----------------------------------------------------------
// Class AerosMessenger is responsible for communicating with
// the AERO-S solver, e.g. for fluid-structure interaction
// simulations
// -----------------------------------------------------------

class AerosMessenger {

  AerosCouplingData &iod_aeros;

  MPI_Comm &m2c_comm; //!< This is the M2C communicator
  MPI_Comm &joint_comm; //!< This is the joint communicator of M2C and AERO-S

  int m2c_rank, m2c_size, joint_rank, joint_size; 

  int numAerosProcs;
  int bufsize; //!< number of DOFs per node (6)

  int nNodes, totalNodes;
  int nElems, totalElems;
  int elemType; //!< 3 or 4
  TriangulatedSurface *surface; //!< the embedded surface
  CrackingSurface *cracking; //!< activated only if cracking is considered in the structure code


public:

  AerosMessenger(IoData &iod_, MPI_Comm &m2c_comm_, MPI_Comm &joint_comm_);
  ~AerosMessenger();

  void Destroy();


};

You should have a reference to an EmbeddedSurface, and make operations such as send / receive forces and xxx


Main should create and own a bunch of embedded surfaces. Time stepper may also get a reference.












#endif
