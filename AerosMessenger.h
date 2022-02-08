#ifndef _AEROS_MESSENGER_H_
#define _AEROS_MESSENGER_H_

#include<mpi.h>
#include<IoData.h>
#include<Vector3D.h>
#include<TriangulatedSurface.h>
#include<CrackingSurface.h>

// -----------------------------------------------------------
// Class AerosMessenger is responsible for communicating with
// the AERO-S solver, e.g. for fluid-structure interaction
// simulations
// -----------------------------------------------------------

class AerosMessenger {

  AerosCouplingData &iod_aeros;

  MPI_Comm &m2c_comm; //!< This is the M2C communicator
  MPI_Comm &joint_comm; //!< This is the joint communicator of M2C and AERO-S

  int m2c_rank, m2c_size;

  int numAerosProcs;
  int (*numStrNodes)[2];  //!< numStrNodes[AEROS-proc-num][0]: the num of nodes; [1]: index of first node
  vector<int> matchNodes; //!< AERO-S node id -> local (surface) node id
  int bufsize; //!< number of DOFs per node (6)

  int nNodes, totalNodes;
  int nElems, totalElems;
  int elemType; //!< 3 or 4
  TriangulatedSurface surface; //!< the embedded surface
  CrackingSurface *cracking; //!< activated only if cracking is considered in the structure code

  double dt, tmax; //dt and tmax in AERO-S
  int algNum; //algo number received from AERO-S
  int structureSubcycling;

public:

  AerosMessenger(IoData &iod_, MPI_Comm &m2c_comm_, MPI_Comm &joint_comm_);
  ~AerosMessenger();

  double GetTimeStep() {return dt;}
  double GetMaxTime() {return tmax;}

  void Destroy();


protected:

  void GetEmbeddedWetSurfaceInfo(int &eType, bool &crack, int &nStNodes, int &nStElems);
  void GetEmbeddedWetSurface(int nNodes, Vec3D *nodes, int nElems, int *elems, int eType);
  void GetInitialCrackingSetup(int &totalStNodes, int &totalStElems);
  int  SplitQuads(int *quads, int nStElems, vector<Int3> &Tria);
  void Negotiate();
  void GetInfo();
};













#endif
