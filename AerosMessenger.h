#ifndef _AEROS_MESSENGER_H_
#define _AEROS_MESSENGER_H_

#include<mpi.h>
#include<IoData.h>
#include<Vector3D.h>
#include<TriangulatedSurface.h>
#include<CrackingSurface.h>

/*************************************************************
* Class AerosMessenger is responsible for communicating with
* the AERO-S solver, e.g. for fluid-structure interaction
* simulations.
* This class is like (but not the same as) a combination of 
* DynamicNodalTransfer, EmbeddedStructure, StructExc, 
* MatchNodes, and EmbeddedMeshMotionHandler in AERO-F.
************************************************************/

class AerosMessenger {

  AerosCouplingData &iod_aeros;

  MPI_Comm &m2c_comm; //!< This is the M2C communicator
  MPI_Comm &joint_comm; //!< This is the joint communicator of M2C and AERO-S

  int m2c_rank, m2c_size;

  int numAerosProcs;
  int (*numStrNodes)[2];  //!< numStrNodes[AEROS-proc-num][0]: the num of nodes; [1]: index of first node
//  std::vector<int> matchNodes; //!< AERO-S node id -> local (surface) node id
  int bufsize; //!< number of DOFs per node (6)

  int nNodes, totalNodes;
  int nElems, totalElems;
  int elemType; //!< 3 or 4
  TriangulatedSurface &surface; //!< the embedded surface
  CrackingSurface *cracking; //!< activated only if cracking is considered in the structure code
  std::vector<Vec3D> &F;

  double dt, tmax; //!< dt and tmax in AERO-S
  int algNum; //!< algo number received from AERO-S
  int structureSubcycling;

  std::vector<double> temp_buffer;

public:

  AerosMessenger(IoData &iod_, MPI_Comm &m2c_comm_, MPI_Comm &joint_comm_, TriangulatedSurface &surf_, 
                 std::vector<Vec3D> &F_);
  ~AerosMessenger();

  double GetTimeStepSize() {return dt;}
  double GetMaxTime() {return tmax;}
  bool   Cracking()   {return cracking==NULL ? false : true;}

  CrackingSurface *GetPointerToCrackingSurface() {return cracking;} 

  int StructSubcycling() {return structureSubcycling;}

  //! Exchange data w/ AERO-S (called at the beginning of the first time step)
  void FirstExchange();

  //! Exchange data w/ AERO-S (called at every time step)
  void Exchange();

  //! Exchange data w/ AERO-S (called at the last time step)
  void FinalExchange();

protected:

  //! functions called by the constructor
  void GetEmbeddedWetSurfaceInfo(int &eType, bool &crack, int &nStNodes, int &nStElems);
  void GetEmbeddedWetSurface(int nNodes, Vec3D *nodes, int nElems, int *elems, int eType);
  void GetInitialCrackingSetup(int &totalStNodes, int &totalStElems);
  int  SplitQuads(int *quads, int nStElems, std::vector<Int3> &Tria);
  void Negotiate();
  void GetInfo();
  void GetInitialCrack();
  bool GetNewCrackingStats(int& numConnUpdate, int& numLSUpdate, int& newNodes);
  void GetInitialPhantomNodes(int newNodes, std::vector<Vec3D>& xyz, int nNodes);
  void GetNewCracking(int numConnUpdate, int numLSUpdate, int newNodes);
  int  GetNewCracking();
  void GetNewCrackingCore(int numConnUpdate, int numLSUpdate, int *phantoms, double *phi, int *phiIndex,
                          int *new2old, int newNodes);
  int  GetStructSubcyclingInfo();

  
  void GetTimeInfo(); //!< get/update dt and tmax from AERO-S
};













#endif
