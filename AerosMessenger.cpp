#include<AerosMessenger.h>

//---------------------------------------------------------------

AerosMessenger::AerosMessenger(IoData &iod_, MPI_Comm &m2c_comm_, MPI_Comm &joint_comm_) 
              : iod_aeros(iod_.concurrent.aeros), m2c_comm(m2c_comm_), joint_comm(joint_comm_)
{

  MPI_Comm_rank(m2c_comm, &m2c_rank);
  MPI_Comm_size(m2c_comm, &m2c_size);
  MPI_Comm_rank(joint_comm, &joint_rank);
  MPI_Comm_size(joint_comm, &joint_size);

  MPI_Comm_remote_size(joint_comm, &numAerosProcs);
  fprintf(stderr,"AERO-S running on %d procs.\n", numAerosProcs);

  // copy from AERO-F/DynamicNodalTransfer.cpp (written by KW in grad school...)
  bufsize = 6; //6 dofs per node (disp and velo)

  bool crack;
  int nStNodes, nStElems, totalStNodes, totalStElems;

  GetEmbeddedWetSurfaceInfo(elemType, crack, nStNodes, nStElems);

  // initialize cracking information
  if(crack) {
    GetInitialCrackingSetup(totalStNodes, totalStElems);
    cracking = new CrackingSurface(elemType, nStElems, totalStElems, nStNodes, totalStNodes);
  }
  else {
    totalStNodes = nStNodes;
    totalStElems = nStElems;
  }

  // allocate memory for the node list
  totalNodes = totalStNodes;
  nNodes     = nStNodes;
  surface.X.resize(totalNodes, Vec3D(0.0));
  surface.X0.resize(totalNodes, Vec3D(0.0));

  // allocate memory for the element topology list
  int tmpTopo[nStElems][4];
  switch(elemType) {
    case 3: // all triangles
      totalElems = totalStElems;
      nElems     = nStElems;
      surface.elems.resize(totalElems, Int3D(0));
      Tria = new int[totalElems][3];
      GetEmbeddedWetSurface(nNodes, X0.data(), nElems, (int *)Tria, elemType);
      break;
    case 4: // quadrangles include triangles represented as degenerated quadrangles.
      structExc->getEmbeddedWetSurface(nNodes, X0.data(), nStElems, (int *)tmpTopo, elemType);
      if(cracking) {
        totalElems = totalStElems * 2;
        Tria = new int[totalElems][3];
        nElems = cracking->splitQuads((int *)tmpTopo, nStElems, Tria);
      }
      else {
        splitQuads((int *)tmpTopo, nStElems); // memory for Tria will be allocated
      }
      break;
    default:
      com.fprintf(stderr, "*** Error: Element type (%d) of the wet surface not recognized! Must be 3 or 4.\n", elemType);
      exit(-1);
  }
  for(int i = 0; i < nNodes; i++) {
    for(int j = 0; j < 3; ++j) {
      X[i][j] = X0[i][j];
    }
  }
  if(com.cpuNum() == 0) {
    mns[0]->autoInit(totalNodes); // in case of cracking, match all the nodes including inactive ones.
    structExc->updateMNS(mns);
  }
} 
