#include<AerosMessenger.h>

#define WET_SURF_TAG1 555
#define WET_SURF_TAG2 666
#define WET_SURF_TAG3 888
#define WET_SURF_TAG4 999

#define NEGO_NUM_TAG 10000
#define NEGO_BUF_TAG 10001

//---------------------------------------------------------------

AerosMessenger::AerosMessenger(IoData &iod_, MPI_Comm &m2c_comm_, MPI_Comm &joint_comm_) 
              : iod_aeros(iod_.concurrent.aeros), m2c_comm(m2c_comm_), joint_comm(joint_comm_)
{

  MPI_Comm_rank(m2c_comm, &m2c_rank);
  MPI_Comm_size(m2c_comm, &m2c_size);

  // verification of MPI / communicators
  int joint_rank, joint_size;
  MPI_Comm_rank(joint_comm, &joint_rank);
  MPI_Comm_size(joint_comm, &joint_size);
  assert(m2c_rank == joint_rank);
  assert(m2c_size == joint_size)

  MPI_Comm_remote_size(joint_comm, &numAerosProcs);
  fprintf(stderr,"I am [%d][%d]. AERO-S running on %d procs.\n", m2c_rank, joint_rank, numAerosProcs);

  //----------------------
  // copied from AERO-F/DynamicNodalTransfer.cpp (written by KW in grad school...)
  //----------------------
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
      GetEmbeddedWetSurface(nNodes, X0.data(), nElems, (int*)surface.elems.data(), elemType);
      break;
    case 4: // quadrangles include triangles represented as degenerated quadrangles.
      GetEmbeddedWetSurface(nNodes, X0.data(), nStElems, (int*)tmpTopo, elemType);
      if(cracking) {
        totalElems = totalStElems * 2;
        surface.elems.resize(totalElems, Int3D(0));
        nElems = cracking->splitQuads((int*)tmpTopo, nStElems, (int*)surface.elems.data());
      }
      else {
        nElems = SplitQuads((int *)tmpTopo, nStElems, surface.elems); 
        totalElems = nElems;
      }
      break;
    default:
      print_error("*** Error: Element type (%d) of the wet surface not recognized! Must be 3 or 4.\n", elemType);
      exit_mpi();
  }
  for(int i = 0; i < nNodes; i++)
    X[i] = X0[i];

  Negotiate(); // Following the AERO-F/S function name, although misleading

} 

//---------------------------------------------------------------

void
AerosMessenger::GetEmbeddedWetSurfaceInfo(int &eType, bool &crack, int &nStNodes, int &nStElems)
{
  int info[4];
  if(m2c_rank==0)
    MPI_Recv(info, 4, MPI_INT, MPI_ANY_SOURCE, WET_SURF_TAG1, joint_comm, MPI_STATUS_IGNORE);

  MPI_Bcast(info, 4, MPI_INT, 0, m2c_comm);

  eType    = info[0];
  crack    = info[1] ? true : false;
  nStNodes = info[2]; 
  nStElems = info[3];
}

//---------------------------------------------------------------

void
AerosMessenger::GetEmbeddedWetSurface(int nNodes, Vec3D *nodes, int nElems, int *elems, int eType)
{
  if(m2c_rank==0)
    MPI_Recv((double*)nodes, nNodes*3, MPI_DOUBLE, MPI_ANY_SOURCE, WET_SURF_TAG2, joint_comm, MPI_STATUS_IGNORE);

  MPI_Bcast((double*)nodes, nNodes*3, MPI_DOUBLE, 0, m2c_comm);

  if(m2c_rank==0)
    MPI_Recv(elems, nElems*eType, MPI_INT, MPI_ANY_SOURCE, WET_SURF_TAG3, joint_comm, MPI_STATUS_IGNORE);

  MPI_Bcast(elems, nElems*eType, MPI_INT, 0, m2c_comm);
}

//---------------------------------------------------------------

void
AerosMessenger::GetInitialCrackingSetup(int &totalStNodes, int &totalStElems)
{
  int info[2];
  if(m2c_rank==0)
    MPI_Recv(info, 2, MPI_INT, MPI_ANY_SOURCE, WET_SURF_TAG4, joint_comm, MPI_STATUS_IGNORE);

  MPI_Bcast(info, 2, MPI_INT, 0, m2c_comm);
  
  totalStNodes = info[0];
  totalStElems = info[1];
}

//---------------------------------------------------------------

int
AerosMessenger::SplitQuads(int *quads, int nStElems, vector<Int3> &Tria)
{
  int nTrias = 0;
  for(int i = 0; i < nStElems; ++i) {
    if(quads[i * 4 + 2] == quads[i * 4 + 3]) {
      nTrias += 1;
    }
    else {
      nTrias += 2;
    }
  }
  Tria.resize(nTrias, Int3(0));
  int count = 0;
  for(int i = 0; i < nStElems; ++i) {
    Tria[count][0] = quads[i * 4];
    Tria[count][1] = quads[i * 4 + 1];
    Tria[count][2] = quads[i * 4 + 2];
    count++;
    if(quads[i * 4 + 2] == quads[i * 4 + 3]) {
      continue;
    }
    Tria[count][0] = quads[i * 4];
    Tria[count][1] = quads[i * 4 + 2];
    Tria[count][2] = quads[i * 4 + 3];
    count++;
  }

  assert(count == nTrias);

  return nTrias;
}

//---------------------------------------------------------------

void
AerosMessenger::Negotiate()
{
  int numCPUMatchedNodes = m2c_rank? 0 : totalNodes;

  for (int proc = 0; proc < numAerosProcs; proc++) {

    MPI_Isend(&numCPUMatchedNodes, 1, MPI_INT, proc, NEGO_NUM_TAG, joint_comm, xxx);

  }


}

//---------------------------------------------------------------









//---------------------------------------------------------------




