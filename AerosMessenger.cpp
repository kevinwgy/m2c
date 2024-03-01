/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<AerosMessenger.h>
#include<cassert>
#include<climits>

#define SUGGEST_DT_TAG 444
#define WET_SURF_TAG1 555
#define WET_SURF_TAG2 666
#define WET_SURF_TAG3 888
#define WET_SURF_TAG4 999
#define SUBCYCLING_TAG 777

#define NEGO_NUM_TAG 10000
#define NEGO_BUF_TAG 10001

#define FORCE_TAG 1000
#define DISP_TAG 2000
#define INFO_TAG 3000

#define CRACK_TAG1 22
#define CRACK_TAG2 33
#define CRACK_TAG3 44
#define CRACK_TAG4 55


using std::vector;

extern int verbose;

//---------------------------------------------------------------
// TODO(KW): When I or someone has time, should re-visit the different coupling algorithms (C0, A6, etc.).
//
AerosMessenger::AerosMessenger(AerosCouplingData &iod_aeros_, MPI_Comm &m2c_comm_, MPI_Comm &joint_comm_, 
                               TriangulatedSurface &surf_, vector<Vec3D> &F_)
              : iod_aeros(iod_aeros_), m2c_comm(m2c_comm_), joint_comm(joint_comm_),
                surface(surf_), F(F_), cracking(NULL), numStrNodes(NULL)
{

  MPI_Comm_rank(m2c_comm, &m2c_rank);
  MPI_Comm_size(m2c_comm, &m2c_size);

  // verification of MPI / communicators
  int joint_rank, joint_size;
  MPI_Comm_rank(joint_comm, &joint_rank);
  MPI_Comm_size(joint_comm, &joint_size);
  assert(m2c_rank == joint_rank);
  assert(m2c_size == joint_size);

  MPI_Comm_remote_size(joint_comm, &numAerosProcs);
  //fprintf(stdout,"I am [%d]. Aero-S running on %d procs.\n", m2c_rank, numAerosProcs);

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
  surface.Udot.resize(totalNodes, Vec3D(0.0));
  F.resize(totalNodes, Vec3D(0.0));
  temp_buffer.resize(2*bufsize*totalNodes, 0.0); //a long temporary buffer, can be used for anything

  // allocate memory for the element topology list
  int tmpTopo[nStElems][4];
  switch(elemType) {
    case 3: // all triangles
      totalElems = totalStElems;
      nElems     = nStElems;
      surface.elems.resize(totalElems, Int3(0));
      GetEmbeddedWetSurface(nNodes, surface.X0.data(), nElems, (int*)surface.elems.data(), elemType);
      break;
    case 4: // quadrangles include triangles represented as degenerated quadrangles.
      GetEmbeddedWetSurface(nNodes, surface.X0.data(), nStElems, (int*)tmpTopo, elemType);
      if(cracking) {
        totalElems = totalStElems * 2;
        surface.elems.resize(totalElems, Int3(0));
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
    surface.X[i] = surface.X0[i];


  surface.BuildConnectivities();
  surface.CalculateNormalsAndAreas();


  surface.active_nodes = nNodes;
  surface.active_elems = nElems;

  Negotiate(); // Following the AERO-F/S function name, although misleading

  GetInfo(); // Get algorithm number, dt, and tmax

  if(cracking) {
    GetInitialCrack();
    cracking->setNewCrackingFlag(false);
  }

  structureSubcycling = (algNum == 22) ? GetStructSubcyclingInfo() : 0; 
  //currently, M2C does not support "structure subcycling".

  if(algNum == 6)
    print("- Coupled with Aero-S (running on %d processors) using the A6 algorithm "
          "(Assuming fixed time-step in Aero-S).\n", numAerosProcs);
  else if(algNum == 22)
    print("- Coupled with Aero-S (running on %d processors) using the C0 algorithm.\n", numAerosProcs);



//debug
/*
  if(m2c_rank==0) {
    vector<Vec3D> &x(surface.X);
    for(int i=0; i<x.size(); i++)
      fprintf(stdout,"%d %e %e %e.\n", i, x[i][0], x[i][1], x[i][2]);
    vector<Int3> &e(surface.elems);
    for(int i=0; i<e.size(); i++)
      fprintf(stdout,"%d %d %d %d.\n", i, e[i][0], e[i][1], e[i][2]);
  }
  MPI_Barrier(m2c_comm);
  exit_mpi();
*/
} 

//---------------------------------------------------------------

AerosMessenger::~AerosMessenger()
{
  if(numStrNodes)
    delete [] numStrNodes;

  if(cracking)
    delete cracking;
}

//---------------------------------------------------------------

void
AerosMessenger::Destroy()
{ }

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

  print("- Embedded Surface from Aero-S: Number of Nodes/Elements: %d/%d, Element Type = %d, Fracture = %d.\n", 
        nStNodes, nStElems, eType, (int)crack);
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

  if(verbose>=1)
    print("- Received nodes and elements of embedded surface from Aero-S.\n");

/*
  if(m2c_rank==0) {
    for(int i=0; i<nNodes; i++) {
      fprintf(stdout,"%d  %e %e %e\n", i, nodes[i][0], nodes[i][1], nodes[i][2]);
    }
  }
*/
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

  print("- Received fracture info from Aero-S: total nodes/elements = %d/%d.\n", totalStNodes, totalStElems);
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
  Tria.assign(nTrias, Int3(0));
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

  int numCPUMatchedNodes = m2c_rank ? 0 : totalNodes;

  vector<int> ibuffer;
  if(numCPUMatchedNodes>0) {
    ibuffer.resize(numCPUMatchedNodes);
    for(int i=0; i<numCPUMatchedNodes; i++)
      ibuffer[i] = i; //trivial for M2C/Embedded boundary. (non-trivial in AERO-F/S w/ ALE)
  }

  // send the matched node numbers of this fluid CPU to all the structure CPUs
  vector<MPI_Request> send_requests;
  for (int proc = 0; proc < numAerosProcs; proc++) {
    send_requests.push_back(MPI_Request());
    MPI_Isend(&numCPUMatchedNodes, 1, MPI_INT, proc, NEGO_NUM_TAG, 
              joint_comm, &(send_requests[send_requests.size()-1]));
    if(numCPUMatchedNodes>0)
      MPI_Isend(ibuffer.data(), numCPUMatchedNodes, MPI_INT, proc, NEGO_BUF_TAG, 
                joint_comm, &(send_requests[send_requests.size()-1]));
  }
  MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);

  // receive the list of matched nodes that each structure CPU contains
  if(numCPUMatchedNodes > 0) {

    local2pack.assign(totalNodes,-INT_MAX);

    numStrNodes = new int[numAerosProcs][2];
    for(int proc = 0; proc < numAerosProcs; ++proc) {
      MPI_Recv(&numStrNodes[proc][0], 1, MPI_INT, proc, NEGO_NUM_TAG, joint_comm, MPI_STATUS_IGNORE);
      //fprintf(stdout,"Got %d from proc %d.\n", numStrNodes[proc][0], proc);

      if(numStrNodes[proc][0] > 0) {
        MPI_Recv(ibuffer.data(), numStrNodes[proc][0], MPI_INT, proc, NEGO_BUF_TAG, joint_comm, MPI_STATUS_IGNORE);

        pack2local.reserve(pack2local.size() + numStrNodes[proc][0]);

        for(int i = 0; i < numStrNodes[proc][0]; ++i) {
          int idx = ibuffer[i];
          assert(idx>=0 && idx<totalNodes);

          pack2local.push_back(idx);
          local2pack[idx] = pack2local.size()-1;
        }
      }
    }

    if((int)pack2local.size() != numCPUMatchedNodes) {
      fprintf(stdout, "\033[0;31m*** Error (proc %d): wrong number of matched nodes (%d instead of %d).\n\033[0m",
              m2c_rank, (int)pack2local.size(), numCPUMatchedNodes);
      exit(-1);
    }

    for(auto it = local2pack.begin(); it != local2pack.end(); it++) {
      if(*it < 0) {
        fprintf(stdout, "\033[0;31m*** Error (proc %d): found unmatched node (local id: %d).\n\033[0m",
              m2c_rank, (int)(it - local2pack.begin()));
        exit(-1);
      }
    }

    numStrNodes[0][1] = 0;
    for(int proc = 1; proc < numAerosProcs; ++proc) {
      numStrNodes[proc][1] = numStrNodes[proc - 1][1] + numStrNodes[proc - 1][0];
    }
  }
  
}

//---------------------------------------------------------------

void
AerosMessenger::GetInfo()
{

  double info[5];
  if(m2c_rank == 0) {
    MPI_Recv(info, 5, MPI_DOUBLE, MPI_ANY_SOURCE, INFO_TAG, joint_comm, MPI_STATUS_IGNORE);
  }
  MPI_Bcast(info, 5, MPI_DOUBLE, 0, m2c_comm);

  algNum = int(info[0]);
  dt     = info[1];
  tmax   = info[2];

  //int rstrt = int(info[3]); //not used
  //int smode = int(info[4]); //not used

  // check for consistency in algorithm number
  if(iod_aeros.fsi_algo == AerosCouplingData::A6)
    assert(algNum == 6);
  if(iod_aeros.fsi_algo == AerosCouplingData::C0)
    assert(algNum == 22);

  if(algNum != 6 && algNum != 22) {
    print_error("*** Error: Detected unsupported FSI algorithm from Aero-S (%d).\n", algNum);
    exit_mpi();
  }

  if(algNum == 6) {
    tmax -= 0.5 * dt;
  }
  if(algNum == 22) {
    tmax += 0.5 * dt;
  }

  if(verbose>1) 
    print("- Received from Aero-S: dt = %e, tmax(+/-0.5dt) = %e.\n", dt, tmax);
}

//---------------------------------------------------------------

void
AerosMessenger::GetInitialCrack()
{
  int numConnUpdate, numLSUpdate, newNodes;
  bool need2update = GetNewCrackingStats(numConnUpdate, numLSUpdate, newNodes); // inputs will be modified
  if(!need2update) {
    return;
  }

  // get initial phantom nodes.
  GetInitialPhantomNodes(newNodes, surface.X, nNodes);

  // NOTE: nNodes will be updated in "getNewCracking"
  // get initial phantom elements (topo change)
  GetNewCracking(numConnUpdate, numLSUpdate, newNodes);
}

//---------------------------------------------------------------

bool
AerosMessenger::GetNewCrackingStats(int& numConnUpdate, int& numLSUpdate, int& newNodes)
{
  int size = 4;
  int nNew[size];

  if(m2c_rank == 0)
    MPI_Recv(nNew, size, MPI_INT, MPI_ANY_SOURCE, CRACK_TAG1, joint_comm, MPI_STATUS_IGNORE);

  MPI_Bcast(nNew, 4, MPI_INT, 0, m2c_comm);

  numConnUpdate = nNew[1];
  numLSUpdate = nNew[2];
  newNodes = nNew[3];

  return nNew[0] > 0;
}

//---------------------------------------------------------------

void
AerosMessenger::GetInitialPhantomNodes(int newNodes, vector<Vec3D>& xyz, int nNodes)
{
  int size = newNodes*3;

  double coords[size];

  if(m2c_rank == 0) {
    // assume the correct ordering: nNodes, nNodes+1, ..., nNodes+newNodes-1
    MPI_Recv(coords, size, MPI_DOUBLE, MPI_ANY_SOURCE, CRACK_TAG4, joint_comm, MPI_STATUS_IGNORE);
  }

  MPI_Bcast(coords, size, MPI_DOUBLE, 0, m2c_comm);

  for(int i = 0; i < newNodes; i++)
    for(int j = 0; j < 3; j++) {
      xyz[nNodes + i][j] = coords[i * 3 + j];
    }
}

//---------------------------------------------------------------

void
AerosMessenger::GetNewCracking(int numConnUpdate, int numLSUpdate, int newNodes)
{

  if(numConnUpdate < 1) {
    return;
  }

  int phantElems[5 * numConnUpdate]; // elem id and node id
  double phi[4 * numLSUpdate];
  int phiIndex[numLSUpdate];
  int new2old[std::max(1, newNodes * 2)]; // in the case of element deletion, newNodes might be 0

  GetNewCrackingCore(numConnUpdate, numLSUpdate, phantElems, phi, phiIndex, new2old, newNodes);

  if(elemType != 4) {
    print_error("*** Error: only support quadrangle elements for cracking!\n");
    exit_mpi();
  }

  nNodes += newNodes;
  nElems += cracking->updateCracking(numConnUpdate, numLSUpdate, phantElems, phi, phiIndex, surface.elems, 
                                     nNodes, new2old, newNodes);

  //update info stored in Triangulated surface
  surface.active_nodes = nNodes;
  surface.active_elems = nElems; 

  if(nElems != cracking->usedTrias()) {
    print_error("*** Error: inconsistency in the number of used triangles. (Software bug.)\n");
    exit_mpi();
  }

  if(newNodes!=0)
    assert(numAerosProcs==1); //I think this is a current limitation

  if(m2c_rank == 0 && newNodes) 
    numStrNodes[0][0] = nNodes; //Assuming only talking to one structure proc?


  if(newNodes!=0 && verbose>=1) {
    print("- Received %d new nodes from Aero-S (due to fracture).\n", newNodes);
  }
}

//---------------------------------------------------------------

int 
AerosMessenger::GetNewCracking()
{
  int newNodes, numConnUpdate, numLSUpdate;
  bool need2update = GetNewCrackingStats(numConnUpdate, numLSUpdate, newNodes); // inputs will be modified
  if(!need2update) {
    assert(numConnUpdate == 0);
    return 0;
  }
  GetNewCracking(numConnUpdate, numLSUpdate, newNodes);
  return numConnUpdate;
}

//---------------------------------------------------------------

void
AerosMessenger::GetNewCrackingCore(int numConnUpdate, int numLSUpdate, int *phantoms, double *phi, int *phiIndex, 
                                   int *new2old, int newNodes)
{

  int integer_pack_size = 5 * numConnUpdate + numLSUpdate + 2 * newNodes;
  int integer_pack[integer_pack_size]; // KW: This should be a short array since there will not be many new cracked elements in
                                       //     one time-step. Therefore it should be OK to create and destroy it repeatedly.

  if(m2c_rank == 0) {
    MPI_Recv(integer_pack, integer_pack_size, MPI_INT, MPI_ANY_SOURCE, CRACK_TAG2, joint_comm, MPI_STATUS_IGNORE);
  }
  MPI_Bcast(integer_pack, integer_pack_size, MPI_INT, 0, m2c_comm);

  for(int i = 0; i < 5 * numConnUpdate; i++) {
    phantoms[i] = integer_pack[i];
    assert(integer_pack[i] >= 0);
  }
  for(int i = 0; i < numLSUpdate; i++) {
    phiIndex[i] = integer_pack[5 * numConnUpdate + i];
  }
  for(int i = 0; i < newNodes; i++) {
    new2old[2 * i] = integer_pack[5 * numConnUpdate + numLSUpdate + 2 * i];
    new2old[2 * i + 1] = integer_pack[5 * numConnUpdate + numLSUpdate + 2 * i + 1];
  }

  if(m2c_rank == 0) {
    MPI_Recv(phi, 4*numLSUpdate, MPI_DOUBLE, MPI_ANY_SOURCE, CRACK_TAG3, joint_comm, MPI_STATUS_IGNORE);
  }
  MPI_Bcast(phi, 4*numLSUpdate, MPI_DOUBLE, 0, m2c_comm);

}

//---------------------------------------------------------------

int
AerosMessenger::GetStructSubcyclingInfo()
{
  double info;
  if(m2c_rank == 0) {
    MPI_Recv(&info, 1, MPI_DOUBLE, MPI_ANY_SOURCE, SUBCYCLING_TAG, joint_comm, MPI_STATUS_IGNORE);
  }
  MPI_Bcast(&info, 1, MPI_DOUBLE, 0, m2c_comm);

  return (int)info;
}

//---------------------------------------------------------------

void
AerosMessenger::SendForce()
{
  //IMPORTANT: Assuming that the force has been assembled on Proc 0
  if(m2c_rank == 0) {

    //TODO: Need to take care of "staggering"
    //

    // prepare package
    for(int i=0; i<(int)pack2local.size(); i++) {
      for(int j=0; j<3; j++)
        temp_buffer[3*i+j] = F[pack2local[i]][j];
    //  fprintf(stdout,"temp_buffer[%d] = %e %e %e.\n", i, temp_buffer[3*i], temp_buffer[3*i+1], temp_buffer[3*i+2]);
    }

    vector<MPI_Request> send_requests;

    for(int proc = 0; proc < numAerosProcs; proc++) {
      if(numStrNodes[proc][0] > 0) {
        int size = 3*numStrNodes[proc][0];

        double *localBuffer = temp_buffer.data() + 3*numStrNodes[proc][1];

        send_requests.push_back(MPI_Request());
        MPI_Isend(localBuffer, size, MPI_DOUBLE, proc, FORCE_TAG, 
              joint_comm, &(send_requests[send_requests.size()-1]));
      }
    }
    MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);

  }

}

//---------------------------------------------------------------

void
AerosMessenger::GetDisplacementAndVelocity()
{
  if(m2c_rank == 0) {

    assert(bufsize==6);

    // get disp and velo to temp_buffer
    for(int proc = 0; proc < numAerosProcs; proc++) {
      if(numStrNodes[proc][0] > 0) {
        int size = bufsize*numStrNodes[proc][0]; 
        int first_index =  bufsize*numStrNodes[proc][1];

        double *localBuffer = temp_buffer.data() + first_index;

        MPI_Recv(localBuffer, size, MPI_DOUBLE, proc, DISP_TAG, joint_comm, MPI_STATUS_IGNORE);
      }
    }

    // apply disp and velo to surface and Udot
    for(int i=0; i<nNodes; i++) {
      int id = local2pack[i];
      for(int j=0; j<3; j++)
        surface.X[i][j] = surface.X0[i][j] + temp_buffer[bufsize*id+j];
      for(int j=0; j<3; j++)
        surface.Udot[i][j] = temp_buffer[bufsize*id+3+j];
      if(algNum==6) {//A6
        for(int j=0; j<3; j++)
          surface.X[i][j] += 0.5*dt*surface.Udot[i][j];
      }
    }
  }

  
  //broadcast surface and Udot
  MPI_Bcast((double*)surface.X.data(), 3*nNodes, MPI_DOUBLE, 0, m2c_comm);
  MPI_Bcast((double*)surface.Udot.data(), 3*nNodes, MPI_DOUBLE, 0, m2c_comm);

}

//---------------------------------------------------------------

void
AerosMessenger::SendM2CSuggestedTimeStep(double dtf0)
{
  if(m2c_rank == 0) {

    vector<MPI_Request> send_requests;

    for(int proc = 0; proc < numAerosProcs; proc++) {
      if(numStrNodes[proc][0] > 0) {
        send_requests.push_back(MPI_Request());
        MPI_Isend(&dtf0, 1, MPI_DOUBLE, proc, SUGGEST_DT_TAG, 
              joint_comm, &(send_requests[send_requests.size()-1]));
      }
    }
    MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);
  }
}

//---------------------------------------------------------------

void
AerosMessenger::CommunicateBeforeTimeStepping()
{
  if(algNum == 6) //A6
    CommunicateBeforeTimeSteppingForA6();
  else if(algNum == 22)
    CommunicateBeforeTimeSteppingForC0();
  else {
    print_error("*** Error: Detected unsupported Aero-S algNum: %d.\n", algNum);
    exit_mpi();
  }
}

//---------------------------------------------------------------

void
AerosMessenger::FirstExchange()
{
  if(algNum == 6) //A6
    FirstExchangeForA6();
  else if(algNum == 22)
    FirstExchangeForC0();
  else {
    print_error("*** Error: Detected unsupported Aero-S algNum: %d.\n", algNum);
    exit_mpi();
  }
}

//---------------------------------------------------------------

void
AerosMessenger::Exchange()
{
  if(algNum == 6) //A6
    ExchangeForA6();
  else if(algNum == 22)
    ExchangeForC0();
  else {
    print_error("*** Error: Detected unsupported Aero-S algNum: %d.\n", algNum);
    exit_mpi();
  }
}

//---------------------------------------------------------------

void
AerosMessenger::FinalExchange()
{
  if(algNum == 6) //A6
    FinalExchangeForA6();
  else if(algNum == 22)
    FinalExchangeForC0();
  else {
    print_error("*** Error: Detected unsupported Aero-S algNum: %d.\n", algNum);
    exit_mpi();
  }
}

//---------------------------------------------------------------

void
AerosMessenger::CommunicateBeforeTimeSteppingForA6()
{
  dt *= 0.5;
  GetDisplacementAndVelocity();
}

//---------------------------------------------------------------

void
AerosMessenger::FirstExchangeForA6()
{
  //nothing special
  ExchangeForA6(); 
}

//---------------------------------------------------------------

void
AerosMessenger::ExchangeForA6()
{
  SendForce();  

  if(cracking)
    GetNewCracking();

  GetDisplacementAndVelocity();
}

//---------------------------------------------------------------

void
AerosMessenger::FinalExchangeForA6()
{
  SendForce();
}

//---------------------------------------------------------------

void
AerosMessenger::CommunicateBeforeTimeSteppingForC0()
{
  GetDisplacementAndVelocity();

  SendForce();

  GetInfo(); //get dt, tmax

  dt *= 0.5;
}

//---------------------------------------------------------------

void
AerosMessenger::FirstExchangeForC0()
{
  GetInfo(); //get dt, tmax

  if(cracking)
    GetNewCracking();

  GetDisplacementAndVelocity();
}

//---------------------------------------------------------------

void
AerosMessenger::ExchangeForC0()
{
  SendForce();  

  GetInfo();

  if(cracking)
    GetNewCracking();

  GetDisplacementAndVelocity();
}

//---------------------------------------------------------------

void
AerosMessenger::FinalExchangeForC0()
{
  SendForce();
  dt = 0.0;
}

//---------------------------------------------------------------















//---------------------------------------------------------------

