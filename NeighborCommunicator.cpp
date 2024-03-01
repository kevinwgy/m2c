/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <NeighborCommunicator.h>
#include <SpaceVariable.h>
#include <cassert>
#include <map>

using std::vector;

//----------------------------------------------------
NeighborCommunicator::NeighborCommunicator(MPI_Comm& comm_, GlobalMeshInfo& global_mesh_)
                    : comm(comm_), global_mesh(global_mesh_)
{
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  assert(rank>=0 && rank<size);

  subD_neighbors_all = &global_mesh.GetAllNeighborsOfSub(rank);
  subD_neighbors_face_edge = &global_mesh.GetFaceEdgeNeighborsOfSub(rank);
  subD_neighbors_face = &global_mesh.GetFaceNeighborsOfSub(rank);
}

//----------------------------------------------------

NeighborCommunicator::~NeighborCommunicator()
{ }

//----------------------------------------------------

void
NeighborCommunicator::Send(int exchange_type, // 0~all, 1~face_edge, 2~face
                           vector<vector<double> > &Export, vector<vector<double> > &Import)
{

  //Sanity check
  assert(exchange_type==2 || exchange_type==1 || exchange_type==0);

  std::vector<int> *neighbors_ptr; //points to the correct neighborhood
  if(exchange_type==2)
    neighbors_ptr = subD_neighbors_face;
  else if(exchange_type==1)
    neighbors_ptr = subD_neighbors_face_edge;
  else
    neighbors_ptr = subD_neighbors_all;

  int Nneigh = neighbors_ptr->size();
  
  assert((int)Export.size()==Nneigh);
  assert((int)Import.size()==Nneigh);


  //Step 1: Get the size of each package
  vector<int> counter(Nneigh, 0);
   
  int proc, count;
  vector<MPI_Request> send_requests;
  vector<MPI_Request> recv_requests;
  for(int p=0; p<Nneigh; p++) {
    send_requests.push_back(MPI_Request());
    proc = (*neighbors_ptr)[p];
    count = Export[p].size();
    MPI_Isend(&count, 1, MPI_INT, proc, rank/*tag*/, comm, &send_requests.back());

    recv_requests.push_back(MPI_Request());
    MPI_Irecv(&counter[p], 1, MPI_INT, proc, proc, comm, &recv_requests.back());
  }

  MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE); //not necessary(?)
  MPI_Waitall(recv_requests.size(), recv_requests.data(), MPI_STATUSES_IGNORE);
  send_requests.clear();
  recv_requests.clear();


  //Step 2: Resize Import
  for(int p=0; p<Nneigh; p++)
    Import[p].resize(counter[p]);


  //Step 3: Get the actual packages
  for(int p=0; p<Nneigh; p++) {
    if(!Export[p].empty()) {
      send_requests.push_back(MPI_Request());
      proc = (*neighbors_ptr)[p];
      MPI_Isend(Export[p].data(), Export[p].size(), MPI_DOUBLE, proc, rank/*tag*/, 
                comm, &send_requests.back());
    }
    if(counter[p]!=0) {
      recv_requests.push_back(MPI_Request());
      proc = (*neighbors_ptr)[p];
      MPI_Irecv(Import[p].data(), counter[p], MPI_DOUBLE, proc, proc, 
                comm, &recv_requests.back());
    }
  }

  MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE); //not necessary(?)
  MPI_Waitall(recv_requests.size(), recv_requests.data(), MPI_STATUSES_IGNORE);

}

//----------------------------------------------------

void
NeighborCommunicator::Request(SpaceVariable3D &V, vector<Int3> &Request, vector<double> &Received,
                              int exchange_type) //!< 0~all, 1~face_edge, 2~face)
{

  int dof = V.NumDOF();
  
  double*** v = V.GetDataPointer();

  //Sanity check
  assert(exchange_type==2 || exchange_type==1 || exchange_type==0);

  std::vector<int> *neighbors_ptr; //points to the correct neighborhood
  if(exchange_type==2)
    neighbors_ptr = subD_neighbors_face;
  else if(exchange_type==1)
    neighbors_ptr = subD_neighbors_face_edge;
  else
    neighbors_ptr = subD_neighbors_all;

  int Nneigh = neighbors_ptr->size();
  

  Received.resize(dof*Request.size());

  //Step 1: Construct array for each neighbor
  vector<vector<Int3> > neighbor_requests(Nneigh); 
  vector<vector<int> > neighbor_received_map(Nneigh); 
  int owner, nei;
  for(int n=0; n<(int)Request.size(); n++) {
    owner = global_mesh.GetOwnerOfNode(Request[n][0], Request[n][1], Request[n][2], true);
    if(owner==rank) {//weird, but fine.
      for(int d=0; d<dof; d++)
        Received[n*dof+d] = v[Request[n][2]][Request[n][1]][dof*Request[n][0]+d];
      continue;
    }

    nei = std::find(neighbors_ptr->begin(), neighbors_ptr->end(), owner)
        - neighbors_ptr->begin();
    assert(nei != Nneigh);

    neighbor_requests[nei].push_back(Request[n]);
    neighbor_received_map[nei].push_back(n);
  }


  vector<vector<double> > neighbor_received(Nneigh);
  for(int p=0; p<Nneigh; p++)
    neighbor_received[p].resize(dof*neighbor_requests[p].size());


  //Step 2: Get the size of each package
  vector<int> counter(Nneigh, 0);
   
  int proc, count;
  vector<MPI_Request> send_requests;
  vector<MPI_Request> recv_requests;
  for(int p=0; p<Nneigh; p++) {
    send_requests.push_back(MPI_Request());
    proc = (*neighbors_ptr)[p];
    count = neighbor_requests[p].size();
    MPI_Isend(&count, 1, MPI_INT, proc, rank/*tag*/, comm, &send_requests.back());

    recv_requests.push_back(MPI_Request());
    MPI_Irecv(&counter[p], 1, MPI_INT, proc, proc, comm, &recv_requests.back());
  }

  MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE); //not necessary(?)
  MPI_Waitall(recv_requests.size(), recv_requests.data(), MPI_STATUSES_IGNORE);
  send_requests.clear();
  recv_requests.clear();


  //Step 3: Transfer the indices (i,j,k) of each package
  std::map<int, vector<Int3> > ExportInds;
  for(int p=0; p<Nneigh; p++) {
    if(!neighbor_requests[p].empty()) {
      send_requests.push_back(MPI_Request());
      proc = (*neighbors_ptr)[p];
      MPI_Isend((int*)neighbor_requests[p].data(), 3*neighbor_requests[p].size(), MPI_INT, proc, rank/*tag*/, 
                comm, &send_requests.back());
    }
    if(counter[p]!=0) {
      ExportInds[p] = vector<Int3>(counter[p]);
      recv_requests.push_back(MPI_Request());
      proc = (*neighbors_ptr)[p];
      MPI_Irecv((int*)ExportInds[p].data(), 3*counter[p], MPI_INT, proc, proc, 
                comm, &recv_requests.back());
    }
  }
  MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE); //not necessary(?)
  MPI_Waitall(recv_requests.size(), recv_requests.data(), MPI_STATUSES_IGNORE);
  send_requests.clear();
  recv_requests.clear();


  //Step 4: Extract the data
  std::map<int, vector<double> > Export;
  for(auto&& ex : ExportInds) {
    Export[ex.first] = vector<double>(dof*ex.second.size());
    for(int n=0; n<(int)ex.second.size(); n++) {
      Int3 &ijk(ex.second[n]);
      assert(V.IsHere(ijk[0],ijk[1],ijk[2],true)); //allow accessing thost layer (needed for outer ghosts)
      for(int d=0; d<dof; d++)
        Export[ex.first][dof*n+d] = v[ijk[2]][ijk[1]][dof*ijk[0]+d];
    }
  }
  V.RestoreDataPointerToLocalVector();
 

  //Step 5: Transfer the data
  for(auto it = Export.begin(); it != Export.end(); it++) { 
    send_requests.push_back(MPI_Request());
    proc = (*neighbors_ptr)[it->first];
    MPI_Isend(it->second.data(), it->second.size(), MPI_DOUBLE, proc, rank/*tag*/, 
              comm, &send_requests.back());
  }
  for(int p=0; p<Nneigh; p++) {
    if(neighbor_received[p].size()!=0) {
      recv_requests.push_back(MPI_Request());
      proc = (*neighbors_ptr)[p];
      MPI_Irecv(neighbor_received[p].data(), neighbor_received[p].size(), MPI_DOUBLE, proc, proc, 
                comm, &recv_requests.back());
    }
  }

  MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE); //not necessary(?)
  MPI_Waitall(recv_requests.size(), recv_requests.data(), MPI_STATUSES_IGNORE);


  //Step 6: Apply data to "Received"
  int ind;
  for(int p=0; p<Nneigh; p++)
    for(int n=0; n<(int)neighbor_received_map[p].size(); n++)
      for(int d=0; d<dof; d++) {
        ind = dof*neighbor_received_map[p][n]+d;
        assert(ind<(int)Received.size());
        Received[ind] = neighbor_received[p][dof*n+d];
      }

}

//----------------------------------------------------

void
NeighborCommunicator::Request(double*** v, int dof,
                              vector<Int3> &Request, vector<double> &Received,
                              int exchange_type) //!< 0~all, 1~face_edge, 2~face)
{

  //Sanity check
  assert(exchange_type==2 || exchange_type==1 || exchange_type==0);

  std::vector<int> *neighbors_ptr; //points to the correct neighborhood
  if(exchange_type==2)
    neighbors_ptr = subD_neighbors_face;
  else if(exchange_type==1)
    neighbors_ptr = subD_neighbors_face_edge;
  else
    neighbors_ptr = subD_neighbors_all;

  int Nneigh = neighbors_ptr->size();
  

  Received.resize(dof*Request.size());

  //Step 1: Construct array for each neighbor
  vector<vector<Int3> > neighbor_requests(Nneigh); 
  vector<vector<int> > neighbor_received_map(Nneigh); 
  int owner, nei;
  for(int n=0; n<(int)Request.size(); n++) {
    owner = global_mesh.GetOwnerOfNode(Request[n][0], Request[n][1], Request[n][2], true);
    if(owner==rank) {//weird, but fine.
      for(int d=0; d<dof; d++)
        Received[n*dof+d] = v[Request[n][2]][Request[n][1]][dof*Request[n][0]+d];
      continue;
    }

    nei = std::find(neighbors_ptr->begin(), neighbors_ptr->end(), owner)
        - neighbors_ptr->begin();
    assert(nei != Nneigh);

    neighbor_requests[nei].push_back(Request[n]);
    neighbor_received_map[nei].push_back(n);
  }


  vector<vector<double> > neighbor_received(Nneigh);
  for(int p=0; p<Nneigh; p++)
    neighbor_received[p].resize(dof*neighbor_requests[p].size());


  //Step 2: Get the size of each package
  vector<int> counter(Nneigh, 0);
   
  int proc, count;
  vector<MPI_Request> send_requests;
  vector<MPI_Request> recv_requests;
  for(int p=0; p<Nneigh; p++) {
    send_requests.push_back(MPI_Request());
    proc = (*neighbors_ptr)[p];
    count = neighbor_requests[p].size();
    MPI_Isend(&count, 1, MPI_INT, proc, rank/*tag*/, comm, &send_requests.back());

    recv_requests.push_back(MPI_Request());
    MPI_Irecv(&counter[p], 1, MPI_INT, proc, proc, comm, &recv_requests.back());
  }

  MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE); //not necessary(?)
  MPI_Waitall(recv_requests.size(), recv_requests.data(), MPI_STATUSES_IGNORE);
  send_requests.clear();
  recv_requests.clear();


  //Step 3: Transfer the indices (i,j,k) of each package
  std::map<int, vector<Int3> > ExportInds;
  for(int p=0; p<Nneigh; p++) {
    if(!neighbor_requests[p].empty()) {
      send_requests.push_back(MPI_Request());
      proc = (*neighbors_ptr)[p];
      MPI_Isend((int*)neighbor_requests[p].data(), 3*neighbor_requests[p].size(), MPI_INT, proc, rank/*tag*/, 
                comm, &send_requests.back());
    }
    if(counter[p]!=0) {
      ExportInds[p] = vector<Int3>(counter[p]);
      recv_requests.push_back(MPI_Request());
      proc = (*neighbors_ptr)[p];
      MPI_Irecv((int*)ExportInds[p].data(), 3*counter[p], MPI_INT, proc, proc, 
                comm, &recv_requests.back());
    }
  }
  MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE); //not necessary(?)
  MPI_Waitall(recv_requests.size(), recv_requests.data(), MPI_STATUSES_IGNORE);
  send_requests.clear();
  recv_requests.clear();


  //Step 4: Extract the data
  Int3& ijk0(global_mesh.subD_ijk_min[rank]);
  Int3& ijkmax(global_mesh.subD_ijk_max[rank]);

  std::map<int, vector<double> > Export;
  for(auto&& ex : ExportInds) {
    Export[ex.first] = vector<double>(dof*ex.second.size());
    for(int n=0; n<(int)ex.second.size(); n++) {
      Int3 &ijk(ex.second[n]);
      assert(ijk[0]>=ijk0[0]-1 && ijk[0]<=ijkmax[0] &&
             ijk[1]>=ijk0[1]-1 && ijk[1]<=ijkmax[1] &&
             ijk[2]>=ijk0[2]-1 && ijk[2]<=ijkmax[2]); //allow accessing thost layer (needed for outer ghosts)

      for(int d=0; d<dof; d++)
        Export[ex.first][dof*n+d] = v[ijk[2]][ijk[1]][dof*ijk[0]+d];
    }
  }
 

  //Step 5: Transfer the data
  for(auto it = Export.begin(); it != Export.end(); it++) { 
    send_requests.push_back(MPI_Request());
    proc = (*neighbors_ptr)[it->first];
    MPI_Isend(it->second.data(), it->second.size(), MPI_DOUBLE, proc, rank/*tag*/, 
              comm, &send_requests.back());
  }
  for(int p=0; p<Nneigh; p++) {
    if(neighbor_received[p].size()!=0) {
      recv_requests.push_back(MPI_Request());
      proc = (*neighbors_ptr)[p];
      MPI_Irecv(neighbor_received[p].data(), neighbor_received[p].size(), MPI_DOUBLE, proc, proc, 
                comm, &recv_requests.back());
    }
  }

  MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE); //not necessary(?)
  MPI_Waitall(recv_requests.size(), recv_requests.data(), MPI_STATUSES_IGNORE);


  //Step 6: Apply data to "Received"
  int ind;
  for(int p=0; p<Nneigh; p++)
    for(int n=0; n<(int)neighbor_received_map[p].size(); n++)
      for(int d=0; d<dof; d++) {
        ind = dof*neighbor_received_map[p][n]+d;
        assert(ind<(int)Received.size());
        Received[ind] = neighbor_received[p][dof*n+d];
      }

}

//----------------------------------------------------




//----------------------------------------------------


