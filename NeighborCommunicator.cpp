#include <NeighborCommunicator.h>
#include <SpaceVariable.h>
#include <cassert>
#include <map>

using std::vector;

//----------------------------------------------------
NeighborCommunicator::NeighborCommunicator(MPI_Comm& comm_,
                                           vector<int> &subD_neighbors_all_,
                                           vector<int> &subD_neighbors_face_edge_,
                                           vector<int> &subD_neighbors_face_)
                    : comm(comm_),
                      subD_neighbors_all(subD_neighbors_all_),
                      subD_neighbors_face_edge(subD_neighbors_face_edge_),
                      subD_neighbors_face(subD_neighbors_face_)
{
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  assert(rank>=0 && rank<size);
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
    neighbors_ptr = &subD_neighbors_face;
  else if(exchange_type==1)
    neighbors_ptr = &subD_neighbors_face_edge;
  else
    neighbors_ptr = &subD_neighbors_all;

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
NeighborCommunicator::Request(int exchange_type, // 0~all, 1~face_edge, 2~face
                              SpaceVariable3D &V,
                              vector<vector<Int3> > &Request, vector<vector<double> > &Received)
{

  int dim = V.NumDOF();
  
  //Sanity check
  assert(exchange_type==2 || exchange_type==1 || exchange_type==0);

  std::vector<int> *neighbors_ptr; //points to the correct neighborhood
  if(exchange_type==2)
    neighbors_ptr = &subD_neighbors_face;
  else if(exchange_type==1)
    neighbors_ptr = &subD_neighbors_face_edge;
  else
    neighbors_ptr = &subD_neighbors_all;

  int Nneigh = neighbors_ptr->size();
  
  assert((int)Request.size()==Nneigh);
  assert((int)Received.size()==Nneigh);


  //Step 1: Resize Received to have the correct sizes
  for(int p=0; p<Nneigh; p++)
    Received[p].resize(dim*Request[p].size());

  //Step 2: Get the size of each package
  vector<int> counter(Nneigh, 0);
   
  int proc, count;
  vector<MPI_Request> send_requests;
  vector<MPI_Request> recv_requests;
  for(int p=0; p<Nneigh; p++) {
    send_requests.push_back(MPI_Request());
    proc = (*neighbors_ptr)[p];
    count = Request[p].size();
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
    if(!Request[p].empty()) {
      send_requests.push_back(MPI_Request());
      proc = (*neighbors_ptr)[p];
      MPI_Isend((int*)Request[p].data(), 3*Request[p].size(), MPI_INT, proc, rank/*tag*/, 
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
  double*** v = V.GetDataPointer();
  std::map<int, vector<double> > Export;
  for(auto&& ex : ExportInds) {
    Export[ex.first] = vector<double>(dim*ex.second.size());
    for(int n=0; n<(int)ex.second.size(); n++) {
      Int3 &ijk(ex.second[n]);
      assert(V.IsHere(ijk[0],ijk[1],ijk[2],false));
      for(int d=0; d<dim; d++)
        Export[ex.first][dim*n+d] = v[ijk[2]][ijk[1]][dim*ijk[0]+d];
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
    if(Received[p].size()!=0) {
      recv_requests.push_back(MPI_Request());
      proc = (*neighbors_ptr)[p];
      MPI_Irecv(Received[p].data(), Received[p].size(), MPI_DOUBLE, proc, proc, 
                comm, &recv_requests.back());
    }
  }

  MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE); //not necessary(?)
  MPI_Waitall(recv_requests.size(), recv_requests.data(), MPI_STATUSES_IGNORE);

}

//----------------------------------------------------



//----------------------------------------------------


