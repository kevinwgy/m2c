#include <NeighborCommunicator.h>
#include <cassert>

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
NeighborCommunicator::Exchange(int exchange_type, // 0~all, 1~face_edge, 2~face
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






//----------------------------------------------------


