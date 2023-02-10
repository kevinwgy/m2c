#include <NeighborCommunicator.h>

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

  buf_all.resize(subD_neighbors_all.size());
  buf_face_edge.resize(subD_neighbors_face_edge.size());
  buf_face.resize(subD_neighbors_face.size());
}

//----------------------------------------------------

NeighborCommunicator::~NeighborCommunicator()
{ }

//----------------------------------------------------

void
NeighborCommunicator::ExchangeAll(vector<vector<double> > &Export, vector<vector<double> > &Import)
{

  //Step 1: Get the size of each package
  vector<int> counter(subD_neighbors_all.size(), 0);
   
  vector<MPI_Request> send_requests;

  int proc, count;
  for(int p=0; p<(int)subD_neighbors_all.size(); p++) {
    send_requests.push_back(MPI_Request());
    proc = subD_neighbors_all[p];
    count = Export[p].size();
    MPI_Isend(&count, 1, MPI_INT, proc, rank/*tag*/, comm, &send_requests.back());
  }
  MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);
    send_requests.clear();

 I AM HERE
}

//----------------------------------------------------


//----------------------------------------------------


