#ifndef _NEIGHBOR_COMMUNICATOR_H_
#define _NEIGHBOR_COMMUNICATOR_H_

#include<mpi.h>
#include<Vector3D.h>
#include<vector>

/*******************************************************
 * Class CustomCommunicator sets up a communication
 * channel between adjacent subdomains. It allows them
 * to exchange data of arbitrary size
 * Note: The data type is assumed to be "double". Integers
 * will be represented as doubles.
 ******************************************************/

class NeighborCommunicator {

  //! MPI info
  MPI_Comm& comm;
  int rank;
  int size;

  std::vector<int>& subD_neighbors_all; //!< all real neighbors, w/o self and non-exist ones
  std::vector<int>& subD_neighbors_face_edge; //!< real neighbors, excluding corners (at most 19)
  std::vector<int>& subD_neighbors_face; //!< only real neighbors with face-contact (at most 6)  

  std::vector<std::vector<double> > buf_all;
  std::vecotr<std::vector<double> > buf_face_edge;
  std::vector<std::vector<double> > buf_face;

public:

  //! Create the communicator. Setup the buffers
  NeighborCommunicator(MPI_Comm& comm_, std::vector<int> &subD_neighbors_all_,
                       std::vector<int> &subD_neighbors_face_edge_,
                       std::vector<int> &subD_neighbors_face_);

  ~NeighborCommunicator();

  //! Get neighbor lists
  std::vector<int> &GetAllNeighbors() {return subD_neighbors_all;}
  std::vector<int> &GetFaceEdgeNeighbors() {return subD_neighbors_face_edge;}
  std::vector<int> &GetFaceNeighbors() {return subD_neighbors_face;}

  //! Data exchange 
  void ExchangeAll(std::vector<std::vector<double> > &Export, std::vector<std::vector<double> > &Import);
  void ExchangeFaceEdge(std::vector<std::vector<double> > &Export, std::vector<std::vector<double> > &Import);
  void ExchangeFace(std::vector<std::vector<double> > &Export, std::vector<std::vector<double> > &Import);
   

};


#endif
