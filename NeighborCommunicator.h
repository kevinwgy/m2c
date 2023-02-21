/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _NEIGHBOR_COMMUNICATOR_H_
#define _NEIGHBOR_COMMUNICATOR_H_

#include<GlobalMeshInfo.h>

class SpaceVariable3D;

/*******************************************************
 * Class NeighborCommunicator sets up a communication
 * channel between adjacent subdomains. It allows them
 * to exchange data of arbitrary size
 * Note: The data type is assumed to be "double". Integers
 * will have to be represented as doubles.
 * Note: Some of the functions can also be implemented
 * using one-sided communication / remote memory access
 * (RMA), which might be more efficient. However, this
 * is not done here. We are still using (non-blocking)
 * two-sided communications --- for clarity (and also
 * robustness?). The communication cost should be trivial
 * compared to the total computation time regardless of
 * one-sided or two-sided comm. 
 ******************************************************/

class NeighborCommunicator {

  //! MPI info
  MPI_Comm& comm;
  int rank;
  int size;

  //! Global mesh info
  GlobalMeshInfo& global_mesh;

  //! Neighbor info pertaining to the current subdomain
  std::vector<int>* subD_neighbors_all; //!< all real neighbors, w/o self and non-exist ones
  std::vector<int>* subD_neighbors_face_edge; //!< real neighbors, excluding corners (at most 19)
  std::vector<int>* subD_neighbors_face; //!< only real neighbors with face-contact (at most 6)  

public:

  //! Create the communicator. Setup the buffers
  NeighborCommunicator(MPI_Comm& comm_, GlobalMeshInfo& global_mesh_);

  ~NeighborCommunicator();

  //! Get neighbor lists
  std::vector<int> &GetAllNeighbors() {return *subD_neighbors_all;}
  std::vector<int> &GetFaceEdgeNeighbors() {return *subD_neighbors_face_edge;}
  std::vector<int> &GetFaceNeighbors() {return *subD_neighbors_face;}

  //! Data exchange (Import[p] will be resized to exactly the size of data passed to it)
  void Send(int exchange_type, //!< 0~all, 1~face_edge, 2~face
            std::vector<std::vector<double> > &Export, std::vector<std::vector<double> > &Import);
   
  void Request(SpaceVariable3D &V, std::vector<Int3> &Request, std::vector<double> &Received,
               int exchange_type = 0); //!< 0~all, 1~face_edge, 2~face);

  void Request(double*** v, int dof,
               std::vector<Int3> &Request, std::vector<double> &Received,
               int exchange_type = 0); //!< 0~all, 1~face_edge, 2~face); //!< V --> v
   
};


#endif
