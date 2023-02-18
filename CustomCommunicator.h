/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _CUSTOM_COMMUNICATOR_H_
#define _CUSTOM_COMMUNICATOR_H_

#include<SpaceVariable.h>
#include<Vector3D.h>

/*******************************************************
 * Class CustomCommunicator sets up a communication
 * channel for SpaceVariable3D, in which only
 * a pre-specified subset of ghost nodes/cells 
 * receive information from the neighbor(s) that own 
 * these nodes/cells. In other words, only a subset of 
 * ghost nodes/cells participate in the communication.
 * Note: The data type is assumed to be "double". Integers
 * will be represented as doubles.
 ******************************************************/

class CustomCommunicator {

  //! MPI info
  MPI_Comm& comm;
  int rank;
  int size;

  //! Variable type
  int dof;
  int ghost_width;

  //! A nested structure that stores the package to/from a neighbor. Everything except data is fixed.
  struct Package {
    enum Type {SEND = 0, RECEIVE = 1} type;
    int rank; //rank of the package sender/receiver
    std::vector<double> buffer;
    std::vector<Int3> index; //location to extract/apply the data
    Package(Type type_, int rank_ = -1) : type(type_), rank(rank_) {}
    ~Package() {}
  };

  std::vector<Package> send_pack;
  std::vector<Package> recv_pack;

public:

  //! Create the communicator. ghost_nodes: list of ghost nodes that participate in the communication
  CustomCommunicator(MPI_Comm& comm_, SpaceVariable3D &V, std::vector<Int3> &ghost_nodes);

  //! Create a new communicator that is the same as the one given, except different "dof" and "ghost_width"
  //! For safety, it is required that the new ghost_width is equal to or larger than that of the given comm.
  CustomCommunicator(const CustomCommunicator &cc, int dof_, int ghost_width_);

  ~CustomCommunicator();

  //! Data exchange (Passing "v" reduces the cost by <20% (i.e. not much) compared to passing "V")
  void ExchangeAndInsert(SpaceVariable3D &V);
  void ExchangeAndInsert(double*** v); //note: this function does not verify compatibility

  //! Can easily add other types of communication (e.g., add, max, etc.)

  //! Get info
  int NumDOF() {return dof;}
  int NumGhostLayers() {return ghost_width;}

};


#endif
