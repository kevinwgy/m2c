/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<CommunicationTools.h>
#include<Utils.h>
#include<Vector2D.h>
#include<Vector3D.h>
#include<Vector5D.h>

//------------------------------------------------------------------------------

template<typename T>
void
CommunicationTools::GatherArray(MPI_Comm& comm, int gatherer, std::vector<T>& my_data,
                                std::vector<T>& all_data)
{
  // figure out the data type
  int multiplier = 1;
  bool double_or_int = true;
  if(std::is_same<T, double>::value) {
    /*nothing*/
  } else if(std::is_same<T, Vec2D>::value) {
    multiplier = 2;
  } else if(std::is_same<T, Vec3D>::value) {
    multiplier = 3;
  } else if(std::is_same<T, Vec5D>::value) {
    multiplier = 5;
  } else if(std::is_same<T, int>::value) {
    double_or_int = false;
  } else if(std::is_same<T, Int2>::value) {
    double_or_int = false;
    multiplier = 2;
  } else if(std::is_same<T, Int3>::value) {
    double_or_int = false;
    multiplier = 3;
  } else {
    print_error(comm, "*** Error: gather_array detected unknown data type.\n");
    exit(-1);
  }

  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  // talk to get package size
  std::vector<int> counts;
  int my_data_count = multiplier*my_data.size();
  if(mpi_rank == gatherer) { //I am the receiver
    counts.resize(mpi_size,-1);
    MPI_Gather(&my_data_count, 1, MPI_INT, counts.data(), 1, MPI_INT, gatherer, comm);
  } else //I am a sender
    MPI_Gather(&my_data_count, 1, MPI_INT, NULL, 1, MPI_INT, gatherer, comm);

  // work
  std::vector<int> displacements;
  int counter = 0;
  if(mpi_rank == gatherer) {
    displacements.resize(mpi_size,-1);
    for(int i=0; i<(int)displacements.size(); i++) {
      displacements[i] = counter;
      counter += counts[i];
    }
  }

  // talk again
  if(mpi_rank == gatherer) { //receiver
    all_data.resize(counter/multiplier);
    if(double_or_int) 
      MPI_Gatherv((double*)my_data.data(), my_data_count, MPI_DOUBLE, (double*)all_data.data(), counts.data(),
                  displacements.data(), MPI_DOUBLE, gatherer, comm);
    else
      MPI_Gatherv((int*)my_data.data(), my_data_count, MPI_INT, (int*)all_data.data(), counts.data(),
                  displacements.data(), MPI_INT, gatherer, comm);
  } else {//sender
    all_data.clear();
    if(double_or_int) 
      MPI_Gatherv((double*)my_data.data(), my_data_count, MPI_DOUBLE, NULL, NULL,
                  NULL, MPI_DOUBLE, gatherer, comm);
    else
      MPI_Gatherv((int*)my_data.data(), my_data_count, MPI_INT, NULL, NULL,
                  NULL, MPI_DOUBLE, gatherer, comm);
  }
  
}

//------------------------------------------------------------------------------

template void
CommunicationTools::GatherArray<double>(MPI_Comm& comm, int gatherer,
                                        std::vector<double>& my_data,
                                        std::vector<double>& all_data);

//------------------------------------------------------------------------------


