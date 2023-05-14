/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<CommunicationTools.h>
#include<Utils.h>
#include<Vector2D.h>
#include<Vector3D.h>
#include<Vector5D.h>
#include<algorithm> //copy()
#include<cassert>

//------------------------------------------------------------------------------

template<typename T>
void
CommunicationTools::GatherVector(MPI_Comm& comm, int gatherer, std::vector<T>& my_data,
                                std::vector<T>* all_data_ptr)
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
    print_error(comm, "*** Error: GatherVector detected unknown data type.\n");
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
    assert(all_data_ptr);
    all_data_ptr->resize(counter/multiplier);
    if(double_or_int) 
      MPI_Gatherv((double*)my_data.data(), my_data_count, MPI_DOUBLE, (double*)all_data_ptr->data(), counts.data(),
                  displacements.data(), MPI_DOUBLE, gatherer, comm);
    else
      MPI_Gatherv((int*)my_data.data(), my_data_count, MPI_INT, (int*)all_data_ptr->data(), counts.data(),
                  displacements.data(), MPI_INT, gatherer, comm);
  } else {//sender
    //if(all_data_ptr)
    //  all_data_ptr->clear();
    if(double_or_int) 
      MPI_Gatherv((double*)my_data.data(), my_data_count, MPI_DOUBLE, NULL, NULL,
                  NULL, MPI_DOUBLE, gatherer, comm);
    else
      MPI_Gatherv((int*)my_data.data(), my_data_count, MPI_INT, NULL, NULL,
                  NULL, MPI_DOUBLE, gatherer, comm);
  }
  
}

//------------------------------------------------------------------------------

template<typename T>
void
CommunicationTools::AllGatherVector(MPI_Comm& comm, std::vector<T>& data)
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
    print_error(comm, "*** Error: AllGatherVector detected unknown data type.\n");
    exit(-1);
  }

  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  // talk to get package size
  std::vector<int> counts(mpi_size, 0);
  counts[mpi_rank] = multiplier*data.size();
  MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, counts.data(), 1, MPI_INT, comm);

  // work
  std::vector<int> displacements(mpi_size, 0);
  int counter = 0;
  for(int i=0; i<(int)displacements.size(); i++) {
    displacements[i] = counter;
    counter += counts[i];
  }

  std::vector<T> my_data = data; //temporarily stores my data
  data.resize(counter/multiplier); //existing data outside my range will be overwritten after "Allgatherv"
  int T0 = displacements[mpi_rank]/multiplier;
  for(int i=0; i<(int)my_data.size(); i++)
    data[T0 + i] = my_data[i]; //add my data to the full data vector
  
  // talk again to exchange data
  if(double_or_int)
    MPI_Allgatherv(MPI_IN_PLACE, counts[mpi_rank], MPI_DOUBLE, (double*)data.data(), counts.data(),
                   displacements.data(), MPI_DOUBLE, comm);
  else
    MPI_Allgatherv(MPI_IN_PLACE, counts[mpi_rank], MPI_INT, (int*)data.data(), counts.data(),
                   displacements.data(), MPI_INT, comm);

}


//------------------------------------------------------------------------------

template void
CommunicationTools::GatherVector<double>(MPI_Comm& comm, int gatherer,
                                        std::vector<double>& my_data,
                                        std::vector<double>* all_data_ptr);

template void
CommunicationTools::GatherVector<Vec3D>(MPI_Comm& comm, int gatherer,
                                       std::vector<Vec3D>& my_data,
                                       std::vector<Vec3D>* all_data_ptr);

template void
CommunicationTools::AllGatherVector<double>(MPI_Comm& comm, std::vector<double>& data);

template void
CommunicationTools::AllGatherVector<Vec3D>(MPI_Comm& comm, std::vector<Vec3D>& data);

//------------------------------------------------------------------------------


