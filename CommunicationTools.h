/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _COMMUNICATION_TOOLS_H_
#define _COMMUNICATION_TOOLS_H_

#include<vector>
#include<mpi.h>

/*********************************************************************
 * class CommunicationTools stores tools for MPI communication
 *********************************************************************
*/

class CommunicationTools
{

public:

  template<typename T>
  static void GatherVector(MPI_Comm& comm, int gatherer, std::vector<T>& my_data,
                          std::vector<T>* all_data_ptr);

  template<typename T>
  static void AllGatherVector(MPI_Comm& comm, std::vector<T>& data); //!< directly update "data" (MPI_IN_PLACE)

  //! This function dynamically allocate and free large chunks of memnory. Should not be called frequently
  template<typename T>
  static void GatherVectorOfVectors(MPI_Comm& comm, int gatherer, std::vector<std::vector<T> >& my_data,
                                    std::vector<std::vector<T> >* all_data_ptr); 

  //! This function dynamically allocate and free large chunks of memnory. Should not be called frequently
  template<typename T>
  static void AllGatherVectorOfVectors(MPI_Comm& comm, std::vector<std::vector<T> >& data); 

};

#endif
