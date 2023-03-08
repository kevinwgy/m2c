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
  static void GatherArray(MPI_Comm& comm, int gatherer, std::vector<T>& my_data,
                          std::vector<T>& all_data);

};

#endif
