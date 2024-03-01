/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _AEROF_MESSENGER_H_
#define _AEROF_MESSENGER_H_

#include<mpi.h>
#include<IoData.h>
#include<Vector3D.h>


/*************************************************************
 * Class AerofMessenger is responsible for communicating with
 * the AERO_F solver in an implementation of the overset grids
 * method. 
 ************************************************************/

class AerofMessenger {

  IoData &iod;

  MPI_Comm &m2c_comm; //!< This is the M2C communicator
  MPI_Comm &joint_comm; //!< This is the joint communicator of M2C and AERO-F

  int m2c_rank, m2c_size; 

  std::vector<double> temp_buffer;

public:

  AerofMessenger(IoData &iod_, MPI_Comm &m2c_comm_, MPI_Comm &joint_comm_);
  ~AerofMessenger();
  void Destroy();

  //! Exchange data w/ AERO-F (called before the first time step)
  void CommunicateBeforeTimeStepping();

  //! Exchange data w/ AERO-F (called at the first time step)
  void FirstExchange();

  //! Exchange data w/ AERO-F (called at every time step except first and last)
  void Exchange();

  //! Exchange data w/ AERO-F (called at the last time step)
  void FinalExchange();

};

#endif
