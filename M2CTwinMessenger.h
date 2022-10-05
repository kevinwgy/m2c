#ifndef _M2C_TWIN_MESSENGER_H_
#define _M2C_TWIN_MESSENGER_H_

#include<mpi.h>
#include<IoData.h>
#include<Vector3D.h>


/*************************************************************
 * Class M2CTwinMessenger is responsible for communicating with
 * the M2C Twin solver in an implementation of the overset grids
 * method. 
 ************************************************************/

class M2CTwinMessenger {

  IoData &iod;

  MPI_Comm &m2c_comm; //!< This is the M2C communicator
  MPI_Comm &joint_comm; //!< This is the joint communicator of M2C and M2C Twin

  int twinning_status; //!< 0~non-existent, 1~I am the ``leader'', 2~I am the ``follower''
  int m2c_rank, m2c_size; 

  std::vector<double> temp_buffer;

public:

  M2CTwinMessenger(IoData &iod_, MPI_Comm &m2c_comm_, MPI_Comm &joint_comm_, int status_);
  ~M2CTwinMessenger();
  void Destroy();

  //! Exchange data w/ M2C Twin (called before the first time step)
  void CommunicateBeforeTimeStepping();

  //! Exchange data w/ M2C Twin (called at the first time step)
  void FirstExchange();

  //! Exchange data w/ M2C Twin (called at every time step except first and last)
  void Exchange();

  //! Exchange data w/ M2C Twin (called at the last time step)
  void FinalExchange();

};

#endif
