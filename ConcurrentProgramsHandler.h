#ifndef _CONCURRENT_PROGRAMS_HANDLER_
#define _CONCURRENT_PROGRAMS_HANDLER_

#include<IoData.h>
#include<SpaceVariable.h>
#include<AerosMessenger.h>
#include<AerofMessenger.h>
#include<M2CTwinMessenger.h>

/***********************************************************************
* Class ConcurrentProgramsHandler is responsible for splitting the MPI
* communicator among multiple programs and sending/receiving data
* to/from other programs that are coupled with M2C.
* Notes:
* - This class is largely a wrapper in the sense that the the
*   actually communications between M2C and other programs are
*   performed in the "messenger" classes (e.g., AerosMessenger).
* - This class and the messengers only ensure proper communications.
*   They are not responsible for preparing the correct data (except for 
*   taking care of "staggering"). This should be done by other classes 
*   who owns or creates the data needed for communication.
***********************************************************************/

//! Concurrent programs handler
class ConcurrentProgramsHandler {

  IoData &iod; //!< user inputs

  bool coupled; //!< whether M2C is coupled to (i.e. running concurrently w/) any other programs

  int m2c_color; //!< the id ("color") of M2C in the MPI split
  int maxcolor; //!< the total number of "colors" (must be the same in all concurrent programs)

  MPI_Comm global_comm; //!< the global communicator
  int global_size, global_rank;

  MPI_Comm m2c_comm; //!< the communicator for M2C
  int m2c_size, m2c_rank;

  //! the communicators between m2c and each of the other programs
  std::vector<MPI_Comm> c;
 
  //! time-step size suggested by other solvers
  double dt;
  double tmax;

  //! other concurrent/coupled programs
  AerosMessenger *aeros; //!< takes care of communications w/ AERO-S
  MPI_Comm aeros_comm;  //!< this is just c[aeros_color], a communicator that includes
                        //!< M2C and AERO-S processes

  AerofMessenger *aerof;
  MPI_Comm aerof_comm;

  M2CTwinMessenger *m2c_twin;
  MPI_Comm m2c_twin_comm;
  int twinning_status; //!< 0~non-existent, 1~I am the ``leader'', 2~I am the ``follower''
  
public:

  //! The constructor calls MPI_Comm_split together with all the concurrent programs
  ConcurrentProgramsHandler(IoData &iod_, MPI_Comm global_comm_, MPI_Comm &comm_); 
  ~ConcurrentProgramsHandler();

  void InitializeMessengers(TriangulatedSurface *surf_, vector<Vec3D> *F_);
  void Destroy();

  bool Coupled() {return coupled;}
  double GetTimeStepSize() {return dt;}
  double GetMaxTime() {return tmax;}

  //! The main functions that handle communications
  void CommunicateBeforeTimeStepping(); //!< called before the 1st time step
  void FirstExchange(); //!< called at the the 1st time step
  void Exchange(); //!< called every time step (except 1st and last)
  void FinalExchange(); //!< at the last time step

private:

  void SetupCommunicators(); //!< called by the constructor

};

#endif
