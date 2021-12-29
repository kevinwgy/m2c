#ifndef _CONCURRENT_PROGRAMS_HANDLER_
#define _CONCURRENT_PROGRAMS_HANDLER_

#include<IoData.h>
#include<SpaceVariable.h>
#include<mpi.h>
//--------------------------------------------------------------
// Class ConcurrentProgramsHandler is responsible for splitting the MPI
// communicator among multiple programs and sending/receiving data
// to/from other programs that are coupled with M2C
//--------------------------------------------------------------

//! Concurrent programs handler
class ConcurrentProgramsHandler {

  ConcurrentProgramsData &iod_concurrent; //!< user inputs

  bool coupled; //!< whether M2C is coupled to (i.e. running concurrently w/) any other programs

  int M2C_Comm_Tag; //!< the id ("color") of M2C in the MPI split

  MPI_Comm* global_comm; //!< the global communicator
  MPI_Comm* m2c_comm; //!< the communicator for M2C

  AerosMessenger *aeros; //!< takes care of communications w/ AERO-S
  MPI_Comm* aeros_comm;

  //! TODO: other software/programs can be added later
  
public:

  //! The constructor calls MPI_Comm_split together with all the concurrent programs
  ConcurrentProgramsHandler(IoData &iod_, MPI_Comm &global_comm_, MPI_Comm &comm_); 
  ~ConcurrentProgramsHandler();

  void Destroy();

  void Init1(double *in, double *out);
  void Init2(double *in, double *out);

  void Update1(double *in, double *out);
  void Update2(double *in, double *out);

  void Finalize1(double *in, double *out);
  void Finalize2(double *in, double *out);

private:

  void Split(int color, int maxcolor, vector<MPI_Comm*> &c);

};

#endif
