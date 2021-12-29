#include <ConcurrentProgramsHandler.h>

//---------------------------------------------------------

ConcurrentProgramsHandler::ConcurrentProgramsHandler(IoData &iod_, MPI_Comm &global_comm_, MPI_Comm &comm_)
                         : iod_concurrent(iod_.concurrent), global_comm(&global_comm_), 
                           m2c_comm(&global_comm_), M2C_Comm_Tag(-INT_MAX),
                           aeros_comm(NULL) 
{

  // check if M2C is coupled with any other programs 
  if(iod_concurrent.aeros.fsi_algo != AerosCouplingData::NONE) {
    coupled = true;
    M2C_Comm_Tag = 1; // "negotiated" with AERO-S. Note that although in AERO-S the color "FLUID_ID" is set 
                      // to 0, the MPI_Comm_split command actually takes "color+1"
  }
  else
    coupled = false;


  // -------------------------------------------
  // simultaneous operations w/ other programs 
  // -------------------------------------------
  if(coupled) {

  }

  comm_ = *comm;
}

//---------------------------------------------------------

void
ConcurrentProgramsHandler::Split(int color, int maxcolor, vector<MPI_Comm*> &c)
{





}

//---------------------------------------------------------

