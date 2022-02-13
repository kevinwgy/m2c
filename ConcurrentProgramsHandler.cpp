#include <ConcurrentProgramsHandler.h>
#include <cassert>

//---------------------------------------------------------

ConcurrentProgramsHandler::ConcurrentProgramsHandler(IoData &iod_, MPI_Comm global_comm_, MPI_Comm &comm_)
                         : iod_concurrent(iod_.concurrent), global_comm(global_comm_), 
                           m2c_comm(global_comm_), aeros_comm(), aeros(NULL)
{

  // check if M2C is coupled with any other programs 
  int aeros_color = -1;
  if(iod_concurrent.aeros.fsi_algo != AerosCouplingData::NONE) {
    coupled = true;
    //The following parameters are the same as "FLUID_ID", "STRUCT_ID" and "MAX_CODES" in AERO-S
    m2c_color = 0; 
    aeros_color = 1;
    maxcolor = 4; 
  }
  else
    coupled = false;


  // simultaneous operations w/ other programs 
  if(coupled)
    SetupCommunicators();

  // create inter-communicators 
  if(iod_concurrent.aeros.fsi_algo != AerosCouplingData::NONE) {
    aeros_comm = c[aeros_color];
  }

  // time-step size suggested by other solvers, will be updated
  dt = 0.0;
  tmax = 0.0;

  // outputs the m2c communicator
  comm_ = m2c_comm;
}

//---------------------------------------------------------

ConcurrentProgramsHandler::~ConcurrentProgramsHandler()
{
  if(aeros) delete aeros;
}

//---------------------------------------------------------

void
ConcurrentProgramsHandler::InitializeMessengers(TriangulatedSurface *surf_, vector<Vec3D> *F_) //for AERO-S messengers
{
  if(iod_concurrent.aeros.fsi_algo != AerosCouplingData::NONE) {

    assert(surf_); //cannot be NULL
    assert(F_); //cannot be NULL
    aeros = new AerosMessenger(iod_concurrent.aeros, m2c_comm, aeros_comm, *surf_, *F_); 

    dt = aeros->GetTimeStepSize();
    tmax = aeros->GetMaxTime();
  }
}

//---------------------------------------------------------

void
ConcurrentProgramsHandler::SetupCommunicators()
{

  MPI_Comm_rank(global_comm, &global_rank);
  MPI_Comm_size(global_comm, &global_size);

  MPI_Comm_split(global_comm, m2c_color + 1, global_rank, &m2c_comm);
  MPI_Comm_rank(m2c_comm, &m2c_rank);
  MPI_Comm_size(m2c_comm, &m2c_size);

  c.resize(maxcolor);

  c[m2c_color] = m2c_comm;

  vector<int> leaders(maxcolor, -1);
  vector<int> newleaders(maxcolor, -1);

  if(m2c_rank == 0) {
    leaders[m2c_color] = global_rank;
  }
  MPI_Allreduce(leaders.data(), newleaders.data(), maxcolor, MPI_INTEGER, MPI_MAX, global_comm);

  for(int i=0; i<maxcolor; i++) {
    if(i != m2c_color && newleaders[i] >= 0) {
      // create a communicator between m2c and program i
      int tag;
      if(m2c_color < i)
        tag = maxcolor * (m2c_color + 1) + i + 1;
      else
        tag = maxcolor * (i + 1) + m2c_color + 1;

      MPI_Intercomm_create(m2c_comm, 0, global_comm, newleaders[i], tag, &c[i]);
    }
  }

}

//---------------------------------------------------------

void
ConcurrentProgramsHandler::Destroy()
{
  for(int i=0; i<c.size(); i++)
    MPI_Comm_free(&c[i]);
}

//---------------------------------------------------------

void
ConcurrentProgramsHandler::CommunicateBeforeTimeStepping()
{
  if(aeros) {
    aeros->CommunicateBeforeTimeStepping();
    dt = aeros->GetTimeStepSize();
    tmax = aeros->GetMaxTime();
  }
}

//---------------------------------------------------------

void
ConcurrentProgramsHandler::FirstExchange()
{
  if(aeros) {
    aeros->FirstExchange();
    dt = aeros->GetTimeStepSize();
    tmax = aeros->GetMaxTime();
  }
}

//---------------------------------------------------------

void
ConcurrentProgramsHandler::Exchange()
{
  if(aeros) {
    aeros->Exchange();
    dt = aeros->GetTimeStepSize();
    tmax = aeros->GetMaxTime();
  }
}

//---------------------------------------------------------

void
ConcurrentProgramsHandler::FinalExchange()
{
  if(aeros) {
    aeros->FinalExchange();
    dt = aeros->GetTimeStepSize();
    tmax = aeros->GetMaxTime();
  }
}

//---------------------------------------------------------


//---------------------------------------------------------



