#include<M2CTwinMessenger.h>
#include<cassert>
#include<climits>

using std::vector;

extern int verbose;

//---------------------------------------------------------------

M2CTwinMessenger::M2CTwinMessenger(IoData &iod_, MPI_Comm &m2c_comm_, MPI_Comm &joint_comm_,
                                   int status_)
              : iod(iod_), m2c_comm(m2c_comm_), joint_comm(joint_comm_),
                twinning_status(status_)
{
//TODO


}

//---------------------------------------------------------------

M2CTwinMessenger::~M2CTwinMessenger()
{ }

//---------------------------------------------------------------

void
M2CTwinMessenger::Destroy()
{ }

//---------------------------------------------------------------

void
M2CTwinMessenger::CommunicateBeforeTimeStepping()
{

}

//---------------------------------------------------------------

void
M2CTwinMessenger::FirstExchange()
{

}

//---------------------------------------------------------------

void
M2CTwinMessenger::Exchange()
{

}

//---------------------------------------------------------------

void
M2CTwinMessenger::FinalExchange()
{

}

//---------------------------------------------------------------


//---------------------------------------------------------------

