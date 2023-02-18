/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<AerofMessenger.h>
#include<cassert>
#include<climits>

using std::vector;

extern int verbose;

//---------------------------------------------------------------

AerofMessenger::AerofMessenger(IoData &iod_, MPI_Comm &m2c_comm_, MPI_Comm &joint_comm_)
              : iod(iod_), m2c_comm(m2c_comm_), joint_comm(joint_comm_)
{
//TODO


}

//---------------------------------------------------------------

AerofMessenger::~AerofMessenger()
{ }

//---------------------------------------------------------------

void
AerofMessenger::Destroy()
{ }

//---------------------------------------------------------------

void
AerofMessenger::CommunicateBeforeTimeStepping()
{

}

//---------------------------------------------------------------

void
AerofMessenger::FirstExchange()
{

}

//---------------------------------------------------------------

void
AerofMessenger::Exchange()
{

}

//---------------------------------------------------------------

void
AerofMessenger::FinalExchange()
{

}

//---------------------------------------------------------------


//---------------------------------------------------------------

