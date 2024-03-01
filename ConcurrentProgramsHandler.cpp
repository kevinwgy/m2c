/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <ConcurrentProgramsHandler.h>
#include <cassert>

//---------------------------------------------------------

ConcurrentProgramsHandler::ConcurrentProgramsHandler(IoData &iod_, MPI_Comm global_comm_, MPI_Comm &comm_)
                         : iod(iod_), global_comm(global_comm_), 
                           m2c_comm(global_comm_), aeros_comm(), aeros(NULL),
                           aerof_comm(), aerof(NULL), m2c_twin_comm(), m2c_twin(NULL)
{
  coupled = false; //default

  // check if M2C is coupled with any other programs 
  int aeros_color    = -1;
  int aerof_color    = -1;
  int m2c_twin_color = -1;

  // If coupled with another instantiation of M2C (through overset grids), figure out my role
  twinning_status = NONE; //0~non-existent, 1~I am the ``leader'', 2~I am the ``follower''
  if(iod.concurrent.m2c_twin.type == M2CTwinningData::OVERSET_GRIDS) {
    // IMPORTANT: If mesh has "OVERSET" boundaries --> "leader"; otherwise --> follower
    if(iod_.mesh.bc_x0 == MeshData::OVERSET || iod_.mesh.bc_xmax == MeshData::OVERSET ||
       iod_.mesh.bc_y0 == MeshData::OVERSET || iod_.mesh.bc_ymax == MeshData::OVERSET ||
       iod_.mesh.bc_z0 == MeshData::OVERSET || iod_.mesh.bc_zmax == MeshData::OVERSET)
      twinning_status = LEADER;
    else
      twinning_status = FOLLOWER;
  }

  if(twinning_status == FOLLOWER) { 
    coupled = true;
    m2c_color = 2; // my color
    m2c_twin_color = 0; // the leader's color
    maxcolor = 4;

    // Because M2C cannot be coupled with M2C-follower and AERO-F at the same time, AERO-F and the
    // M2C-follower share the same ``color'' (2).

    //A follower to another M2C instantiation cannot be coupled with others
    if(iod.concurrent.aeros.fsi_algo != AerosCouplingData::NONE ||
       iod.concurrent.aerof.type     != AerofCouplingData::NONE) {
      fprintf(stdout,"\033[0;31m*** Error: A follower M2C instantiation cannot be coupled with "
                     "other software.\033[0m\n");
      exit(-1);
    }
  }
  // Either not coupled with anything, or if coupled, this is the leader M2C instantiation
  else if(iod.concurrent.aeros.fsi_algo != AerosCouplingData::NONE ||
          iod.concurrent.aerof.type     != AerofCouplingData::NONE ||
          iod.concurrent.m2c_twin.type  != M2CTwinningData::NONE) {
    // Within the family... Common codes allow multiple (more than 2) solvers coupled together
    coupled = true;
    //The following parameters are the same as "FLUID_ID" and "MAX_CODES" in AERO-S and AERO-F
    m2c_color = 0; // my color
    maxcolor = 4; 

    if(iod.concurrent.aeros.fsi_algo != AerosCouplingData::NONE)
      aeros_color = 1; //"STRUCT_ID" in AERO-S
    if(iod.concurrent.aerof.type     != AerofCouplingData::NONE)
      aerof_color = 2; //"AEROF_ID_FOR_M2C" in AERO-F
    if(iod.concurrent.m2c_twin.type  != M2CTwinningData::NONE) {
      assert(twinning_status = LEADER);
      m2c_twin_color = 2;
    }

    if(iod.concurrent.m2c_twin.type != M2CTwinningData::NONE &&
       iod.concurrent.aerof.type    != AerofCouplingData::NONE) {
      fprintf(stdout,"\033[0;31m*** Error: Cannot be coupled with AERO-F and another instantiation of M2C"
                     "at the same time.\033[0m\n");
      exit(-1);
    }
  }

  // simultaneous operations w/ other programs 
  if(coupled)
    SetupCommunicators();

  // create inter-communicators 
  if(iod.concurrent.aeros.fsi_algo != AerosCouplingData::NONE) {
    aeros_comm = c[aeros_color];
    int aeros_size(-1);
    MPI_Comm_size(aeros_comm, &aeros_size);
    assert(aeros_size>0);
  }
  if(iod.concurrent.aerof.type != AerofCouplingData::NONE) {
    aerof_comm = c[aerof_color];
    int aerof_size(-1);
    MPI_Comm_size(aerof_comm, &aerof_size);
    assert(aerof_size>0);
  }
  if(iod.concurrent.m2c_twin.type != M2CTwinningData::NONE) {
    m2c_twin_comm = c[m2c_twin_color];
    int m2c_twin_size(-1);
    MPI_Comm_size(m2c_twin_comm, &m2c_twin_size);
    assert(m2c_twin_size>0);
  }

  // time-step size suggested by other solvers, will be updated
  dt = -1.0;
  tmax = -1.0;

  // outputs the m2c communicator
  comm_ = m2c_comm;
}

//---------------------------------------------------------

ConcurrentProgramsHandler::~ConcurrentProgramsHandler()
{
  if(aeros) delete aeros;
  if(aerof) delete aerof;
  if(m2c_twin) delete m2c_twin;
}

//---------------------------------------------------------

void
ConcurrentProgramsHandler::InitializeMessengers(TriangulatedSurface *surf_, vector<Vec3D> *F_) //for AERO-S messengers
{
  if(iod.concurrent.aeros.fsi_algo != AerosCouplingData::NONE) {

    assert(surf_); //cannot be NULL
    assert(F_); //cannot be NULL
    aeros = new AerosMessenger(iod.concurrent.aeros, m2c_comm, aeros_comm, *surf_, *F_); 

    dt = aeros->GetTimeStepSize();
    tmax = aeros->GetMaxTime();
  }

  if(iod.concurrent.aerof.type != AerofCouplingData::NONE) {
    aerof = new AerofMessenger(iod, m2c_comm, aerof_comm);
  }

  if(iod.concurrent.m2c_twin.type != M2CTwinningData::NONE) {
    m2c_twin = new M2CTwinMessenger(iod, m2c_comm, m2c_twin_comm, twinning_status);
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
  assert(m2c_rank<m2c_size); //rank must be 0 -- (size-1)

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
  if(aeros)
    aeros->Destroy();

  if(aerof)
    aerof->Destroy();

  if(m2c_twin)
    m2c_twin->Destroy();
/*
  for(int i=0; i<(int)c.size(); i++)
    MPI_Comm_free(&c[i]);
*/
}

//---------------------------------------------------------

void
ConcurrentProgramsHandler::CommunicateBeforeTimeStepping(SpaceVariable3D *coordinates_,
                               DataManagers3D *dms_,
                               std::vector<GhostPoint> *ghost_nodes_inner_,
                               std::vector<GhostPoint> *ghost_nodes_outer_,
                               GlobalMeshInfo *global_mesh_, SpaceVariable3D *V,
                               SpaceVariable3D *ID, std::set<Int3> *spo_frozen_nodes)
{
  if(aeros) {
    aeros->CommunicateBeforeTimeStepping();
    dt = aeros->GetTimeStepSize();
    tmax = aeros->GetMaxTime();
  }

  if(aerof) {
    assert(coordinates_);
    assert(dms_);
    assert(ghost_nodes_inner_);
    assert(ghost_nodes_outer_);
    assert(global_mesh_);
    assert(ID);
    assert(spo_frozen_nodes);
    aerof->CommunicateBeforeTimeStepping();
  }

  if(m2c_twin) {
    assert(coordinates_);
    assert(dms_);
    assert(ghost_nodes_inner_);
    assert(ghost_nodes_outer_);
    assert(global_mesh_);
    assert(V);
    assert(ID);
    assert(spo_frozen_nodes);
    m2c_twin->CommunicateBeforeTimeStepping(*coordinates_, *dms_, *ghost_nodes_inner_,
                                            *ghost_nodes_outer_, *global_mesh_, 
                                            *V, *ID, *spo_frozen_nodes);
    if(twinning_status==FOLLOWER) { 
      // in terms of time-stepping, the follower is one-step behind the leader --> O(dt) error
      dt = m2c_twin->GetTimeStepSize();
      tmax = m2c_twin->GetMaxTime();
    }
  }

}

//---------------------------------------------------------

void
ConcurrentProgramsHandler::FirstExchange(SpaceVariable3D *V, double dt0, double tmax0)
{
  if(aeros) {
    aeros->FirstExchange();
    dt = aeros->GetTimeStepSize();
    tmax = aeros->GetMaxTime();
  }

  if(aerof) {
    assert(V);
    assert(dt0>=0.0);
    assert(tmax0>=0.0);
    aerof->FirstExchange();
  }

  if(m2c_twin) {
    if(twinning_status==LEADER) {
      assert(dt0>=0.0);
      assert(tmax0>=0.0);
    }
    assert(V);
    m2c_twin->FirstExchange(*V, dt0, tmax0);
    if(twinning_status==FOLLOWER) {
      dt = m2c_twin->GetTimeStepSize();
      tmax = m2c_twin->GetMaxTime();
    }
  }

}

//---------------------------------------------------------

void
ConcurrentProgramsHandler::Exchange(SpaceVariable3D *V, double dt0, double tmax0)
{
  if(aeros) {
    aeros->Exchange();
    dt = aeros->GetTimeStepSize();
    tmax = aeros->GetMaxTime();
  }

  if(aerof) {
    assert(V);
    assert(dt0>=0.0);
    assert(tmax0>=0.0);
    aerof->Exchange();
  }

  if(m2c_twin) {
    if(twinning_status==LEADER) {
      assert(dt0>=0.0);
      assert(tmax0>=0.0);
    }
    assert(V);
    m2c_twin->Exchange(*V, dt0, tmax0);
    if(twinning_status==FOLLOWER) {
      dt = m2c_twin->GetTimeStepSize();
      tmax = m2c_twin->GetMaxTime();
    }
  }

}

//---------------------------------------------------------

void
ConcurrentProgramsHandler::FinalExchange(SpaceVariable3D *V)
{
  if(aeros) {
    aeros->FinalExchange();
    dt = aeros->GetTimeStepSize();
    tmax = aeros->GetMaxTime();
  }

  if(aerof) {
    aerof->FinalExchange();
  }

  if(m2c_twin) {
    assert(V);
    m2c_twin->FinalExchange(*V);
  }

}

//---------------------------------------------------------


//---------------------------------------------------------



