/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _M2C_TWIN_MESSENGER_H_
#define _M2C_TWIN_MESSENGER_H_

#include<IoData.h>
#include<GlobalMeshInfo.h>
#include<FloodFill.h>

/*************************************************************
 * Class M2CTwinMessenger is responsible for communicating with
 * the M2C Twin solver in an implementation of the overset grids
 * method. Both of the twins will activate this class.
 ************************************************************/

class M2CTwinMessenger {

  IoData &iod;

  MPI_Comm &m2c_comm; //!< This is the M2C communicator
  MPI_Comm &joint_comm; //!< This is the joint communicator of M2C and M2C Twin
  int m2c_rank, m2c_size; 

  enum TwinningStatus {LEADER = 1, FOLLOWER = 2} twinning_status; 

  SpaceVariable3D *coordinates;
  std::vector<GhostPoint> *ghost_nodes_inner;
  std::vector<GhostPoint> *ghost_nodes_outer;
  GlobalMeshInfo *global_mesh;

  std::vector<std::vector<double> > import_buffer;
  std::vector<std::vector<double> > export_buffer;

  //! For both the ``leader'' and the ``follower'', one package for each remote processor
  std::vector<std::vector<Int3> > import_nodes;
  std::vector<std::vector<GhostPoint> > export_points;
  GlobalMeshInfo global_mesh_twin;

  //! For the ``follower''
  SpaceVariable3D *TMP; //!< dim = 1
  SpaceVariable3D *TMP3; //!< dim = 3
  SpaceVariable3D *Color; //!< dim = 1
  FloodFill *floodfiller;

  double dt, tmax;

public:

  M2CTwinMessenger(IoData &iod_, MPI_Comm &m2c_comm_, MPI_Comm &joint_comm_, int status_);
  ~M2CTwinMessenger();
  void Destroy();

  //! Exchange data w/ M2C Twin (called before the first time step)
  void CommunicateBeforeTimeStepping(SpaceVariable3D &coordinates_, DataManagers3D &dms_,
                                     vector<GhostPoint> &ghost_nodes_inner_,
                                     vector<GhostPoint> &ghost_nodes_outer_,
                                     GlobalMeshInfo &global_mesh_,
                                     SpaceVariable3D &V, SpaceVariable3D &ID,
                                     std::set<Int3> &spo_frozen_nodes);

  //! Exchange data w/ M2C Twin (called at the end of the first time step)
  //! dt0 is the size of the last completed time-step
  void FirstExchange(SpaceVariable3D &V, double dt0, double tmax0);

  //! Exchange data w/ M2C Twin (called at every time step except first and last)
  //! dt0 is the size of the last completed time-step
  void Exchange(SpaceVariable3D &V, double dt0, double tmax0);

  //! Exchange data w/ M2C Twin (called at the last time step)
  void FinalExchange(SpaceVariable3D &V);

  //! Get time step size
  double GetTimeStepSize() {assert(twinning_status==FOLLOWER); return dt;}

  //! Get max time 
  double GetMaxTime() {assert(twinning_status==FOLLOWER); return tmax;}

private:

  //! Exchange data between M2C twins
  void ExchangeData(SpaceVariable3D &V);

  //! Interpolate and transfer
  void InterpolateDataAndTransfer(SpaceVariable3D &V, TwinningStatus sender);

};

#endif
