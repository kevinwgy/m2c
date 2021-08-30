#ifndef _LASER_ABSORPTION_SOLVER_H_
#define _LASER_ABSORPTION_SOLVER_H_

#include<SpaceVariable.h>
#include<GhostPoint.h>
#include<VarFcnBase.h>
#include<tuple>

class LaserAbsorptionSolver {

  //! Variable fcn
  vector<VarFcnBase*>& varFcn;

  //! IoData
  IoData &iod;

  //! Mesh info
  SpaceVariable3D &coordinates;
  SpaceVariable3D &delta_xyz;
  SpaceVariable3D &volume;

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain

  //! Pointer to ghost node lists (outside domain boundary)
  vector<GhostPoint> &ghost_nodes_inner;
  vector<GhostPoint> &ghost_nodes_outer;

  //! Internal variables
  SpaceVariable3D L0; //!< ``old'' L
  SpaceVariable3D FluxIn;
  SpaceVariable3D FluxOut;
  SpaceVariable3D Phi; //!< distance from each node to source plane; -1 if node is out of scope
  SpaceVariable3D Tag;

  //! Ordering of nodes
  std::vector<std::tuple<Int3,double,int> > sortedNodes;  //!< node ind, phi, queue-level, order-on-the-level
  std::vector<int> queueCounter;  //!< number of nodes on each level of the queue


  //! Absorption coefficients
  std::vector<std::pair<double,double> > abs;   // !< pair<slope,alpha0> for each fluid model

  //! stores time history of source power (if provided by user)
  std::vector<std::pair<double,double> > source_power_timehistory;


public:

  LaserAbsorptionSolver(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_, std::vector<VarFcnBase*> &varFcn_,
                        SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_, SpaceVariable3D &volume_,
                        std::vector<GhostPoint> &ghost_nodes_inner_, std::vector<GhostPoint> &ghost_nodes_outer_);

  ~LaserAbsorptionSolver();

  void Destroy();

  void ComputeLaserRadiance(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &L);

};

#endif
