#ifndef _LASER_ABSORPTION_SOLVER_H_
#define _LASER_ABSORPTION_SOLVER_H_

#include<SpaceVariable.h>
#include<CustomCommunicator.h>
#include<GhostPoint.h>
#include<VarFcnBase.h>
#include<tuple>

//------------------------------------------------------------
// Class LaserAbsorptionSolver is responsible for solving the
// laser radiation equation and coupling it with the Navier-
// Stokes equations.
//------------------------------------------------------------

//! Structure of nodal status
struct NodalLaserInfo {
  int i,j,k;
  double phi;
  int level;
  NodalLaserInfo(int i_, int j_, int k_, double phi_, int ql_)
    : i(i_),j(j_),k(k_),phi(phi_),level(ql_) {}
  ~NodalLaserInfo() {}
};


class LaserAbsorptionSolver {

  MPI_Comm& comm;

  //! Communication operator for each level
  vector<CustomCommunicator*> levelcomm;

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
  SpaceVariable3D Level; //!< -1 if out of scope
  SpaceVariable3D Tag;
  //Phi, Level, and Tag are undefined in the ghost layer outside the physical domain

  //! Ordering of nodes
  std::vector<NodalLaserInfo> sortedNodes;
  std::vector<int> queueCounter;  //!< number of nodes on each level of the queue
  
  //! Nodes/cells in laser domain 
  int numNodesInScope; //!< within the subdomain interior

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

private:

  void CheckForInputErrors();

  void ReadUserSpecifiedPowerFile(const char *source_power_timehistory_file);

  void CalculateDistanceToSource(SpaceVariable3D &Phi);

  static bool LaserDistCompare(NodalLaserInfo& node1, NodalLaserInfo& node2); 

  void BuildSortedNodeList();

  void BuildCustomizedCommunicators();

  void VerifyCustomizedCommunicators();

  void ResetTag();

};

#endif
