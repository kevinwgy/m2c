#ifndef _LASER_ABSORPTION_SOLVER_H_
#define _LASER_ABSORPTION_SOLVER_H_

#include<CustomCommunicator.h>
#include<GhostPoint.h>
#include<VarFcnBase.h>
#include<EmbeddedBoundaryFormula.h>

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

//! Laser source info
struct SourceGeometry {
  double radius;//!< radius on the source "plane" (user-specified)
  double angle; //!< radians
  Vec3D  x0;    //!< center of source "plane"
  Vec3D  dir;   //!< central axis 
  double R;     //!< radius of curvature (infinite for parallel laser)
  Vec3D  O;     //!< origin of the sector (relevant only for angle != 0)
};

//! Embedded boundary method
struct LaserEBM{

  std::vector<std::pair<Int3, EmbeddedBoundaryFormula> > ghostNodes1; //ghost nodes whose image is inside this subdomain

  std::vector<std::vector<Int3> > ghostNodes2; //ghost nodes whose image is inside another subdomain
  std::vector<int> ghostNodes2_sender;
  
  std::vector<std::vector<std::pair<Int3, EmbeddedBoundaryFormula> > > friendsGhostNodes;
  std::vector<int> friendsGhostNodes_receiver;

};


//! Laser absorption solver
class LaserAbsorptionSolver {

  MPI_Comm& comm;
  int mpi_rank;
  int mpi_size;

  //! Communication operator for each level
  vector<CustomCommunicator*> levelcomm;

  //! Variable fcn
  vector<VarFcnBase*>& varFcn;

  //! IoData
  IoData &iod;

  //! Laser info
  SourceGeometry source;
  double laser_range; //!< range

  //! Mesh info
  SpaceVariable3D &coordinates;
  SpaceVariable3D &delta_xyz;
  SpaceVariable3D &volume;

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain

  Vec3D dxmin_glob;
  Vec3D xmin_glob, xmax_glob; //coords of the first and last nodes inside the physical domain
  Vec3D xmin_glob_ghost, xmax_glob_ghost; //coords of the first and last nodes in the ghost boundary layer

  //! Pointer to ghost node lists (outside domain boundary)
  vector<GhostPoint> &ghost_nodes_inner;
  vector<GhostPoint> &ghost_nodes_outer;

  //! Internal variables
  SpaceVariable3D L0; //!< ``old'' L
  SpaceVariable3D FluxIn;
  SpaceVariable3D FluxOut;
  SpaceVariable3D Phi; //!< distance from each node to source; -1 if node is out of scope
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

  //! embedded boundary method (for applying embedded and non-embedded boundary conditions)
  LaserEBM ebm;

public:

  LaserAbsorptionSolver(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_, std::vector<VarFcnBase*> &varFcn_,
                        SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_, SpaceVariable3D &volume_,
                        std::vector<GhostPoint> &ghost_nodes_inner_, std::vector<GhostPoint> &ghost_nodes_outer_);

  ~LaserAbsorptionSolver();

  void Destroy();

  void SetSourceRadiance(SpaceVariable3D &L, double t);

  void ComputeLaserRadiance(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &L);

private:

  void CheckForInputErrors();

  void CalculateGlobalMeshInfo();

  void ReadUserSpecifiedPowerFile(const char *source_power_timehistory_file);

  void CalculateLaserInfo();

  double DistanceToSource(Vec3D &x, int &loc, double *image = NULL);

  int FindImagePoint(Vec3D &x, Vec3D &image);

  void CalculateDistanceToSource(SpaceVariable3D &Phi);

  static bool LaserDistCompare(NodalLaserInfo& node1, NodalLaserInfo& node2); 

  void BuildSortedNodeList();

  void BuildCustomizedCommunicators();
  //! for debug purpose only
  void VerifySortedNodesAndCommunicators();

  void SetupLaserGhostNodes();

  void UpdateGhostNodes(double ***l);

  void ResetTag();

};

#endif
