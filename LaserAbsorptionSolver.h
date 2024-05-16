/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _LASER_ABSORPTION_SOLVER_H_
#define _LASER_ABSORPTION_SOLVER_H_

#include<CustomCommunicator.h>
#include<GhostPoint.h>
#include<Vector5D.h>
#include<VarFcnBase.h>
#include<EmbeddedBoundaryFormula.h>
#include<SpaceOperator.h>
#include<MeshMatcher.h>
#include<tuple>

/**************************************************************
* Class LaserAbsorptionSolver is responsible for solving the
* laser radiation equation and coupling it with the Navier-
* Stokes equations.
*
* Note: The term "radiance" in this class is actually irradiance
*
**************************************************************/

//! Structure of nodal status
struct NodalLaserInfo {
  int i,j,k;
  double phi;
  int level;
  double fin, fout; //fluxes
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
  double depth; //!< depth of the layer in front of source that will get source radiance
  inline void GetDirection(Vec3D& x, Vec3D& d) { //!< calc. laser dir at x (assuming x is inside/near laser domain)
    if(angle==0) {d = dir; return;}
    if(angle>0)  {d = O - x; d /= d.norm(); return;}
    if(angle<0)  {d = x - O; d /= d.norm(); return;}
  }
};

//! Embedded boundary method
struct LaserEBM{

  std::vector<std::pair<Int3, EmbeddedBoundaryFormula> > ghostNodes1; //ghost nodes whose image is inside this subdomain

  std::vector<std::vector<Int3> > ghostNodes2; //ghost nodes whose image is inside another subdomain
  std::vector<int> ghostNodes2_sender;
  
  std::vector<std::vector<std::pair<Int3, EmbeddedBoundaryFormula> > > friendsGhostNodes;
  std::vector<int> friendsGhostNodes_receiver;

  // laser radiance at ghost nodes
  std::vector<double> l1; //ghostNodes1
  std::vector<std::vector<double> > l2; //ghostNodes2
};


//! Laser absorption solver
class LaserAbsorptionSolver {

  //! Laser MPI (may or may not be the same as the N-S MPI)
  bool active_core;  //!< whether this core is used for laser computation
  MPI_Comm comm; //!< creating a new copy (a small object anyway)
  int mpi_rank;
  int mpi_size;

  //! N-S MPI (can be different from comm only if iod.laser.parallel == BALANCED)
  MPI_Comm& nscomm;

  //! Transfer data between the two mesh partitions (used when iod.laser.parallel==BALANCED)
  MeshMatcher* NS2Laser;
  MeshMatcher* Laser2NS;
  DataManagers3D *dms;
  SpaceOperator *spo;
  
  //! Cutoff radiance (L must be non-negative)
  double lmin;

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
  SpaceVariable3D *coordinates;
  SpaceVariable3D *delta_xyz;
  SpaceVariable3D *volume;

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain

  Vec3D dxmin_glob;
  Vec3D xmin_glob, xmax_glob; //coords of the first and last nodes inside the physical domain
  Vec3D xmin_glob_ghost, xmax_glob_ghost; //coords of the first and last nodes in the ghost boundary layer

  //! Pointer to ghost node lists (outside domain boundary)
  vector<GhostPoint> *ghost_nodes_inner;
  vector<GhostPoint> *ghost_nodes_outer;

  //! Internal variable on the NS mesh
  SpaceVariable3D TemperatureNS;

  //! Internal variables (on the Laser mesh)
  SpaceVariable3D Temperature;
  SpaceVariable3D *ID; //!< material id 
  SpaceVariable3D *L;  //!< latest L
  bool L_initialized;  //!< whether L has been initialized
  SpaceVariable3D L0; //!< ``old'' L
  SpaceVariable3D Lbk; //!< backup of input L (for "failsafe")
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
  std::vector<std::tuple<double,double,double> > absorption;   // !< <slope,T0(Kelvin),eta0> for each fluid model

  //! stores time history of source power (if provided by user)
  std::vector<std::pair<double,double> > source_power_timehistory;

  //! embedded boundary method (for applying embedded and non-embedded boundary conditions)
  LaserEBM ebm;

  //! parameters in the mean flux method and SOR (Set to input values, adjusted only in "failsafe")
  double mfm_alpha;
  double sor_relax;

public:

  LaserAbsorptionSolver(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_, std::vector<VarFcnBase*> &varFcn_,
                        SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_, SpaceVariable3D &volume_,
                        std::vector<GhostPoint> &ghost_nodes_inner_, std::vector<GhostPoint> &ghost_nodes_outer_);

  //! This constructor is used when load balancing is requested (i.e. a sub-mesh is constructed/partitioned for the laser domain)
  LaserAbsorptionSolver(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_, std::vector<VarFcnBase*> &varFcn_,
                        SpaceVariable3D &coordinates_, FluxFcnBase &fluxFcn_, ExactRiemannSolverBase &riemann_, 
                        vector<double> &x, vector<double> &y, vector<double> &z,
                        vector<double> &dx, vector<double> &dy, vector<double> &dz);
                        
  ~LaserAbsorptionSolver();

  void Destroy();

  //! Compute raser radiance L 
  void ComputeLaserRadiance(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &L, 
                            const double t, int time_step); 

  //! Compute eta*L and add to the 5th entry of R. (Not multiplying cell volume)
  void AddHeatToNavierStokesResidual(SpaceVariable3D &R, SpaceVariable3D &L, SpaceVariable3D &ID, 
                                     SpaceVariable3D *V = NULL); //if NULL, use stored temperature

  inline double GetAbsorptionCoefficient(double T, int id) { //T must be in Kelvin
    return id<(int)absorption.size() ?
             std::get<0>(absorption[id])*(T - std::get<1>(absorption[id])) + std::get<2>(absorption[id])
           : 0.0; //could get here if there are "inactive" nodes
  }
  
private:

  void CheckForInputErrors();

  void CalculateGlobalMeshInfo();

  void CalculateLaserInfo();

  double SpecifySourceDepth(vector<double> *x_ptr = NULL, vector<double> *y_ptr = NULL,
                            vector<double> *dx_ptr = NULL, vector<double> *dy_ptr = NULL);

  //! Create/partion a new sub-mesh for laser calculation
  void SetupLoadBalancing(SpaceVariable3D &coordinates_, FluxFcnBase &fluxFcn_, ExactRiemannSolverBase &riemann_,
                          vector<double> &x, vector<double> &y, vector<double> &z,
                          vector<double> &dx, vector<double> &dy, vector<double> &dz);

  void PopulateLaserMesh(SpaceVariable3D &V_, SpaceVariable3D &ID_, SpaceVariable3D &L_);

  void PopulateRadianceOnNavierStokesMesh(SpaceVariable3D &L_);

  void ReadUserSpecifiedPowerFile(const char *source_power_timehistory_file);

  double DistanceToSource(Vec3D &x, int &loc, double *image = NULL);

  int FindImagePoint(Vec3D &x, Vec3D &image);

  void CalculateDistanceToSource(SpaceVariable3D &Phi);

  static bool LaserDistCompare(NodalLaserInfo& node1, NodalLaserInfo& node2); 

  void BuildSortedNodeList();

  void BuildCustomizedCommunicators();
  void VerifySortedNodesAndCommunicators(); //!< for debug purpose only

  void SetupLaserGhostNodes();

  void ResetTag();

  void CopyValues(double*** from, double*** to);

  bool ComputeLaserRadianceMeanFluxMethod(const double t, double alpha, double relax_coeff);

  //-------------------------------------------------------------
  // functions called by ComputeLaserRadianceMeanFluxMethod
  void ComputeTemperatureInLaserDomain(Vec5D*** v, double*** id, double*** T);

  void UpdateGhostNodes(double ***l, int custom_max_iter = 0);

  void UpdateGhostNodesOneIteration(double ***l);

  double GetSourcePower(double t);

  void SetSourceRadiance(double*** l, Vec3D*** coords, double t);

  void AdjustRadianceToPowerChange(double*** l);

  void InitializeLaserDomainAndGhostNodes(double*** l, double*** level);

  void ComputeErrorsInLaserDomain(double*** lold, double*** lnew, double &max_error, double &avg_error);

  void RunMeanFluxMethodOneIteration(double*** l, double*** T, Vec3D*** coords, Vec3D*** 
                                     dxyz, double*** vol, double*** id, double alpha, double relax);

  void ComputeLaserHeating(double*** l, double*** T, double*** id, double*** s);

  void CleanUpGhostNodes(double*** l, bool backup); //!< remove radiance from ghost nodes (for outputting)

  void ApplyStoredGhostNodesRadiance(double*** l);
  //-------------------------------------------------------------

  void AddHeatToNavierStokesResidualSingleMesh(SpaceVariable3D &R, SpaceVariable3D &L, SpaceVariable3D &ID, 
                                               SpaceVariable3D *V = NULL); //if NULL, use stored temperature
  
/*
  inline double GetAbsorptionCoefficient(double T, int id) { //T must be in Kelvin
    return id<(int)absorption.size() ? 
             std::get<0>(absorption[id])*(T - std::get<1>(absorption[id])) + std::get<2>(absorption[id])
           : 0.0; //could get here if there are "inactive" nodes
  }
*/

};

#endif
