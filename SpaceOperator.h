/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _SPACEOPERATOR_H_
#define _SPACEOPERATOR_H_
#include <ExactRiemannSolverBase.h>
#include <GlobalMeshInfo.h>
#include <SymmetryOperator.h>
#include <ViscosityOperator.h>
#include <HeatDiffusionOperator.h>
#include <HyperelasticityOperator.h>
#include <SmoothingOperator.h>
#include <FluxFcnBase.h>
#include <Reconstructor.h>
#include <RiemannSolutions.h>

class EmbeddedBoundaryDataSet;
class TriangulatedSurface;
class IntersectionPoint;

/*******************************************
 * class SpaceOperator drives computations
 * that require domain/mesh information
 ******************************************/
class SpaceOperator
{
  MPI_Comm&                 comm;
  DataManagers3D&           dm_all;
  IoData&                   iod;
  FluxFcnBase&              fluxFcn;

  FluxFcnBase*              interfluxFcn; //!< used only when Multi-Material Flux = LocalLaxFriedrichs

  vector<VarFcnBase*>& varFcn; //!< each material has a varFcn

  //! Exact Riemann problem solver (multi-phase)
  ExactRiemannSolverBase &riemann;

  //! Mesh info
  SpaceVariable3D coordinates;
  SpaceVariable3D delta_xyz;
  SpaceVariable3D volume; //!< volume of node-centered control volumes
  
  vector<GhostPoint> ghost_nodes_inner; //!< ghost nodes inside the physical domain (shared with other subd)
  vector<GhostPoint> ghost_nodes_outer; //!< ghost nodes outside the physical domain

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain
  int NX, NY, NZ; //!< global size

  GlobalMeshInfo &global_mesh;

  //! Nodes/cells where residual should be re-set to 0. (Assigned by other modules)
  std::set<Int3> *frozen_nodes_ptr;

  //! Class for spatial reconstruction
  Reconstructor rec;

  //! Class for imposing spherical or cylindrical symmetry (sink terms placed on the left-hand-side!)
  SymmetryOperator* symm;

  //! Class for computing viscous fluxes
  ViscosityOperator* visco;

  //! Class for calculating heat diffusion fluxes
  HeatDiffusionOperator* heat_diffusion;

  //! Class for smoothing the solution
  HyperelasticityOperator* heo;

  //! Class for smoothing the solution
  SmoothingOperator* smooth;

  //! Reconstructed primitive state variables at cell boundaries
  SpaceVariable3D Vl, Vr, Vb, Vt, Vk, Vf;

  //! For temporary variable (5D)
  SpaceVariable3D Utmp;

  //! internal variable for temporary use (1D)
  SpaceVariable3D Tag;

  //! State variables at overset ghost/boundary nodes (only type == FACE)
  bool domain_has_overset; //!< whether the entire domain has overset boundaries
  vector<std::pair<Int3, Vec5D> > ghost_overset; //!< overset ghost nodes outside the physical domain


public:
  SpaceOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_,
                vector<VarFcnBase*> &varFcn_, FluxFcnBase &fluxFcn_,
                ExactRiemannSolverBase &riemann_,
                GlobalMeshInfo &global_mesh_,
                bool screenout = true); 
  ~SpaceOperator();

  //! Reset the coords of ghost layer nodes  (a NULL pointer means that value does not need to be reset)
  void ResetGhostLayer(double* xminus, double* xplus, double* yminus,  double* yplus,
                       double* zminus, double* zplus, double* dxminus, double* dxplus, double* dyminus,
                       double* dyplus, double* dzminus, double* dzplus);

  void ConservativeToPrimitive(SpaceVariable3D &U, SpaceVariable3D &ID, SpaceVariable3D &V,
                               bool workOnGhost = false);
  void PrimitiveToConservative(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &U,
                               bool workOnGhost = false);
  int  ClipDensityAndPressure(SpaceVariable3D &V, SpaceVariable3D &ID, 
                              bool workOnGhost = false, bool checkState = true);

  void SetupViscosityOperator(InterpolatorBase *interpolator_, GradientCalculatorBase *grad_,
                              bool with_embedded_boundary = false);

  ViscosityOperator* GetPointerToViscosityOperator() {return visco;} //!< can be NULL!

  void SetupHeatDiffusionOperator(InterpolatorBase *interpolator_, GradientCalculatorBase *grad_);

  void SetHyperelasticityOperatorPointer(HyperelasticityOperator *heo_) {heo = heo_;}

  void ApplyBoundaryConditions(SpaceVariable3D &V);

  void ApplySmoothingFilter(double time, double dt, int time_step, SpaceVariable3D &V, SpaceVariable3D &ID);

  void FindExtremeValuesOfFlowVariables(SpaceVariable3D &V, SpaceVariable3D &ID,
                                        double *Vmin, double *Vmax, double &cmin, 
                                        double &cmax, double &Machmax, double &char_speed_max,
                                        double &dx_over_char_speed_min);

  void ComputeTimeStepSize(SpaceVariable3D &V, SpaceVariable3D &ID, double &dt, double &cfl,
                           SpaceVariable3D *LocalDt = NULL);

  void ComputeLocalTimeStepSizes(SpaceVariable3D &V, SpaceVariable3D &ID, double &dt, double &cfl,
                                 SpaceVariable3D &LocalDt);

  double ComputeTimeStepSizeSurfaceTension(SpaceVariable3D &V, SpaceVariable3D &ID);

  //! Compute the RHS of the ODE system (Only for cells inside the physical domain)
  void ComputeResidual(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &R, double time,
                       RiemannSolutions *riemann_solutions = NULL,
                       vector<int> *ls_mat_id = NULL, vector<SpaceVariable3D*> *Phi = NULL,
                       vector<SpaceVariable3D*> *KappaPhi = NULL,
                       vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS = nullptr,
                       SpaceVariable3D *Xi = NULL, bool run_heat = true);

  SpaceVariable3D& GetMeshCoordinates() {return coordinates;}
  SpaceVariable3D& GetMeshDeltaXYZ()    {return delta_xyz;}
  SpaceVariable3D& GetMeshCellVolumes() {return volume;}

  vector<GhostPoint>* GetPointerToInnerGhostNodes() {return &ghost_nodes_inner;}
  vector<GhostPoint>* GetPointerToOuterGhostNodes() {return &ghost_nodes_outer;}

  GlobalMeshInfo& GetGlobalMeshInfo() {return global_mesh;}

  void SetPointerToFrozenNodes(std::set<Int3>* fnodes_) {frozen_nodes_ptr = fnodes_;}

  void UpdateOversetGhostNodes(SpaceVariable3D &V);

  void Destroy();


private:

  void SetupMesh(vector<double> &x, vector<double> &y, vector<double> &z,
                 vector<double> &dx, vector<double> &dy, vector<double> &dz);
  void SetupMeshUniformRectangularDomain();
  void PopulateGhostBoundaryCoordinates();

  void CreateGhostNodeLists(bool screenout);

  void ApplyBoundaryConditionsGeometricEntities(Vec5D*** v);

  void CheckReconstructedStates(SpaceVariable3D &V,
                                SpaceVariable3D &Vl, SpaceVariable3D &Vr, SpaceVariable3D &Vb,
                                SpaceVariable3D &Vt, SpaceVariable3D &Vk, SpaceVariable3D &Vf,
                                SpaceVariable3D &ID);

  void ComputeAdvectionFluxes(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &F,
                              RiemannSolutions *riemann_solutions = NULL,
                              vector<int> *ls_mat_id = NULL, vector<SpaceVariable3D*> *Phi = NULL, vector<SpaceVariable3D*> *KappaPhi = NULL,
                              vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS = nullptr);

  Vec3D GetNormalForOneSidedRiemann(int d,/*0,1,2*/
                                    int forward_or_backward,/*1~wall is in the +x/y/z dir of material, -1~-x/y/z*/
                                    Vec3D& nwall);

  Vec3D GetNormalForBimaterialRiemann(int d/*0,1,2*/, int i, int j, int k, Vec3D*** coords, Vec3D*** dxyz,
                                      int myid, int neighborid, vector<int> *ls_mat_id,
                                      vector<double***> *phi);

  Vec3D CalculateGradPhiAtCellInterface(int d/*0,1,2*/, int i, int j, int k, Vec3D*** coords, Vec3D*** dxyz,
                                        int myid, int neighborid, vector<int> *ls_mat_id,
                                        vector<double***> *phi);

  //! Calculate the gradient of any variable phi at cell interface. 
  Vec3D CalculateGradientAtCellInterface(int d/*0,1,2*/, int i, int j, int k, Vec3D*** coords,
                                         Vec3D*** dxyz, double*** phi);

  double CalculateCurvatureAtCellInterface(int d/*0,1,2*/, double*** phi, double*** kappaPhi, int i, int j, int k);  

  bool TagNodesOutsideConRecDepth(vector<SpaceVariable3D*> *Phi, 
                                  vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS,
                                  SpaceVariable3D &Tag0);

  //! Find intersections (forward and backward) between an edge and embedded surface(s)
  bool FindEdgeSurfaceIntersections(int dir/*0~x,1~y,2~z*/, int i, int j, int k,
                                    vector<TriangulatedSurface*>& surfaces,
                                    vector<vector<IntersectionPoint>*>& intersections,
                                    vector<Vec3D***>& xf, vector<Vec3D***>& xb,
                                    Vec3D& vwallf, Vec3D& vwallb, Vec3D& nwallf, Vec3D& nwallb);

  //! Utility
  inline double CentralDifferenceLocal(double phi0, double phi1, double phi2, double x0, double x1, double x2) {
    double c0 = -(x2-x1)/((x1-x0)*(x2-x0));
    double c1 = 1.0/(x1-x0) - 1.0/(x2-x1);
    double c2 = (x1-x0)/((x2-x0)*(x2-x1));
    return c0*phi0 + c1*phi1 + c2*phi2;
  }
};


#endif
