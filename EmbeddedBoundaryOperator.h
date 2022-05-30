#ifndef _EMBEDDED_BOUNDARY_OPERATOR_H_
#define _EMBEDDED_BOUNDARY_OPERATOR_H_

#include<Intersector.h>
#include<UserDefinedDynamics.h>
#include<LagrangianOutput.h>
#include<cassert>
#include<map>
#include<tuple>

struct Vec5D;

/******************************************************************
 * Class EmbeddedBoundaryOperator stores data related to one or
 * multiple embedded surfaces, which can be provided either by a
 * structural dynamics solver (e.g., AERO-S) or directly by the
 * user in a file. The class also contains functions for tracking
 * the embedded surfaces and enforcing interface/boundary conditions
 *****************************************************************/

class EmbeddedBoundaryOperator {

  MPI_Comm &comm;
  IoData &iod;

  bool hasSurfFromOtherSolver; //!< currently, at most 1 surface can come from another solver (to be generalized)

  //! data and tools for each surface
  vector<EmbeddedSurfaceData*> iod_embedded_surfaces; //!< iodata
  vector<TriangulatedSurface> surfaces; //!< embedded surfaces, usually updated at each time step
  vector<TriangulatedSurface> surfaces_prev; //!< saves the topology and coords in the previous time step
  vector<vector<Vec3D> > F; //!< forces
  vector<vector<Vec3D> > F_prev; //!< forces at the previous time step
  vector<EmbeddedSurfaceData::Type> surface_type;
  vector<Intersector*> intersector; //!< one intersector for each embedded surface (initialized to NULL)
 
  vector<LagrangianOutput> lagout; //!< output displacement and *nodal load* on embedded surfaces

  //! inactive closures: pair of <surface number, color>, not including color = 0 (occluded)
  std::set<std::pair<int,int> > inactive_colors;

  //! for each surface (i), inactive_elem_status[i][j] (j: 0 -- surfaces[i].elems.size()) shows weather one or both
  //! sides of triangle element j is part of the inward-facing side of any inactive region. 
  //! Needed for force computation
  std::vector<std::vector<int> > inactive_elem_status;

  vector<std::tuple<UserDefinedDynamics*, void*, DestroyUDD*> > dynamics_calculator; //!< the 1st one is the calculator

  //! Mesh info (Not used when the class is used for special purposes, e.g., DynamicLoadCalculator)
  //! These information are generally needed when the surface needs to be "tracked" within the M2C mesh
  DataManagers3D* dms_ptr;
  SpaceVariable3D* coordinates_ptr;
  std::vector<GhostPoint>* ghost_nodes_inner_ptr;
  std::vector<GhostPoint>* ghost_nodes_outer_ptr;
  GlobalMeshInfo *global_mesh_ptr;

public:
   
  //! Constructor: Without intersector (pointers are set to NULL) or M2C mesh info
  EmbeddedBoundaryOperator(MPI_Comm &comm_, IoData &iod_, bool surface_from_other_solver = false);
  ~EmbeddedBoundaryOperator();

  void Destroy();

  int  NumberOfSurfaces() {return surfaces.size();}
  bool HasSurfaceFromOtherSolver() {return hasSurfFromOtherSolver;} //!< index of surface from other solver is 0

  vector<TriangulatedSurface> *GetPointerToSurfaces() {return &surfaces;}
  TriangulatedSurface         *GetPointerToSurface(int i) {assert(i>=0 && i<surfaces.size()); return &surfaces[i];}
  vector<vector<Vec3D> >      *GetPointerToForces() {return &F;}
  vector<Vec3D>               *GetPointerToForcesOnSurface(int i) {assert(i>=0 && i<F.size()); return &F[i];}
  vector<Intersector*>        *GetPointerToIntersectors() {return &intersector;}
  Intersector                 *GetPointerToIntersector(int i) {assert(i>=0 && i<intersector.size()); return intersector[i];}

  std::unique_ptr<std::vector<std::unique_ptr<EmbeddedBoundaryDataSet> > > GetPointerToEmbeddedBoundaryData();
  std::unique_ptr<EmbeddedBoundaryDataSet> GetPointerToEmbeddedBoundaryData(int i); 


  void SetCommAndMeshInfo(DataManagers3D &dms_, SpaceVariable3D &coordinates_, 
                          std::vector<GhostPoint> &ghost_nodes_inner_, std::vector<GhostPoint> &ghost_nodes_outer_,
                          GlobalMeshInfo &global_mesh_);
  void SetupIntersectors();

  //build inactive_colors, and inactive_elem_status
  void FindSolidBodies(std::multimap<int, std::pair<int,int> > &id2closure); 

  void ComputeForces(SpaceVariable3D &V, SpaceVariable3D &ID);

  void TrackSurfaces();
  void TrackUpdatedSurfaces();

  void ApplyUserDefinedSurfaceDynamics(double t, double dt);

  void UpdateSurfacesPrevAndFPrev(bool partial_copy=true); //!< copy surfaces.nodes/elements to surfaces_prev; also copy F to F_prev

  //! Output embedded surfaces to files (meshes)
  void OutputSurfaces();

  //! Output displacements and nodal forces on embedded surfaces
  void OutputResults(double time, double dt, int time_step, bool force_write=false);

private:

  void ReadMeshFile(const char *filename, EmbeddedSurfaceData::Type& surface_type,
                    vector<Vec3D> &Xs, vector<Int3> &Es);

  void SetupUserDefinedDynamicsCalculator(); //!< setup dynamics_calculator

  double CalculateLoftingHeight(Vec3D &p, double factor);

  Vec3D CalculateTractionAtPoint(Vec3D &p, int side/*0 or 1*/, double area, Vec3D &normal/*towards the "side"*/, 
                                 Int3 &tnodes, vector<Vec3D> &Xs, Vec5D*** v, double*** id);

};




#endif
