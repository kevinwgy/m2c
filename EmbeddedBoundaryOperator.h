#ifndef _EMBEDDED_BOUNDARY_OPERATOR_H_
#define _EMBEDDED_BOUNDARY_OPERATOR_H_

#include<Intersector.h>
#include<cassert>
#include<map>

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

  vector<Intersector*> intersector; //one intersector for each embedded surface (initialized to NULL)
 
  //! Mesh info (Not used when the class is used for special purposes, e.g., DynamicLoadCalculator)
  //! These information are generally needed when the surface needs to be "tracked" within the M2C mesh
  DataManagers3D* dms_ptr;
  SpaceVariable3D* coordinates_ptr;
  std::vector<GhostPoint>* ghost_nodes_inner_ptr;
  std::vector<GhostPoint>* ghost_nodes_outer_ptr;
  std::vector<double> *x_glob_ptr, *y_glob_ptr, *z_glob_ptr;
  std::vector<double> *dx_glob_ptr, *dy_glob_ptr, *dz_glob_ptr;

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



  void SetCommAndMeshInfo(DataManagers3D &dms_, SpaceVariable3D &coordinates_, 
                          std::vector<GhostPoint> &ghost_nodes_inner_, std::vector<GhostPoint> &ghost_nodes_outer_,
                          std::vector<double> &x_, std::vector<double> &y_, std::vector<double> &z_,
                          std::vector<double> &dx_, std::vector<double> &dy_, std::vector<double> &dz_);
  void SetupIntersectors();

  void ComputeForces(SpaceVariable3D &V, SpaceVariable3D &ID);

  void TrackSurfaces();
  void TrackUpdatedSurfaceFromOtherSolver();

private:

  void ReadMeshFile(const char *filename, EmbeddedSurfaceData::Type& surface_type,
                    vector<Vec3D> &Xs, vector<Int3> &Es);

  void UpdateSurfacesPrevAndFPrev(bool partial_copy=true); //!< copy surfaces.nodes/elements to surfaces_prev; also copy F to F_prev
};




#endif
