/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _MULTI_SURFACE_INTERSECTOR_H_
#define _MULTI_SURFACE_INTERSECTOR_H_

#include<Intersector.h>

class MultiSurfaceIntersector {

  MPI_Comm &comm;

  //! Mesh info
  SpaceVariable3D& coordinates;

  std::vector<GhostPoint> &ghost_nodes_inner; //!< ghost nodes inside the physical domain (shared with other subd)
  std::vector<GhostPoint> &ghost_nodes_outer; //!< ghost nodes outside the physical domain

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain
  int ii0_in, jj0_in, kk0_in, iimax_in, jjmax_in, kkmax_in; //!< corners of the ghosted subdomain, excluding external ghosts
  int NX, NY, NZ; //!< global size (number of cells in the real domain)

  GlobalMeshInfo &global_mesh;

  int numSurfaces; //!< number of surfaces considered (has to be 1 or 2 at the moment)

  std::vector<TriangulatedSurface*> surface; //!< size = numSurfaces;
  std::vector<Intersector*> intersector; //!< size = numSurfaces
  std::vector<std::vector<bool> > elems_active;

  Intersector* joint_intersector; //!< constructed even for self-intersection (i.e., numSurfaces==1)
  TriangulatedSurface joint_surface; //!< surface provided to joint_intersector. 
  EmbeddedSurfaceData iod_surface_dummy; //!< a dummy (except for thickness) provided to joint_intersector. 

  int ruling_surface_id; //!< the index (starting at 0) in the `intersector' vector (-1: inactive)

  std::vector<int> new_enclosure_color; //!< color(s) from joint_intersector that represent new enclosure(s) 
  std::vector<std::vector<int> > elem_new_status; //!< for each new enclosure, status (0,1,2,3) marks its boundary

public:

  MultiSurfaceIntersector(MPI_Comm &comm_, DataManagers3D &dms_, SurfaceIntersectionData &iod_surfX_,
                          SpaceVariable3D &coordinates_, std::vector<TriangulatedSurface> &surface_,
                          std::vector<Intersector*> &intersector_, std::vector<GhostPoint> &ghost_nodes_inner_,
                          std::vector<GhostPoint> &ghost_nodes_outer_, GlobalMeshInfo &global_mesh_);

  ~MultiSurfaceIntersector();

  void Destroy();

  int GetNumberOfSurfaces() {return numSurfaces;}

  void UpdateJointSurface();

  bool CheckSurfaceIntersections();

  int FindNewEnclosures();

  void UpdateIntersectors();


private:

  bool CheckSurfaceIntersectionsOneWay(bool surf1_surf2);

  int FindNewEnclosuresByFloodFill(); //!< use joint_intersector; returns num. of new enclosures
  int DetectNewEnclosures(int nPossiblePositiveColors,
                          int nRegions1, int nRegions2, int nRegions_jnt, double*** color1, double*** color2,
                          double*** color_jnt, std::vector<int>& new_enclosure_color);


  void FindNewEnclosureBoundary(); //!< fills elem_new_status

  void MoidfyIntersectionsAndOccludedNodes(int id, std::vector<bool> elem_drop_status,
                                           std::set<int> elem_to_drop); //!< modifies intersector[id]

  void ModifyShortestDistance(int id, std::vector<bool> elem_drop_status,
                              std::set<int> elem_to_drop); //!< modifies intersector[id]);



};





#endif
