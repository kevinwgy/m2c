/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _EMBEDDED_BOUNDARY_DATA_SET_H_
#define _EMBEDDED_BOUNDARY_DATA_SET_H_

#include <TriangulatedSurface.h>
#include <SpaceVariable.h>

/*****************************************************************************
 * class IntersectionPoint is a utility class that stores information about an
 * intersection point between an edge an a triangle
 *
 * class ClosestPoint is a utility class that stores information about the
 * closest point on a triangle to another point (e.g., a node in the mesh)
 *
 * class EmbeddedBoundaryDataSet stores pointers to the results obtained
 * by an "intersector", including edge-surface intersections (always computed),
 * occluded nodes (always computed), "color" (always computed), 
 * unsigned distance (not always computed), and nodes swept by the surface 
 * in one step (not always computed).
 * 
 * Note: 
 *   (1) This class does not check which data is computed (i.e. usable), which
 *       is not.
 *   (2) If there are multiple embedded boundaries, each one of them should
 *       have an instantiation of this class.
 ****************************************************************************/


//! stores information about an intersection point between an edge and a triangle
struct IntersectionPoint {
  Int3 n0; //!< the first (i.e., left, bottom, or back) node
  int dir; //!< the direction of the edge (0~x, 1~y, 2~z)
  double dist; //!< dist from n0 to the intersection point along dir
  int tid; //!< id of the triangle that it intersects
  double xi[3]; //!< barycentric coords of the intersection point within the triangle it intersects

  /** NOTE: In most cases, the point obtained using tid and xi is IDENTICAL to the point
   *        obtained using n0, dir, and dist.
   *        However, the two points may not be the same when the intersection is IMPOSED to an edge because
   *        (1) one (or two) of the vertices of the edge is occluded and
   *        (2) an intersection cannot be identified normally (i.e. w/o imposing thickness).
   *        In this case, the second point (obtained using n0, dir, and dist) would be the occluded vertex,
   *        and the first point is the closest point on the triangle to the occluded vertex.
   *        If both vertices of the edge are occluded. The two vertices are considered as two intersection
   *        points. */

  IntersectionPoint() : n0(-1,-1,-1), dir(-1), dist(-1), tid(-1) {xi[0] = xi[1] = xi[2] = -1;}

  IntersectionPoint(int i, int j, int k, int dir_, double dist_, int tid_, double* xi_)
    : n0(Int3(i,j,k)), dir(dir_), dist(dist_), tid(tid_) {xi[0] = xi_[0]; xi[1] = xi_[1]; xi[2] = xi_[2];}

  IntersectionPoint &operator=(const IntersectionPoint& p2) {
    n0 = p2.n0;  dir = p2.dir;  dist = p2.dist;  tid = p2.tid;  for(int i=0; i<3; i++) xi[i] = p2.xi[i]; return *this;}
};


//! Utility class to store the closest point on the surface to a node (an arbitrary point)
struct ClosestPoint {
  int tid;
  double dist; //!< dist to the node in question
  double xi[3]; //!< closest point is at xi[0]*node[0] + xi[1]*node[1] + xi[2]*node[2]

  ClosestPoint(int tid_, double dist_, double xi_[3]) : tid(tid_), dist(dist_) {
    for(int i=0; i<3; i++) xi[i] = xi_[i];}

  ClosestPoint &operator=(const ClosestPoint& p2) {
    tid = p2.tid;  dist = p2.dist;
    for(int i=0; i<3; i++) {xi[i] = p2.xi[i];} return *this;}
};



class EmbeddedBoundaryDataSet {

public:

  TriangulatedSurface *surface_ptr;
  double half_thickness;

  SpaceVariable3D *XForward_ptr;
  SpaceVariable3D *XBackward_ptr;
  SpaceVariable3D *Phi_ptr;
  int Phi_nLayer;
  SpaceVariable3D *Color_ptr; 
  std::vector<int> *ColorReachesBoundary_ptr;
  bool hasInlet, hasInlet2, hasOutlet;
  int nRegions;
  SpaceVariable3D *ClosestPointIndex_ptr;
  std::vector<std::pair<Int3, ClosestPoint> > *closest_points_ptr;
  std::vector<IntersectionPoint> *intersections_ptr;

  std::set<Int3> *occluded_ptr;
  std::set<Int3> *firstLayer_ptr;
  std::set<Int3> *imposed_occluded_ptr;

  std::set<Int3> *swept_ptr;


public:
  EmbeddedBoundaryDataSet() : surface_ptr(nullptr), half_thickness(0.0),
                              XForward_ptr(nullptr), XBackward_ptr(nullptr), Phi_ptr(nullptr),
                              Phi_nLayer(0), Color_ptr(nullptr), ColorReachesBoundary_ptr(nullptr), 
                              hasInlet(false), hasInlet2(false), hasOutlet(false), nRegions(0),
                              ClosestPointIndex_ptr(nullptr), closest_points_ptr(nullptr), 
                              occluded_ptr(nullptr), firstLayer_ptr(nullptr),
                              imposed_occluded_ptr(nullptr), swept_ptr(nullptr)
  { }

  ~EmbeddedBoundaryDataSet() {}

};

#endif
