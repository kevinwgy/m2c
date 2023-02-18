/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _TRIANGULATED_SURFACE_H_
#define _TRIANGULATED_SURFACE_H_
#include <Vector3D.h>
#include <set>
#include <vector>

/*****************************************************************************
 * A utility class to store a triangulated surface (including the degenerate 
 * scenario of a set of line segments in the x-y plane; but cannot have both
 * triangles and line segments)
 *****************************************************************************/
struct TriangulatedSurface {

  bool degenerate; //!< line segments in 2D (x-y)

  std::vector<Vec3D> X0;   //!< Original config. (for points in 2D: z-coord = 0)
  std::vector<Vec3D> X;    //!< Current config. (for point in 2D: z-coord = 0)
  std::vector<Vec3D> Udot; //!< Velocity vector

  std::vector<Int3> elems; //!< line segment is recognized as a triangle with node2 = node3

  //! number of active nodes and elements (not equal to the size of X0, elems in case of fracture)
  int active_nodes;
  int active_elems;

  std::vector<Vec3D> elemNorm;
  std::vector<double> elemArea;
  std::vector<std::set<int> > node2node;
  std::vector<std::set<int> > node2elem;
  std::vector<std::set<int> > elem2elem;

  TriangulatedSurface(bool degen_ = false) : degenerate(degen_), active_nodes(0), active_elems(0) { }

  TriangulatedSurface(std::vector<Vec3D> &X_, std::vector<Int3>& e_, bool degen_ = false) 
      : degenerate(degen_), X(X_), elems(e_) {
    active_nodes = X.size();
    active_elems = elems.size();
    BuildConnectivities();
    CalculateNormalsAndAreas();
  }

  ~TriangulatedSurface() {} 

  void BuildConnectivities();
  void CalculateNormalsAndAreas(); //!< calculate the normal and area of each element

  bool CheckSurfaceOrientation(); //!<  check whether all the elements have consistent normal directions
  bool CheckSurfaceClosedness(); //!< check whether this is a closed surface
  bool CheckSurfaceOrientationAndClosedness(); //!< more efficient that calling the above two functions seperately


};

#endif
