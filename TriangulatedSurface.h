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

  std::vector<Int3> elems; //line segment is recognized as a triangle with node2 = node3

  std::vector<Vec3D> elemNorm;
  std::vector<std::set<int> > node2node;
  std::vector<std::set<int> > node2elem;

  TriangulatedSurface(bool degen_ = false) : degenerate(degen_) { }

  TriangulatedSurface(std::vector<Vec3D> &X_, std::vector<Int3>& e_, bool degen_ = false) 
      : degenerate(degen_), X(X_), elems(e_) {
    BuildConnectivity();
    BuildElementNormals();
  }

  ~TriangulatedSurface() {} 

  void BuildConnectivity();
  void BuildElementNormals();

};

#endif
