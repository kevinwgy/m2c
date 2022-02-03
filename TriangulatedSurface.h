#ifndef _TRIANGULATED_SURFACE_H_
#define _TRIANGULATED_SURFACE_H_
#include <Vector3D.h>
#include <set>
#include <vector>
using std::set;
using std::vector;

/*****************************************************************************
 * A utility class to store a triangulated surface (including the degenerate 
 * scenario of a set of line segments in the x-y plane; but cannot have both
 * triangles and line segments)
 *****************************************************************************/
struct TriangulatedSurface {

  bool degenerate; //!< line segments in 2D (x-y)

  vector<Vec3D> X0;    //point in 2D: z-coord = 0
  vector<Vec3D> X;    //point in 2D: z-coord = 0
  vector<Int3> elems; //line segment is recognized as a triangle with node2 = node3

  vector<Vec3D> elemNorm;
  vector<set<int> > node2node;
  vector<set<int> > node2elem;

  TriangulatedSurface(vector<Vec3D> &X_, vector<Int3>& e_) : X(X_), elems(e_) {
    degenerate = false;
    BuildConnectivity();
    BuildElementNormals();
  }

  ~TriangulatedSurface() {} 

private:
  void BuildConnectivity();
  void BuildElementNormals();
};

#endif
