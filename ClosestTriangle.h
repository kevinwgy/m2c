/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _CLOSEST_TRIANGLE_H_
#define _CLOSEST_TRIANGLE_H_

#include <Vector3D.h>
#include <map>
#include <set>

/***********************************************************
 * class ClosestTriangle is a utility class that finds the 
 * signed distance from a point to a set of triangles that
 * form a CLOSED SURFACE (e.g., an embedded structure). This class was
 * originally written by KW and ML around 2009. See
 * Wang et al., IJNMF, 2012 (The projection-based method).
 * Usage: It is usually not efficient nor necessary to run
 * "checkTriangle" over all the triangles, as some of triangles
 * may be very far from the point. It is usually a good idea
 * to create a tree structure (e.g., KDTree) to store the
 * triangles, and use it to pre-select the triangles that
 * are relatively close to the point (e.g., using a bounding
 * box). Then, run "checkTriange" over these selected triangles.
 **********************************************************/
class ClosestTriangle {

  int (*triNodes)[3];
  Vec3D *structX;
  Vec3D *structNorm;
  std::set<int> *node2node;
  std::set<int> *node2elem;

protected:
  bool fail;
  bool isFirst;
  int bestTrId;
  int n1, n2; //!< if both ns are non negative, the best point is on an edge
  Vec3D n;
  double minDist; //!< Signed distance to the surface
  Vec3D x;
  std::map<int,int> nd2tri; //needed for the vertex case (only if isConsistent = false)
  std::map<int,int> nd2tester;
  double dist2othernode;

  static const int maxNPairs = 100;

  bool isConsistent, isPositive, hasBeenFixed;
  int nPairs;
  int periTri[maxNPairs];

  bool checkEdge(int trId, int p1, int p2, int p3, double trDist);
  void checkVertex(int vn, int trId, double trDist);
  int registerNodes(int ip1, int trId, int& repeated1, int& repeated2);
  double project(Vec3D x0, int tria, double& xi1, double& xi2) const;
  double edgeProject(Vec3D x0, int n1, int n2, double &alpha) const;
  double getSignedVertexDistance() const;
  double findSignedVertexDistance();
  int findTesterStatus(Vec3D) const;
  void checkEdgeForTester(Vec3D xt, int trId, int ip1, int ip2, int p3, double trDist, int &nn1, int &nn2,
                          double &mindist, int &myMode, int &bestTriangle, const double eps) const;

public:

  ClosestTriangle(int (*triNodes)[3], Vec3D *structX, Vec3D *sN, std::set<int> *n2n, std::set<int> *n2e);
  ~ClosestTriangle() {}

  void start(Vec3D x);
  void checkTriangle(int trId);
  double signedDistance() {
    if(n1 < 0 || n2 >= 0) return minDist;
    else
      return findSignedVertexDistance();
      //return getSignedVertexDistance();
  }

  int bestTriangle() const { return bestTrId; }

  int mode;

};



#endif













