/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <Intersections.h>
#include <GeoTools.h>
#include <set>
#include <utility> //std::pair
#include <cassert>
#include <cfloat> //DBL_MAX
#include <algorithm> //std::sort

using std::vector;
using std::set;

namespace GeoTools {

// ------------------------------------------------------------------------------------------------------------
//! Moller-Trumbore intersection (Moller and Trumbore, 1997)
bool
RayIntersectsTriangle(Vec3D O, Vec3D D, //!< origina and direction of ray
                      Vec3D& V0, Vec3D& V1, Vec3D& V2, //!< vertices of triangle
                      double* d, //!< optional output: dist from origin to intersection
                      Vec3D* xp, //!< optional output: intersection point
                      Vec3D* baryCoords, //!< optional output: barycentric coords of intersection point
                      bool D_normalized) //!< input: whether D is normalized
{
  if(!D_normalized) {
    double Dnorm = D.norm();
    assert(Dnorm>0.0);
    D /= Dnorm;
  }

  Vec3D E1, E2, P, T, Q;
  double denom,f,u,v;

  E1 = V1 - V0;
  E2 = V2 - V0;
  P = D^E2;
  denom = P*E1;

  if (denom==0) {
    if(d)  *d  = DBL_MAX;
    if(xp) *xp = DBL_MAX;
    if(baryCoords) *baryCoords = DBL_MAX;
    return false;    // This ray is parallel to this triangle.
  }

  f = 1.0/denom;
  T = O - V0;

  // Calculate u
  u = f*(P*T);
  if (u < -INTERSECTIONS_EPSILON || u > 1.0+INTERSECTIONS_EPSILON) {
    if(d || xp || baryCoords) { // calculate optional outputs
      Q = T^E1;
      v = f*(Q*D);
      double t = f*(Q*E2);
      if(d)  *d  = t; 
      if(xp) *xp = O + t*D;
      if(baryCoords) *baryCoords = Vec3D(1.0-u-v, u, v);
    }
    return false;
  }

  // Calculate v
  Q = T^E1;
  v = f*(Q*D);
  if (v < -INTERSECTIONS_EPSILON || u + v > 1.0+INTERSECTIONS_EPSILON) {
    if(d || xp || baryCoords) { // calculate optional outputs
      double t = f*(Q*E2);
      if(d)  *d  = t; 
      if(xp) *xp = O + t*D;
      if(baryCoords) *baryCoords = Vec3D(1.0-u-v, u, v);
    }
    return false;
  }

  // Calculate t to find out where the intersection point is on the line.
  double t = f*(Q*E2);
  if(d)  *d  = t;
  if(xp) *xp = O + D * t;
  if(baryCoords) *baryCoords = Vec3D(1.0-u-v, u, v);
  if (t > -INTERSECTIONS_EPSILON) {
    return true;
  } else 
    return false;
}

// ------------------------------------------------------------------------------------------------------------

// for axis-aligned directions
bool
AxisIntersectsTriangle(Vec3D O, int dir, //!< dir = 0 (x-axis), 1 (y-axis), or 2 (z-axis)
                       Vec3D& V0, Vec3D& V1, Vec3D& V2,
                       double* d, Vec3D* xp, Vec3D* baryCoords)
{

  Vec3D P = V1 - V0;
  Vec3D Q = V2 - V0;
  Vec3D R = O  - V0;
  double p1(0.0), p2(0.0), q1(0.0), q2(0.0), r1(0.0), r2(0.0);

  if(dir == 0) { //x-axis
    p1 = P[1];
    p2 = P[2];
    q1 = Q[1];
    q2 = Q[2];
  } 
  else if(dir == 1) { //y-axis
    p1 = P[0];
    p2 = P[2];
    q1 = Q[0];
    q2 = Q[2];
  }
  else if(dir == 2) { //z-axis
    p1 = P[0];
    p2 = P[1];
    q1 = Q[0];
    q2 = Q[1];
  }
  else {
    fprintf(stdout,"*** Error: In AxisIntersectsTriangle, dir must be 0, 1, or 2. Found %d.\n", dir);
    exit(-1);
  }

  double denom = p1*q2 - p2*q1;
  if(denom==0) {
    if(d)  *d  = DBL_MAX;
    if(xp) *xp = DBL_MAX;
    if(baryCoords) *baryCoords = DBL_MAX;
    return false;    // This ray is parallel to this triangle.
  }

  // solve the 2x2 linear system
  if(dir == 0) {
    r1 = R[1];
    r2 = R[2];
  } else if (dir == 1) {
    r1 = R[0];
    r2 = R[2];
  } else {
    r1 = R[0];
    r2 = R[1];
  }

  double f = 1.0/denom;

  double u = f*(r1*q2 - r2*q1);
  if(u < -INTERSECTIONS_EPSILON || u > 1.0+INTERSECTIONS_EPSILON) {
    if(d || xp || baryCoords) { // calculate optional outputs
      double v = f*(r2*p1 - r1*p2);
      if(d)  *d  = P[dir]*u + Q[dir]*v - R[dir]; 
      if(xp) *xp = (1.0-u-v)*V0 + u*V1 + v*V2;
      if(baryCoords) *baryCoords = Vec3D(1.0-u-v, u, v);
    }
    return false;
  }

  double v = f*(r2*p1 - r1*p2);
  if (v < -INTERSECTIONS_EPSILON || u + v > 1.0+INTERSECTIONS_EPSILON) {
    if(d || xp || baryCoords) { // calculate optional outputs
      if(d)  *d  = P[dir]*u + Q[dir]*v - R[dir]; 
      if(xp) *xp = (1.0-u-v)*V0 + u*V1 + v*V2;
      if(baryCoords) *baryCoords = Vec3D(1.0-u-v, u, v);
    }
    return false;
  }

  // Calculate t to find out where the intersection point is on the line.
  double t = P[dir]*u + Q[dir]*v - R[dir];
  if(d)  *d  = t;
  if(xp) *xp = (1.0-u-v)*V0 + u*V1 + v*V2;
  if(baryCoords) *baryCoords = Vec3D(1.0-u-v, u, v);
  if (t > -INTERSECTIONS_EPSILON) {
    return true;
  } else 
    return false;
    
}

// ------------------------------------------------------------------------------------------------------------
//! Moller-Trumbore intersection algorithm (Moller and Trumbore, 1997)
bool
LineSegmentIntersectsTriangle(Vec3D X0, Vec3D X1, //!< vertices of line segment
                              Vec3D& V0, Vec3D& V1, Vec3D& V2, //!< vertices of triangle
                              double* d, //!< optional output: dist from X0 to intersection
                              Vec3D* xp, //!< optional output: intersection point
                              Vec3D* baryCoords) //!< optional output: barycentric coords of intersection point
{

  Vec3D D = X1 - X0;
  double len = D.norm();
  assert(len>0);
  D /= len;

  double t(0.0);
  bool intersect = RayIntersectsTriangle(X0, D, V0, V1, V2, &t, xp, baryCoords, true);

  if(d) *d = t;

  if(!intersect || t>len + INTERSECTIONS_EPSILON)
    return false;
  else
    return true;
}
    
// ------------------------------------------------------------------------------------------------------------
// for axis-aligned directions
bool
LineSegmentIntersectsTriangle(Vec3D O, int dir, //!< dir = 0 (x-axis), 1 (y-axis), or 2 (z-axis)
                              double len, //!< length of the line segment (along dir)
                              Vec3D& V0, Vec3D& V1, Vec3D& V2,
                              double* d, Vec3D* xp, Vec3D* baryCoords)
{
  assert(len>0);

  double t(0.0);
  bool intersect = AxisIntersectsTriangle(O, dir, V0, V1, V2, &t, xp, baryCoords);

  if(d) *d = t;

  if(!intersect || t>len + INTERSECTIONS_EPSILON)
    return false;
  else
    return true;
}

// ------------------------------------------------------------------------------------------------------------

bool
LineSegmentIntersectsPlane(Vec3D X0, Vec3D X1, //!< vertices of line segment
                           Vec3D V0, Vec3D dir, //!< a point on the plane, and its normal
                           double* d, //!< optional output: dist from X0 to intersection
                           Vec3D* xp, //!< optional output: intersection point
                           bool N_normalized) //!< input: whether dir is normalized
{
  if(!N_normalized) {
    double Nnorm = dir.norm();
    assert(Nnorm>0.0);
    dir /= Nnorm;
  }

  Vec3D X01 = X1 - X0;
  double denom = X01*dir;

  if(denom==0) {
    if(d)  *d  = DBL_MAX;
    if(xp) *xp = DBL_MAX;
    return false;    // This line segment is parallel to the plane
  }

  double d0 = ProjectPointToPlane(X0, V0, dir, true);
  double d1 = ProjectPointToPlane(X1, V0, dir, true);

  if(d0*d1<0) {
    if(d || xp) {
      double alpha = fabs(d0)/(fabs(d0)+fabs(d1));
      if(d) 
        *d = alpha*X01.norm();
      if(xp)
        *xp = X0 + alpha*X01;
    }
    return true;
  }
  else if(fabs(d0)<INTERSECTIONS_EPSILON) {
    if(d)
      *d = 0.0;
    if(xp)
      *xp = X0;
    return true;
  }
  else if(fabs(d1)<INTERSECTIONS_EPSILON) {
    if(d)
      *d = X01.norm();
    if(xp)
      *xp = X1;
    return true;
  }
  else {//on the same side
    if(d)
      *d = DBL_MAX;
    if(xp)
      *xp = DBL_MAX;
    return false;
  }

  return false; //will never get here
}

// ------------------------------------------------------------------------------------------------------------
// Connectivities of 8-noded box element  
//    7-------------6
//  4/_|________5_/|
//   | |         | |
//   | |         | |
//   | |3________+_|2
//  0|/_________1|/
//
int E2N[12][2] = //edge to node connectivity [min-node, max-node]
    {{0,1}, {1,2}, {3,2}, {0,3}, {0,4}, {1,5}, {2,6}, {3,7}, {4,5}, {5,6}, {7,6}, {4,7}};
int N2N[8][3] = //node to node connectivity
      {{1,3,4}, {0,2,5}, {1,3,6}, {0,2,7}, {0,5,7}, {1,4,6}, {2,5,7}, {3,4,6}};
int N2E[8][3] = //node to edge connectivity
      {{0,3,4}, {0,1,5}, {1,2,6}, {2,3,7}, {4,8,11}, {5,8,9}, {6,9,10}, {7,10,11}};
int ECoPlane[12][6] = //edge to co-planar edge connectivity
      {{1,2,3,4,5,8}, {0,2,3,5,6,9}, {0,1,3,6,7,10}, {0,1,2,4,7,11}, {0,5,8,3,7,11},
       {0,4,8,1,6,9}, {1,5,9,2,7,10}, {2,6,10,3,4,11}, {0,4,5,9,10,11}, {1,5,6,8,10,11},
       {2,6,7,8,9,11}, {3,4,7,8,9,10}};
int NOppo[8] = {6, 7, 4, 5, 2, 3, 0, 1};//node to opposite corner node 

// ------------------------------------------------------------------------------------------------------------

int
PlaneCuttingAxisAlignedBox(Vec3D& O, Vec3D& N, Vec3D& Vmin, Vec3D& Vmax, 
                           vector<Vec3D> *intersections, double tolerance)
{

  if(N.norm()==0)
    return 0;

  N /= N.norm();

  Vec3D dV = Vmax - Vmin;
  assert(dV[0]>=0 && dV[1]>=0 && dV[2]>=0); 

  assert(tolerance>=0.0);

  // fast check
  double upper_bound = 1.01*dV.norm() + tolerance;
  if(fabs((Vmax-O)*N)>upper_bound || fabs((Vmin-O)*N)>upper_bound)
    return 0;


  // order the 8 vertices based on rule specified above
  Vec3D V[8] = {Vmin, Vec3D(Vmax[0],Vmin[1],Vmin[2]), Vec3D(Vmax[0],Vmax[1],Vmin[2]),
                Vec3D(Vmin[0],Vmax[1],Vmin[2]),
                Vec3D(Vmin[0],Vmin[1],Vmax[2]), Vec3D(Vmax[0],Vmin[1],Vmax[2]),
                Vmax, Vec3D(Vmin[0],Vmax[1],Vmax[2])};
                       
  int status[8];
  double dist[8];
  set<int> xnodes;
  int pos_count(0), zer_count(0), neg_count(0);
  for(int i=0; i<8; i++) {
    dist[i] = (V[i]-O)*N;
    if(dist[i]>tolerance) {
      status[i] = 1;
      pos_count++;
    } else if(dist[i]<-tolerance) {
      status[i] = -1;
      neg_count++;
    } else {
      status[i] = 0;
      zer_count++;
      xnodes.insert(i);
    }
  }

  if(zer_count==0 && (pos_count==0 || neg_count==0))
    return 0;

  // Now, there must be intersections!

  set<int> xedges;
  for(int i=0; i<12; i++) {
    if(status[E2N[i][0]]==0 || status[E2N[i][1]]==0)
      continue;
    if(status[E2N[i][0]] != status[E2N[i][1]]) 
      xedges.insert(i);
  }
  
  int nPoints = xedges.size()+xnodes.size();
  assert(nPoints>0);

  if(!intersections)
    return nPoints; // no need to compute anything else

  // fill intersections
  if((int)intersections->size()<nPoints)
    intersections->resize(nPoints);
 
  int counter = 0;
  std::pair<bool, int> tip; //<node?, node or edge number>

  // find the starter
  if(!xnodes.empty()) {
    (*intersections)[counter++] = V[*xnodes.begin()];
    tip.first = true;
    tip.second = *xnodes.begin();
    xnodes.erase(xnodes.begin());
  } else {
    double d1 = fabs(dist[E2N[*xedges.begin()][0]]);
    double d2 = fabs(dist[E2N[*xedges.begin()][1]]);
    double sum = d1 + d2;
    assert(sum>0);
    (*intersections)[counter++] = (d2*V[E2N[*xedges.begin()][0]] + 
                                   d1*V[E2N[*xedges.begin()][1]])/sum;
    tip.first = false;
    tip.second = *xedges.begin();
    xedges.erase(xedges.begin());
  }

  // complete the loop
  int iter = 0;
  while(!xnodes.empty() || !xedges.empty()) {

    int me = tip.second;

    for(auto it = xnodes.begin(); it != xnodes.end(); it++) {
      if(tip.first) { //tip is a node
        if(NOppo[me] != *it) { //co-planar, acceptable
          (*intersections)[counter++] = V[*it];
          tip.second = *it;
          me = tip.second;
          xnodes.erase(it); 
          break;
        }
      } else { //tip is an edge
        if(E2N[me][0] == *it || E2N[me][1] == *it) {
          (*intersections)[counter++] = V[*it];
          tip.first = true; //new tip is a node
          tip.second = *it;
          me = tip.second;
          xnodes.erase(it);
          break;
        }
      }
    }

    for(auto it = xedges.begin(); it != xedges.end(); it++) {
      if(tip.first) { //tip is a node
        if(N2N[me][0] == E2N[*it][0] || N2N[me][0] == E2N[*it][1] ||
           N2N[me][1] == E2N[*it][0] || N2N[me][1] == E2N[*it][1] ||
           N2N[me][2] == E2N[*it][0] || N2N[me][2] == E2N[*it][1]) {
          double d1 = fabs(dist[E2N[*it][0]]);
          double d2 = fabs(dist[E2N[*it][1]]);
          double sum = d1 + d2;
          assert(sum>0);
          (*intersections)[counter++] = (d2*V[E2N[*it][0]] + d1*V[E2N[*it][1]])/sum;
          tip.first = false;  //new tip is an edge
          tip.second = *it;
          xedges.erase(it);
          break;
        }
      } else { // tip is an edge
        bool found = false;
        for(int j=0; j<6; j++)
          if(ECoPlane[me][j] == *it) {
            found = true;
            break;
          }
        if(found) {
          double d1 = fabs(dist[E2N[*it][0]]);
          double d2 = fabs(dist[E2N[*it][1]]);
          double sum = d1 + d2;
          assert(sum>0);
          (*intersections)[counter++] = (d2*V[E2N[*it][0]] + d1*V[E2N[*it][1]])/sum;
          tip.second = *it;
          xedges.erase(it);
          break;
        }
      }
    } 

    if(++iter>=10) //no way. something is wrong
      break; //avoid infinite loop (should not occur anyway)
  }

  assert(iter<10);
  assert(counter==nPoints);    

  return nPoints;

}

// ------------------------------------------------------------------------------------------------------------


}
