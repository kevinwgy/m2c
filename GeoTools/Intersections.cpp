#include <Intersections.h>
#include <cassert>
#include <cfloat> //DBL_MAX

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
    D /= D.norm();
  }

  Vec3D E1, E2, P, T, Q;
  double denom,f,u,v;

  E1 = V1 - V0;
  E2 = V2 - V0;
  P = D^E2;
  denom = P*E1;

  if (denom > -INTERSECTIONS_EPSILON && denom < INTERSECTIONS_EPSILON) {
    if(d)  *d  = DBL_MAX;
    if(xp) *xp = DBL_MAX;
    if(baryCoords) *baryCoords = DBL_MAX;
    return false;    // This ray is parallel to this triangle.
  }

  f = 1.0/denom;
  T = O - V0;

  // Calculate u
  u = f*(P*T);
  if (u < 0.0 || u > 1.0) {
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
  if (v < 0.0 || u + v > 1.0) {
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
  if (t > INTERSECTIONS_EPSILON) {
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
    fprintf(stderr,"*** Error: In AxisIntersectsTriangle, dir must be 0, 1, or 2. Found %d.\n", dir);
    exit(-1);
  }

  double denom = p1*q2 - p2*q1;
  if(denom > -INTERSECTIONS_EPSILON && denom < INTERSECTIONS_EPSILON) {
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
  if(u < 0.0 || u > 1.0) {
    if(d || xp || baryCoords) { // calculate optional outputs
      double v = f*(r2*p1 - r1*p2);
      if(d)  *d  = P[dir]*u + Q[dir]*v - R[dir]; 
      if(xp) *xp = (1.0-u-v)*V0 + u*V1 + v*V2;
      if(baryCoords) *baryCoords = Vec3D(1.0-u-v, u, v);
    }
    return false;
  }

  double v = f*(r2*p1 - r1*p2);
  if (v < 0.0 || u + v > 1.0) {
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
  if (t > INTERSECTIONS_EPSILON) {
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
  double t(0.0);
  bool intersect = RayIntersectsTriangle(X0, D, V0, V1, V2, &t, xp, baryCoords, false);

  if(d) *d = t;

  if(!intersect || t>1.0 - INTERSECTIONS_EPSILON)
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

  if(!intersect || t>len - INTERSECTIONS_EPSILON*len)
    return false;
  else
    return true;
}

// ------------------------------------------------------------------------------------------------------------


}
