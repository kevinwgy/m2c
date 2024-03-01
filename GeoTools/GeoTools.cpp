/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <GeoTools.h>
#include <cfloat> //DBL_MAX
#include <polynomial_equations.h>

namespace GeoTools {

/**************************************************************************
 * Project a point onto a plane defined by a triangle in 3D, specified by nodal coordinates
 *   Inputs:
 *     x0 -- the point
 *     xA, xB, xC -- coords of the three vertices of the triangle (order matters!)
 *     area (optional) -- area of the triangle
 *     dir (optional) -- unit normal direction of the triangle
 *   Outputs:
 *     xi[3] -- barycentric coordinates (sum = 1) of the projection point, i.e.
 *              the projection point is at xi[0]*xA + xi[1]*xB + x[2]*xC
 *     return value -- SIGNED distance from the point TO THE PLANE defined by the triangle
 *   Note:
 *     If you do not care about the sign of the returned value, the order of A, B, C, and
 *     the orientation of "dir" do not matter. But be careful, the returned value may
 *     be negative.
 *     If you care about the sign of the returned value, xA, xB, xC have to be ordered
 *     such that the normal of the triangle can be determined based on the right-hand rule.
 *     This is the case even if you explicitly specify "dir".
 */
double ProjectPointToPlane(Vec3D& x0, Vec3D& xA, Vec3D& xB, Vec3D& xC, double xi[3],
                           double* area, Vec3D* dir)
{
  double dist, areaPBC, areaPCA, areaABC;
  Vec3D xp;

/*
  fprintf(stdout,"x0 = %e %e %e : xA = %e %e %e, xB = %e %e %e, xC = %e %e %e. area = %e, dir = %e %e %e.\n",
                  x0[0], x0[1], x0[2], xA[0], xA[1], xA[2], xB[0], xB[1], xB[2], xC[0], xC[1], xC[2], 
                  *area, (*dir)[0], (*dir)[1], (*dir)[2]);
*/
  if(area && dir) {

    //calculate the projection.
    dist = (x0-xA)*(*dir);
    xp = x0 - dist*(*dir);

    //calculate barycentric coords.
    areaPBC = (0.5*(xB-xp)^(xC-xp))*(*dir);
    areaPCA = (0.5*(xC-xp)^(xA-xp))*(*dir);
    areaABC = *area;

  } else {
    Vec3D mydir;
    areaABC = GetNormalAndAreaOfTriangle(xA, xB, xC, mydir);

    //calculate the projection.
    dist = (x0-xA)*mydir;
    xp = x0 - dist*mydir;

    areaPBC = (0.5*(xB-xp)^(xC-xp))*mydir;
    areaPCA = (0.5*(xC-xp)^(xA-xp))*mydir;
  }

  //calculate barycentric coords.
  xi[0] = areaPBC/areaABC;
  xi[1] = areaPCA/areaABC;
  xi[2] = 1.0 - xi[0] - xi[1];

  return dist;
}

/**************************************************************************
 * Project a point onto a triangle, that is, find the closest point on the triangle
 * to the point.
 *   Inputs:
 *     x0 -- the point
 *     xA, xB, xC -- coords of the three vertices of the triangle (order matters!)
 *     area (optional) -- area of the triangle
 *     dir (optional) -- unit normal direction of the triangle
 *     return_signed_distance -- whether the returned distance is the signed distance or
 *                               just the magnitude (usually it should be the latter)
 *   Outputs:
 *     xi[3] -- barycentric coordinates (sum = 1) of the closest point, i.e.
 *              the point is at xi[0]*xA + xi[1]*xB + x[2]*xC. All the coords are in [0, 1]
 *     return value -- Distance from the point to the triangle, unsigned by default
 *   Note:
 *     - Usually, you should avoid getting/using the "sign" of the returned value, as it has
 *        a sharp discontinuity when a point moves across the plane at a point outside the triange.
 *     - If you do not care about the sign of the returned value, the order of A, B, C, and
 *       the orientation of "dir" do not matter. But be careful, if "return_signed_distance" is
 *       mistakenly turned on, the returned value may be negative.
 *     - If you care about the sign of the returned value, xA, xB, xC have to be ordered
 *       such that the normal of the triangle can be determined based on the right-hand rule.
 *       This is the case even if you explicitly specify "dir".
 */
double ProjectPointToTriangle(Vec3D& x0, Vec3D& xA, Vec3D& xB, Vec3D& xC, double xi[3],
                              double* area, Vec3D* dir, 
                              bool return_signed_distance)
{
  
  double dist = ProjectPointToPlane(x0, xA, xB, xC, xi, area, dir);
  int sign = dist>=0 ? 1 : -1; //NOTE: if the point is exactly on the plane, sign = 1


  dist = fabs(dist); // from now all, dist is unsigned distance

  if(xi[0] >= 0.0 && xi[1] >= 0.0 && xi[2] >= 0.0) // projection point is within triangle
    return return_signed_distance ? sign*dist : dist;

  dist = DBL_MAX;
  double d2p[3] = {-1.0, -1.0, -1.0}; //dist to xA, xB, xC
  Vec3D *nodes_ptr[3] = {&xA, &xB, &xC};
  double alpha(0.0), d2l(0.0);
  int p1, p2;
  for(int i=0; i<3; i++) { //check the edges
    if(xi[i]<0) {
      p1 = (i+1)%3;
      p2 = (i+2)%3;
      d2l = ProjectPointToLine(x0, *nodes_ptr[p1], *nodes_ptr[p2], alpha);
      if(alpha >= 0.0 && alpha <= 1.0) { //along edge
        dist   = d2l;
        xi[i]  = 0.0;
        xi[p1] = 1.0-alpha;
        xi[p2] = alpha;
        return return_signed_distance ? sign*dist : dist;
      } 
      else if(alpha < 0.0) {
        if(d2p[p1]<0)
          d2p[p1] = (x0 - *nodes_ptr[p1]).norm(); //dist to point
        if(d2p[p1] < dist) {
          dist   = d2p[p1];
          xi[i]  = 0.0;
          xi[p1] = 1.0;
          xi[p2] = 0.0;
        }
      }
      else { //alpha > 1.0
        if(d2p[p2]<0)
          d2p[p2] = (x0 - *nodes_ptr[p2]).norm(); //dist to point
        if(d2p[p2] < dist) {
          dist   = d2p[p2];
          xi[i]  = 0.0;
          xi[p1] = 0.0;
          xi[p2] = 1.0;
        }  
      }
    }
  }

  return return_signed_distance ? sign*dist : dist;

}



/**************************************************************************
 * Project a point onto a parallelogram. Find the closest point to the point
 *   Inputs:
 *     x0 -- the point
 *     xA -- coords of one vertex of the parallelogram
 *     AB, AC -- the two edges that have xA as a vertex (vector)
 *     area (optional) -- area of the parallelogram 
 *     dir (optional) -- unit normal direction of the parallelogram 
 *     return_signed_distance -- whether the returned distance is the signed distance or
 *                               just the magnitude (usually it should be the latter)
 *   Outputs:
 *     xi[2] -- coordinates of the closest point. Specifically,
 *              the point is at xA + xi[0]*AB + xi[1]*AC. Both coords are in [0, 1]
 *     return value -- Distance from the point to the parallelogram, unsigned by default
 *   Note:
 *     - Usually, you should avoid getting/using the "sign" of the returned value, as it has
 *        a sharp discontinuity when a point moves across the plane at a point outside the parallelogram
 *     - If you do not care about the sign of the returned value, the order of AB and AC and
 *       the orientation of "dir" do not matter. But be careful, if "return_signed_distance" is
 *       mistakenly turned on, the returned value may be negative.
 *     - If you care about the sign of the returned value, AB and AC must be ordered
 *       such that the normal of the triangle can be determined based on the right-hand rule.
 *       This is the case even if you explicitly specify "dir".
 */

double ProjectPointToParallelogram(Vec3D& x0, Vec3D& xA, Vec3D& AB, Vec3D& AC, double xi[2],
                                   double* area, Vec3D* dir, bool return_signed_distance)
{
  double bary[3], dist;

  Vec3D xB = xA + AB;
  Vec3D xC = xA + AC;

  if(area && dir) {
    assert(*area>0.0);
    double tri_area = 0.5*(*area);
    dist = ProjectPointToPlane(x0, xA, xB, xC, bary, &tri_area, dir); //pass the area of triangle ABC
  }
  else {
    Vec3D mydir = AB^AC; //cross product
    double myarea = mydir.norm(); //area of ABCD
    assert(myarea!=0.0);
    mydir = 1.0/myarea*mydir;
    myarea /= 2.0;
    dist = ProjectPointToPlane(x0, xA, xB, xC, bary, &myarea, &mydir); //pass the area of triangle ABC
  }

  int sign = dist>=0 ? 1 : -1; //NOTE: if the point is exactly on the plane, sign = 1
  dist = fabs(dist); // from now all, dist is unsigned distance

  xi[0] = bary[1];
  xi[1] = bary[2];
  if(xi[0] >= 0.0 && xi[0] <= 1.0 && xi[1] >= 0.0 && xi[1] <= 1.0) // projection point is within parallelogram
    return return_signed_distance ? sign*dist : dist;


  Vec3D xD = xB + AC;
  dist = DBL_MAX;
  double d2p[4] = {-1.0, -1.0, -1.0, -1.0}; //dist to xA, xB, xD, xC
  Vec3D *nodes_ptr[4] = {&xA, &xB, &xD, &xC};
  double alpha(0.0), d2l(0.0);
  int p1, p2, p;
  for(int i=0; i<2; i++) { //check the edges
    if(xi[i]<0) {
      p1 = (i+3)%4;
      p2 = i;
      d2l = ProjectPointToLine(x0, *nodes_ptr[p1], *nodes_ptr[p2], alpha);
      if(alpha >= 0.0 && alpha <= 1.0) { //along edge
        dist        = d2l;
        xi[i]       = 0.0;
        xi[(i+1)%2] = (i==1) ? alpha : 1.0-alpha;
        return return_signed_distance ? sign*dist : dist;
      } 
      else {
        p = alpha<0.0 ? p1 : p2;
        if(d2p[p]<0)
          d2p[p] = (x0 - *nodes_ptr[p]).norm(); //dist to point
        if(d2p[p] < dist) {
          dist        = d2p[p];
          xi[i]       = 0.0;
          xi[(i+1)%2] = alpha<0.0 ? 1-i : i;
        }
      }
    }
    else if (xi[i]>1.0) {
      p1 = i+1;
      p2 = i+2;
      d2l = ProjectPointToLine(x0, *nodes_ptr[p1], *nodes_ptr[p2], alpha);
      if(alpha >= 0.0 && alpha <= 1.0) { //along edge
        dist        = d2l;
        xi[i]       = 1.0;
        xi[(i+1)%2] = (i==0) ? alpha : 1.0-alpha;
        return return_signed_distance ? sign*dist : dist;
      } 
      else {
        p = alpha<0.0 ? p1 : p2;
        if(d2p[p]<0)
          d2p[p] = (x0 - *nodes_ptr[p]).norm(); //dist to point
        if(d2p[p] < dist) {
          dist        = d2p[p];
          xi[i]       = 1.0;
          xi[(i+1)%2] = alpha<0.0 ? i : 1-i;
        }
      }
    }
  }

  return return_signed_distance ? sign*dist : dist;
}




/**************************************************************************
 * Find if the distance from a point to a triangle is less than "half_thickness"
 *   Inputs:
 *     x0 -- the point
 *     xA, xB, xC -- coords of the three vertices of the triangle (order matters!)
 *     half_thickness -- half the thickness of the (thickened & widened) triangle
 *     area (optional) -- area of the triangle
 *     dir (optional) -- unit normal direction of the triangle
 */
bool IsPointInsideTriangle(Vec3D& x0, Vec3D& xA, Vec3D& xB, Vec3D& xC, double half_thickness,
                           double* area, Vec3D* dir, double* xi_out)
{
  assert(half_thickness>=0.0);

  double areaABC;
  Vec3D mydir;
  if(area && dir) {
    areaABC = *area;
    mydir = *dir;
  } else
    areaABC = GetNormalAndAreaOfTriangle(xA, xB, xC, mydir);

  double dist = (x0-xA)*mydir;

  if(fabs(dist)>half_thickness)
    return false; //distance from point to plane is already greater than half_thickness
  
  Vec3D xp = x0 - dist*mydir;

  double areaPBC = (0.5*(xB-xp)^(xC-xp))*mydir;
  double areaPCA = (0.5*(xC-xp)^(xA-xp))*mydir;

  //calculate barycentric coords.
  double xi[3];
  xi[0] = areaPBC/areaABC;
  xi[1] = areaPCA/areaABC;
  xi[2] = 1.0 - xi[0] - xi[1];

  if(xi[0]>=0.0 && xi[1]>=0.0 && xi[2]>=0.0) {
    if(xi_out) 
      for(int i=0; i<3; i++) xi_out[i] = xi[i];
    return true;
  }

  // copying part of the code in ProjectPointToTriangle
  double d2p[3] = {-1.0, -1.0, -1.0}; //dist to xA, xB, xC
  Vec3D *nodes_ptr[3] = {&xA, &xB, &xC};
  double alpha(0.0), d2l(0.0);
  int p1, p2;
  for(int i=0; i<3; i++) { //check the edges
    if(xi[i]<0) {  //one or two xi's may be negative
      p1 = (i+1)%3;
      p2 = (i+2)%3;
      d2l = ProjectPointToLine(x0, *nodes_ptr[p1], *nodes_ptr[p2], alpha);
      if(d2l>half_thickness)
        return false;
      if(alpha >= 0.0 && alpha <= 1.0) { //along edge
        if(xi_out) {
          xi_out[i]  = 0.0;
          xi_out[p1] = 1.0-alpha;
          xi_out[p2] = alpha;
        }
        return true;
      }
      else if(alpha < 0.0) {
        if(d2p[p1]<0) { //has not tested this vertex
          d2p[p1] = (x0 - *nodes_ptr[p1]).norm(); //dist to point
          if(d2p[p1] <= half_thickness) {
            if(xi_out) {
              xi[i]  = 0.0;
              xi[p1] = 1.0;
              xi[p2] = 0.0;
            }
            return true;
          }
        } else //tested this vertex and the distance is longer than half_thickness
          return false;
      }
      else { //alpha > 1.0
        if(d2p[p2]<0) { //has not tested this vertex
          d2p[p2] = (x0 - *nodes_ptr[p2]).norm(); //dist to point
          if(d2p[p2] <= half_thickness) {
            if(xi_out) {
              xi[i]  = 0.0;
              xi[p1] = 0.0;
              xi[p2] = 1.0;
            }
            return true;
          }
        } else //tested this vertex and the distance is longer than half_thickness
          return false;
      }
    }
  }

  // I don't think it would ever get here. But I could be wrong.
  return false;
}


/**************************************************************************
 * Check if a ray collides with a moving triangle (continuous collision 
 * detection (CCD)). Vertices of the triangle are assumed to move along 
 * straight lines at constant speeds. 
 *   Inputs:
 *     x0 -- initial position of the point
 *     x  -- current position of the point
 *     A0,B0,C0 -- initial positions of the three vertices of the triangle
 *     A,B,C -- current positions of the three vertices (same order!)
 *     half_thickness (optional) -- half of the thickness of the triangle. Default: 0
 *     area0, dir0 (optional) -- area and normal of the original triangle (A0-B0-C0)
 *     area, dir (optional) -- area and normal of the current triangle (A-B-C)
 *   Output:
 *     time (optional): time of collision, [0, 1]
 *   Ref:
 *     See M2C notes.
 */
bool ContinuousRayTriangleCollision(Vec3D& x0, Vec3D& x, Vec3D& A0, Vec3D& B0, Vec3D& C0, Vec3D& A, Vec3D& B, Vec3D& C,
                                    double* time, double half_thickness, double* area0, Vec3D* dir0,
                                    double* area, Vec3D* dir)
{

  // Step 1. Check if x0 is occluded by A0-B0-C0
  double areaABC0;
  Vec3D mydir0;
  if(area0 && dir0) {
    areaABC0 = *area0;
    mydir0 = *dir0;
  } else
    areaABC0 = GetNormalAndAreaOfTriangle(A0, B0, C0, mydir0);

  if(IsPointInsideTriangle(x0,A0,B0,C0,half_thickness,&areaABC0,&mydir0)) {
    if(*time) *time = 0;
    return true;
  }

  // Step 2. Check if x is occluded by A-B-C
  double areaABC;
  Vec3D mydir;
  if(area && dir) {
    areaABC = *area;
    mydir = *dir;
  } else
    areaABC = GetNormalAndAreaOfTriangle(A, B, C, mydir);

  if(IsPointInsideTriangle(x,A,B,C,half_thickness,&areaABC,&mydir)) {
    if(*time) *time = 1.0;
    return true;
  }

  // Step 3. Check for point-plane collision at 0<t<1
  Vec3D Dx=x-x0, DA=A-A0, DB=B-B0, DC=C-C0;
  Vec3D B0_A0 = B0 - A0;
  Vec3D C0_A0 = C0 - A0;
  Vec3D DB_DA = DB - DA;
  Vec3D DC_DA = DC - DA;

  Vec3D r0 = x0 - A0;
  Vec3D r1 = Dx - DA;
  Vec3D v0 = B0_A0^C0_A0;
  Vec3D v1 = (DB_DA^C0_A0) + (B0_A0^DC_DA);
  Vec3D v2 = DB_DA^DC_DA;

  // form cubic equation at^3 + bt^2 + ct + d = 0;
  double a = r1*v2;
  double b = r0*v2 + r1*v1;
  double c = r0*v1 + r1*v0;
  double d = r0*v0;

  // solve the equation analytically.
  double t[3] = {-1,-1,-1};
  int nReal = MathTools::cubic_equation_solver(a,b,c,d,t[0],t[1],t[2]);
  if(nReal==0)
    return false;

  // Step 4. Check for point-triangle collision at possible time instants
  for(int i=0; i<nReal; i++) {
    if(t[i]<0 || t[i]>1.0)
      continue;
    Vec3D xt = x0+t[i]*Dx;
    Vec3D At = A0+t[i]*DA;
    Vec3D Bt = B0+t[i]*DB;
    Vec3D Ct = C0+t[i]*DC;
    if(IsPointInsideTriangle(xt, At, Bt, Ct, half_thickness)) {
      if(*time) *time = t[i];
      return true;
    }
  }
  return false;  
}


/**************************************************************************
 * For a given vector, find two unit vectors such that the three form an
 * orthonormal basis, satisfying the right hand rule (U0-U1-U2). (If the 
 * given vector is NOT normalized, this function does not change it.
 * The VALUE of U0 is passed in, not a reference.)
 *   Inputs:
 *     U0 -- a given vector
 *     U0_normalized -- (T/F) whether U0 is normalized.
 *   Outputs:
 *     U1, U2: two unit vectors that are orthogonal to each other, and to U0.
 */
void GetOrthonormalVectors(Vec3D U0, Vec3D &U1, Vec3D &U2, bool U0_normalized)
{

  if(!U0_normalized) {
    double norm = U0.norm();
    assert(norm != 0.0);
    U0 /= norm;
  }

  U1 = U2 = 0.0;
  bool done = false;
  for(int i=0; i<3; i++) {
    if(U0[i]==0) {
      U1[i] = 1.0; //got U1
      bool gotU2 = false;
      for(int j=i+1; j<3; j++) {
        if(U0[j]==0) {
          U2[j] = 1.0; //got U2;
          gotU2 = true;
          break;
        }
      }
      if(!gotU2) {
        int i1 = (i+1) % 3;
        int i2 = (i+2) % 3;
        U2[i1] = -U0[i2];
        U2[i2] = U0[i1];
        U2 /= U2.norm();
      }
      done = true;
      break;
    }
  }
  if(!done) { //!< all the three components of U0 are nonzero
    U1[0] = 1.0;
    U1[1] = 0.0;
    U1[2] = -U0[0]/U0[2];
    U1 /= U1.norm();
    U2[0] = 1.0;
    U2[1] = -(U0[2]*U0[2]/U0[0] + U0[0])/U0[1];
    U2[2] = U0[2]/U0[0];
    U2 /= U2.norm();
  }

  if((U0^U1)*U2<0.0) {//swap U1 and U2 to satisfy the right-hand rule
    Vec3D Utmp = U1;
    U1 = U2;
    U2 = Utmp;
  }

}



} //end of namespace
