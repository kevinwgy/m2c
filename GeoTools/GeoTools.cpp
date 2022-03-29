#include <GeoTools.h>
#include <cfloat> //DBL_MAX

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

  if(area && dir) {

    //calculate the projection.
    dist = (x0-xA)*(*dir);
    xp = x0 - dist*(*dir);

    //calculate barycentric coords.
    areaPBC = (((xB-xp)^(xC-xp))*(*dir));
    areaPCA = (((xC-xp)^(xA-xp))*(*dir));
    areaABC = *area;

  } else {
    Vec3D mydir;
    areaABC = GetNormalAndAreaOfTriangle(xA, xB, xC, mydir);

    //calculate the projection.
    dist = (x0-xA)*mydir;
    xp = x0 - dist*mydir;

    areaPBC = (((xB-xp)^(xC-xp))*mydir);
    areaPCA = (((xC-xp)^(xA-xp))*mydir);
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
 *       the orientation of "dir" do not matter. But be careful, the returned value may
 *       be negative, if "return_signed_distance" is (mistakenly) turned on.
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


  dist = abs(dist); // from now all, dist is unsigned distance

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
 * Find if the distance from a point to a triangle is less than "half_thickness"
 *   Inputs:
 *     x0 -- the point
 *     xA, xB, xC -- coords of the three vertices of the triangle (order matters!)
 *     half_thickness -- half the thickness of the (thickened & widened) triangle
 *     area (optional) -- area of the triangle
 *     dir (optional) -- unit normal direction of the triangle
 */
bool IsPointInThickenedTriangle(Vec3D& x0, Vec3D& xA, Vec3D& xB, Vec3D& xC, double half_thickness,
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

  double areaPBC = (((xB-xp)^(xC-xp))*mydir);
  double areaPCA = (((xC-xp)^(xA-xp))*mydir);

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

} //end of namespace
