#pragma once
#include <Vector3D.h>
#include <cassert>

/**************************************************************************
 * This file declares some basic geometry tools
 *************************************************************************/
namespace GeoTools {

/**************************************************************************
 * Check if a point in 2D is within a disk
 *   Inputs:
 *     x,y          -- coords of the point
 *     cen_x, cen_y -- coords of the center of the disk
 *     r            -- radius of the disk
 *   Outputs: true or false
 */
inline bool IsPointInDisk(double x, double y, double cen_x, double cen_y, double r)
{
  return ( (x-cen_x)*(x-cen_x) + (y-cen_y)*(y-cen_y) <= r*r);
}

/**************************************************************************
 * Check if a point in 2D is within a rectangle
 *   Inputs:
 *     x,y          -- coords of the point
 *     cen_x, cen_y -- coords of the center of the rectangle
 *     lx, ly       -- dimensions of the rectangle in x and y directions
 *   Outputs: true or false
 */
inline bool IsPointInRectangle(double x, double y, double cen_x, double cen_y, double lx, double ly)
{
  return x >= cen_x - 0.5*lx  &&  x <= cen_x + 0.5*lx  &&  y >= cen_y - 0.5*ly  &&  y <= cen_y + 0.5*ly;
}

/**************************************************************************
 * Project a point onto a line in 3D, specified using an edge (line segment) 
 *   Inputs:
 *     x0 -- coords of the point
 *     xA -- coords of the start point of the edge
 *     xB -- coords of the end point of the edge
 *   Outputs:
 *     alpha -- affine coordinate (i.e. xA + alpha*AB = projection point)
 *     return value -- distance from the point to the line
 *   Note: If the "edge" is actually a point (i.e. xA = xB), alpha will be 0,
 *         and the distance will be the distance to that point
 */
inline double ProjectPointToLine(Vec3D& x0, Vec3D& xA, Vec3D& xB, double &alpha)
{
  Vec3D AB= xB-xA;
  Vec3D AX = x0-xA;
  double length2 = AB*AB;

  alpha = (length2 != 0) ? AB*AX/length2 : 0.0;
  Vec3D P = xA + alpha*AB;
  return (P-x0).norm();
}

/**************************************************************************
 * Calculate the shortest distance from a point to a line segement in 3D 
 *   (this is unsigned distance, i.e. always positive)
 *   Inputs:
 *     x0 -- coords of the point
 *     xA -- coords of the start point of the edge
 *     xB -- coords of the end point of the edge
 *   Outputs:
 *     alpha -- affine coordinate of the closest point (between 0 and 1)
 *     return value -- shortest distance from x0 to the line segment AB
 *   Note: This function can handle the degenerate case of a point (i.e.
 *         xA = xB)
 */
inline double GetShortestDistanceFromPointToLineSegment(Vec3D& x0, Vec3D& xA, Vec3D& xB,
                                                        double &alpha)
{
  double dist = ProjectPointToLine(x0, xA, xB, alpha);
  if(alpha>1.0) {
    dist = (x0-xB).norm();
    alpha = 1.0; 
  } else if (alpha<0.0 || !std::isfinite(alpha)/*xA=xB*/) {
    dist = (x0-xA).norm();
    alpha = 0.0;
  }
  return dist;
}

/**************************************************************************
 * Calculate the normal direction and area of a triangle (by cross product)
 *   Inputs:
 *     xA, xB, xC -- coords of the three vertices of the triangle
 *                    (the order matters!)
 *   Outputs:
 *     dir -- unit normal dir (xB-xA)^(xC-xA)
 *     return value -- area of the triangle
 */
inline double GetNormalAndAreaOfTriangle(Vec3D& xA, Vec3D& xB, Vec3D& xC, 
                                         Vec3D& dir)
{
  Vec3D ABC = (xB-xA)^(xC-xA); //cross product
  double area= ABC.norm();
  assert(area != 0.0);
  dir = 1.0/area*ABC;
  return area;
}                                  

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
                           double* area = NULL, Vec3D* dir = NULL);

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
 *       be negative.
 *     - If you care about the sign of the returned value, xA, xB, xC have to be ordered
 *       such that the normal of the triangle can be determined based on the right-hand rule.
 *       This is the case even if you explicitly specify "dir".
 */
double ProjectPointToTriangle(Vec3D& x0, Vec3D& xA, Vec3D& xB, Vec3D& xC, double xi[3],
                              double* area = NULL, Vec3D* dir = NULL, 
                              bool return_signed_distance = false);

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
                                double* area = NULL, Vec3D* dir = NULL);

/**************************************************************************
 * Trilinear interpolation
 *   Inputs:
 *     xi = (xi1, xi2, xi3) -- local coordinates the interpolation point 
 *     c000, c100, c010, c110, c001, c101, c011, c111 --- 8 nodal values (i,j,k)
 *   Outputs:
 *     return value: interpolated value.
 */
inline double TrilinearInterpolation(Vec3D &xi, double c000, double c100, double c010, double c110,
                                     double c001, double c101, double c011, double c111)
{
  double c00 = c000*(1.0 - xi[0]) + c100*xi[0];
  double c01 = c001*(1.0 - xi[0]) + c101*xi[0];
  double c10 = c010*(1.0 - xi[0]) + c110*xi[0];
  double c11 = c011*(1.0 - xi[0]) + c111*xi[0];

  double c0  = c00*(1.0 - xi[1]) + c10*xi[1];
  double c1  = c01*(1.0 - xi[1]) + c11*xi[1];

  return c0*(1.0 - xi[2]) + c1*xi[2];
}

/**************************************************************************
 * Find the smallest 3D axis-aligned boundingbox of a circle 
 *   Inputs:
 *     x0 -- center of the circle
 *     dir -- normal direction
 *     r -- radius of circle
 *   Outputs:
 *     xyzmin, xyzmax -- the bounding box, defined by two points: (xmin,ymin,zmin) and (xmax,ymax,zmax)
 */
inline void BoundingBoxOfCircle3D(Vec3D &x0, Vec3D &dir, double r, Vec3D &xyzmin, Vec3D &xyzmax)
{
  Vec3D dir0 = dir / dir.norm();
  double r_axis;
  for(int i=0; i<3; i++) {
    r_axis = sqrt(1.0 - dir0[i]*dir0[i])*r;
    xyzmin[i] = x0[i] - r_axis;
    xyzmax[i] = x0[i] + r_axis;
  }
}

/**************************************************************************
 * Find the smallest 3D axis-aligned boundingbox of a cylinder or truncated cone
 *   Inputs:
 *     x0 -- center of the base
 *     r0 -- radius of the base
 *     x1 -- center of the cap 
 *     r1 -- radius of the cap 
 *     dir -- normal direction
 *   Outputs:
 *     bbmin, bbmax -- the bounding box, defined by two points: (xmin,ymin,zmin) and (xmax,ymax,zmax)
 */
inline void BoundingBoxOfCylinder(Vec3D &x0, double r0, Vec3D &x1, double r1, Vec3D &dir, Vec3D &bbmin, Vec3D &bbmax)
{
  BoundingBoxOfCircle3D(x0, dir, r0, bbmin, bbmax);
  Vec3D xyzmin(0.0), xyzmax(0.0);
  BoundingBoxOfCircle3D(x1, dir, r1, xyzmin, xyzmax);
  for(int i=0; i<3; i++) {
    if(bbmin[i] > xyzmin[i])  bbmin[i] = xyzmin[i];
    if(bbmax[i] < xyzmax[i])  bbmax[i] = xyzmax[i];
  }
}








} //end of namespace
