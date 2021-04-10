#pragma once
#include <Vector3D.h>

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
 *   Note: If the "edge" is actually a point (i.e. xA = xB), alpha and 
 *         the distance will both be nan.
 */
inline double ProjectPointToLine(Vec3D& x0, Vec3D& xA, Vec3D& xB, double &alpha)
{
  Vec3D AB= xB-xA;
  Vec3D AX = x0-xA;
  alpha = AB*AX/(AB*AB);
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
  } else if (alpha<0.0 || isnan(alpha)/*xA=xB*/) {
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
  dir = 1.0/area*ABC;
  return area;
}                                  

/**************************************************************************
 * Project a point onto a triangle in 3D, specified by nodal coordinates
 *   Inputs:
 *     x0 -- the point
 *     xA, xB, xC -- coords of the three vertices of the triangle (order matters!)
 *     area (optional) -- area of the triangle
 *     dir (optional) -- unit normal direction of the triangle
 *   Outputs:
 *     xi1, xi2 -- barycentric coordinates (the third one is 1 - xi1 - xi2) of 
 *                 the projection point
 *     return value -- distance from the point TO THE PLANE defined by the triangle
 */
inline double ProjectPointToTriangle(Vec3D& x0, Vec3D& xA, Vec3D& xB, Vec3D& xC, 
                                     double& xi1, double& xi2,
                                     double* area = NULL, Vec3D* dir = NULL)
{
  double dist, areaPBC, areaPCA;
  Vec3D xp;

  if(area && dir) {

    //calculate the projection.
    dist = (x0-xA)*(*dir);
    xp = x0 - dist*(*dir);

    //calculate barycentric coords.
    areaPBC = (((xB-xp)^(xC-xp))*(*dir));
    areaPCA = (((xC-xp)^(xA-xp))*(*dir));
    xi1 = areaPBC/(*area);
    xi2 = areaPCA/(*area);

  } else {

    Vec3D mydir;
    double areaABC = GetNormalAndAreaOfTriangle(xA, xB, xC, mydir);

    //calculate the projection.
    dist = (x0-xA)*mydir;
    xp = x0 - dist*mydir;

    //calculate barycentric coords.
    areaPBC = (((xB-xp)^(xC-xp))*mydir);
    areaPCA = (((xC-xp)^(xA-xp))*mydir);
    xi1 = areaPBC/areaABC;
    xi2 = areaPCA/areaABC;
  }

  return dist;
}



} //end of namespace
