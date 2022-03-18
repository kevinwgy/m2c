#pragma once

#include<Vec3D>

namespace GeoTools {

// Moller-Trumbore intersection
// copied from Wikepedia, with modifications

bool RayIntersectsTriangleMT(Vec3D rayOrigin, Vec3D rayVector, 
                             Vec3D& vertex0, Vec3D& vertex1, Vec3D& vertex2,
                             double* d, Vec3D* outIntersectionPoint, 
                             bool rayVector_normalized = false)
{
    const double EPSILON = 1.0e-14;

    if(!rayVector_normalized)
      rayVector /= rayVector.norm();

    Vec3D edge1, edge2, h, s, q;
    double a,f,u,v;
    edge1 = vertex1 - vertex0;
    edge2 = vertex2 - vertex0;
    h = rayVector^edge2;
    a = edge1*h;
    if (a > -EPSILON && a < EPSILON)
        return false;    // This ray is parallel to this triangle.
    f = 1.0/a;
    s = rayOrigin - vertex0;
    u = f * (s*h);
    if (u < 0.0 || u > 1.0)
        return false;
    q = s^edge1;
    v = f * (rayVector*q);
    if (v < 0.0 || u + v > 1.0)
        return false;
    // At this stage we can compute t to find out where the intersection point is on the line.
    double t = f * (edge2*q);
    if (t > EPSILON) {// ray intersection
      if(d)
        *d = t;
      if(outIntersectionPoint)
        *outIntersectionPoint = rayOrigin + rayVector * t;
      return true;
    } else // This means that there is a line intersection but not a ray intersection.
      return false;
}

// for axis-aligned directions
bool RayIntersectsTriangleMT(Vec3D rayOrigin, Vec3D rayVector, 
                             Vec3D& vertex0, Vec3D& vertex1, Vec3D& vertex2,
                             Vec3D* outIntersectionPoint, bool rayVector_normalized = false)
{
  I AM HERE
} //end of namespace
