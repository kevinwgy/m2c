#pragma once

#include<Vector3D.h>

namespace GeoTools {

const double INTERSECTIONS_EPSILON = 1.0e-20;

//! Moller-Trumbore intersection algorithm (Moller and Trumbore, 1997)
bool RayIntersectsTriangle(Vec3D O, Vec3D D, //!< origina and direction of ray
                           Vec3D& V0, Vec3D& V1, Vec3D& V2, //!< vertices of triangle
                           double* d = NULL, //!< optional output: dist from origin to intersection
                           Vec3D* xp = NULL, //!< optional output: intersection point
                           Vec3D* baryCoords = NULL, //!< optional output: barycentric coords of intersection point
                           bool D_normalized = false); //!< input: whether D is normalized)

// for axis-aligned directions
bool AxisIntersectsTriangle(Vec3D O, int dir, //!< dir = 0 (x-axis), 1 (y-axis), or 2 (z-axis)
                            Vec3D& V0, Vec3D& V1, Vec3D& V2,
                            double* d = NULL, Vec3D* xp = NULL, Vec3D* baryCoords = NULL);

//! Moller-Trumbore intersection algorithm (Moller and Trumbore, 1997)
bool LineSegmentIntersectsTriangle(Vec3D X0, Vec3D X1, //!< vertices of line segment
                                   Vec3D& V0, Vec3D& V1, Vec3D& V2, //!< vertices of triangle
                                   double* d = NULL, //!< optional output: dist from X0 to intersection
                                   Vec3D* xp = NULL, //!< optional output: intersection point
                                   Vec3D* baryCoords = NULL); //!< optional output: barycentric coords of intersection point

// for axis-aligned directions
bool LineSegmentIntersectsTriangle(Vec3D O, int dir, //!< dir = 0 (x-axis), 1 (y-axis), or 2 (z-axis)
                                   double len, //!< length of the line segment (along dir)
                                   Vec3D& V0, Vec3D& V1, Vec3D& V2,
                                   double* d = NULL, Vec3D* xp = NULL, Vec3D* baryCoords = NULL);

} //end of namespace
