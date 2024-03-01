/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#pragma once

#include<Vector3D.h>
#include<vector>

namespace GeoTools {

const double INTERSECTIONS_EPSILON = 1.0e-14; //this should be a tolerance smaller than "half_thickness"

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

//! for axis-aligned directions
bool LineSegmentIntersectsTriangle(Vec3D O, int dir, //!< dir = 0 (x-axis), 1 (y-axis), or 2 (z-axis)
                                   double len, //!< length of the line segment (along dir)
                                   Vec3D& V0, Vec3D& V1, Vec3D& V2,
                                   double* d = NULL, Vec3D* xp = NULL, Vec3D* baryCoords = NULL);

//! Checks whether a line segment intersects a plane
bool LineSegmentIntersectsPlane(Vec3D X0, Vec3D X1, //!< vertices of line segment
                                Vec3D V0, Vec3D dir, //!< a point on the plane, and its normal
                                double* d = NULL, //!< optional output: dist from X0 to intersection
                                Vec3D* xp = NULL, //!< optional output: intersection point
                                bool N_normalized = false); //!< input: whether dir is normalized

/** Checks whether a plane cuts an axis-aligned box. If yes, find edge-plane intersections.
 *  "intersections" are ordered such that the the points form the intersection polygon.
 *  Returns the number of intersection points. */
int PlaneCuttingAxisAlignedBox(Vec3D& O, Vec3D& N, //!< point and normal that define the plane
                               Vec3D& Vmin, Vec3D& Vmax, //!< two corners of the axis-aligned box
                               std::vector<Vec3D> *intersections = NULL, //!< output: intersections (dim might change!!)
                               double tolerance = 0.0); //!< tolerance for distance (abs. value)
                               

} //end of namespace
