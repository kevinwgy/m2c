/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#pragma once
#include <Vector3D.h>
#include <cassert>

/**************************************************************************
 * Functions for finding bounding boxes for different shapes. The axes of
 * the bounding box can be specified by the user, and do not have to be
 * orthogonal to each other. The norm of each basis vector may not be 1.
 * Note: The result may NOT be the smallest bounding box.
 *************************************************************************/

namespace GeoTools {

/**************************************************************************
 * Finds a bounding box for a cylinder-cone in user-specified axes
 *   Inputs:
 *     p0: center of the base of the cylinder (in x,y,z coords)
 *     n: axis of the cylinder (in x,y,z coords)
 *     r, L: radius and height of the cylinder
 *     tan_alpha: tangent of the opening angle
 *     H: height of the cone (on top of the cylinder)
 *     O, dir0, dir1, dir2: origin and basis vectors of the user-specified axes.
 *                          the basis vectors may not have length = 1
 *     scaling_factor: scale up/down the b.b. Default is 1.0. 
 *   Output:
 *     lmin, lmax: the min and max of the bounding box coordinates (in the
 *     user-specified axes)
 *   Note: The M2C & A2C solvers allow the user to specify cone height AND/OR
 *         opening angle. When calling this function, the actually cone height
 *         should be passed in. 
 **************************************************************************/
void
GetBoundingBoxOfCylinderCone(Vec3D &p0, Vec3D &n, double r, double L, 
                             [[maybe_unused]] double tan_alpha,
                             double H, Vec3D &O, Vec3D &dir0, Vec3D &dir1, Vec3D &dir2,
                             Vec3D &lmin, Vec3D &lmax, //outputs
                             double scaling_factor = 1.0);


/**************************************************************************
 * Finds a bounding box for a cylinder-sphere in user-specified axes
 *   Inputs:
 *     p0: center of the base of the cylinder (in x,y,z coords)
 *     n: axis of the cylinder (in x,y,z coords)
 *     r, L: radius and height of the cylinder
 *     front_cap: (T/F) whether the front/top cap is on
 *     back_cap: (T/F) whether the back/bottom cap is on
 *     O, dir0, dir1, dir2: origin and basis vectors of the user-specified axes.
 *                          the basis vectors may not have length = 1
 *     scaling_factor: scale up/down the b.b. Default is 1.0. 
 *   Output:
 *     lmin, lmax: the min and max of the bounding box coordinates (in the
 *     user-specified axes)
 **************************************************************************/
void
GetBoundingBoxOfCylinderSphere(Vec3D &p0, Vec3D &n, double r, double L, bool front_cap,
                               bool back_cap, Vec3D &O, Vec3D &dir0, Vec3D &dir1, Vec3D &dir2,
                               Vec3D &lmin, Vec3D &lmax, //outputs
                               double scaling_factor = 1.0);


/**************************************************************************
 * Finds a bounding box for a sphere in user-specified axes
 *   Inputs:
 *     p0: center of the sphere
 *     r: radius of the sphere
 *     O, dir0, dir1, dir2: origin and basis vectors of the user-specified axes.
 *                          the basis vectors may not have length = 1
 *     scaling_factor: scale up/down the b.b. Default is 1.0. 
 *   Output:
 *     lmin, lmax: the min and max of the bounding box coordinates (in the
 *     user-specified axes)
 **************************************************************************/
void
GetBoundingBoxOfSphere(Vec3D &p0, double r, Vec3D &O, Vec3D &dir0, Vec3D &dir1, Vec3D &dir2,
                       Vec3D &lmin, Vec3D &lmax, //outputs
                       double scaling_factor = 1.0);


/**************************************************************************
 * Finds a bounding box for a parallelepiped in user-specified axes
 *   Inputs:
 *     p0: the ``origin'' (one of the 8 vertices).
 *     pa, pb, pc: vectors from the origin to three edges, which defines the parallelepiped
 *     O, dir0, dir1, dir2: origin and basis vectors of the user-specified axes.
 *                          the basis vectors may not have length = 1
 *     scaling_factor: scale up/down the b.b. Default is 1.0. 
 *   Output:
 *     lmin, lmax: the min and max of the bounding box coordinates (in the
 *     user-specified axes)
 **************************************************************************/
void
GetBoundingBoxOfParallelepiped(Vec3D &p0, Vec3D &pa, Vec3D &pb, Vec3D &pc, 
                               Vec3D &O, Vec3D &dir0, Vec3D &dir1, Vec3D &dir2,
                               Vec3D &lmin, Vec3D &lmax, //outputs
                               double scaling_factor = 1.0);


/**************************************************************************
 * Finds a bounding box for a spheroid in user-specified axes
 *   Inputs:
 *     p0: the center of the spheroid
 *     n: revolutionary axis of the spheroid
 *     semi_length, r: half length along axis & max radius on transverse plane
 *     O, dir0, dir1, dir2: origin and basis vectors of the user-specified axes.
 *                          the basis vectors may not have length = 1
 *     scaling_factor: scale up/down the b.b. Default is 1.0. 
 *   Output:
 *     lmin, lmax: the min and max of the bounding box coordinates (in the
 *     user-specified axes)
 **************************************************************************/
void
GetBoundingBoxOfSpheroid(Vec3D &p0, Vec3D &n, double semi_length, double r,
                         Vec3D &O, Vec3D &dir0, Vec3D &dir1, Vec3D &dir2,
                         Vec3D &lmin, Vec3D &lmax, //outputs
                         double scaling_factor = 1.0);


}; //end of namespace

