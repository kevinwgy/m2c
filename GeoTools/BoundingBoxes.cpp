/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/
#include <BoundingBoxes.h>
#include <GeoTools.h>
#include <cfloat> //DBL_MAX

namespace GeoTools {

//NOTE: dir0, dir1, and dir2 may not be orthogonal to each other. They may not have norm = 1 either.

//---------------------------------------------------------------------

void
GetBoundingBoxOfCylinderCone(Vec3D &p0, Vec3D &n, double r, double L,
                             [[maybe_unused]] double tan_alpha,
                             double H, Vec3D &O, Vec3D &dir0, Vec3D &dir1, Vec3D &dir2,
                             Vec3D &lmin, Vec3D &lmax, //outputs
                             double scaling_factor)
{
  
  // ----------------------------------------------------
  // Step 1: Get three *orthonormal* basis vectors for the object
  // ----------------------------------------------------
  Vec3D U0(0.0), U1(0.0), U2(0.0);
  double h = n.norm();
  assert(h>0.0);
  U0 = n/h;
  GetOrthonormalVectors(U0, U1, U2, true);
 
  // ----------------------------------------------------
  // Step 2: Get the intermediate bounding box in (U0,U1,U2)
  // ----------------------------------------------------
  Vec3D P  = p0 - r*U1 - r*U2;
  Vec3D PA = 2.0*r*U1;
  Vec3D PB = 2.0*r*U2;
  Vec3D PC = (L+H)*U0;

  // ----------------------------------------------------
  // Step 3: Get the bounding box in dir0, dir1, dir2 axes that contains the 8 vertices
  // ----------------------------------------------------
  GetBoundingBoxOfParallelepiped(P, PA, PB, PC, O, dir0, dir1, dir2, lmin, lmax,
                                 scaling_factor);

}

//---------------------------------------------------------------------

void
GetBoundingBoxOfCylinderSphere(Vec3D &p0, Vec3D &n, double r, double L, bool front_cap,
                               bool back_cap, Vec3D &O, Vec3D &dir0, Vec3D &dir1, Vec3D &dir2,
                               Vec3D &lmin, Vec3D &lmax, //outputs
                               double scaling_factor)
{

  // ----------------------------------------------------
  // Step 1: Get three *orthonormal* basis vectors for the object
  // ----------------------------------------------------
  Vec3D U0(0.0), U1(0.0), U2(0.0);
  double h = n.norm();
  assert(h>0.0);
  U0 = n/h;
  GetOrthonormalVectors(U0, U1, U2, true);

  // ----------------------------------------------------
  // Step 2: Get the intermediate bounding box in (U0,U1,U2)
  Vec3D P = p0 - r*U1 - r*U2;
  if(back_cap)
    P -= r*U0;
  Vec3D PA = 2.0*r*U1;
  Vec3D PB = 2.0*r*U2;
  Vec3D PC = (L + (back_cap ? r : 0.0) + (front_cap ? r : 0.0))*U0;

  // ----------------------------------------------------
  // Step 3: Get the bounding box in dir0, dir1, dir2 axes that contains the 8 vertices
  // ----------------------------------------------------
  GetBoundingBoxOfParallelepiped(P, PA, PB, PC, O, dir0, dir1, dir2, lmin, lmax,
                                 scaling_factor);

}

//---------------------------------------------------------------------

void
GetBoundingBoxOfSphere(Vec3D &p0, double r, Vec3D &O, Vec3D &dir0, Vec3D &dir1, Vec3D &dir2,
                       Vec3D &lmin, Vec3D &lmax, //outputs
                       double scaling_factor)
{
  Vec3D dir[3];
  dir[0] = dir0;
  dir[1] = dir1;
  dir[2] = dir2;

  Vec3D dir_norm;
  for(int i=0; i<3; i++) {
    dir_norm[i] = dir[i].norm();
    assert(dir_norm[i]!=0.0);
    dir[i] /= dir_norm[i]; //now, dir becomes unit vectors
  }    

  Vec3D Op0 = p0 - O;
  double scaled = scaling_factor*r, projection;
  for(int i=0; i<3; i++) {
    projection = Op0*dir[i];
    lmin[i] = (projection - scaled)/dir_norm[i];
    lmax[i] = (projection + scaled)/dir_norm[i];
  }
}

//---------------------------------------------------------------------

void
GetBoundingBoxOfParallelepiped(Vec3D &p0, Vec3D &pa, Vec3D &pb, Vec3D &pc,
                               Vec3D &O, Vec3D &dir0, Vec3D &dir1, Vec3D &dir2,
                               Vec3D &lmin, Vec3D &lmax, //outputs
                               double scaling_factor)
{

//         F________________
//         /\              /\G
//       C/__\____________/E \.
//        \   \...........\...\D
//         \  /B           \  /
//        P0\/______________\/A
//
  // ----------------------------------------------------
  // Step 1: Get the 8 vertices of the parallelepiped
  // ----------------------------------------------------
  Vec3D v[8];
  v[0] = p0; //P0
  v[1] = p0 + pa; //A
  v[2] = p0 + pb; //B
  v[3] = p0 + pc; //C
  v[4] = v[1] + pb; //D
  v[5] = v[1] + pc; //E
  v[6] = v[2] + pc; //F
  v[7] = v[4] + pc; //G 

  // ----------------------------------------------------
  // Step 2: Get the bounding box in dir0, dir1, dir2 axes that contains the 8 vertices
  // ----------------------------------------------------
  Vec3D dir[3], dir_norm;
  dir[0] = dir0;
  dir[1] = dir1;
  dir[2] = dir2;
  for(int i=0; i<3; i++) {
    dir_norm[i] = dir[i].norm();
    assert(dir_norm[i]!=0.0);
    dir[i] /= dir_norm[i]; //now, dirs have unit norm.
  }

  lmin = DBL_MAX;
  lmax = -DBL_MAX;
  double coord;
  Vec3D Ov;
  for(int i=0; i<8; i++) {
    Ov = v[i]-O;
    for(int j=0; j<3; j++) {
      coord = Ov*dir[j];
      if(coord<lmin[j])
        lmin[j] = coord;
      if(coord>lmax[j])
        lmax[j] = coord;
    }
  }

  for(int i=0; i<3; i++)
    assert(lmax[i]>lmin[i]);
  
  // ----------------------------------------------------
  // Step 3: Apply expansion if needed
  // ----------------------------------------------------
  if(scaling_factor != 1.0) {
    assert(scaling_factor>0.0);
    double mid, wid;
    for(int i=0; i<3; i++) {
      mid = 0.5*(lmin[i] + lmax[i]);
      wid = scaling_factor*(mid - lmin[i]);
      lmin[i] = mid - wid;
      lmax[i] = mid + wid;
    }
  }

  // ----------------------------------------------------
  // Step 4: Get the coordinates in terms of the original dirs
  // ----------------------------------------------------
  for(int i=0; i<3; i++) {
    lmin[i] /= dir_norm[i];
    lmax[i] /= dir_norm[i];
  }

}

//---------------------------------------------------------------------

void
GetBoundingBoxOfSpheroid(Vec3D &p0, Vec3D &n, double semi_length, double r,
                         Vec3D &O, Vec3D &dir0, Vec3D &dir1, Vec3D &dir2,
                         Vec3D &lmin, Vec3D &lmax, //outputs
                         double scaling_factor)
{

  // ----------------------------------------------------
  // Step 1: Get three *orthonormal* basis vectors for the object
  // ----------------------------------------------------
  Vec3D U0(0.0), U1(0.0), U2(0.0);
  double h = n.norm();
  assert(h>0.0);
  U0 = n/h;
  GetOrthonormalVectors(U0, U1, U2, true);
  
  // ----------------------------------------------------
  // Step 2: Get the intermediate bounding box in (U0,U1,U2)
  // ----------------------------------------------------
  Vec3D P  = p0 - semi_length*U0 - r*U1 - r*U2;
  Vec3D PA = 2.0*r*U1;
  Vec3D PB = 2.0*r*U2;
  Vec3D PC = 2.0*semi_length*U0;

  // ----------------------------------------------------
  // Step 3: Get the bounding box in dir0, dir1, dir2 axes that contains the 8 vertices
  // ----------------------------------------------------
  GetBoundingBoxOfParallelepiped(P, PA, PB, PC, O, dir0, dir1, dir2, lmin, lmax,
                                 scaling_factor);

}

//---------------------------------------------------------------------




//---------------------------------------------------------------------




}; //end of namespace
