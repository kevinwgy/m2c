/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/
#include <BoundingBoxes.h>
#include <cfloat> //DBL_MAX

namespace GeoTools {

//---------------------------------------------------------------------

void
GetBoundingBoxOfCylinderCone(Vec3D &p0, Vec3D &n, double r, double L, double cone_height,
                             Vec3D &O, Vec3D &dir0, Vec3D &dir1, Vec3D &dir2,
                             Vec3D &lmin, Vec3D &lmax, //outputs
                             bool dir_normalized, double scaling_factor)
{
  
  // ----------------------------------------------------
  // Step 1: Get three *orthogonal* axes for the object
  // ----------------------------------------------------
  Vec3D U0(0.0), U1(0.0), U2(0.0);
  double h = n.norm();
  assert(h>=0.0);
  U0 = n/h;
  // calculate U1 and U2
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
  
  // ----------------------------------------------------
  // Step 2: Get the 8 vertices of the intermediate bounding box in (U0,U1,U2)
  // ----------------------------------------------------
  Vec3D v[8];
  v[0] = p0 - r*U1 - r*U2;
  v[1] = p0 + r*U1 - r*U2;
  v[2] = p0 - r*U1 + r*U2;
  v[3] = p0 + r*U1 + r*U2;
  Vec3D disp = (L + cone_height)*U0;
  for(int i=4; i<8; i++)
    v[i] = v[i-4] + disp;

  // ----------------------------------------------------
  // Step 3: Get the bounding box in dir0, dir1, dir2 axes that contains the 8 vertices
  // ----------------------------------------------------
  Vec3D dir[3];
  if(dir_normalized) {
    dir[0] = dir0;
    dir[1] = dir1;
    dir[2] = dir2;
  }
  else {
    h = dir0.norm(); assert(h!=0.0);
    dir[0] = dir0/h;
    h = dir1.norm(); assert(h!=0.0);
    dir[1] = dir1/h;
    h = dir2.norm(); assert(h!=0.0);
    dir[2] = dir2/h;
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
  // Step 4: Apply expansion if needed
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

}

//---------------------------------------------------------------------

void
GetBoundingBoxOfCylinderSphere(Vec3D &p0, Vec3D &n, double r, double L, bool front_cap,
                               bool back_cap, Vec3D &O, Vec3D &dir0, Vec3D &dir1, Vec3D &dir2,
                               Vec3D &lmin, Vec3D &lmax, //outputs
                               bool dir_normalized, double scaling_factor)
{
  
  // ----------------------------------------------------
  // Step 1: Get three *orthogonal* axes for the object
  // ----------------------------------------------------
  Vec3D U0(0.0), U1(0.0), U2(0.0);
  double h = n.norm();
  assert(h>=0.0);
  U0 = n/h;
  // calculate U1 and U2
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
  
  // ----------------------------------------------------
  // Step 2: Get the 8 vertices of the intermediate bounding box in (U0,U1,U2)
  // ----------------------------------------------------
  Vec3D v[8];
  v[0] = p0 - r*U1 - r*U2;
  v[1] = p0 + r*U1 - r*U2;
  v[2] = p0 - r*U1 + r*U2;
  v[3] = p0 + r*U1 + r*U2;
  if(back_cap) {
    for(int i=0; i<3; i++)
      v[i] -= r*U0;
  }

  h = L + (back_cap ? r : 0.0) + (front_cap ? r : 0.0);
  Vec3D disp = h*U0;
  for(int i=4; i<8; i++)
    v[i] = v[i-4] + disp;

  // ----------------------------------------------------
  // Step 3: Get the bounding box in dir0, dir1, dir2 axes that contains the 8 vertices
  // ----------------------------------------------------
  Vec3D dir[3];
  if(dir_normalized) {
    dir[0] = dir0;
    dir[1] = dir1;
    dir[2] = dir2;
  }
  else {
    h = dir0.norm(); assert(h!=0.0);
    dir[0] = dir0/h;
    h = dir1.norm(); assert(h!=0.0);
    dir[1] = dir1/h;
    h = dir2.norm(); assert(h!=0.0);
    dir[2] = dir2/h;
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
  // Step 4: Apply expansion if needed
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

}

//---------------------------------------------------------------------

void
GetBoundingBoxOfSphere(Vec3D &p0, double r, Vec3D &O, Vec3D &dir0, Vec3D &dir1, Vec3D &dir2,
                       Vec3D &lmin, Vec3D &lmax, //outputs
                       bool dir_normalized, double scaling_factor)
{
  Vec3D dir[3];
  if(dir_normalized) {
    dir[0] = dir0;
    dir[1] = dir1;
    dir[2] = dir2;
  }
  else {
    double h = dir0.norm(); assert(h!=0.0);
    dir[0] = dir0/h;
    h = dir1.norm(); assert(h!=0.0);
    dir[1] = dir1/h;
    h = dir2.norm(); assert(h!=0.0);
    dir[2] = dir2/h;
  }

  Vec3D Op0 = p0 - O;
  double scaled = scaling_factor*r, projection;
  for(int i=0; i<3; i++) {
    projection = Op0*dir[i];
    lmin[i] = projection - scaled;
    lmax[i] = projection + scaled;
  }
}

//---------------------------------------------------------------------

void
GetBoundingBoxOfParallelepiped(Vec3D &p0, Vec3D &pa, Vec3D &pb, Vec3D &pc,
                               Vec3D &O, Vec3D &dir0, Vec3D &dir1, Vec3D &dir2,
                               Vec3D &lmin, Vec3D &lmax, //outputs
                               bool dir_normalized, double scaling_factor)
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
  Vec3D dir[3];
  if(dir_normalized) {
    dir[0] = dir0;
    dir[1] = dir1;
    dir[2] = dir2;
  }
  else {
    double h = dir0.norm(); assert(h!=0.0);
    dir[0] = dir0/h;
    h = dir1.norm(); assert(h!=0.0);
    dir[1] = dir1/h;
    h = dir2.norm(); assert(h!=0.0);
    dir[2] = dir2/h;
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

}

//---------------------------------------------------------------------

void
GetBoundingBoxOfSpheroid(Vec3D &p0, Vec3D &n, double semi_length, double r,
                         Vec3D &O, Vec3D &dir0, Vec3D &dir1, Vec3D &dir2,
                         Vec3D &lmin, Vec3D &lmax, //outputs
                         bool dir_normalized, double scaling_factor)
{
  
  // ----------------------------------------------------
  // Step 1: Get three *orthogonal* axes for the spheroid
  // ----------------------------------------------------
  Vec3D U0(0.0), U1(0.0), U2(0.0);
  double h = n.norm();
  assert(h>=0.0);
  U0 = n/h;
  // calculate U1 and U2
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
  
  // ----------------------------------------------------
  // Step 2: Get the 8 vertices of the intermediate bounding box in (U0,U1,U2)
  // ----------------------------------------------------
  Vec3D v[8];
  Vec3D disp = semi_length*U0;
  v[0] = p0 - disp - r*U1 - r*U2;
  v[1] = p0 - disp + r*U1 - r*U2;
  v[2] = p0 - disp - r*U1 + r*U2;
  v[3] = p0 - disp + r*U1 + r*U2;
  disp *= 2.0;
  for(int i=4; i<8; i++)
    v[i] = v[i-4] + disp;

  // ----------------------------------------------------
  // Step 3: Get the bounding box in dir0, dir1, dir2 axes that contains the 8 vertices
  // ----------------------------------------------------
  Vec3D dir[3];
  if(dir_normalized) {
    dir[0] = dir0;
    dir[1] = dir1;
    dir[2] = dir2;
  }
  else {
    h = dir0.norm(); assert(h!=0.0);
    dir[0] = dir0/h;
    h = dir1.norm(); assert(h!=0.0);
    dir[1] = dir1/h;
    h = dir2.norm(); assert(h!=0.0);
    dir[2] = dir2/h;
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
  // Step 4: Apply expansion if needed
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

}

//---------------------------------------------------------------------




//---------------------------------------------------------------------




}; //end of namespace
