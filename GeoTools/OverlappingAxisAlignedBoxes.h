/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#pragma once
#include <Vector2D.h>
#include <Vector3D.h>
#include <cassert>

/**************************************************************************
 * Functions for checking whether two axis-aligned boxes (1, 2, or 3D) overlap
 *************************************************************************/

namespace GeoTools {

inline bool BoxesOverlapping1D(double b1min, double b1max, double b2min, double b2max) {
  return b1max <= b2min || b2max <= b1min;
}

bool BoxesOverlapping2D(Vec2D &b1min, Vec2D &b1max, Vec2D &b2min, Vec2D &b2max) {
  for(int i=0; i<2; i++)
    if(b1max[i] < b2min[i] || b2max[i] < b1min[i])
      return false;
  return true;
}

bool BoxesOverlapping3D(Vec3D &b1min, Vec3D &b1max, Vec3D &b2min, Vec3D &b2max) {
  for(int i=0; i<3; i++)
    if(b1max[i] < b2min[i] || b2max[i] < b1min[i])
      return false;
  return true;
} 

};
