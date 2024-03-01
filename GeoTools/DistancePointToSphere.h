/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/
#pragma once

#include<Vector3D.h>

namespace GeoTools {

/************************************************************
 * Calculate signed distance from an arbitrary point in 3D to a
 * sphere.
 ***********************************************************/

class DistanceFromPointToSphere{

  Vec3D x0; //!< center of the sphere 
  double R; //!< radius of the sphere

public:

  //! Constructor
  DistanceFromPointToSphere(double *x0_, double R_) : R(R_) {
    for(int i=0; i<3; i++)
      x0[i]  = x0_[i];
  }

  ~DistanceFromPointToSphere() {}

  //! Calculates the signed distance, including the projection point
  double Calculate(double *Q_, double *P = NULL) {
    Vec3D vec = *(Vec3D *)Q_ - x0;
    double dist = vec.norm();
    if(P)
      *(Vec3D *)P = x0 + (dist==0.0 ? Vec3D(1.0,0.0,0.0) // arbitrary dir
                                    : R*vec/dist);
    return dist - R;
  }

};



}; //end of namespace
