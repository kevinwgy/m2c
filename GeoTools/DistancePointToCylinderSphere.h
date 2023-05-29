/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/
#pragma once

#include<GeoTools.h>
#include<vector>
#include<utility> //pair, make_pair

namespace GeoTools {

/************************************************************
 * Calculate signed distance from an arbitrary point in 3D to a
 * cylinder possibly with spherical end cap(s)
 ***********************************************************/

class DistanceFromPointToCylinderSphere {

  Vec3D xcen; //!< center of the entire cylinder (not the base disk)
  Vec3D dir; //!< axis of the cylinder-cone (from base to the other end)
  double L, R; //!< cylinder height and radius;
  bool front_cap, back_cap; //!< whether front and back caps are on

  double Lhalf; //!< L/2.0
  Vec3D xf, xb; //!< x0 +/- Lhalf*dir (center of the front/back face)

  std::vector<std::pair<Vec3D, Vec3D> > lineSegments; //!< boundary of the cylinder


public:

  //! Constructor (x0_base is the center of the BASE DISK)
  DistanceFromPointToCylinderSphere(double *x0_base, double *dir_, double L_, double R_,
                                    bool front_cap_, bool back_cap_)
      : L(L_), R(R_), front_cap(front_cap_), back_cap(back_cap_)
  {       
    for(int i=0; i<3; i++) {
      xcen[i] = x0_base[i];
      dir[i] = dir_[i];
    }
    double norm = dir.norm();
    assert(norm!=0.0);
    dir /= norm;

    Lhalf = L/2.0;

    xcen += Lhalf*dir; //xcen is the center of the cylinder

    xf = xcen + Lhalf*dir; //center of the front face
    xb = xcen - Lhalf*dir; //center of the back face

    // define the geometry
    Vec3D p0(-Lhalf, 0.0, 0.0);
    Vec3D p1(-Lhalf, R, 0.0);
    Vec3D p2(Lhalf, R, 0.0);
    Vec3D p3(Lhalf, 0.0, 0.0);
    lineSegments.push_back(std::make_pair(p1,p2));
    if(!back_cap)
      lineSegments.push_back(std::make_pair(p0,p1));
    if(!front_cap)
      lineSegments.push_back(std::make_pair(p2,p3));
  }

  ~DistanceFromPointToCylinderSphere() {}


  //! Calculates the signed distance, including a projection point (may not be unique) 
  double Calculate(double *Q_, double *P = NULL) {

    double     x = (*(Vec3D *)Q_ - xcen)*dir;
    Vec3D   rdir = *(Vec3D *)Q_ - xcen - x*dir;
    double     r = rdir.norm();
    Vec3D     xp = Vec3D(x,r,0.0);
    if(r!=0.0)
      rdir /= r;

    //calculate unsigned distance from node to the boundary of the cylinder
    double dist, rb(0.0), rf(0.0);
    if(back_cap && x<=-Lhalf) {
      Vec3D vec = *(Vec3D *)Q_ - xb;
      rb = vec.norm();
      dist = fabs(rb - R);
      if(P)
        *(Vec3D *)P = xb + R*(rb==0 ? -dir : vec/rb);
    }
    else if(front_cap && x>= Lhalf) {
      Vec3D vec = *(Vec3D *)Q_ - xf;
      rf = vec.norm();
      dist = fabs(rf - R);
      if(P)
        *(Vec3D *)P = xf + R*(rf==0 ? dir : vec/rf);
    }
    else {
      double dist_tmp, coeff;
      Vec3D xc(0.0);
      dist = DBL_MAX;
      for(int l=0; l<(int)lineSegments.size(); l++) {
        dist_tmp = GeoTools::GetShortestDistanceFromPointToLineSegment(xp,
                                 lineSegments[l].first, lineSegments[l].second, coeff);
        if(dist_tmp<dist) {
          dist = dist_tmp;
          xc = (1.0-coeff)*lineSegments[l].first + coeff*lineSegments[l].second;
        }
      }

      if(P) 
        *(Vec3D *)P = xcen + xc[0]*dir + xc[1]*(r==0 ? 0.0 : rdir);
    }

    //figure out the sign, and update dist
    if( (x>-Lhalf && x<Lhalf && r<R) || (back_cap && x<=-Lhalf && rb<R) ||
       (front_cap && x>=Lhalf && rf<R) ) //inside
      dist = -dist;

    return dist;

  }

};





}; //end of namespace
