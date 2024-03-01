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
 * cylinder-cone.
 ***********************************************************/

class DistanceFromPointToCylinderCone {

/***********************************
  __________________ 
 |          L      |\
 |__dir            |_\cone height
 |x0      |R       | /_ opening angle
 |________|________|/_/_ _
***********************************/

  Vec3D x0; //!< center of the base disk of the cylinder
  Vec3D dir; //!< axis of the cylinder-cone (from base to the other end)
  double L, R; //!< cylinder height and radius;
  double tan_alpha; //!< tangent of the opening angle
  double H; //!< cone height (Note: the cone might have a truncated head)

  double Hmax; //!< = R/tan_alpha (greater than or equal to H)

  std::vector<std::pair<Vec3D, Vec3D> > lineSegments; //!< boundary of the cylinder-cone


public:

  //! Constructor
  DistanceFromPointToCylinderCone(double *x0_, double *dir_, double L_, double R_,
                                  double tan_alpha_, double H_) 
      : L(L_), R(R_), tan_alpha(tan_alpha_), H(H_)
  {       
    for(int i=0; i<3; i++) {
      x0[i] = x0_[i];
      dir[i] = dir_[i];
    }
    double norm = dir.norm();
    assert(norm!=0.0);
    dir /= norm;

    Hmax = R/tan_alpha;

    // define the geometry
    Vec3D p0(0.0, 0.0, 0.0); //x0
    Vec3D p1(0.0, R, 0.0);
    Vec3D p2(L, R, 0.0);
    Vec3D p3(L+H, (Hmax-H)*tan_alpha, 0.0);
    Vec3D p4(L+H, 0.0, 0.0);
    lineSegments.push_back(std::make_pair(p0,p1));
    lineSegments.push_back(std::make_pair(p1,p2));
    lineSegments.push_back(std::make_pair(p2,p3));
    lineSegments.push_back(std::make_pair(p3,p4));
  }

  ~DistanceFromPointToCylinderCone() {}

  //! Calculates the signed distance, including a projection point (may not be unique) 
  double Calculate(double *Q_, double *P = NULL) {

    double     x = (*(Vec3D *)Q_ - x0)*dir;
    Vec3D   rdir = *(Vec3D *)Q_ - x0 - x*dir;
    double     r = rdir.norm();
    Vec3D     xp = Vec3D(x,r,0.0);
    if(r!=0.0)
      rdir /= r;

    //calculate unsigned distance from node to the boundary of the cylinder-cone
    double dist = DBL_MAX, dist_tmp, coeff;
    Vec3D xc(0.0);
    for(int l=0; l<(int)lineSegments.size(); l++) {
      dist_tmp = GeoTools::GetShortestDistanceFromPointToLineSegment(xp,
                                lineSegments[l].first, lineSegments[l].second, coeff); //unsigned dist!
      if(dist_tmp<dist) {
        dist = dist_tmp;
        xc = (1.0-coeff)*lineSegments[l].first + coeff*lineSegments[l].second;
      }
    }

    if( (x>0 && x<L && r<R) || (x>=L && x<L+H && r<(L+Hmax-x)*tan_alpha) ) //inside
      dist = -dist; //getting the signed distance

    if(P) 
      *(Vec3D *)P = x0 + xc[0]*dir + xc[1]*(r==0 ? 0.0 : rdir);

    return dist;
  }

};





}; //end of namespace
