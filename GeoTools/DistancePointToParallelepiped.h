/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#pragma once

#include <GeoTools.h>
#include <float.h> //DBL_MAX

namespace GeoTools {

/************************************************************
 * Calculate distance from an arbitrary point in 3D to a
 * parallelepiped.
 ***********************************************************/

class DistanceFromPointToParallelepiped {

//         F________________ 
//         /\              /\G
//       C/__\____________/E \.     
//        \   \...........\...\D
//         \  /B           \  /
//         O\/______________\/A
//

  //! Parallelpiped info
  Vec3D  face_o[6], face_e1[6], face_e2[6]; 
  double face_area[6];
  Vec3D  face_normal[6]; //unit outward normal

  Vec3D O, OA, OB, OC; //!< points and edges 
   
public:

  //! Constructor
  DistanceFromPointToParallelepiped(double *O_, double *OA_, double *OB_, double *OC_) {
    for(int i=0; i<3; i++) {
      O[i] = O_[i];
      OA[i] = OA_[i];
      OB[i] = OB_[i];
      OC[i] = OC_[i];
    }  
    assert((OA^OB)*OC>0.0); //!< otherwise, orientation is wrong

    Vec3D A = O + OA;
    Vec3D B = O + OB;
    Vec3D C = O + OC;
    Vec3D D = A + OB;
    Vec3D E = A + OC;
    Vec3D F = B + OC;  

    // Face 0: OBA
    face_o[0] = O; face_e1[0] = OB;  face_e2[0] = OA;
    // Face 1: OAC
    face_o[1] = O; face_e1[1] = OA;  face_e2[1] = OC;
    // Face 2: ADE
    face_o[2] = A; face_e1[2] = D-A; face_e2[2] = E-A;
    // Face 3: BFD 
    face_o[3] = B; face_e1[3] = F-B; face_e2[3] = D-B;
    // Face 4: OCB 
    face_o[4] = O; face_e1[4] = OC;  face_e2[4] = OB;
    // Face 5: CEF 
    face_o[5] = C; face_e1[5] = E-C; face_e2[5] = F-C;

    for(int i=0; i<6; i++) {
      face_normal[i] = face_e1[i]^face_e2[i];
      face_area[i] = face_normal[i].norm();
      assert(face_area[i]!=0.0);
      face_normal[i] /= face_area[i];
    }

  }

  ~DistanceFromPointToParallelepiped() {}


  /****************************************************************************
   * Calculate signed distance (>0 outside) from an arbitrary point Q to the 
   * parallelepiped, including a closest point (P) (may not be unique)
   ***************************************************************************/
  double Calculate(double *Q_, double *P = NULL) {

    Vec3D Q(Q_);

    double dist, min_dist = DBL_MAX;
    int sign = -1; //default: inside
    double xi[2];

    // Loop through the 6 faces
    // Because parallelepipes are convex, an interior point must be on the negative side
    // of all the 6 faces.
    for(int i=0; i<6; i++) {
      dist = ProjectPointToParallelogram(Q, face_o[i], face_e1[i], face_e2[i], xi,
                                         &face_area[i], &face_normal[i], true);
      if(fabs(dist)<min_dist) {
        min_dist = fabs(dist);
        if(P) {
          for(int k=0; k<3; k++)
            P[k] = face_o[i][k] + xi[0]*face_e1[i][k] + xi[1]*face_e2[i][k];                            
        }
      }
      if(dist>0.0)
        sign = 1.0;
    }

    return sign*min_dist;
  }

};

}//end of namespace
