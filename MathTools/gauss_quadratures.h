/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#pragma once
#include<Vector3D.h>

namespace MathTools {

/**********************************************
 * class GaussQuadraturesTriangle stores parameters
 * of Gauss quadratures from Dunavant, IJNME,
 * 1985 (Appendix II).
 *********************************************/
class GaussQuadraturesTriangle {

  static constexpr double one_third  = 1.0/3.0;
  static constexpr double two_thirds = 2.0/3.0;
  static constexpr double one_sixth  = 1.0/6.0;
  static constexpr double w41        = -0.5625;
  static constexpr double w42        = 1.5625/3.0; 
  static constexpr double w61        = 0.223381589678011;
  static constexpr double p611       = 0.108103018168070;
  static constexpr double p612       = (1.0 - p611)/2.0; //0.445948490915965;
  static constexpr double w62        = (1.0 - 3.0*w61)/3.0; //0.109951743655322;
  static constexpr double p621       = 0.816847572980459;
  static constexpr double p622       = (1.0 - p621)/2.0; //0.091576213509771;

public:
  static void GetParameters(int nPoint, double* w, Vec3D* p) {
    switch (nPoint) {       
      case 1: //degree-1
        w[0] = 1.0;
        p[0][0] = p[0][1] = p[0][2] = one_third;
        break;
      case 3: //degree-2
        w[0] = w[1] = w[2] = one_third;
        p[0][0] = two_thirds;  p[0][1] = p[0][2] = one_sixth;
        p[1][1] = two_thirds;  p[1][0] = p[1][2] = one_sixth;
        p[2][2] = two_thirds;  p[2][0] = p[2][1] = one_sixth;
        break;
      case 4: //degree-3
        w[0] = w41;
        w[1] = w[2] = w[3] = w42;
        p[0][0] = p[0][1] = p[0][2] = one_third;
        p[1][0] = 0.6;  p[1][1] = p[1][2] = 0.2;
        p[2][1] = 0.6;  p[2][0] = p[2][2] = 0.2;
        p[3][2] = 0.6;  p[3][0] = p[3][1] = 0.2;
        break;
      case 6: //degree-4
        w[0] = w[1] = w[2] = w61;
        w[3] = w[4] = w[5] = w62;
        for(int i=0; i<3; i++) {
          p[i][i] = p611;
          p[i][(i+1)%3] = p612;
          p[i][(i+2)%3] = p612;
          p[i+3][i] = p621;
          p[i+3][(i+1)%3] = p622;
          p[i+3][(i+2)%3] = p622;
        }
        break;
      default:
        fprintf(stdout,"*** Error: The %d-point Gauss quadrature rule has not been implemented.\n", nPoint);
        exit(-1);
    }
  }

};


}
