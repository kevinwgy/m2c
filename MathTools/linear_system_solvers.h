#pragma once

namespace MathTools {

inline bool SolveLinearSystem2x2(double a11, double a12, double a21, double a22, double b1, double b2,
                                 double &x1, double &x2)
{
  double det = a11*a22 - a12*a21;

  if(det==0) {
    x1 = x2 = DBL_MAX;
    return false;
  }   

  x1 = b1*a22 - b2*a12;
  x2 = a11*b2 - b1*a21;

  return true;
}








}
