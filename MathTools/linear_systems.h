#pragma once
#include <Vector3D.h>

namespace MathTools {

struct LinearSystems {

static bool 
SolveLinearSystem2x2(double a11, double a12, double a21, double a22, double b1, double b2,
                     double &x1, double &x2)
{
  double det = a11*a22 - a12*a21;

  if(det==0) 
    return false;

  double x1_tmp = (b1*a22 - b2*a12)/det;
  double x2_tmp = (a11*b2 - b1*a21)/det;

  if(!std::isfinite(x1_tmp) || !std::isfinite(x2_tmp))
    return false;

  x1 = x1_tmp;
  x2 = x2_tmp;
  return true;
}


static double
CalculateDeterminant3x3(double a11, double a12, double a13, double a21, double a22, double a23,
                        double a31, double a32, double a33)
{
  return a11*a22*a33 + a12*a23*a31 + a13*a21*a32 - a13*a22*a31 - a12*a21*a33 - a11*a23*a32;
}


static double
CalculateDeterminant3x3(Vec3D &a, Vec3D &b, Vec3D &c) //a, b, c are the three columns
{
  return a[0]*b[1]*c[2] + b[0]*c[1]*a[2] + c[0]*a[1]*b[2] - c[0]*b[1]*a[2] - b[0]*a[1]*c[2] - a[0]*c[1]*b[2];
}


static bool
SolveLinearSystem3x3(Vec3D &a, Vec3D &b, Vec3D &c, Vec3D &d, Vec3D &x) //[a b c][x] = [d]
{
  double det = CalculateDeterminant3x3(a,b,c);
  if(det==0)
    return false;

  double one_over_det = 1.0/det;

  double xtmp[3];
  xtmp[0] = one_over_det*CalculateDeterminant3x3(d,b,c);
  xtmp[1] = one_over_det*CalculateDeterminant3x3(a,d,c);
  xtmp[2] = one_over_det*CalculateDeterminant3x3(a,b,d);

  if(!std::isfinite(xtmp[0]) || !std::isfinite(xtmp[1]) || !std::isfinite(xtmp[2]))
    return false;

  for(int i=0; i<3; i++)
    x[i] = xtmp[i];

  return true;
}


static bool
SolveLinearSystem3x3(double a11, double a12, double a13, double a21, double a22, double a23,
                     double a31, double a32, double a33, double b1, double b2, double b3,
                     double &x1, double &x2, double &x3)
{
  Vec3D a(a11,a21,a31), b(a12,a22,a32), c(a13,a23,a33), d(b1,b2,b3); 
  Vec3D x(0,0,0);
  bool success = SolveLinearSystem3x3(a,b,c,d,x);
  if(success) {
    x1 = x[0];
    x2 = x[1];
    x3 = x[2];
  }
  return success;
}

}; //end of class

}
