#pragma once

namespace MathTools {

template<typename T>
T trilinear_interpolation(T c000, T c001, T c010, T c011, 
                          T c100, T c101, T c110, T c111, double *xi)
{
  T c00 = (1-xi[0])*c000 + xi[0]*c001;
  T c10 = (1-xi[0])*c100 + xi[0]*c101;
  T c01 = (1-xi[0])*c010 + xi[0]*c011;
  T c11 = (1-xi[0])*c110 + xi[0]*c111; 

  T c0 = (1-xi[1])*c00 + xi[1]*c01;
  T c1 = (1-xi[1])*c10 + xi[1]*c11;

  return (1-xi[2])*c0 + xi[2]*c1; 
}

template<typename T>
T trilinear_interpolation(T c000, T c001, T c010, T c011, 
                          T c100, T c101, T c110, T c111, double *x0, double *dx, double *x)
{
  double xi[3];
  for(int i=0; i<3; i++)
    xi[i] = (x[i]-x0[i])/dx[i]; 

  return trilinear_interpolation(c000, c001, c010, c011, c100, c101, c110, c111, xi);
}

}
