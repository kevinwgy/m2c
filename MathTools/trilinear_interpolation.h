/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#pragma once

namespace MathTools {

//----------------------------------------------------------------------------------------

template<typename T>
T trilinear_interpolation(T c000, T c001, T c010, T c011, 
                          T c100, T c101, T c110, T c111, T *xi)
{
  T c00 = (1-xi[0])*c000 + xi[0]*c001;
  T c10 = (1-xi[0])*c100 + xi[0]*c101;
  T c01 = (1-xi[0])*c010 + xi[0]*c011;
  T c11 = (1-xi[0])*c110 + xi[0]*c111; 

  T c0 = (1-xi[1])*c00 + xi[1]*c01;
  T c1 = (1-xi[1])*c10 + xi[1]*c11;

  return (1-xi[2])*c0 + xi[2]*c1; 
}

//----------------------------------------------------------------------------------------

template<typename T>
T trilinear_interpolation(T c000, T c001, T c010, T c011, 
                          T c100, T c101, T c110, T c111, T *x0, T *dx, T *x)
{
  T xi[3];
  for(int i=0; i<3; i++)
    xi[i] = (x[i]-x0[i])/dx[i]; 

  return trilinear_interpolation(c000, c001, c010, c011, c100, c101, c110, c111, xi);
}

//----------------------------------------------------------------------------------------
/**************************************************************************
 * Trilinear interpolation
 *   Inputs:
 *     xi = (xi1, xi2, xi3) -- local coordinates the interpolation point
 *     c000, c100, c010, c110, c001, c101, c011, c111 --- 8 nodal values (i,j,k)
 *   Outputs:
 *     return value: interpolated value.
 */
template<typename T>
T trilinear_interpolation(Vec3D &xi, T c000, T c100, T c010, T c110,
                          T c001, T c101, T c011, T c111)
{
  T c00 = c000*(1.0 - xi[0]) + c100*xi[0];
  T c01 = c001*(1.0 - xi[0]) + c101*xi[0];
  T c10 = c010*(1.0 - xi[0]) + c110*xi[0];
  T c11 = c011*(1.0 - xi[0]) + c111*xi[0];

  T c0  = c00*(1.0 - xi[1]) + c10*xi[1];
  T c1  = c01*(1.0 - xi[1]) + c11*xi[1];

  return c0*(1.0 - xi[2]) + c1*xi[2];
}

//----------------------------------------------------------------------------------------

} //end of namespace 
