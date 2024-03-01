/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#pragma once

namespace MathTools {

//----------------------------------------------------------------------------------------

template<typename T>
T bilinear_interpolation(T c00, T c01, T c10, T c11, T xi0, T xi1);
{
  T c0 = (1-xi0)*c00 + xi0*c01;
  T c1 = (1-xi0)*c10 + xi0*c11;
  return (1-xi1)*c0 + xi1*c1; 
}

//----------------------------------------------------------------------------------------

template<typename T>
T bilinear_interpolation(T c00, T c01, T c10, T c11, T *x0, T *dx, T *x)
{
  T xi[2];
  for(int i=0; i<2; i++)
    xi[i] = (x[i]-x0[i])/dx[i]; 

  return bilinear_interpolation(c00, c01, c10, c11, xi[0], xi[1]);
}

//----------------------------------------------------------------------------------------

} //end of namespace 
