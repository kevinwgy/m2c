/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#pragma once

namespace MathTools {

//----------------------------------------------------------------------------------------

template<typename T>
T quadratic_interpolation(double x, double x0, double x1, double x2, T y0, T y1, T y2)
{
  double x01 = x1 - x0, x02 = x2 - x0, x12 = x2 - x1;
  return (x-x1)*(x-x2)/(x01*x02)*y0 - (x-x0)*(x-x2)/(x01*x12)*y1 + (x-x0)*(x-x1)/(x02*x12)*y2;
}

//----------------------------------------------------------------------------------------

template<typename T>
T quadratic_derivative(double x, double x0, double x1, double x2, T y0, T y1, T y2)
{
  double x01 = x1 - x0, x02 = x2 - x0, x12 = x2 - x1;
  return (2.0*x-x1-x2)/(x01*x02)*y0 - (2.0*x-x0-x2)/(x01*x12)*y1 + (2.0*x-x0-x1)/(x02*x12)*y2;
}

//----------------------------------------------------------------------------------------

template<typename T>
T quadratic_derivative(double x0, double x1, double x2, T y0, T y1, T y2) //calc. deriv at x1
{
  double x01 = x1 - x0, x02 = x2 - x0, x12 = x2 - x1;
  return -x12/(x01*x02)*y0 - (x01-x12)/(x01*x12)*y1 + x01/(x02*x12)*y2;
}

//----------------------------------------------------------------------------------------

template<typename T>
T quadratic_second_derivative(double x0, double x1, double x2, T y0, T y1, T y2)
{
  double x01 = x1 - x0, x02 = x2 - x0, x12 = x2 - x1;
  return 2.0*(y0/(x01*x02) - y1/(x01*x12) + y2/(x02*x12));
}

//----------------------------------------------------------------------------------------

} //end of namespace 
