/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#pragma once
#include<cmath>

namespace MathTools {

/*************************************************************************
 * Solves ax + b = 0. Returns the number of (real) roots: 0 or 1
 ************************************************************************/
inline int linear_equation_solver(double a, double b, double &x)
{
  if(a==0) 
    return 0;

  double x0 = -b/a;
  if(!std::isfinite(x0)) 
    return 0;

  x = x0;
  return 1;
}  
   
/*************************************************************************
 * Solves ax^2 + bx + c = 0. Returns the number of distinct real roots: 0 or 1 or 2.
 ************************************************************************/
int quadratic_equation_solver(double a, double b, double c, double &x1, double &x2,
                              int *num_complex = NULL, double *x1_im = NULL, double *x2_im = NULL);

/*************************************************************************
 * Solves ax^3 + bx^2 + c*x + d = 0. Returns the number of distinct real roots: 0,1,2,3.
 ************************************************************************/
int cubic_equation_solver(double a, double b, double c, double d, double &x1, double &x2, double &x3,
                          int *num_complex = NULL, double *x1_im = NULL, double *x2_im = NULL,
                          double *x3_im = NULL);

//----------------------------------------------------------------------------------------
}
