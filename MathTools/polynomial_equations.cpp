/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <polynomial_equations.h>

namespace MathTools {

/*************************************************************************
 * Solves ax^2 + bx + c = 0. Returns the number of distinct real roots: 0 or 1 or 2.
 ************************************************************************/
int quadratic_equation_solver(double a, double b, double c, double &x1, double &x2,
                              int *num_complex, double *x1_im, double *x2_im)
{
  if(a==0) {
    if(num_complex) *num_complex = 0;
    return linear_equation_solver(b, c, x1);
  }
  double b1 = b/a;
  double c1 = c/a;
  if(!std::isfinite(b1) || !std::isfinite(c1)) {//a is too small, like 0
    if(num_complex) *num_complex = 0;
    return linear_equation_solver(b, c, x1);
  }

  double coeff2 = b1*b1 - 4.0*c1;

  if(!std::isfinite(coeff2)) {
    if(num_complex) *num_complex = 0;
    return linear_equation_solver(b, c, x1);
  }
    
  if(coeff2>0) {
    if(num_complex) *num_complex = 0;
    double coeff = sqrt(coeff2); 
    x1 = (-b1 + coeff)/2.0;
    double x2_tmp = (-b1 - coeff)/2.0;
    if(x2_tmp == x1)
      return 1; 
    x2 = x2_tmp;
    return 2; 
  }
  else if(coeff2==0) {
    x1 = -b1/2.0;
    if(num_complex) *num_complex = 0;
    return 1;
  }
  else { //two complex roots
    if(!num_complex) //the user does not care about complex roots
      return 0;
    *num_complex = 2;
    if(!x1_im || !x2_im) //incomplete inputs, just return 0
      return 0;
    // Now, do the real calculation
    x1 = x2 = -b1/2;
    *x1_im = sqrt(-coeff2)/2.0;
    *x2_im = -(*x1_im);
    return 0;
  }

  return 0; //will not get here.
}

/*************************************************************************
 * Solves ax^3 + bx^2 + c*x + d = 0. Returns the number of distinct real roots: 0,1,2,3.
 ************************************************************************/
int cubic_equation_solver(double a, double b, double c, double d, double &x1, double &x2, double &x3,
                          int *num_complex, double *x1_im, double *x2_im, double *x3_im)
{
  if(a==0) 
    return quadratic_equation_solver(b,c,d,x1,x2,num_complex,x1_im,x2_im);

  double b1 = b/a, c1 = c/a, d1 = d/a;
  if(!std::isfinite(b1) || !std::isfinite(c1) || !std::isfinite(d1))
    return quadratic_equation_solver(b,c,d,x1,x2,num_complex,x1_im,x2_im);


  double disc, q, r, dum1, s, t, term1, r13;
  q = (3.0*c1 - (b1*b1))/9.0;
  r = -(27.0*d1) + b1*(9.0*c1 - 2.0*(b1*b1));
  r /= 54.0;
  disc = q*q*q + r*r;
  term1 = (b1/3.0);
    
  double x1_real, x2_real, x3_real;

  if (!std::isfinite(disc))
    return quadratic_equation_solver(b,c,d,x1,x2,num_complex,x1_im,x2_im);

  if (disc>0) {// One root real, two are complex
    
    s = r + sqrt(disc);
    s = s<0 ? -cbrt(-s) : cbrt(s);
    t = r - sqrt(disc);
    t = t<0 ? -cbrt(-t) : cbrt(t);
    x1_real = -term1 + s + t;

    x1 = x1_real;

    if(!num_complex) //The user does not care about complex roots
      return 1;
    
    *num_complex = 2;

    if(!x1_im || !x2_im || !x3_im)
      return 1;

    term1 += (s + t)/2.0;
    x3 = x2 = -term1;
    term1 = sqrt(3.0)*(-t + s)/2;
    *x2_im = term1;
    *x3_im = -term1;
    return 1;
  } 
  // The remaining options are all real
  else if (disc == 0) {// All roots real, at least two are equal.
    if(num_complex)
      *num_complex = 0;
    r13 = r<0 ? -cbrt(-r) : cbrt(r);
    x1_real = -term1 + 2.0*r13;
    x2_real = -(r13 + term1);
    if(x1_real == x2_real) {
      x1 = x1_real;
      return 1;
    } else {
      x1 = x1_real;
      x2 = x2_real;
      return 2;
    }
  }
  // Only option left is that all roots are real and unequal (to get here, q < 0)
  else {
    if(num_complex)
      *num_complex = 0;

    double PI = acos(-1.0);
    q = -q;
    dum1 = q*q*q;
    dum1 = acos(r/sqrt(dum1));
    r13 = 2.0*sqrt(q);
    x1_real = -term1 + r13*cos(dum1/3.0);
    x2_real = -term1 + r13*cos((dum1 + 2.0*PI)/3.0);
    x3_real = -term1 + r13*cos((dum1 + 4.0*PI)/3.0);

    int count = 1;
    x1 = x1_real;
    if(x2_real != x1) {
      count++;
      x2 = x2_real;
    }
    if(count==2) {
      if(x3_real != x1 && x3_real != x2) {
        x3 = x3_real;
        count++;
      }
    } else {//count==1
      if(x3_real != x1) {
        x2 = x3_real;
        count++;
      }
    }
    return count;
  }
    
  return 0; //will not get here.
}

//----------------------------------------------------------------------------------------
}
