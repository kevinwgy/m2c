/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#pragma once
#include<cmath>

namespace MathTools {

/****************************************************************************
 * Evaluates the polylogarithm function 
 *   Li_s(z) = \sum_{k=1}^{kmax} ( z^k / k^s )
 ***************************************************************************/
double polylogarithm_function(int s, double z, int kmax, double rel_tol = 0.0)
{
  // closed form expressions for s = 1, 0, -1, -2, -3, -4
  if(s==1) {
    assert(z<1.0);
    return -log(1.0-z);
  }
  else if(s==0) {
    assert(z!=1.0);
    return z/(1.0-z); 
  } 
  else if(s==-1) {
    assert(z!=1.0);
    return z/((1.0-z)*(1.0-z));
  }
  else if(s==-2) {
    assert(z!=1.0);
    return z*(1.0+z)/pow(1.0-z,3);
  }
  else if(s==-3) {
    assert(z!=1.0);
    return z*(1.0+4.0*z+z*z)/pow(1.0-z,4);
  }
  else if(s==-4) {
    assert(z!=1.0);
    return z*(1.0+z)*(1.0+10.0*z+z*z)/pow(1.0-z,5);
  }

  // Now, implementing the general summation formula
  assert(kmax>=1);

  double numerator = z;
  double res = z; //the first term
  double term(0.0);
  for(int k=2; k<=kmax; k++) {
    numerator *= z;
    term = numerator/pow(k,s);
    res += term;
    if(res!=0.0 && fabs(term/res)<=rel_tol)
      break;
  }

  return res;
}

/****************************************************************************
 * Evaluates the derivative of the polylogarithm function 
 *  d(Li_s(z))/dz = (1/z)*Li_(s-1)
 ***************************************************************************/
double polylogarithm_derivative(int s, double z, int kmax, double rel_tol = 0.0)
{
  if(z==0.0)
    return 1.0; 

  return polylogarithm_function(s-1,z,kmax,rel_tol)/z;
}

}
