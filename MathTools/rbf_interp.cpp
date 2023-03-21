#include <rbf_interp.hpp>
#include <iostream>
#include <cmath>
#include <ctime>

namespace MathTools
{

double vec_dot_product(int n, double a1[], double a2[])
{
  double value=0.0;
  for (int i=0; i<n; ++i)
    value += a1[i]*a2[i];

  return value;
}

void phi1 ( int n, double r[], double r0, double v[] )

//****************************************************************************80
//
//  Purpose:
//
//    PHI1 evaluates the multiquadric radial basis function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
//    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
//    Third Edition,
//    Cambridge University Press, 2007,
//    ISBN13: 978-0-521-88068-8,
//    LC: QA297.N866.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double R[N], the radial separation.
//    0 < R.
//
//    Input, double R0, a scale factor.
//
//    Output, double V[N], the value of the radial basis function.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    v[i] = sqrt ( r[i] * r[i] + r0 * r0 );
  }
  return;
}
//****************************************************************************80

void phi2 ( int n, double r[], double r0, double v[] )

//****************************************************************************80
//
//  Purpose:
//
//    PHI2 evaluates the inverse multiquadric radial basis function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
//    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
//    Third Edition,
//    Cambridge University Press, 2007,
//    ISBN13: 978-0-521-88068-8,
//    LC: QA297.N866.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double R[N], the radial separation.
//    0 < R.
//
//    Input, double R0, a scale factor.
//
//    Output, double V[N], the value of the radial basis function.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    v[i] = 1.0 / sqrt ( r[i] * r[i] + r0 * r0 );
  }
  return;
}
//****************************************************************************80

void phi3 ( int n, double r[], double r0, double v[] )

//****************************************************************************80
//
//  Purpose:
//
//    PHI3 evaluates the thin-plate spline radial basis function.
//
//  Discussion:
//
//    Note that PHI3(R,R0) is negative if R < R0.  Thus, for this basis function,
//    it may be desirable to choose a value of R0 smaller than any possible R.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
//    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
//    Third Edition,
//    Cambridge University Press, 2007,
//    ISBN13: 978-0-521-88068-8,
//    LC: QA297.N866.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double R[N], the radial separation.
//    0 < R.
//
//    Input, double R0, a scale factor.
//
//    Output, double V[N], the value of the radial basis function.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    if ( r[i] <= 0.0 )
    {
      v[i] = 0.0;
    }
    else
    {
      v[i] = r[i] * r[i] * log ( r[i] / r0 );
    }
  }
  return;
}
//****************************************************************************80

void phi4 ( int n, double r[], double r0, double v[] )

//****************************************************************************80
//
//  Purpose:
//
//    PHI4 evaluates the gaussian radial basis function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
//    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
//    Third Edition,
//    Cambridge University Press, 2007,
//    ISBN13: 978-0-521-88068-8,
//    LC: QA297.N866.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double R[N], the radial separation.
//    0 < R.
//
//    Input, double R0, a scale factor.
//
//    Output, double V[N], the value of the radial basis function.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    v[i] = exp ( - 0.5 * r[i] * r[i] / r0 / r0 );
  }
  return;
}

void rbf_weight ( int m, int nd, double xd[], double r0,
  void phi ( int n, double r[], double r0, double v[] ),
  double fd[], double w[] )
{

  Eigen::MatrixXd A(nd, nd);
  Eigen::VectorXd b(nd,  1);

  double *r = new double[nd];
  double *v = new double[nd];

  // populate linear system to calculate weights
  for (int i=0; i<nd; ++i) {
    for (int j=0; j<nd; ++j) {
   
      // calculate radial distance   
      r[j] = 0.0;
      for (int dim=0; dim<m; ++dim)
        r[j] += pow(xd[m*i+dim] - xd[m*j+dim], 2); 
   
      r[j] = sqrt(r[j]);
     
    }

    // basis function matrix
    phi(nd, r, r0, v);
    for (int j=0; j<nd; ++j)
      A(i,j) = v[j];

    // target values
    b(i) = fd[i];
  }

  delete [] r;
  delete [] v;

  Eigen::VectorXd x(nd, 1);
  x = A.colPivHouseholderQr().solve(b);

  for (int j=0; j<nd; ++j)
    w[j] = x(j); 

}

void rbf_interp ( int m, int nd, double xd[], double r0,
  void phi ( int n, double r[], double r0, double v[] ), double w[],
  int ni, double xi[], double fi[] )
{

  double *r = new double[nd];
  double *v = new double[nd];

  for (int i=0; i<ni; ++i) {
    for (int j=0; j<nd; ++j) {
      r[j] = 0.0;
      for (int dim=0; dim<m; ++dim)
        r[j] += pow( xi[m*i+dim] - xd[m*j+dim], 2);

      r[j] = sqrt(r[j]);
    }

    phi(nd, r, r0, v);

    fi[i] = vec_dot_product(nd, v, w);
  }

  delete [] r;
  delete [] v;

}

} // end of namespace
