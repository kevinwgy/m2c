#include <linear_algebra.h>
#include <Eigen/Dense>

namespace MathTools {

//---------------------------------------------------------------------------------

bool 
LinearAlgebra::SolveLinearSystem2x2(double a11, double a12, double a21, double a22, 
                                    double b1, double b2,
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

//---------------------------------------------------------------------------------

double
LinearAlgebra::CalculateDeterminant3x3(double a11, double a12, double a13, 
                                       double a21, double a22, double a23,
                                       double a31, double a32, double a33)
{
  return a11*a22*a33 + a12*a23*a31 + a13*a21*a32 - a13*a22*a31 - a12*a21*a33 - a11*a23*a32;
}

//---------------------------------------------------------------------------------

double
LinearAlgebra::CalculateDeterminant3x3(Vec3D &a, Vec3D &b, Vec3D &c) //a, b, c are the three columns
{
  return a[0]*b[1]*c[2] + b[0]*c[1]*a[2] + c[0]*a[1]*b[2] - c[0]*b[1]*a[2] - b[0]*a[1]*c[2] - a[0]*c[1]*b[2];
}

//---------------------------------------------------------------------------------

double
LinearAlgebra::CalculateDeterminant3x3(double *A) //column first, i.e. A[0] is A(1,1), A[1] = A(2,1)
{
  return   A[0]*A[4]*A[8] + A[3]*A[7]*A[2] + A[6]*A[1]*A[5]
         - A[6]*A[4]*A[2] - A[3]*A[1]*A[8] - A[0]*A[7]*A[5];
}

//---------------------------------------------------------------------------------

double
LinearAlgebra::CalculateFirstPrincipalInvariant3x3(double *A) //column first, i.e. A[0] is A(1,1), A[1] = A(2,1), ...
{
  return A[0] + A[4] + A[8];
}

//---------------------------------------------------------------------------------

double
LinearAlgebra::CalculateSecondPrincipalInvariant3x3(double *A) //column first, i.e. A[0] is A(1,1), A[1] = A(2,1), ...
{
  return A[0]*A[4] + A[4]*A[8] + A[0]*A[8] - A[3]*A[1] - A[7]*A[5] - A[6]*A[2];
}

//---------------------------------------------------------------------------------

double
LinearAlgebra::CalculateThirdPrincipalInvariant3x3(double *A) //column first, i.e. A[0] is A(1,1), A[1] = A(2,1), ...
{
  return CalculateDeterminant3x3(A);
}

//---------------------------------------------------------------------------------

bool
LinearAlgebra::SolveLinearSystem3x3(Vec3D &a, Vec3D &b, Vec3D &c, Vec3D &d, Vec3D &x) //[a b c][x] = [d]
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

//---------------------------------------------------------------------------------

bool
LinearAlgebra::SolveLinearSystem3x3(double a11, double a12, double a13, double a21, double a22, double a23,
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

//---------------------------------------------------------------------------------

/** A and vectors are both column-first. For example A[2] is A(3,1).
 *  vectors are the right eigenvectors --- the i-th column (i=1,2,3) is
 *  the unit eigenvector corresponding to the i-th eigenvalue.
 *  Note: The eigenvalues are SORTED FROM LOWEST TO HIGHEST. 
 *        They are sorted using ther actual values, NOT THE ABSOLUTE VALUES! */
bool
LinearAlgebra::ComputeEigenSymmetricMatrix3x3(double *A, double *values, double *vectors)
{

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(*(Eigen::Matrix3d *)A); //computation done here.

  if(eigensolver.info() != Eigen::Success)
    return false;

  if(values) {
    Eigen::Vector3d *lam = (Eigen::Vector3d *)values;
    *lam = eigensolver.eigenvalues();
  }

  if(vectors) {
    Eigen::Matrix3d *R = (Eigen::Matrix3d *)vectors;
    *R = eigensolver.eigenvectors();
  }

  return true;
}

//---------------------------------------------------------------------------------

/** A, vectors_real, vectors_imag are all column-first. For example A[2] is A(3,1).
 *  vectors_{real,imag} are the right eigenvectors --- the i-th column (i=1,2,3) is
 *  the unit eigenvector corresponding to the i-th eigenvalue 
 *  NOTE: The eigenvalues are NOT sorted! 
 */
bool
LinearAlgebra::ComputeEigen3x3(double *A, double *values_real, double *values_imag,
                               double *vectors_real, double *vectors_imag)
{

  Eigen::EigenSolver<Eigen::Matrix3d> eigensolver(*(Eigen::Matrix3d *)A); //computation done here.

  if(eigensolver.info() != Eigen::Success)
    return false;

  if(values_real) {
    assert(values_imag);
    Eigen::Vector3cd lam = eigensolver.eigenvalues();
    for(int i=0; i<3; i++) {
      values_real[i] = lam(i).real();
      values_imag[i] = lam(i).imag();
    }
  }

  if(vectors_real) {
    assert(vectors_imag);
    Eigen::Matrix3cd R = eigensolver.eigenvectors();
    for(int j=0; j<3; j++)
      for(int i=0; i<3; i++) {
        vectors_real[j*3+i] = R(i,j).real();
        vectors_imag[j*3+i] = R(i,j).imag();
      }
  }

  return true;
}

//---------------------------------------------------------------------------------


}
