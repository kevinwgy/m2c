/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <linear_algebra.h>
#include <Eigen/Dense>
#include <cassert>

//---------------------------------------------------------------------------------
// Note: Some of the functions below convert between "double*" and "Eigen::Matrix2d*"
//       or "Eigen::Vector2d*", assuming that (1) Eigen uses a fixed and continuous array
//       to store the matrix/vector, and (2) The Eigen structures do not store anything
//       else. This assumption is valid for Eigen 3, see:
//       https://eigen.tuxfamily.org/dox/classEigen_1_1Matrix.html
//---------------------------------------------------------------------------------

namespace MathTools {

//---------------------------------------------------------------------------------

double
LinearAlgebra::CalculateDeterminant2x2(double a11, double a12, double a21, double a22)
{
  return a11*a22 - a12*a21;
}

//---------------------------------------------------------------------------------

double
LinearAlgebra::CalculateDeterminant2x2(Vec2D &a, Vec2D &b) //a, b are the two columns
{
  return a[0]*b[1] - a[1]*b[0];
}

//---------------------------------------------------------------------------------

double
LinearAlgebra::CalculateDeterminant2x2(double *A) //column first, i.e. A[0] is A(1,1), A[1] = A(2,1)
{
  return   A[0]*A[3] - A[1]*A[2];
}

//---------------------------------------------------------------------------------

double
LinearAlgebra::CalculateFirstPrincipalInvariant2x2(double *A) //column first, i.e. A[0] is A(1,1), A[1] = A(2,1), ...
{
  return CalculateMatrixTrace2x2(A);
}

//---------------------------------------------------------------------------------

double
LinearAlgebra::CalculateSecondPrincipalInvariant2x2(double *A) //column first, i.e. A[0] is A(1,1), A[1] = A(2,1), ...
{
  return CalculateDeterminant2x2(A);
}

//---------------------------------------------------------------------------------

bool
LinearAlgebra::SolveLinearSystem2x2(Vec2D &a, Vec2D &b, Vec2D &d, Vec2D &x) //[a b][x] = [d]
{
  double det = CalculateDeterminant2x2(a,b);
  if(det==0)
    return false;

  double one_over_det = 1.0/det;

  double xtmp[2];
  xtmp[0] = one_over_det*CalculateDeterminant2x2(d,b);
  xtmp[1] = one_over_det*CalculateDeterminant2x2(a,d);

  if(!std::isfinite(xtmp[0]) || !std::isfinite(xtmp[1]))
    return false;

  for(int i=0; i<2; i++)
    x[i] = xtmp[i];

  return true;
}

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

/** A and vectors are both column-first. For example A[1] is A(2,1).
 *  vectors are the right eigenvectors --- the i-th column (i=1,2) is
 *  the unit eigenvector corresponding to the i-th eigenvalue.
 *  Returns whether eigen decomposition is successful.
 *  Note: The eigenvalues are SORTED FROM LOWEST TO HIGHEST. 
 *        They are sorted using ther actual values, NOT THE ABSOLUTE VALUES!
 *        The input A is assumed to be double[4], storing a symmetric matrix.
 *        The Eigen library only uses the lower-triangular part of A.*/
bool
LinearAlgebra::CalculateEigenSymmetricMatrix2x2(double *A, double *values, double *vectors)
{
  int option = vectors ? Eigen::ComputeEigenvectors : Eigen::EigenvaluesOnly;

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver(*(Eigen::Matrix2d *)A, option);

  if(eigensolver.info() != Eigen::Success)
    return false;

  assert(values);
  Eigen::Vector2d *lam = (Eigen::Vector2d *)values;
  *lam = eigensolver.eigenvalues();

  if(vectors) {
    Eigen::Matrix2d *R = (Eigen::Matrix2d *)vectors;
    *R = eigensolver.eigenvectors();
  }

  return true;
}

//---------------------------------------------------------------------------------

/** A, vectors_real, vectors_imag are all column-first. For example A[1] is A(2,1).
 *  vectors_{real,imag} are the right eigenvectors --- the i-th column (i=1,2) is
 *  the unit eigenvector corresponding to the i-th eigenvalue 
 *  Returns whether eigen decomposition is successful.
 *  NOTE: The eigenvalues are NOT sorted! 
 */
bool
LinearAlgebra::CalculateEigen2x2(double *A, double *values_real, double *values_imag,
                               double *vectors_real, double *vectors_imag)
{

  Eigen::EigenSolver<Eigen::Matrix2d> eigensolver(*(Eigen::Matrix2d *)A,
                                                  vectors_real || vectors_imag);

  if(eigensolver.info() != Eigen::Success)
    return false;

  assert(values_real && values_imag);
  Eigen::Vector2cd lam = eigensolver.eigenvalues();
  for(int i=0; i<2; i++) {
    values_real[i] = lam(i).real();
    values_imag[i] = lam(i).imag();
  }

  if(vectors_real) {
    assert(vectors_imag);
    Eigen::Matrix2cd R = eigensolver.eigenvectors();
    for(int j=0; j<2; j++)
      for(int i=0; i<2; i++) {
        vectors_real[j*2+i] = R(i,j).real();
        vectors_imag[j*2+i] = R(i,j).imag();
      }
  }

  return true;
}

//---------------------------------------------------------------------------------

/** A and Ainv are both column-first. For example A[1] is A(2,1).
 *  Returns whether matrix is invertible (up to a default tolerance in Eigen, which
 *  can be changed is needed. */
bool
LinearAlgebra::CalculateMatrixInverseAndDeterminant2x2(double *A, double *Ainv, double *det)
{

  bool invertible = false;
  double det_ = 0.0;
  Eigen::Matrix2d *A_    = (Eigen::Matrix2d *)A;
  Eigen::Matrix2d *Ainv_ = (Eigen::Matrix2d *)Ainv;

  A_->computeInverseAndDetWithCheck(*Ainv_, det_, invertible);

  if(det)
    *det = det_;

  return invertible;

}

//---------------------------------------------------------------------------------
// A, U, and V are all column-first. For example, A[1] is A(2,1).
void
LinearAlgebra::CalculateSVD2x2(double *A, double *svalues, double *U, double *V)
{
  int option = (U ? Eigen::ComputeFullU : 0) | (V ? Eigen::ComputeFullV : 0);

  Eigen::JacobiSVD<Eigen::Matrix2d> svd(*(Eigen::Matrix2d*)A, option);

  Eigen::Vector2d *svalues_ = (Eigen::Vector2d *)svalues;
  *svalues_ = svd.singularValues();

  if(U) {
    Eigen::Matrix2d *U_ = (Eigen::Matrix2d *)U;
    *U_ = svd.matrixU();
  }

  if(V) {
    Eigen::Matrix2d *V_ = (Eigen::Matrix2d *)V;
    *V_ = svd.matrixV();
  }
}

//---------------------------------------------------------------------------------

// A = R*U: R is orthogonal, U is symmetric, positive semi-definite
void
LinearAlgebra::CalculateRightPolarDecomposition2x2(double *A, double *R, double *U)
{
  double S[2];
  CalculateSVD2x2(A, S, R, U);

  Eigen::DiagonalMatrix<double, 2> S_(S[0],S[1]);
  Eigen::Matrix2d *R_ = (Eigen::Matrix2d*)R;
  Eigen::Matrix2d *U_ = (Eigen::Matrix2d*)U;

  *R_ = (*R_)*(U_->transpose());
  *U_ = (*U_)*S_*(U_->transpose());
}

//---------------------------------------------------------------------------------

// A = U*R: R is orthogonal, U is symmetric, positive semi-definite
void
LinearAlgebra::CalculateLeftPolarDecomposition2x2(double *A, double *U, double *R)
{
  double S[2];
  CalculateSVD2x2(A, S, U, R);

  Eigen::DiagonalMatrix<double, 2> S_(S[0],S[1]);
  Eigen::Matrix2d *U_ = (Eigen::Matrix2d*)U;
  Eigen::Matrix2d *R_ = (Eigen::Matrix2d*)R;

  *R_ = (*U_)*(R_->transpose());
  *U_ = (*U_)*S_*(U_->transpose());
}

//---------------------------------------------------------------------------------

void
LinearAlgebra::CalculateATransposeA2x2(double *A, double *ATA)
{
  ATA[0] = A[0]*A[0] + A[1]*A[1];
  ATA[1] = A[2]*A[0] + A[3]*A[1];
  ATA[2] = ATA[1];
  ATA[3] = A[2]*A[2] + A[3]*A[3];
}

//---------------------------------------------------------------------------------

void
LinearAlgebra::CalculateATransposeB2x2(double *A, double *B, double *ATB)
{
  ATB[0] = A[0]*B[0] + A[1]*B[1];
  ATB[1] = A[2]*B[0] + A[3]*B[1];
  ATB[2] = A[0]*B[2] + A[1]*B[3];
  ATB[3] = A[2]*B[2] + A[3]*B[3];
}

//---------------------------------------------------------------------------------

void
LinearAlgebra::CalculateAATranspose2x2(double *A, double *AAT)
{
  AAT[0] = A[0]*A[0] + A[2]*A[2];
  AAT[1] = A[1]*A[0] + A[3]*A[2];
  AAT[2] = AAT[1];
  AAT[3] = A[1]*A[1] + A[3]*A[3];
}

//---------------------------------------------------------------------------------

void
LinearAlgebra::CalculateABTranspose2x2(double *A, double *B, double *ABT)
{
  ABT[0] = A[0]*B[0] + A[2]*B[2];
  ABT[1] = A[1]*B[0] + A[3]*B[2];
  ABT[2] = A[0]*B[1] + A[2]*B[3];
  ABT[3] = A[1]*B[1] + A[3]*B[3];
}

//---------------------------------------------------------------------------------

void
LinearAlgebra::CalculateTranspose2x2(double *A, double *AT)
{
  AT[0] = A[0];
  AT[1] = A[2];
  AT[2] = A[1];
  AT[3] = A[3];
}

//---------------------------------------------------------------------------------

void
LinearAlgebra::CalculateMatrixMatrixProduct2x2(double *A, double *B, double *AB)
{
  for(int i=0; i<2; i++) 
    for(int j=0; j<2; j++)
      AB[2*i+j] = A[j]*B[2*i] + A[2+j]*B[2*i+1];
}

//---------------------------------------------------------------------------------

void
LinearAlgebra::CalculateMatrixVectorProduct2x2(double *A, double *x, double *Ax)
{
  for(int i=0; i<2; i++) 
    Ax[i] = A[i]*x[0] + A[2+i]*x[1];
}

//---------------------------------------------------------------------------------

void
LinearAlgebra::CalculateMatrixAPlusB2x2(double *A, double *B, double *Sum)
{
  for(int i=0; i<4; i++)
    Sum[i] = A[i] + B[i];
}

//---------------------------------------------------------------------------------

void
LinearAlgebra::CalculateMatrixC1APlusC2B2x2(double C1, double *A, 
                                            double C2, double *B, double *Sum)
{
  for(int i=0; i<4; i++)
    Sum[i] = C1*A[i] + C2*B[i];
}

//---------------------------------------------------------------------------------

void
LinearAlgebra::CalculateCTimesMatrixA2x2(double C, double *A, double *CA)
{
  for(int i=0; i<4; i++)
    CA[i] = C*A[i];
}

//---------------------------------------------------------------------------------

}
