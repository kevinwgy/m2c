/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <linear_algebra.h>
#include <Eigen/Dense>
#include <cassert>

//---------------------------------------------------------------------------------
// Note: Some of the functions below convert between "double*" and "Eigen::Matrix3d*"
//       or "Eigen::Vector3d*", assuming that (1) Eigen uses a fixed and continuous array
//       to store the matrix/vector, and (2) The Eigen structures do not store anything
//       else. This assumption is valid for Eigen 3, see:
//       https://eigen.tuxfamily.org/dox/classEigen_1_1Matrix.html
//---------------------------------------------------------------------------------

namespace MathTools {

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
  return CalculateMatrixTrace3x3(A);
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
 *  Returns whether eigen decomposition is successful.
 *  Note: The eigenvalues are SORTED FROM LOWEST TO HIGHEST. 
 *        They are sorted using ther actual values, NOT THE ABSOLUTE VALUES!
 *        The input A is supposed to be double[9], storing a symmetric matrix.
 *        The Eigen library only uses the lower-triangular part of A. */
bool
LinearAlgebra::CalculateEigenSymmetricMatrix3x3(double *A, double *values, double *vectors)
{
  int option = vectors ? Eigen::ComputeEigenvectors : Eigen::EigenvaluesOnly;

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(*(Eigen::Matrix3d *)A, option);

  if(eigensolver.info() != Eigen::Success)
    return false;

  assert(values);
  Eigen::Vector3d *lam = (Eigen::Vector3d *)values;
  *lam = eigensolver.eigenvalues();

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
 *  Returns whether eigen decomposition is successful.
 *  NOTE: The eigenvalues are NOT sorted! 
 */
bool
LinearAlgebra::CalculateEigen3x3(double *A, double *values_real, double *values_imag,
                               double *vectors_real, double *vectors_imag)
{

  Eigen::EigenSolver<Eigen::Matrix3d> eigensolver(*(Eigen::Matrix3d *)A,
                                                  vectors_real || vectors_imag);

  if(eigensolver.info() != Eigen::Success)
    return false;

  assert(values_real && values_imag);
  Eigen::Vector3cd lam = eigensolver.eigenvalues();
  for(int i=0; i<3; i++) {
    values_real[i] = lam(i).real();
    values_imag[i] = lam(i).imag();
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

/** A and Ainv are both column-first. For example A[2] is A(3,1).
 *  Returns whether matrix is invertible (up to a default tolerance in Eigen, which
 *  can be changed is needed. */
bool
LinearAlgebra::CalculateMatrixInverseAndDeterminant3x3(double *A, double *Ainv, double *det)
{

  bool invertible = false;
  double det_ = 0.0;
  Eigen::Matrix3d *A_    = (Eigen::Matrix3d *)A;
  Eigen::Matrix3d *Ainv_ = (Eigen::Matrix3d *)Ainv;

  A_->computeInverseAndDetWithCheck(*Ainv_, det_, invertible);

  if(det)
    *det = det_;

  return invertible;

}

//---------------------------------------------------------------------------------
// A, U, and V are all column-first. For example, A[2] is A(3,1).
void
LinearAlgebra::CalculateSVD3x3(double *A, double *svalues, double *U, double *V)
{
  int option = (U ? Eigen::ComputeFullU : 0) | (V ? Eigen::ComputeFullV : 0);

  Eigen::JacobiSVD<Eigen::Matrix3d> svd(*(Eigen::Matrix3d*)A, option);

  Eigen::Vector3d *svalues_ = (Eigen::Vector3d *)svalues;
  *svalues_ = svd.singularValues();

  if(U) {
    Eigen::Matrix3d *U_ = (Eigen::Matrix3d *)U;
    *U_ = svd.matrixU();
  }

  if(V) {
    Eigen::Matrix3d *V_ = (Eigen::Matrix3d *)V;
    *V_ = svd.matrixV();
  }
}

//---------------------------------------------------------------------------------

// A = R*U: R is orthogonal, U is symmetric, positive semi-definite
void
LinearAlgebra::CalculateRightPolarDecomposition3x3(double *A, double *R, double *U)
{
  double S[3];
  CalculateSVD3x3(A, S, R, U);

  Eigen::DiagonalMatrix<double, 3> S_(S[0],S[1],S[2]);
  Eigen::Matrix3d *R_ = (Eigen::Matrix3d*)R;
  Eigen::Matrix3d *U_ = (Eigen::Matrix3d*)U;

  *R_ = (*R_)*(U_->transpose());
  *U_ = (*U_)*S_*(U_->transpose());
}

//---------------------------------------------------------------------------------

// A = U*R: R is orthogonal, U is symmetric, positive semi-definite
void
LinearAlgebra::CalculateLeftPolarDecomposition3x3(double *A, double *U, double *R)
{
  double S[3];
  CalculateSVD3x3(A, S, U, R);

  Eigen::DiagonalMatrix<double, 3> S_(S[0],S[1],S[2]);
  Eigen::Matrix3d *U_ = (Eigen::Matrix3d*)U;
  Eigen::Matrix3d *R_ = (Eigen::Matrix3d*)R;

  *R_ = (*U_)*(R_->transpose());
  *U_ = (*U_)*S_*(U_->transpose());
}

//---------------------------------------------------------------------------------

void
LinearAlgebra::CalculateATransposeA3x3(double *A, double *ATA)
{
  ATA[0] = A[0]*A[0] + A[1]*A[1] + A[2]*A[2];
  ATA[1] = A[3]*A[0] + A[4]*A[1] + A[5]*A[2];
  ATA[2] = A[6]*A[0] + A[7]*A[1] + A[8]*A[2];
  ATA[3] = ATA[1];
  ATA[4] = A[3]*A[3] + A[4]*A[4] + A[5]*A[5];
  ATA[5] = A[6]*A[3] + A[7]*A[4] + A[8]*A[5];
  ATA[6] = ATA[2];
  ATA[7] = ATA[5];
  ATA[8] = A[6]*A[6] + A[7]*A[7] + A[8]*A[8];
}

//---------------------------------------------------------------------------------

void
LinearAlgebra::CalculateATransposeB3x3(double *A, double *B, double *ATB)
{
  ATB[0] = A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
  ATB[1] = A[3]*B[0] + A[4]*B[1] + A[5]*B[2];
  ATB[2] = A[6]*B[0] + A[7]*B[1] + A[8]*B[2];
  ATB[3] = A[0]*B[3] + A[1]*B[4] + A[2]*B[5];
  ATB[4] = A[3]*B[3] + A[4]*B[4] + A[5]*B[5];
  ATB[5] = A[6]*B[3] + A[7]*B[4] + A[8]*B[5];
  ATB[6] = A[0]*B[6] + A[1]*B[7] + A[2]*B[8];
  ATB[7] = A[3]*B[6] + A[4]*B[7] + A[5]*B[8];
  ATB[8] = A[6]*B[6] + A[7]*B[7] + A[8]*B[8];
}

//---------------------------------------------------------------------------------

void
LinearAlgebra::CalculateAATranspose3x3(double *A, double *AAT)
{
  AAT[0] = A[0]*A[0] + A[3]*A[3] + A[6]*A[6];
  AAT[1] = A[1]*A[0] + A[4]*A[3] + A[7]*A[6];
  AAT[2] = A[2]*A[0] + A[5]*A[3] + A[8]*A[6];
  AAT[3] = AAT[1];
  AAT[4] = A[1]*A[1] + A[4]*A[4] + A[7]*A[7];
  AAT[5] = A[2]*A[1] + A[5]*A[4] + A[8]*A[7];
  AAT[6] = AAT[2];
  AAT[7] = AAT[5];
  AAT[8] = A[2]*A[2] + A[5]*A[5] + A[8]*A[8];
}

//---------------------------------------------------------------------------------

void
LinearAlgebra::CalculateABTranspose3x3(double *A, double *B, double *ABT)
{
  ABT[0] = A[0]*B[0] + A[3]*B[3] + A[6]*B[6];
  ABT[1] = A[1]*B[0] + A[4]*B[3] + A[7]*B[6];
  ABT[2] = A[2]*B[0] + A[5]*B[3] + A[8]*B[6];
  ABT[3] = A[0]*B[1] + A[3]*B[4] + A[6]*B[7];
  ABT[4] = A[1]*B[1] + A[4]*B[4] + A[7]*B[7];
  ABT[5] = A[2]*B[1] + A[5]*B[4] + A[8]*B[7];
  ABT[6] = A[0]*B[2] + A[3]*B[5] + A[6]*B[8];
  ABT[7] = A[1]*B[2] + A[4]*B[5] + A[7]*B[8];
  ABT[8] = A[2]*B[2] + A[5]*B[5] + A[8]*B[8];
}

//---------------------------------------------------------------------------------

void
LinearAlgebra::CalculateTranspose3x3(double *A, double *AT)
{
  AT[0] = A[0];
  AT[1] = A[3];
  AT[2] = A[6];
  AT[3] = A[1];
  AT[4] = A[4];
  AT[5] = A[7];
  AT[6] = A[2];
  AT[7] = A[5];
  AT[8] = A[8];
}

//---------------------------------------------------------------------------------

void
LinearAlgebra::CalculateMatrixMatrixProduct3x3(double *A, double *B, double *AB)
{
  for(int i=0; i<3; i++) 
    for(int j=0; j<3; j++)
      AB[3*i+j] = A[j]*B[3*i] + A[3+j]*B[3*i+1] + A[6+j]*B[3*i+2];
}

//---------------------------------------------------------------------------------

void
LinearAlgebra::CalculateMatrixVectorProduct3x3(double *A, double *x, double *Ax)
{
  for(int i=0; i<3; i++) 
    Ax[i] = A[i]*x[0] + A[3+i]*x[1] + A[6+i]*x[2];
}

//---------------------------------------------------------------------------------

void
LinearAlgebra::CalculateMatrixAPlusB3x3(double *A, double *B, double *Sum)
{
  for(int i=0; i<9; i++)
    Sum[i] = A[i] + B[i];
}

//---------------------------------------------------------------------------------

void
LinearAlgebra::CalculateMatrixC1APlusC2B3x3(double C1, double *A, 
                                            double C2, double *B, double *Sum)
{
  for(int i=0; i<9; i++)
    Sum[i] = C1*A[i] + C2*B[i];
}

//---------------------------------------------------------------------------------

void
LinearAlgebra::CalculateCTimesMatrixA3x3(double C, double *A, double *CA)
{
  for(int i=0; i<9; i++)
    CA[i] = C*A[i];
}

//---------------------------------------------------------------------------------

}
