/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#pragma once
#include <Vector2D.h>
#include <Vector3D.h>

namespace MathTools {

struct LinearAlgebra {

  //! ----------------------------
  //! 3D vectors and matrices
  //! ----------------------------

  static double
  CalculateDeterminant3x3(double a11, double a12, double a13, double a21, double a22, double a23,
                          double a31, double a32, double a33);

  static double
  CalculateDeterminant3x3(Vec3D &a, Vec3D &b, Vec3D &c); //a, b, c are the three columns

  static double
  CalculateDeterminant3x3(double *A); //column first, i.e. A[0] is A(1,1), A[1] = A(2,1)

  static void 
  CalculateTranspose3x3(double *A, double *AT);
  
  static void
  CalculateMatrixMatrixProduct3x3(double *A, double *B, double *AB);

  static void
  CalculateMatrixVectorProduct3x3(double *A, double *x, double *Ax);

  static void
  CalculateMatrixAPlusB3x3(double *A, double *B, double *Sum);

  static void
  CalculateMatrixC1APlusC2B3x3(double C1, double *A, double C2, double *B, double *Sum);

  static void
  CalculateCTimesMatrixA3x3(double C, double *A, double *CA);

  static double
  CalculateFirstPrincipalInvariant3x3(double *A); //column first, i.e. A[0] is A(1,1), A[1] = A(2,1), ...

  static double
  CalculateSecondPrincipalInvariant3x3(double *A); //column first, i.e. A[0] is A(1,1), A[1] = A(2,1), ...

  static double
  CalculateThirdPrincipalInvariant3x3(double *A); //column first, i.e. A[0] is A(1,1), A[1] = A(2,1), ...

  static bool
  SolveLinearSystem3x3(Vec3D &a, Vec3D &b, Vec3D &c, Vec3D &d, Vec3D &x); //[a b c][x] = [d]

  static bool
  SolveLinearSystem3x3(double a11, double a12, double a13, double a21, double a22, double a23,
                       double a31, double a32, double a33, double b1, double b2, double b3,
                       double &x1, double &x2, double &x3);


  /** A and vectors are both column-first. For example A[2] is A(3,1).
   *  vectors are the right eigenvectors --- the i-th column (i=1,2,3) is
   *  the unit eigenvector corresponding to the i-th eigenvalue.
   *  Note: The eigenvalues are SORTED FROM LOWEST TO HIGHEST.
   *        They are sorted using ther actual values, NOT THE ABSOLUTE VALUES!
   *        Matrix A is assumed to be double[9], storing a symmetric matrix.
   *        The Eigen library only uses the lower-triangular part of A. */
  static bool
  CalculateEigenSymmetricMatrix3x3(double *A, double *values, double *vectors = NULL);


  /** A, vectors_real, vectors_imag are all column-first. For example A[2] is A(3,1).
   *  vectors_{real,imag} are the right eigenvectors --- the i-th column (i=1,2,3) is
   *  the unit eigenvector corresponding to the i-th eigenvalue
   *  NOTE: The eigenvalues are NOT sorted!
   */
  static bool
  CalculateEigen3x3(double *A, double *values_real, double *values_imag,
                  double *vectors_real = NULL, double *vectors_imag = NULL);

  /** A and Ainv are both column-first. For example A[2] is A(3,1).
   *  Returns whether matrix is invertible (up to a default tolerance in Eigen, which
   *  can be changed is needed. */
  static bool
  CalculateMatrixInverseAndDeterminant3x3(double *A, double *Ainv, double *det = NULL);

  //! A, U, and V are all column-first. For example, A[2] is A(3,1). 
  static void
  CalculateSVD3x3(double *A, double *svalues, double *U = NULL, double *V = NULL);

  //! A = R*U: R is orthogonal, U is symmetric, positive semi-definite
  static void
  CalculateRightPolarDecomposition3x3(double *A, double *R, double *U);

  //! A = U*R: R is orthogonal, U is symmetric, positive semi-definite
  static void
  CalculateLeftPolarDecomposition3x3(double *A, double *U, double *R);

  static void
  CalculateATransposeA3x3(double *A, double *ATA);

  static void
  CalculateATransposeB3x3(double *A, double *B, double *ATB);

  static void
  CalculateAATranspose3x3(double *A, double *AAT);

  static void
  CalculateABTranspose3x3(double *A, double *B, double *ABT);

  static double 
  CalculateMatrixTrace3x3(double *A) {return A[0] + A[4] + A[8];}


  //! ----------------------------
  //! 2D vectors and matrices
  //! ----------------------------

  static double
  CalculateDeterminant2x2(double a11, double a12, double a21, double a22);

  static double
  CalculateDeterminant2x2(Vec2D &a, Vec2D &b); //a, b are the two columns

  static double
  CalculateDeterminant2x2(double *A); //column first, i.e. A[0] is A(1,1), A[1] = A(2,1)

  static void 
  CalculateTranspose2x2(double *A, double *AT);
  
  static void
  CalculateMatrixMatrixProduct2x2(double *A, double *B, double *AB);

  static void
  CalculateMatrixVectorProduct2x2(double *A, double *x, double *Ax);

  static void
  CalculateMatrixAPlusB2x2(double *A, double *B, double *Sum);

  static void
  CalculateMatrixC1APlusC2B2x2(double C1, double *A, double C2, double *B, double *Sum);

  static void
  CalculateCTimesMatrixA2x2(double C, double *A, double *CA);

  static double
  CalculateFirstPrincipalInvariant2x2(double *A); //column first, i.e. A[0] is A(1,1), A[1] = A(2,1), ...

  static double
  CalculateSecondPrincipalInvariant2x2(double *A); //column first, i.e. A[0] is A(1,1), A[1] = A(2,1), ...

  static bool
  SolveLinearSystem2x2(Vec2D &a, Vec2D &b, Vec2D &d, Vec2D &x); //[a b][x] = [d]

  static bool 
  SolveLinearSystem2x2(double a11, double a12, double a21, double a22, double b1, double b2,
                       double &x1, double &x2);


  /** A and vectors are both column-first. For example A[1] is A(2,1).
   *  vectors are the right eigenvectors --- the i-th column (i=1,2) is
   *  the unit eigenvector corresponding to the i-th eigenvalue.
   *  Note: The eigenvalues are SORTED FROM LOWEST TO HIGHEST.
   *        They are sorted using ther actual values, NOT THE ABSOLUTE VALUES!
   *        The input A is assumed to be double[4], storing a symmetric matrix.
   *        The Eigen library only uses the lower-triangular part of A.*/
  static bool
  CalculateEigenSymmetricMatrix2x2(double *A, double *values, double *vectors = NULL);


  /** A, vectors_real, vectors_imag are all column-first. For example A[1] is A(2,1).
   *  vectors_{real,imag} are the right eigenvectors --- the i-th column (i=1,2) is
   *  the unit eigenvector corresponding to the i-th eigenvalue
   *  NOTE: The eigenvalues are NOT sorted!
   */
  static bool
  CalculateEigen2x2(double *A, double *values_real, double *values_imag,
                  double *vectors_real = NULL, double *vectors_imag = NULL);

  /** A and Ainv are both column-first. For example A[1] is A(2,1).
   *  Returns whether matrix is invertible (up to a default tolerance in Eigen, which
   *  can be changed is needed. */
  static bool
  CalculateMatrixInverseAndDeterminant2x2(double *A, double *Ainv, double *det = NULL);

  //! A, U, and V are all column-first. For example, A[1] is A(2,1). 
  static void
  CalculateSVD2x2(double *A, double *svalues, double *U = NULL, double *V = NULL);

  //! A = R*U: R is orthogonal, U is symmetric, positive semi-definite
  static void
  CalculateRightPolarDecomposition2x2(double *A, double *R, double *U);

  //! A = U*R: R is orthogonal, U is symmetric, positive semi-definite
  static void
  CalculateLeftPolarDecomposition2x2(double *A, double *U, double *R);

  static void
  CalculateATransposeA2x2(double *A, double *ATA);

  static void
  CalculateATransposeB2x2(double *A, double *B, double *ATB);

  static void
  CalculateAATranspose2x2(double *A, double *AAT);

  static void
  CalculateABTranspose2x2(double *A, double *B, double *ABT);

  static double 
  CalculateMatrixTrace2x2(double *A) {return A[0] + A[3];}


}; //end of class

}
