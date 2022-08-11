#pragma once
#include <Vector3D.h>

namespace MathTools {

struct LinearAlgebra {

static bool 
SolveLinearSystem2x2(double a11, double a12, double a21, double a22, double b1, double b2,
                     double &x1, double &x2);

static double
CalculateDeterminant3x3(double a11, double a12, double a13, double a21, double a22, double a23,
                        double a31, double a32, double a33);

static double
CalculateDeterminant3x3(Vec3D &a, Vec3D &b, Vec3D &c); //a, b, c are the three columns

static double
CalculateDeterminant3x3(double *A); //column first, i.e. A[0] is A(1,1), A[1] = A(2,1)

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
 *        They are sorted using ther actual values, NOT THE ABSOLUTE VALUES! */
static bool
ComputeEigenSymmetricMatrix3x3(double *A, double *values, double *vectors);


/** A, vectors_real, vectors_imag are all column-first. For example A[2] is A(3,1).
 *  vectors_{real,imag} are the right eigenvectors --- the i-th column (i=1,2,3) is
 *  the unit eigenvector corresponding to the i-th eigenvalue
 *  NOTE: The eigenvalues are NOT sorted!
 */
static bool
ComputeEigen3x3(double *A, double *values_real, double *values_imag,
                double *vectors_real, double *vectors_imag);


}; //end of class

}
