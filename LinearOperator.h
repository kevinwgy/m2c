/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _LINEAR_OPERATOR_H_
#define _LINEAR_OPERATOR_H_
#include<petscmat.h>
#include<SpaceVariable.h>
#include<IoData.h>

/********************************************************************
 * class RowEntries stores information about one or multiple entries 
 * in a row of a matrix operating on SpaceVariable3D.
 *
 * typedef struct {
 *   PetscInt k, j, i, c;
 * } MatStencil;
 *
 * The i,j, and k represent the logical coordinates over the entire grid 
 * (for 2 and 1 dimensional problems the k and j entries are ignored). 
 * The c represents the the degrees of freedom at each grid point (the 
 * dof argument to DMDASetDOF()). If dof is 1 then this entry is ignored.
 *********************************************************************/

struct RowEntries {
  MatStencil row; //!< row number, defined by i,j,k,c
  std::vector<MatStencil> cols; //!< col numbers
  std::vector<double> vals; //entries (one per col)

  RowEntries() {}
  RowEntries(int nEntries) {
    cols.reserve(nEntries);
    vals.reserve(nEntries);
  }
};


/*********************************************************************
 * class LinearOperator handles matrix definition and operations.
 * The actual work is done by PETSc.
 *********************************************************************
*/

class LinearOperator {

protected:

  MPI_Comm &comm;

  DM dm; /**< This is a new DM object "cloned" from the one given in the constructor.\n
              According to PETsc documentation, one DM should be constructed for solving each system.*/
  Mat A; //!< coefficient matrix. type: MATAIJ.

  int i0, j0, k0, ii0, jj0, kk0; //!< same as in SpaceVariable3D
  int imax, jmax, kmax, iimax, jjmax, kkmax;

  int dof; //!< same as in SpaceVariable3D

public:

  LinearOperator(MPI_Comm &comm_, DM &dm_);
  virtual ~LinearOperator(); 
  virtual void Destroy();

  virtual void SetLinearOperator(std::vector<RowEntries>& row_entries);

  void ApplyLinearOperator(SpaceVariable3D &x, SpaceVariable3D &y); //!< calculates Ax -> y (x =/= y!)

  double CalculateMatrixOneNorm();
  double CalculateMatrixInfNorm();
  double CalculateMatrixFrobeniusNorm();

protected:

};

#endif
