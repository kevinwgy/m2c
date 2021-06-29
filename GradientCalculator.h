#ifndef _GRADIENT_CALCULATOR_H_
#define _GRADIENT_CALCULATOR_H_

#include <IoData.h>
#include <SpaceVariable.h>

/****************************************************
 * class GradientCalculatorBase is the base class
 * that handles the calculation of the spatial 
 * gradients of any variable at cell interfaces
 ***************************************************/

class GradientCalculatorBase
{
  //! Mesh info
  SpaceVariable3D &coordinates;
  SpaceVariable3D &delta_xyz;

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain

public:

  GradientCalculator(SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_);

  void CalculateFirstDerivativesAtCellInterfaces(int dir/*0~x,1~y,2~z*/, SpaceVariable3D &V, int dof,
                                                 SpaceVariable3D &DV) {
    print_error("*** Error: CalculateFirstDerivativesAtCellInterface not defined\n");
    exit_mpi();
  }

};


/****************************************************
 * class GradientCalculatorCentral calculates spatial 
 * gradients by central differencing (truly ``central''
 * only for uniform grids). The stencil width
 * is either 2 cells (e.g., when calculating dudx at
 * i +/- 1/2) or 3 cells (e.g., when calculating dudx
 * at j +/- 1/2)
 ***************************************************/

class GradientCalculatorCentral : public GradientCalculatorBase
{
  //! coefficients for a 3-cell stencil
  SpaceVariable3D Cx0, Cx1, Cx2;
  SpaceVariable3D Cy0, Cy1, Cy2;
  SpaceVariable3D Cz0, Cz1, Cz2;

public:

  GradientCalculatorCentral(MPI_Comm &comm_, DataManagers3D &dm_all_,
                            SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_);




};


#endif
