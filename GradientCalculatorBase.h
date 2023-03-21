/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _GRADIENT_CALCULATOR_BASE_H_
#define _GRADIENT_CALCULATOR_BASE_H_

#include <SpaceVariable.h>
#include <Utils.h>

struct Int3;

/****************************************************
 * class GradientCalculatorBase is the base class
 * that handles the calculation of the spatial 
 * gradients of any variable at cell interfaces or
 * cell centers
 * Note: The treatments of the ghost layer outside the
 * physical domain, and of the first layer of cells
 * inside of domain, are NOT unified among the
 * derived classes.
 ***************************************************/

class GradientCalculatorBase
{

protected:
  //! Mesh info
  SpaceVariable3D &coordinates;
  SpaceVariable3D &delta_xyz;

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain

public:

  GradientCalculatorBase(SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_)
                        : coordinates(coordinates_), delta_xyz(delta_xyz_)
  {
    coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
    coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);
  }

  virtual ~GradientCalculatorBase() {}

  //! calculates x-, y-, or j-derivative at nodes 
  virtual void CalculateFirstDerivativeAtNodes([[maybe_unused]] int dir/*0~d/dx,1~d/dy,2~d/dz*/, 
                                               [[maybe_unused]] SpaceVariable3D &V, [[maybe_unused]] std::vector<int> &input_dof,
                                               [[maybe_unused]] SpaceVariable3D &DV, [[maybe_unused]] std::vector<int> &output_dof) {
    print_error("*** Error: CalculateFirstDerivativeAtNodes not defined\n");
    exit_mpi();
  }

  //! calculates x-, y-, or j-derivative at nodes 
  virtual void CalculateFirstDerivativeAtSelectedNodes([[maybe_unused]] int dir/*0~d/dx,1~d/dy,2~d/dz*/,
                                                       [[maybe_unused]] std::vector<Int3> &nodes,
                                                       [[maybe_unused]] SpaceVariable3D &V, [[maybe_unused]] std::vector<int> &input_dof,
                                                       [[maybe_unused]] SpaceVariable3D &DV, [[maybe_unused]] std::vector<int> &output_dof) {
    print_error("*** Error: CalculateFirstDerivativeAtSelectedNodes not defined\n");
    exit_mpi();
  }

  //! calculates x-,y-,or z-derivative at i +/- 1/2, j +/- 1/2, or k +/- 1/2
  virtual void CalculateFirstDerivativeAtCellInterfaces([[maybe_unused]] int dir/*0~d/dx,1~d/dy,2~d/dz*/, 
                                                        [[maybe_unused]] int inter/*0~i +/- 1/2, 1~j +/- 1/2, 2~k +/- 1/2*/,
                                                        [[maybe_unused]] SpaceVariable3D &V, [[maybe_unused]] std::vector<int> &input_dof,
                                                        [[maybe_unused]] SpaceVariable3D &DV, [[maybe_unused]] std::vector<int> &output_dof) {
    print_error("*** Error: CalculateFirstDerivativeAtCellInterfaces not defined\n");
    exit_mpi();
  }

  virtual void Destroy() {}

};

#endif
