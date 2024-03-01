/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _GRADIENT_CALCULATOR_FD3_H_
#define _GRADIENT_CALCULATOR_FD3_H_

#include <GradientCalculatorBase.h>
#include <Vector3D.h>

/****************************************************
 * class GradientCalculatorFD3 calculates spatial 
 * gradients by a 3rd order finite difference scheme,
 * with a stencil specified by user (biased)
 ***************************************************/

class GradientCalculatorFD3 : public GradientCalculatorBase
{
  //! stencil type: -1: {-2, -1, 0, 1},  1: {-1, 0, 1, 2}
  int shift;

  //! coefficients for the 3-cell stencil
  SpaceVariable3D Cx, Cy, Cz;

  //! temporary variable (dim=5) for internal use
  SpaceVariable3D Var;

  //! local copy of coordinates with 2 ghost layers
  SpaceVariable3D coordinatesG2;

public:

  GradientCalculatorFD3(MPI_Comm &comm_, DataManagers3D &dm_all_,
                        SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_, int shift_);

  ~GradientCalculatorFD3() {}

  //! calculates x-, y-, or j-derivative at nodes (domain inteior)
  void CalculateFirstDerivativeAtNodes(int dir/*0~d/dx,1~d/dy,2~d/dz*/, 
                                       SpaceVariable3D &V, std::vector<int> &input_dof,
                                       SpaceVariable3D &DV, std::vector<int> &output_dof);

  //! calculates x-, y-, or j-derivative at selected nodes (domain interior)
  void CalculateFirstDerivativeAtSelectedNodes(int dir/*0~d/dx,1~d/dy,2~d/dz*/, std::vector<Int3> &nodes,
                                       SpaceVariable3D &V, std::vector<int> &input_dof,
                                       SpaceVariable3D &DV, std::vector<int> &output_dof);

  void Destroy() {
    Cx.Destroy();  Cy.Destroy();  Cz.Destroy();  Var.Destroy();  coordinatesG2.Destroy();}

private:

  void CalculateCoefficients();

  void CalculateCoefficientsLeftBiased();
  void CalculateCoefficientsRightBiased();

  void CalculateCoefficientsByTaylerSeries(double dx0, double dx1, double dx2, double dx3, double *coeffs);
};


#endif
