#ifndef _GRADIENT_CALCULATOR_H_
#define _GRADIENT_CALCULATOR_H_

#include <Interpolator.h>

/****************************************************
 * class GradientCalculatorBase is the base class
 * that handles the calculation of the spatial 
 * gradients of any variable at cell interfaces or
 * cell centers
 ***************************************************/

class GradientCalculatorBase
{
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
  virtual void CalculateFirstDerivativeAtNodes(int dir/*0~d/dx,1~d/dy,2~d/dz*/, 
                                             SpaceVariable3D &V, std::vector<int> &input_dof,
                                             SpaceVariable3D &DV, std::vector<int> &output_dof) {
    print_error("*** Error: CalculateFirstDerivativeAtNodes not defined\n");
    exit_mpi();
  }

  //! calculates x-,y-,or z-derivative at i +/- 1/2, j +/- 1/2, or k +/- 1/2
  virtual void CalculateFirstDerivativeAtCellInterfaces(int dir/*0~d/dx,1~d/dy,2~d/dz*/, 
                                                        int inter/*0~i +/- 1/2, 1~j +/- 1/2, 2~k +/- 1/2*/,
                                                        SpaceVariable3D &V, std::vector<int> &input_dof,
                                                        SpaceVariable3D &DV, std::vector<int> &output_dof);
    print_error("*** Error: CalculateFirstDerivativeAtCellInterfaces not defined\n");
    exit_mpi();
  }

  virtual void Destroy();

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
  //! coefficients for the 3-cell stencil
  SpaceVariable3D Cx, Cy, Cz;

  //! temporary variable (dim=5) for internal use
  SpaceVariable3D Var;

  //! interpolator
  Interpolator &interpolator;

public:

  GradientCalculatorCentral(MPI_Comm &comm_, DataManagers3D &dm_all_,
                            SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_,
                            Interpolator &interpolator_);

  ~GradientCalculatorCentral() {}

  //! calculates x-, y-, or j-derivative at nodes 
  void CalculateFirstDerivativeAtNodes(int dir/*0~d/dx,1~d/dy,2~d/dz*/, 
                                       SpaceVariable3D &V, std::vector<int> &input_dof,
                                       SpaceVariable3D &DV, std::vector<int> &output_dof);

  //! calculates x-,y-,or z-derivative at i +/- 1/2, j +/- 1/2, or k +/- 1/2
  void CalculateFirstDerivativeAtCellInterfaces(int dir/*0~d/dx,1~d/dy,2~d/dz*/, 
                                                int inter/*0~i +/- 1/2, 1~j +/- 1/2, 2~k +/- 1/2*/,
                                                SpaceVariable3D &V, std::vector<int> &input_dof,
                                                SpaceVariable3D &DV, std::vector<int> &output_dof);

  void Destroy() {
    Cx.Destroy();  Cy.Destroy();  Cz.Destroy();  Var.Destroy();}

private:

  void CalculateCoefficients();

  //! calculate dV/dx at i +/- 1/2, dV/dy at j +/- 1/2, or dV/dz at k +/- 1/2
  void CentralDifferencingAtCellInterfaces(int dir/*0~d/dx,1~d/dy,2~d/dz*/,
                                           SpaceVariable3D &V, std::vector<int> &input_dof,
                                           SpaceVariable3D &DV, std::vector<int> &output_dof);

};


#endif
