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
  virtual void CalculateFirstDerivativeAtNodes(int dir/*0~d/dx,1~d/dy,2~d/dz*/, 
                                             SpaceVariable3D &V, std::vector<int> &input_dof,
                                             SpaceVariable3D &DV, std::vector<int> &output_dof) {
    print_error("*** Error: CalculateFirstDerivativeAtNodes not defined\n");
    exit_mpi();
  }

  //! calculates x-, y-, or j-derivative at nodes 
  virtual void CalculateFirstDerivativeAtSelectedNodes(int dir/*0~d/dx,1~d/dy,2~d/dz*/, std::vector<Int3> &nodes,
                                                       SpaceVariable3D &V, std::vector<int> &input_dof,
                                                       SpaceVariable3D &DV, std::vector<int> &output_dof) {
    print_error("*** Error: CalculateFirstDerivativeAtSelectedNodes not defined\n");
    exit_mpi();
  }

  //! calculates x-,y-,or z-derivative at i +/- 1/2, j +/- 1/2, or k +/- 1/2
  virtual void CalculateFirstDerivativeAtCellInterfaces(int dir/*0~d/dx,1~d/dy,2~d/dz*/, 
                                                        int inter/*0~i +/- 1/2, 1~j +/- 1/2, 2~k +/- 1/2*/,
                                                        SpaceVariable3D &V, std::vector<int> &input_dof,
                                                        SpaceVariable3D &DV, std::vector<int> &output_dof) {
    print_error("*** Error: CalculateFirstDerivativeAtCellInterfaces not defined\n");
    exit_mpi();
  }

  virtual void Destroy() {}

};

#endif
