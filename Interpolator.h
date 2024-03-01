/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _INTERPOLATOR_H__
#define _INTERPOLATOR_H__

#include <IoData.h>
#include <SpaceVariable.h>

/****************************************************
 * class InterpolatorBase is the base class for 
 * interpolation
 ***************************************************/
class InterpolatorBase 
{

protected:
  //! Mesh info
  SpaceVariable3D &coordinates;
  SpaceVariable3D &delta_xyz;

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain

public:

  InterpolatorBase(SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_)
                 : coordinates(coordinates_), delta_xyz(delta_xyz_)
  {
    coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
    coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);
  }

  virtual ~InterpolatorBase() { }

  //! Interpolate at cell interfaces
  virtual void InterpolateAtCellInterfaces([[maybe_unused]] int dir/*0~x,1~y,2~z*/,
                                           [[maybe_unused]] SpaceVariable3D &Vin, [[maybe_unused]] std::vector<int> &input_dof, 
                                           [[maybe_unused]] SpaceVariable3D &Vout, [[maybe_unused]] std::vector<int> &output_dof) {
    print_error("*** Error: InterpolateAtCellInterfaces not defined\n");
    exit_mpi();
  }
  
  virtual void Destroy() {}

};

/****************************************************
 * class InterpolatorLinear is the class for 
 * calculating linear interpolation
 ***************************************************/
class InterpolatorLinear : public InterpolatorBase
{
  //! interpolation coefficients (left, right, top, bottom, back, front)
  //! For example, Cl[k][j][i] is the coefficient to be multiplied to the variable at
  //! (i,j,k) for interpolating the value at (i+1/2,j,k) --- i.e. when (i,j,k) is the
  //! LEFT cell. In other words, the interpolation formula should be like
  //! v(i+1/2,j,k) = Cl[k][j][i]*v[k][j][i] + Cr[k][j][i+1]*v[k][j][i+1];
  SpaceVariable3D Cl, Cr, Ct, Cb, Ck, Cf; 
                                      

public:

  InterpolatorLinear(MPI_Comm &comm_, DataManagers3D &dm_all_,
                     SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_);

  ~InterpolatorLinear() {}

  void InterpolateAtCellInterfaces(int dir/*0~x,1~y,2~z*/, 
                                   SpaceVariable3D &Vin, std::vector<int> &input_dof, 
                                   SpaceVariable3D &Vout, std::vector<int> &output_dof);

  void Destroy() {
    Cl.Destroy(); Cr.Destroy(); Ct.Destroy(); Cb.Destroy(); Ck.Destroy(); Cf.Destroy();}

private:

  void CalculateInterpolationCoefficients();

};

#endif
