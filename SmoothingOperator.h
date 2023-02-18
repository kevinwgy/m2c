/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _SMOOTHING_Operator_H_
#define _SMOOTHING_Operator_H_
#include <IoData.h>
#include <SpaceVariable.h>

/********************************************************************************
 * Class SmoothingOperator applies smoothing filters (of different types) to a spatial
 * variable of any dimension. By default it does not go across material interfaces.
 * i.e. the stencil for each node will contain only nodes within the same material
 *******************************************************************************/
class SmoothingOperator
{
  MPI_Comm &comm;

  SmoothingData& iod_smooth;

  //! Mesh info
  SpaceVariable3D &coordinates;
  SpaceVariable3D &delta_xyz;
  SpaceVariable3D &volume;

  //! Interval variable to store data temporarily
  SpaceVariable3D V0;

public:

  SmoothingOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, SmoothingData &iod_smooth_,
                    SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_,
                    SpaceVariable3D &volume_);
  ~SmoothingOperator();

  //apply a smoothing filter within the physical domain
  void ApplySmoothingFilter(SpaceVariable3D &V, SpaceVariable3D *ID = NULL);

  void Destroy();

private:

  void ApplyBoxFilter(SpaceVariable3D &V, SpaceVariable3D *ID);

  void ApplyGausianFilter(SpaceVariable3D &V, SpaceVariable3D *ID); 

  void EnforceLocalConservation(SpaceVariable3D &U0, SpaceVariable3D &U, SpaceVariable3D *ID);
};

#endif
