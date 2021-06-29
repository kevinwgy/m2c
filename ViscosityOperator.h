#ifndef VISCOSITY_OPERATOR_H_
#define VISCOSITY_OPERATOR_H_

#include <IoData.h>
#include <SpaceVariable3D.h>
#include <Interpolator.h>

/***************************************************
 * class ViscosityOperator calculates the viscous
 * fluxes (on the LEFT hand side of the Navier-Stokes
 * equations) at cell interfaces, 
 * and adds it to the total flux function
 **************************************************/

class ViscosityOperator
{
  ViscosityModelData& iod_viscosity;
  
  //! Mesh info
  SpaceVariable3D &coordinates;
  SpaceVariable3D &delta_xyz;

  //! interpolator
  InterpolatorBase *interpolator;

  //! internal variables -- storing data at cell interfaces
  SpaceVariable3D G;             //!< viscous fluxes (dim = 5)
  SpaceVariable3D V_i_minus_half;//!< velocity (dim = 3) 
  SpaceVariable3D V_j_minus_half;//!< velocity (dim = 3) 
  SpaceVariable3D V_k_minus_half;//!< velocity (dim = 3) 
  SpaceVariable3D DU;            //!< partial derivatives of x-velocity (dim = 3)  
  SpaceVariable3D DV;            //!< partial derivatives of x-velocity (dim = 3)  
  SpaceVariable3D DW;            //!< partial derivatives of x-velocity (dim = 3)  

  //! internal variables -- storing data at cell center 
  SpaceVariable3D Var;           //!< temporary variable (dim = 1)
 

public:

  ViscosityOperator(DataManagers3D &dm_all_, ViscosityModelData &iod_viscosity_,
                    SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_);
  ~ViscosityOperator();

  //! assign an interpolator (set to a NULL pointer in the constructor)
  void AssignInterpolator(InterpolatorBase *interpolator_) {
    interpolator = interpolator_;
  }

  //! calculate diffusion fluxes (on the left-hand-side of the N-S equations; add them to R) 
  void AddDiffusionFluxes(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &R);

};

#endif
