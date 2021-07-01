#ifndef VISCOSITY_OPERATOR_H_
#define VISCOSITY_OPERATOR_H_

#include <GradientCalculator.h>
#include <VarFcnBase.h>

/***************************************************
 * class ViscosityOperator calculates the viscous
 * fluxes (on the LEFT hand side of the Navier-Stokes
 * equations) at cell interfaces, 
 * and adds it to the total flux function
 **************************************************/

class ViscosityOperator
{
  ViscosityModelData& iod_viscosity;
  
  //! Variable function
  vector<VarFcnBase*>& varFcn; //!< each material has a varFcn

  //! Mesh info
  SpaceVariable3D &coordinates;
  SpaceVariable3D &delta_xyz;

  //! interpolator
  InterpolatorBase *interpolator;

  //! gradient calculator
  GradientCalculatorBase *grad;

  //! internal variables -- storing data at cell interfaces
  SpaceVariable3D G;             //!< viscous fluxes (dim = 5)
  SpaceVariable3D V_i_minus_half;//!< velocity (dim = 3) 
  SpaceVariable3D V_j_minus_half;//!< velocity (dim = 3) 
  SpaceVariable3D V_k_minus_half;//!< velocity (dim = 3) 
  SpaceVariable3D DU;            //!< partial derivatives of x-velocity (dim = 3)  
  SpaceVariable3D DV;            //!< partial derivatives of x-velocity (dim = 3)  
  SpaceVariable3D DW;            //!< partial derivatives of x-velocity (dim = 3)  

  //! internal variables -- storing data at cell center 
  SpaceVariable3D Var;           //!< temporary variable (dim = 5)
 

public:

  ViscosityOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, ViscosityModelData &iod_viscosity_,
                    SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_);
  ~ViscosityOperator();

  //! assign an interpolator (set to a NULL pointer in the constructor)
  void Setup(InterpolatorBase *interpolator_, GradientCalculatorBase *grad_) {
    interpolator = interpolator_;
    grad = grad_;
  }

  //! calculate diffusion fluxes (on the left-hand-side of the N-S equations; add them to R) 
  void AddDiffusionFluxes(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &R);

  //! destroy internal variables
  void Destroy() {
    G.Destroy();  V_i_minus_half.Destroy();  V_j_minus_half.Destroy();  V_k_minus_half.Destroy();
    DU.Destroy();  DV.Destroy();  DW.Destroy();}

};

#endif
