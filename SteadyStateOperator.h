/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _STEADY_STATE_OPERATOR_H_
#define _STEADY_STATE_OPERATOR_H_
#include<SpaceVariable.h>

struct TsData;
class GlobalMeshInfo;

/*********************************************************
 * class SteadyStateOperator handles some operations that
 * are only used in steady-state simulations. In general,
 * it is a "tools class", which does not own a lot of
 * data. This class itself is owned by TimeIntegrator.
 *********************************************************/ 

class SteadyStateOperator
{
  MPI_Comm& comm;

  GlobalMeshInfo &global_mesh;

  /** Currently, only monitors the N-S residual for convergence. Can be easily extended
   *  to consider state variable perturbation or some QoI. */   
  std::vector<double> Rref; //!< ref. residual (mass, momentum and energy fluxes) for non-dimensionalization
  double Rtol; //!< residual tol. (user input)
  bool ref_calculated; //!< whether references has been computed.

  double R1_init, R2_init, Rinf_init; //! initial residual in 1-, 2-, and inf-norm.
  double R1, R2, Rinf;
  bool converged;
  

public:

  SteadyStateOperator(MPI_Comm &comm_, TsData &iod_ts_, GlobalMeshInfo &global_mesh_);
  ~SteadyStateOperator();

  void Destroy();

  bool Converged() {return converged;} 

  void MonitorConvergence(SpaceVariable3D &R, SpaceVariable3D &ID); //!< calculates R1, R2, and Rinf

  double GetResidual1Norm() {return R1;}
  double GetRelativeResidual1Norm() {return R1/R1_init;}
 
  double GetResidual2Norm() {return R2;}
  double GetRelativeResidual2Norm() {return R2/R2_init;}
 
  double GetResidualInfNorm() {return Rinf;}
  double GetRelativeResidualInfNorm() {return Rinf/Rinf_init;}
 
protected:

  void CalculateReferenceResidual(SpaceVariable3D &R, SpaceVariable3D &ID);

};

#endif
