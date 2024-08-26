/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _INCOMPRESSIBLE_OPERATOR_H_
#define _INCOMPRESSIBLE_OPERATOR_H_

#include<Interpolator.h>
#include<GhostFluidOperator.h>
#include<SpaceOperator.h>

struct RowEntries;

/***************************************************
 * class IncompressibleOperator calculates the fluxes
 * of the incompressible Navier-Stokes equations.
 **************************************************/

class IncompressibleOperator
{
  
  MPI_Comm& comm;
  IoData& iod;

  //! Variable function
  std::vector<VarFcnBase*>& vf; //!< each material has a varFcn

  //! Space operator
  SpaceOperator &spo;

  //! dynamic viscosity coefficient, one for each material (INACTIVE has mu = 0)
  vector<double> Mu; //!< 0 if inviscid

  //! Mesh info
  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain
  int NX, NY, NZ;

  //! interpolator
  InterpolatorBase &interpolator;

  //! Ghost fluid method
  GhostFluidOperator *gfo;

  //! Internal space variables
  SpaceVariable3D V3;

  //! "deltas" used by the SIMPLE family for velocities on staggered grids (for SUBDOMAIN only)
  std::vector<std::vector<double> > Dx, Dy, Dz, dx_l, dx_r, dy_b, dy_t, dz_k, dz_f; //!< (5.63) Patankar

public:

  IncompressibleOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_,
                         std::vector<VarFcnBase*> &varFcn_, SpaceOperator &spo_,
                         InterpolatorBase &interp_);
  ~IncompressibleOperator();

  void Destroy();

  void FinalizeInitialCondition(SpaceVariable3D &V, SpaceVariable3D &ID);

  void ComputeTimeStepSize(SpaceVariable3D &V, SpaceVariable3D &ID, double &dt, double &cfl,
                           SpaceVariable3D *LocalDt = NULL);

  void ApplyBoundaryConditions(SpaceVariable3D &V); //!< Also modify non-ghost entries (due to MAC grid)

  void BuildVelocityEquationSIMPLE(int dir, //!< dir: 0,1,2 for x,y,z
                                   Vec5D*** v0, Vec5D*** v, double*** id,
                                   double*** vturb0, //turbulent working term
                                   double*** homo, //wheter each node is in a homogeneous region
                                   std::vector<RowEntries> &vlin_rows, SpaceVariable3D &B, SpaceVariable3D &Ddiag,
                                   bool SIMPLEC, //!< for SIMPLEC, generates a different Ddiag; otherwise the same
                                   double Efactor, double dt,
                                   SpaceVariable3D *LocalDt = NULL); 

  void BuildPressureEquationSIMPLE(Vec5D*** v, double*** homo, SpaceVariable3D &VXstar,
                                   SpaceVariable3D &VYstar, SpaceVariable3D &VZstar,
                                   SpaceVariable3D &DX, SpaceVariable3D &DY, SpaceVariable3D &DZ,
                                   std::vector<RowEntries> &plin_rows, SpaceVariable3D &B,
                                   Int3 *ijk_zero_p = NULL); //!< allows p to be fixed at one node

  //! For specified momentum equation (dir=0,1,2), find coefficients for pressure and velocity equations
  void CalculateCoefficientsSIMPLER(int dir, Vec5D*** v0, Vec5D*** v, double*** id, double*** homo,
                                    std::vector<RowEntries> &vlin_rows, SpaceVariable3D &Bv,
                                    SpaceVariable3D &Vhat, SpaceVariable3D &Ddiag, double Efactor,
                                    double dt, SpaceVariable3D *LocalDt = NULL);

  void UpdateVelocityEquationRHS_SIMPLER(int dir, SpaceVariable3D &P, SpaceVariable3D &B);

  void BuildPressureEquationRHS_SIMPLER(Vec5D*** v, double*** homo,
                                        SpaceVariable3D &VXstar, SpaceVariable3D &VYstar,
                                        SpaceVariable3D &VZstar, SpaceVariable3D &B, Int3 *ijk_zero_p = NULL);

  //! For specified momentum equation (dir=0,1,2), calculate "velocity tilde prime"
  void CalculateVelocityTildePISO(int dir, Vec5D*** v0, Vec5D*** v, double*** id, double*** homo,
                                  SpaceVariable3D &Vprime, SpaceVariable3D &Vtilde, double Efactor,
                                  double dt, SpaceVariable3D *LocalDt = NULL);
 
  void CalculateMomentumChanges(Vec5D*** v0, SpaceVariable3D &V, double*** id, SpaceVariable3D &R3);




  //! For turbulent flows (RANS for now)
  void InitializeTurbulenceVariables(SpaceVariable3D &Vturb);

  void ApplyBoundaryConditionsTurbulenceVariables(SpaceVariable3D &Vturb);

  void BuildSATurbulenceEquationSIMPLE(Vec5D*** v, double*** id,
                                       double*** vturb0, double*** vturb,//added SA eddy viscosity working term
                                       std::vector<RowEntries> &vlin_rows, SpaceVariable3D &B, 
                                       double Efactor, double cw1_reduction,
                                       double dt, SpaceVariable3D *LocalDt = NULL); 

  void ComputeKinematicEddyViscosity(SpaceVariable3D &Vturb, SpaceVariable3D &V, SpaceVariable3D &ID,
                                     SpaceVariable3D &NuT); //compute kin.eddy.vis ==> NuT


private:

  void CheckInputs(IoData &iod); //!< Check input file. Quit if error is found.

  void ComputeLocalTimeStepSizes(SpaceVariable3D &V, SpaceVariable3D &ID, double &dt, double &cfl,
                                 SpaceVariable3D &LocalDt);

  void ApplyBoundaryConditionsGeometricEntities(Vec5D*** v);

  //! Function "A" -- Eq.(5.64) in Patankar's book
  inline double PowerLaw(double pc) {double pp=1.0-0.1*fabs(pc); return pp>0.0 ? pow(pp,5) : 0.0;}
  //inline double PowerLaw([[maybe_unused]] double pc) {return 1.0;} //degenerates to upwinding (useful for debugging)


  //! For turbulent flows (RANS for now) TODO: These functions should be moved to dedicated classes later
  double GetKinematicEddyViscosity(double rho, double mu, double nu_tilde); //!< Spalart-Allmaras
  double GetDynamicEddyViscosity(double rho, double mu, double nu_tilde); //!< Spalart-Allmaras
  double GetDistanceToWall(Vec3D x);
  double ApplyVelocityBoundaryConditionLocal(int dir, GhostPoint::Side side, double vim,
                                             bool flat_plate_wall = false); //!< for computing vorticity

};

#endif
