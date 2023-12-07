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
                           SpaceVariable3D *localDt = NULL);

  void ApplyBoundaryConditions(SpaceVariable3D &V); //!< Also modify non-ghost entries (due to MAC grid)

  void BuildVelocityEquationSIMPLE(int dir, Vec5D*** v, double*** id, std::vector<RowEntries> &vlin_rows,
                                   SpaceVariable3D &B, SemiImplicitTsData &iod_semi); //!< dir: 0,1,2

private:

  void CheckInputs(IoData &iod); //!< Check input file. Quit if error is found.

  void ComputeLocalTimeStepSizes(SpaceVariable3D &V, SpaceVariable3D &ID, double &dt, double &cfl,
                                 SpaceVariable3D &LocalDt);

  void ApplyBoundaryConditionsGeometricEntities(Vec5D*** v);

  //! Function "A" -- Eq.(5.64) in Patankar's book
  inline double PowerLaw(double pc) {double pp=1.0-0.1*fabs(pc); return pp>0.0 ? pow(pp,5) : 0.0;}
};

#endif
