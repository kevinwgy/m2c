/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _TIME_INTEGRATOR_SEMIIMP_H_
#define _TIME_INTEGRATOR_SEMIIMP_H_

#include <TimeIntegrator.h>
#include <IncompressibleOperator.h>
#include <LinearSystemSolver.h>

/********************************************************************
 * SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) 
 * Refs: Patankar's book, Van Doormaal and Raithby, 1984; also see KW's notes
 *******************************************************************/
class TimeIntegratorSIMPLE : public TimeIntegratorBase
{

protected:

  IncompressibleOperator &inco;

  SpaceVariable3D V0; //!< solution at the previous time step
  SpaceVariable3D *Vturb0_ptr; //!< solution of turbulence working variables at previous time step

  SpaceVariable3D VXstar, VYstar, VZstar, Pprime;
  SpaceVariable3D B; //!< generally used as the right-hand-side of linear systems
  SpaceVariable3D Homo; //!< whether each node is within a homogeneous neighborhood (const. rho, mu)

  SpaceVariable3D DX, DY, DZ; //!< coefficients from momentum equations, later used in pressure-correction
  vector<RowEntries> vlin_rows;
  vector<RowEntries> plin_rows;
  vector<RowEntries> vturb_lin_rows; //!< for turbulence closure equation(s), if involved

  LinearSystemSolver vlin_solver; //!< solves the velocity linear systems
  LinearSystemSolver plin_solver; //!< solves the pressure linear system
  LinearSystemSolver vturb_lin_solver; //!< solves the turbulence closure equation linear system(s)

  SpaceVariable3D *R3_ptr; //!< momentum changes (dof = 3), for steady state only

  SpaceVariable3D Tmp1; //!< for temporary use (dof = 1)
  SpaceVariable3D Tmp5; //!< for temporary use (dof = 5)

  //! Relaxation coefficients
  double Efactor;
  double alphaP;

  //! The pressure system is rank-deficient. Can be fixed by fixing pressure to 0 at a point
  bool fix_pressure_at_one_corner;

  //! A corner where pressure is fixed to 0
  Int3 ijk_zero_p; //!< set to [NZ-1][NY-1][NX-1]

public:

  TimeIntegratorSIMPLE(MPI_Comm &comm_, IoData& iod_, DataManagers3D& dms_, SpaceOperator& spo_,
                       IncompressibleOperator &inco_, vector<LevelSetOperator*>& lso_,
                       MultiPhaseOperator &mpo_, LaserAbsorptionSolver* laser_,
                       EmbeddedBoundaryOperator* embed_, HyperelasticityOperator* heo_,
                       PrescribedMotionOperator* pmo_);
  ~TimeIntegratorSIMPLE();

  //! Note: virtual functions in the base class are automatically virtual in derived classes.
  void AdvanceOneTimeStep(SpaceVariable3D &V, SpaceVariable3D &ID,
                          vector<SpaceVariable3D*>& Phi, vector<SpaceVariable3D*> &NPhi,
                          vector<SpaceVariable3D*> &KappaPhi,
                          SpaceVariable3D *L, SpaceVariable3D *Xi, SpaceVariable3D *Vturb,
                          SpaceVariable3D *LocalDt,
                          double time, double dt, int time_step, int subcycle, double dts);

  void Destroy();

protected:

  Int3 FindCornerFixedPressure();

  //! v --> VX, VY, VZ, P
  void ExtractVariableComponents(Vec5D*** v, SpaceVariable3D *VX_ptr, SpaceVariable3D *VY_ptr,
                                 SpaceVariable3D *VZ_ptr, SpaceVariable3D *P_ptr);

  //! Update v and p (for next iteration)
  virtual double UpdateStates(Vec5D*** v, SpaceVariable3D &Pprime, SpaceVariable3D &DX,
                              SpaceVariable3D &DY, SpaceVariable3D &DZ, SpaceVariable3D &VX,
                              SpaceVariable3D &VY, SpaceVariable3D &VZ, double prelax);

};  



/********************************************************************
 * SIMPLER (R: Revised)
 * Refs: Patankar's book, Van Doormaal and Raithby, 1984; also see KW's notes
 *******************************************************************/
class TimeIntegratorSIMPLER : public TimeIntegratorSIMPLE
{

protected:

  vector<RowEntries> ulin_rows;
  vector<RowEntries> wlin_rows;

  LinearSystemSolver ulin_solver;
  LinearSystemSolver wlin_solver;

  SpaceVariable3D Bu, Bv, Bw, P;

public:

  TimeIntegratorSIMPLER(MPI_Comm &comm_, IoData& iod_, DataManagers3D& dms_, SpaceOperator& spo_,
                        IncompressibleOperator &inco_, vector<LevelSetOperator*>& lso_,
                        MultiPhaseOperator &mpo_, LaserAbsorptionSolver* laser_,
                        EmbeddedBoundaryOperator* embed_, HyperelasticityOperator* heo_,
                        PrescribedMotionOperator* pmo_);
  ~TimeIntegratorSIMPLER();

  void Destroy();

  //! Note: virtual functions in the base class are automatically virtual in derived classes.
  void AdvanceOneTimeStep(SpaceVariable3D &V, SpaceVariable3D &ID,
                          vector<SpaceVariable3D*>& Phi, vector<SpaceVariable3D*> &NPhi,
                          vector<SpaceVariable3D*> &KappaPhi,
                          SpaceVariable3D *L, SpaceVariable3D *Xi, SpaceVariable3D *Vturb,
                          SpaceVariable3D *LocalDt,
                          double time, double dt, int time_step, int subcycle, double dts);

protected:

  //! Update v and p (for next iteration)
  double UpdateStates(Vec5D*** v, SpaceVariable3D &P, SpaceVariable3D &Pprime, SpaceVariable3D &DX,
                      SpaceVariable3D &DY, SpaceVariable3D &DZ, SpaceVariable3D &VX,
                      SpaceVariable3D &VY, SpaceVariable3D &VZ);
};  


/********************************************************************
 * SIMPLEC (C: Consistent)
 * Refs: Van Doormaal and Raithby, 1984; also see KW's notes
 *******************************************************************/
class TimeIntegratorSIMPLEC : public TimeIntegratorSIMPLE
{

public:

  TimeIntegratorSIMPLEC(MPI_Comm &comm_, IoData& iod_, DataManagers3D& dms_, SpaceOperator& spo_,
                        IncompressibleOperator &inco_, vector<LevelSetOperator*>& lso_,
                        MultiPhaseOperator &mpo_, LaserAbsorptionSolver* laser_,
                        EmbeddedBoundaryOperator* embed_, HyperelasticityOperator* heo_,
                        PrescribedMotionOperator* pmo_);
  ~TimeIntegratorSIMPLEC();

};  


/********************************************************************
 * PISO (Pressure-Implicit with Split Operator)
 * Refs: Issa, JCP, 1986; also see KW's notes
 *******************************************************************/
class TimeIntegratorPISO : public TimeIntegratorSIMPLE
{

protected:

  SpaceVariable3D VXprime, VYprime, VZprime, VXtildep, VYtildep, VZtildep, Pstar;

public:

  TimeIntegratorPISO(MPI_Comm &comm_, IoData& iod_, DataManagers3D& dms_, SpaceOperator& spo_,
                     IncompressibleOperator &inco_, vector<LevelSetOperator*>& lso_,
                     MultiPhaseOperator &mpo_, LaserAbsorptionSolver* laser_,
                     EmbeddedBoundaryOperator* embed_, HyperelasticityOperator* heo_,
                     PrescribedMotionOperator* pmo_);
  ~TimeIntegratorPISO();

  //! Note: virtual functions in the base class are automatically virtual in derived classes.
  void AdvanceOneTimeStep(SpaceVariable3D &V, SpaceVariable3D &ID,
                          vector<SpaceVariable3D*>& Phi, vector<SpaceVariable3D*> &NPhi,
                          vector<SpaceVariable3D*> &KappaPhi,
                          SpaceVariable3D *L, SpaceVariable3D *Xi, SpaceVariable3D *Vturb,
                          SpaceVariable3D *LocalDt,
                          double time, double dt, int time_step, int subcycle, double dts);

  void Destroy();

protected:

  double UpdateStatesCorrector(SpaceVariable3D &Pstar, SpaceVariable3D &Pprime, SpaceVariable3D &DX,
                               SpaceVariable3D &DY, SpaceVariable3D &DZ, SpaceVariable3D &VXstar,
                               SpaceVariable3D &VYstar, SpaceVariable3D &VZstar,
                               SpaceVariable3D &VXprime, SpaceVariable3D &VYprime,
                               SpaceVariable3D &VZprime, SpaceVariable3D *VXtildep = NULL,
                               SpaceVariable3D *VYtildep = NULL, SpaceVariable3D *VZtildep = NULL);

  void UpdateStatesFinal(Vec5D*** v, SpaceVariable3D &P, SpaceVariable3D &VX,
                         SpaceVariable3D &VY, SpaceVariable3D &VZ);

};  












#endif
