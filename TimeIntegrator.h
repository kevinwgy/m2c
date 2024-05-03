/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _TIME_INTEGRATOR_H_
#define _TIME_INTEGRATOR_H_

#include <SpaceOperator.h>
#include <LevelSetOperator.h>
#include <MultiPhaseOperator.h>
#include <LaserAbsorptionSolver.h>
#include <RiemannSolutions.h>
#include <EmbeddedBoundaryOperator.h>
#include <HyperelasticityOperator.h>
#include <PrescribedMotionOperator.h>
#include <SteadyStateOperator.h>
using std::vector;

/********************************************************************
 * Numerical time-integrator: The base class
 *******************************************************************/
class TimeIntegratorBase
{

protected:

  enum Type {NONE = 0, FORWARD_EULER = 1, RUNGE_KUTTA_2 = 2, RUNGE_KUTTA_3 = 3,
             SIMPLE = 4, SIMPLER = 5, SIMPLEC = 6, PISO = 7} type;


  MPI_Comm&       comm;
  IoData&         iod;
  SpaceOperator&  spo;

  vector<LevelSetOperator*>& lso;
  MultiPhaseOperator& mpo;  

  //! Laser absorption solver (NULL if laser is not activated)
  LaserAbsorptionSolver* laser;

  //! Embedded boundary method (NULL if not activated)
  EmbeddedBoundaryOperator* embed;

  //! Hyperelaticity operator (NULL if not activated)
  HyperelasticityOperator* heo;

  //! Prescribed motion operator (NULL if not activated)
  PrescribedMotionOperator* pmo;
 
  //! Internal variable to temporarily store old ID
  SpaceVariable3D IDn;

  //! Internal variable to temporarily store Phi (e.g., for material ID updates)
  vector<SpaceVariable3D*> Phi_tmp;

  //! Internal variable to store the mat. ID tracked by each level set
  vector<int> ls_mat_id;

  //! Solutions of exact Riemann problems
  RiemannSolutions riemann_solutions;

  //! Variables for steady-state analysis
  SteadyStateOperator *sso;
  bool local_time_stepping;

public:
  TimeIntegratorBase(MPI_Comm &comm_, IoData& iod_, DataManagers3D& dms_, SpaceOperator& spo_, 
                     vector<LevelSetOperator*>& lso_, MultiPhaseOperator& mpo_,
                     LaserAbsorptionSolver* laser_, EmbeddedBoundaryOperator* embed_,
                     HyperelasticityOperator* heo_, PrescribedMotionOperator* pmo_);

  virtual ~TimeIntegratorBase();

  //! Integrate the ODE system for one time-step. Implemented in derived classes
  virtual void AdvanceOneTimeStep([[maybe_unused]] SpaceVariable3D &V, [[maybe_unused]] SpaceVariable3D &ID, 
                                  [[maybe_unused]] vector<SpaceVariable3D*> &Phi,
                                  [[maybe_unused]] vector<SpaceVariable3D*> &NPhi,
                                  [[maybe_unused]] vector<SpaceVariable3D*> &KappaPhi,
                                  [[maybe_unused]] SpaceVariable3D *L, [[maybe_unused]] SpaceVariable3D *Xi,
                                  [[maybe_unused]] SpaceVariable3D *Vturb,
                                  [[maybe_unused]] SpaceVariable3D *LocalDt,
                                  [[maybe_unused]] double time, [[maybe_unused]] double dt,
                                  [[maybe_unused]] int time_step,
                                  [[maybe_unused]] int subcycle, [[maybe_unused]] double dts) {
    print_error("*** Error: AdvanceOneTimeStep function not defined.\n");
    exit_mpi();}

  virtual void Destroy();

  //! All the tasks that are done at the end of a time-step, independent of time integrator
  void UpdateSolutionAfterTimeStepping(SpaceVariable3D &V, SpaceVariable3D &ID,
                                       vector<SpaceVariable3D*> &Phi,
                                       vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS,
                                       SpaceVariable3D *L,
                                       double time, int time_step, int subcycle, double dts);
                                        
  //! Functions pertaining to steady-state computation
  bool Converged() {return sso ? sso->Converged() : false;} //only for steady-state simulations
  double GetRelativeResidual1Norm() {assert(sso); return sso->GetRelativeResidual1Norm();} //function L1 norm
  double GetRelativeResidual2Norm() {assert(sso); return sso->GetRelativeResidual2Norm();} //function L2 norm
  double GetRelativeResidualInfNorm() {assert(sso); return sso->GetRelativeResidualInfNorm();}

protected:

  //! compute U += a*dt*R, where dt can be different for different cells (for steady-state computation)
  void AddFluxWithLocalTimeStep(SpaceVariable3D &U, double a, SpaceVariable3D *Dt, SpaceVariable3D &R);

};

/********************************************************************
 * Numerical time-integrator: Forward Euler (TVD)
 *******************************************************************/
class TimeIntegratorFE : public TimeIntegratorBase
{
  //! conservative state variable at time n
  SpaceVariable3D Un;

  //! "residual", i.e. the right-hand-side of the ODE
  SpaceVariable3D Rn;  
  vector<SpaceVariable3D*> Rn_ls;
  SpaceVariable3D *Rn_xi;

public:
  TimeIntegratorFE(MPI_Comm &comm_, IoData& iod_, DataManagers3D& dms_, SpaceOperator& spo_,
                   vector<LevelSetOperator*>& lso_, MultiPhaseOperator &mpo_,
                   LaserAbsorptionSolver* laser_, EmbeddedBoundaryOperator* embed_,
                   HyperelasticityOperator* heo_, PrescribedMotionOperator* pmo_);
  ~TimeIntegratorFE();

  void AdvanceOneTimeStep(SpaceVariable3D &V, SpaceVariable3D &ID, 
                          vector<SpaceVariable3D*>& Phi, vector<SpaceVariable3D*> &NPhi,
                          vector<SpaceVariable3D*> &KappaPhi,
                          SpaceVariable3D *L, SpaceVariable3D *Xi, SpaceVariable3D *Vturb,
                          SpaceVariable3D *LocalDt,
                          double time, double dt, int time_step, int subcycle, double dts);

  void Destroy();

};

/********************************************************************
 * Numerical time-integrator: 2nd-order Runge-Kutta (Heun's method, TVD)
 *******************************************************************/
class TimeIntegratorRK2 : public TimeIntegratorBase
{
  //! conservative state variable at time n
  SpaceVariable3D Un;
  //! intermediate state
  SpaceVariable3D U1;
  SpaceVariable3D V1;
  //! "residual", i.e. the right-hand-side of the ODE
  SpaceVariable3D R;  

  //! level set variables
  vector<SpaceVariable3D*> Phi1; 
  vector<SpaceVariable3D*> Rls; 

  //! reference map equation
  SpaceVariable3D* Xi1;
  SpaceVariable3D* Rxi;

public:
  TimeIntegratorRK2(MPI_Comm &comm_, IoData& iod_, DataManagers3D& dms_, SpaceOperator& spo_,
                    vector<LevelSetOperator*>& lso_, MultiPhaseOperator &mpo_,
                    LaserAbsorptionSolver* laser_, EmbeddedBoundaryOperator* embed_,
                    HyperelasticityOperator* heo_, PrescribedMotionOperator* pmo_);
  ~TimeIntegratorRK2();

  void AdvanceOneTimeStep(SpaceVariable3D &V, SpaceVariable3D &ID,
                          vector<SpaceVariable3D*>& Phi, vector<SpaceVariable3D*> &NPhi,
                          vector<SpaceVariable3D*> &KappaPhi,
                          SpaceVariable3D *L, SpaceVariable3D *Xi, SpaceVariable3D *Vturb,
                          SpaceVariable3D *LocalDt,
                          double time, double dt, int time_step, int subcycle, double dts);

  void Destroy(); 

};

/********************************************************************
 * Numerical time-integrator: 3rd-order Runge-Kutta (Gottlieb & Shu, TVD)
 *******************************************************************/
class TimeIntegratorRK3 : public TimeIntegratorBase
{
  //! conservative state variable at time n
  SpaceVariable3D Un;
  //! intermediate state
  SpaceVariable3D U1;
  SpaceVariable3D V1;
  SpaceVariable3D V2;
  //! "residual", i.e. the right-hand-side of the ODE
  SpaceVariable3D R;  

  //!level set variables
  vector<SpaceVariable3D*> Phi1;
  vector<SpaceVariable3D*> Rls;

  //! reference map equation
  SpaceVariable3D* Xi1;
  SpaceVariable3D* Rxi;

public:
  TimeIntegratorRK3(MPI_Comm &comm_, IoData& iod_, DataManagers3D& dms_, SpaceOperator& spo_,
                    vector<LevelSetOperator*>& lso_, MultiPhaseOperator &mpo_,
                    LaserAbsorptionSolver* laser_, EmbeddedBoundaryOperator* embed_,
                    HyperelasticityOperator* heo_, PrescribedMotionOperator* pmo_);
  ~TimeIntegratorRK3();

  void AdvanceOneTimeStep(SpaceVariable3D &V, SpaceVariable3D &ID,
                          vector<SpaceVariable3D*>& Phi, vector<SpaceVariable3D*> &NPhi,
                          vector<SpaceVariable3D*> &KappaPhi,
                          SpaceVariable3D *L, SpaceVariable3D *Xi, SpaceVariable3D *Vturb,
                          SpaceVariable3D *LocalDt,
                          double time, double dt, int time_step, int subcycle, double dts);

  void Destroy(); 

};

//----------------------------------------------------------------------

#endif
