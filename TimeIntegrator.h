#ifndef _TIME_INTEGRATOR_H_
#define _TIME_INTEGRATOR_H_

#include <SpaceOperator.h>
#include <LevelSetOperator.h>
#include <MultiPhaseOperator.h>
#include <LaserAbsorptionSolver.h>
#include <RiemannSolutions.h>
using std::vector;

/********************************************************************
 * Numerical time-integrator: The base class
 *******************************************************************/
class TimeIntegratorBase
{
protected:
  MPI_Comm&       comm;
  IoData&         iod;
  SpaceOperator&  spo;

  vector<LevelSetOperator*>& lso;
  MultiPhaseOperator& mpo;  

  //! Laser absorption solver (NULL if laser is not activated)
  LaserAbsorptionSolver* laser;
  SpaceVariable3D* LH; //Laser-induced heating
  bool laser_initialized;

  //!< Internal variable to temporarily store old ID
  SpaceVariable3D IDn;

  //!< Solutions of exact Riemann problems
  RiemannSolutions riemann_solutions;

public:
  TimeIntegratorBase(MPI_Comm &comm_, IoData& iod_, DataManagers3D& dms_, SpaceOperator& spo_, 
                     vector<LevelSetOperator*>& lso_, MultiPhaseOperator& mpo_,
                     LaserAbsorptionSolver* laser_) : 
      comm(comm_), iod(iod_), spo(spo_), lso(lso_), mpo(mpo_), laser(laser_),
      laser_initialized(false), IDn(comm_, &(dms_.ghosted1_1dof)),
      LH(laser_ ? new SpaceVariable3D(comm_, &(dms_.ghosted1_1dof)) : NULL) { }

  virtual ~TimeIntegratorBase() {if(LH) delete LH;}

  // Integrate the ODE system for one time-step. Implemented in derived classes
  virtual void AdvanceOneTimeStep(SpaceVariable3D &V, SpaceVariable3D &ID, 
                                  vector<SpaceVariable3D*> &Phi, SpaceVariable3D *L, double time,
                                  double dt, int time_step) {
    print_error("*** Error: AdvanceOneTimeStep function not defined.\n");
    exit_mpi();}

  virtual void Destroy() {
    IDn.Destroy();
    if(LH) LH->Destroy();
  }

  // all the tasks that are done at the end of a time-step, independent of 
  // time integrator
  void UpdateSolutionAfterTimeStepping(SpaceVariable3D &V, SpaceVariable3D &ID,
                                       vector<SpaceVariable3D*> &Phi, SpaceVariable3D *L,
                                       double time, double dt, int time_step);
                                        
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

public:
  TimeIntegratorFE(MPI_Comm &comm_, IoData& iod_, DataManagers3D& dms_, SpaceOperator& spo_,
                   vector<LevelSetOperator*>& lso_, MultiPhaseOperator &mpo_,
                   LaserAbsorptionSolver* laser_);
  ~TimeIntegratorFE() {}

  void AdvanceOneTimeStep(SpaceVariable3D &V, SpaceVariable3D &ID, 
                          vector<SpaceVariable3D*>& Phi, SpaceVariable3D *L, double time,
                          double dt, int time_step);

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

public:
  TimeIntegratorRK2(MPI_Comm &comm_, IoData& iod_, DataManagers3D& dms_, SpaceOperator& spo_,
                    vector<LevelSetOperator*>& lso_, MultiPhaseOperator &mpo_,
                    LaserAbsorptionSolver* laser_);
  ~TimeIntegratorRK2() {}

  void AdvanceOneTimeStep(SpaceVariable3D &V, SpaceVariable3D &ID,
                          vector<SpaceVariable3D*>& Phi, SpaceVariable3D *L, double time,
                          double dt, int time_step);

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
  //! "residual", i.e. the right-hand-side of the ODE
  SpaceVariable3D R;  

  //!level set variables
  vector<SpaceVariable3D*> Phi1;
  vector<SpaceVariable3D*> Rls;

public:
  TimeIntegratorRK3(MPI_Comm &comm_, IoData& iod_, DataManagers3D& dms_, SpaceOperator& spo_,
                    vector<LevelSetOperator*>& lso_, MultiPhaseOperator &mpo_,
                    LaserAbsorptionSolver* laser_);
  ~TimeIntegratorRK3() {}

  void AdvanceOneTimeStep(SpaceVariable3D &V, SpaceVariable3D &ID,
                          vector<SpaceVariable3D*>& Phi, SpaceVariable3D *L, double time,
                          double dt, int time_step);

  void Destroy();

};

//----------------------------------------------------------------------

#endif
