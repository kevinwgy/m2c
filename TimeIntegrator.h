#ifndef _TIME_INTEGRATOR_H_
#define _TIME_INTEGRATOR_H_

#include <SpaceOperator.h>
#include <LevelSetOperator.h>
#include <vector>
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

public:
  TimeIntegratorBase(MPI_Comm &comm_, IoData& iod_, SpaceOperator& spo_, 
                     vector<LevelSetOperator*>& lso_, MultiPhaseOperator& mpo_) : 
      comm(comm_), iod(iod_), spo(spo_), lso(lso_), mpo(mpo_) { }

  virtual ~TimeIntegratorBase() {}

  // Integrate the ODE system for one time-step. Implemented in derived classes
  virtual void AdvanceOneTimeStep(SpaceVariable3D &V, SpaceVariable3D &ID, 
                                  vector<SpaceVariable3D*> &Phi, double dt) {
    print_error("*** Error: AdvanceOneTimeStep function not defined.\n");
    exit_mpi();}

  virtual void Destroy() {
    print_error("*** Error: TimeIntegratorBase::Destroy() is not defined.\n");
    exit_mpi();}

};

/********************************************************************
 * Numerical time-integrator: Forward Euler (TVD)
 *******************************************************************/
class TimeIntegratorFE : public TimeIntegratorBase
{
  //! conservative state variable at time n
  SpaceVariable3D Un;

  //! material id at time n
  SpaceVariable3D IDn;

  //! "residual", i.e. the right-hand-side of the ODE
  SpaceVariable3D Rn;  
  vector<SpaceVariable3D*> Rn_ls;

public:
  TimeIntegratorFE(MPI_Comm &comm_, IoData& iod_, DataManagers3D& dms_, SpaceOperator& spo_,
                   vector<LevelSetOperator*>& lso_, MultiPhaseOperator &mpo_);
  ~TimeIntegratorFE() {}

  void AdvanceOneTimeStep(SpaceVariable3D &V, SpaceVariable3D &ID, 
                          vector<SpaceVariable3D*>& Phi, double dt);

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

  //! material id at time n
  SpaceVariable3D IDn;

public:
  TimeIntegratorRK2(MPI_Comm &comm_, IoData& iod_, DataManagers3D& dms_, SpaceOperator& spo_,
                    vector<LevelSetOperator*>& lso_, MultiPhaseOperator &mpo_);
  ~TimeIntegratorRK2() {}

  void AdvanceOneTimeStep(SpaceVariable3D &V, SpaceVariable3D &ID,
                          vector<SpaceVariable3D*>& Phi, double dt);

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

  //! material id at time n
  SpaceVariable3D IDn;

public:
  TimeIntegratorRK3(MPI_Comm &comm_, IoData& iod_, DataManagers3D& dms_, SpaceOperator& spo_,
                    vector<LevelSetOperator*>& lso_, MultiPhaseOperator &mpo_);
  ~TimeIntegratorRK3() {}

  void AdvanceOneTimeStep(SpaceVariable3D &V, SpaceVariable3D &ID,
                          vector<SpaceVariable3D*>& Phi, double dt);

  void Destroy();

};

//----------------------------------------------------------------------

#endif
