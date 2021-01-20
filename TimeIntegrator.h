#ifndef _TIME_INTEGRATOR_H_
#define _TIME_INTEGRATOR_H_

#include <SpaceOperator.h>

/********************************************************************
 * Numerical time-integrator: The base class
 *******************************************************************/
class TimeIntegratorBase
{
protected:
  MPI_Comm&       comm;
  IoData&         iod;
  SpaceOperator&  spo;

public:
  TimeIntegratorBase(MPI_Comm &comm_, IoData& iod_, SpaceOperator& spo_) : 
      comm(comm_), iod(iod_), spo(spo_) {/*nothing else to be done*/}

  virtual ~TimeIntegratorBase() {}

  // Integrate the ODE system for one time-step. Implemented in derived classes
  virtual void AdvanceOneTimeStep(SpaceVariable2D &V, double dt) {
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
  SpaceVariable2D Un;
  //! "residual", i.e. the right-hand-side of the ODE
  SpaceVariable2D Rn;  

public:
  TimeIntegratorFE(MPI_Comm &comm_, IoData& iod_, DataManagers2D& dms_, SpaceOperator& spo_) :
    TimeIntegratorBase(comm_, iod_, spo_), 
    Un(comm_, &(dms_.ghosted1_5dof)), Rn(comm_, &(dms_.ghosted1_5dof)) {/*nothing else to be done*/}

  ~TimeIntegratorFE() {}

  void AdvanceOneTimeStep(SpaceVariable2D &V, double dt);

  void Destroy() {Un.Destroy(); Rn.Destroy();}

};

/********************************************************************
 * Numerical time-integrator: 2nd-order Runge-Kutta (Huen's method, TVD)
 *******************************************************************/
class TimeIntegratorRK2 : public TimeIntegratorBase
{
  //! conservative state variable at time n
  SpaceVariable2D Un;
  //! intermediate state
  SpaceVariable2D U1;
  SpaceVariable2D V1;
  //! "residual", i.e. the right-hand-side of the ODE
  SpaceVariable2D R;  

public:
  TimeIntegratorRK2(MPI_Comm &comm_, IoData& iod_, DataManagers2D& dms_, SpaceOperator& spo_) :
    TimeIntegratorBase(comm_, iod_, spo_), 
    Un(comm_, &(dms_.ghosted1_5dof)), U1(comm_, &(dms_.ghosted1_5dof)), 
    V1(comm_, &(dms_.ghosted1_5dof)), R(comm_, &(dms_.ghosted1_5dof)) {/*nothing else to be done*/}

  ~TimeIntegratorRK2() {}

  void AdvanceOneTimeStep(SpaceVariable2D &V, double dt);

  void Destroy() {Un.Destroy(); U1.Destroy(); V1.Destroy(); R.Destroy();}

};

/********************************************************************
 * Numerical time-integrator: 3rd-order Runge-Kutta (Gottlieb & Shu, TVD)
 *******************************************************************/
class TimeIntegratorRK3 : public TimeIntegratorBase
{
  //! conservative state variable at time n
  SpaceVariable2D Un;
  //! intermediate state
  SpaceVariable2D U1;
  SpaceVariable2D V1;
  //! "residual", i.e. the right-hand-side of the ODE
  SpaceVariable2D R;  

public:
  TimeIntegratorRK3(MPI_Comm &comm_, IoData& iod_, DataManagers2D& dms_, SpaceOperator& spo_) :
    TimeIntegratorBase(comm_, iod_, spo_), 
    Un(comm_, &(dms_.ghosted1_5dof)), U1(comm_, &(dms_.ghosted1_5dof)), 
    V1(comm_, &(dms_.ghosted1_5dof)), R(comm_, &(dms_.ghosted1_5dof)) {/*nothing else to be done*/}

  ~TimeIntegratorRK3() {}

  void AdvanceOneTimeStep(SpaceVariable2D &V, double dt);

  void Destroy() {Un.Destroy(); U1.Destroy(); V1.Destroy(); R.Destroy();}

};

//----------------------------------------------------------------------

#endif
