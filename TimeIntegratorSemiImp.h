/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _TIME_INTEGRATOR_SEMIIMP_H_
#define _TIME_INTEGRATOR_SEMIIMP_H_

#include <TimeIntegrator.h>
#include <IncompressibleOperator.h>

/********************************************************************
 * SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) 
 * Refs: Patankar's book, Van Doormaal and Raithby, 1984; also see KW's notes
 *******************************************************************/
class TimeIntegratorSIMPLE : public TimeIntegratorBase
{

protected:

  IncompressibleOperator &inco;

  SpaceVariable3D Vstar, Pprime;
  SpaceVariable3D B; //!< generally used as the right-hand-side of linear systems

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
                          SpaceVariable3D *L, SpaceVariable3D *Xi, SpaceVariable3D *LocalDt,
                          double time, double dt, int time_step, int subcycle, double dts);

  void Destroy();

};  



/********************************************************************
 * SIMPLER (R: Revised)
 * Refs: Patankar's book, Van Doormaal and Raithby, 1984; also see KW's notes
 *******************************************************************/
class TimeIntegratorSIMPLER : public TimeIntegratorSIMPLE
{

public:

  TimeIntegratorSIMPLER(MPI_Comm &comm_, IoData& iod_, DataManagers3D& dms_, SpaceOperator& spo_,
                        IncompressibleOperator &inco_, vector<LevelSetOperator*>& lso_,
                        MultiPhaseOperator &mpo_, LaserAbsorptionSolver* laser_,
                        EmbeddedBoundaryOperator* embed_, HyperelasticityOperator* heo_,
                        PrescribedMotionOperator* pmo_);
  ~TimeIntegratorSIMPLER();

  //! Note: virtual functions in the base class are automatically virtual in derived classes.
  void AdvanceOneTimeStep(SpaceVariable3D &V, SpaceVariable3D &ID,
                          vector<SpaceVariable3D*>& Phi, vector<SpaceVariable3D*> &NPhi,
                          vector<SpaceVariable3D*> &KappaPhi,
                          SpaceVariable3D *L, SpaceVariable3D *Xi, SpaceVariable3D *LocalDt,
                          double time, double dt, int time_step, int subcycle, double dts);

  void Destroy();

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

  //! Note: virtual functions in the base class are automatically virtual in derived classes.
  void AdvanceOneTimeStep(SpaceVariable3D &V, SpaceVariable3D &ID,
                          vector<SpaceVariable3D*>& Phi, vector<SpaceVariable3D*> &NPhi,
                          vector<SpaceVariable3D*> &KappaPhi,
                          SpaceVariable3D *L, SpaceVariable3D *Xi, SpaceVariable3D *LocalDt,
                          double time, double dt, int time_step, int subcycle, double dts);

  void Destroy();

};  


/********************************************************************
 * PISO (Pressure-Implicit with Split Operator)
 * Refs: Issa, JCP, 1986; also see KW's notes
 *******************************************************************/
class TimeIntegratorPISO : public TimeIntegratorSIMPLE
{

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
                          SpaceVariable3D *L, SpaceVariable3D *Xi, SpaceVariable3D *LocalDt,
                          double time, double dt, int time_step, int subcycle, double dts);

  void Destroy();

};  












#endif
