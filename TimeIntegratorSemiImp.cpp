/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<TimeIntegratorSemiImp.h>

//----------------------------------------------------------------------------
// SIMPLE 
//----------------------------------------------------------------------------

TimeIntegratorSIMPLE::TimeIntegratorSIMPLE(MPI_Comm &comm_, IoData& iod_, DataManagers3D& dms_,
                                           SpaceOperator& spo_, IncompressibleOperator &inco_,
                                           vector<LevelSetOperator*>& lso_, MultiPhaseOperator &mpo_,
                                           LaserAbsorptionSolver* laser_, EmbeddedBoundaryOperator* embed_,
                                           HyperelasticityOperator* heo_, PrescribedMotionOperator* pmo_)
                    : TimeIntegratorBase(comm_, iod_, dms_, spo_, lso_, mpo_, laser_, embed_, heo_, pmo_),
                      inco(inco_), Vstar(comm_, &(dms_.ghosted1_1dof)),
                      Pprime(comm_, &(dms_.ghosted1_1dof)), B(comm_, &(dms_.ghosted1_1dof))
{


}

//----------------------------------------------------------------------------

TimeIntegratorSIMPLE::~TimeIntegratorSIMPLE()
{ }

//----------------------------------------------------------------------------

void
TimeIntegratorSIMPLE::Destroy()
{
  Vstar.Destroy();
  Pprime.Destroy();
  B.Destroy();

  TimeIntegratorBase::Destroy();
}

//----------------------------------------------------------------------------

void
TimeIntegratorSIMPLE::AdvanceOneTimeStep(SpaceVariable3D &V, SpaceVariable3D &ID,
                                         vector<SpaceVariable3D*>& Phi, vector<SpaceVariable3D*> &NPhi,
                                         vector<SpaceVariable3D*> &KappaPhi,
                                         SpaceVariable3D *L, SpaceVariable3D *Xi, SpaceVariable3D *LocalDt,
                                         double time, double dt, int time_step, int subcycle, double dts)
{
  print_error("*** Error: TimeIntegratorSIMPLE::AdvanceOneTimeStep has not been implemented yet.\n");

  int maxIter = time_step == 1 ? 10*iod.ts.semi_impl.maxIts : iod.ts.semi_impl.maxIts;
  for(int iter = 0; iter < maxIter; iter++) {

    // build velocity equation

    //  




  }

}

//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
// SIMPLER
//----------------------------------------------------------------------------

TimeIntegratorSIMPLER::TimeIntegratorSIMPLER(MPI_Comm &comm_, IoData& iod_, DataManagers3D& dms_,
                                             SpaceOperator& spo_, IncompressibleOperator &inco_,
                                             vector<LevelSetOperator*>& lso_, MultiPhaseOperator &mpo_,
                                             LaserAbsorptionSolver* laser_, EmbeddedBoundaryOperator* embed_,
                                             HyperelasticityOperator* heo_, PrescribedMotionOperator* pmo_)
                     : TimeIntegratorSIMPLE(comm_, iod_, dms_, spo_, inco_, lso_, mpo_, laser_, embed_, 
                                            heo_, pmo_)
{


}

//----------------------------------------------------------------------------

TimeIntegratorSIMPLER::~TimeIntegratorSIMPLER()
{ }

//----------------------------------------------------------------------------

void
TimeIntegratorSIMPLER::Destroy()
{


  TimeIntegratorBase::Destroy();
}

//----------------------------------------------------------------------------

void
TimeIntegratorSIMPLER::AdvanceOneTimeStep(SpaceVariable3D &V, SpaceVariable3D &ID,
                                          vector<SpaceVariable3D*>& Phi, vector<SpaceVariable3D*> &NPhi,
                                          vector<SpaceVariable3D*> &KappaPhi,
                                          SpaceVariable3D *L, SpaceVariable3D *Xi, SpaceVariable3D *LocalDt,
                                          double time, double dt, int time_step, int subcycle, double dts)
{
  print_error("*** Error: TimeIntegratorSIMPLER::AdvanceOneTimeStep has not been implemented yet.\n");
}

//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
// SIMPLEC
//----------------------------------------------------------------------------

TimeIntegratorSIMPLEC::TimeIntegratorSIMPLEC(MPI_Comm &comm_, IoData& iod_, DataManagers3D& dms_,
                                             SpaceOperator& spo_, IncompressibleOperator &inco_,
                                             vector<LevelSetOperator*>& lso_, MultiPhaseOperator &mpo_,
                                             LaserAbsorptionSolver* laser_, EmbeddedBoundaryOperator* embed_,
                                             HyperelasticityOperator* heo_, PrescribedMotionOperator* pmo_)
                     : TimeIntegratorSIMPLE(comm_, iod_, dms_, spo_, inco_, lso_, mpo_, laser_, embed_, 
                                            heo_, pmo_)
{


}

//----------------------------------------------------------------------------

TimeIntegratorSIMPLEC::~TimeIntegratorSIMPLEC()
{ }

//----------------------------------------------------------------------------

void
TimeIntegratorSIMPLEC::Destroy()
{ 


  TimeIntegratorBase::Destroy();
}

//----------------------------------------------------------------------------

void
TimeIntegratorSIMPLEC::AdvanceOneTimeStep(SpaceVariable3D &V, SpaceVariable3D &ID,
                                          vector<SpaceVariable3D*>& Phi, vector<SpaceVariable3D*> &NPhi,
                                          vector<SpaceVariable3D*> &KappaPhi,
                                          SpaceVariable3D *L, SpaceVariable3D *Xi, SpaceVariable3D *LocalDt,
                                          double time, double dt, int time_step, int subcycle, double dts)
{
  print_error("*** Error: TimeIntegratorSIMPLEC::AdvanceOneTimeStep has not been implemented yet.\n");
}

//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
// PISO 
//----------------------------------------------------------------------------

TimeIntegratorPISO::TimeIntegratorPISO(MPI_Comm &comm_, IoData& iod_, DataManagers3D& dms_,
                                       SpaceOperator& spo_, IncompressibleOperator &inco_,
                                       vector<LevelSetOperator*>& lso_, MultiPhaseOperator &mpo_,
                                       LaserAbsorptionSolver* laser_, EmbeddedBoundaryOperator* embed_,
                                       HyperelasticityOperator* heo_, PrescribedMotionOperator* pmo_)
                     : TimeIntegratorSIMPLE(comm_, iod_, dms_, spo_, inco_, lso_, mpo_, laser_, embed_, 
                                            heo_, pmo_)
{


}

//----------------------------------------------------------------------------

TimeIntegratorPISO::~TimeIntegratorPISO()
{ }

//----------------------------------------------------------------------------

void
TimeIntegratorPISO::Destroy()
{ 


  TimeIntegratorBase::Destroy();
}

//----------------------------------------------------------------------------

void
TimeIntegratorPISO::AdvanceOneTimeStep(SpaceVariable3D &V, SpaceVariable3D &ID,
                                       vector<SpaceVariable3D*>& Phi, vector<SpaceVariable3D*> &NPhi,
                                       vector<SpaceVariable3D*> &KappaPhi,
                                       SpaceVariable3D *L, SpaceVariable3D *Xi, SpaceVariable3D *LocalDt,
                                       double time, double dt, int time_step, int subcycle, double dts)
{
  print_error("*** Error: TimeIntegratorPISO::AdvanceOneTimeStep has not been implemented yet.\n");
}

//----------------------------------------------------------------------------








//----------------------------------------------------------------------------




//----------------------------------------------------------------------------







