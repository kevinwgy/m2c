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
                      inco(inco_), VXstar(comm_, &(dms_.ghosted1_1dof)),
                      VYstar(comm_, &(dms_.ghosted1_1dof)), VZstar(comm_, &(dms_.ghosted1_1dof)),
                      Pprime(comm_, &(dms_.ghosted1_1dof)), B(comm_, &(dms_.ghosted1_1dof)),
                      vlin_solver(comm_, dms_.ghosted1_1dof, iod.ts.semi_impl.velocity_linear_solver),
                      plin_solver(comm_, dms_.ghosted1_1dof, iod.ts.semi_impl.pressure_linear_solver)
{


}

//----------------------------------------------------------------------------

TimeIntegratorSIMPLE::~TimeIntegratorSIMPLE()
{ }

//----------------------------------------------------------------------------

void
TimeIntegratorSIMPLE::Destroy()
{
  VXstar.Destroy();
  VYstar.Destroy();
  VZstar.Destroy();
  Pprime.Destroy();
  B.Destroy();

  vlin_solver.Destroy();
  plin_solver.Destroy();

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

  GlobalMeshInfo &global_mesh(spo.GetGlobalMeshInfo());

  int maxIter = time_step == 1 ? 10*iod.ts.semi_impl.maxIts : iod.ts.semi_impl.maxIts;
  for(int iter = 0; iter < maxIter; iter++) {

    Vec5D***   v = (Vec5D***)V.GetDataPointer();
    double*** id = ID.GetDataPointer();

    ExtractVariableComponents(v, VXstar, VYstar, VZstar, Pprime);

    //-----------------------------------------------------
    // Step 1: Solve the momentum equations for u*, v*, w*
    //-----------------------------------------------------

    // Solve the x-momentum equation
    inco.BuildVelocityEquationSIMPLE(0, v, id, vlin_rows, B, iod.ts.semi_impl);
    vlin_solver.SetLinearOperator(vlin_rows);
    vlin_solver.Solve(B, VXstar);

    // Solve the y-momentum equation
    if(!global_mesh.IsMesh1D()) {
      inco.BuildVelocityEquationSIMPLE(1, v, id, vlin_rows, B, iod.ts.semi_impl);
      vlin_solver.SetLinearOperator(vlin_rows);
      vlin_solver.Solve(B, VYstar);
    }

    // Solve the z-momentum equation
    if(!global_mesh.IsMesh1D() && !global_mesh.IsMesh2D()) {
      inco.BuildVelocityEquationSIMPLE(2, v, id, vlin_rows, B, iod.ts.semi_impl);
      vlin_solver.SetLinearOperator(vlin_rows);
      vlin_solver.Solve(B, VZstar);
    }

    
    //-----------------------------------------------------
    // Step 2: Solve the p' equation
    //-----------------------------------------------------
    inco.BuildPressureEquationSIMPLE(v, id, VXstar, VYstar, VZstar, plin_rows, B, iod.ts.semi_impl);
    plin_solver.SetLinearOperator(plin_rows);
    plin_solver.Solve(B, Pprime);


    //-----------------------------------------------------
    // Step 3: Update p, u, v, w, and compute relative error
    //-----------------------------------------------------
    double err = UpdateStates(v, id, VXstar, VYstar, VZstar, Pprime); 

    I AM HERE!

    if(err<tol)
      break; 
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







