/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<TimeIntegrator.h>
#include<EmbeddedBoundaryDataSet.h>
#include<memory> //unique_ptr
using std::cout;
using std::endl;
using std::unique_ptr;

//extern std::ofstream lam_file;
//extern std::ofstream de_file;
//extern std::ofstream residual_file;

//----------------------------------------------------------------------------
// BASE
//----------------------------------------------------------------------------

TimeIntegratorBase::TimeIntegratorBase(MPI_Comm &comm_, IoData& iod_, DataManagers3D& dms_, 
                        SpaceOperator& spo_, vector<LevelSetOperator*>& lso_, MultiPhaseOperator& mpo_,
                        LaserAbsorptionSolver* laser_, EmbeddedBoundaryOperator* embed_,
                        HyperelasticityOperator* heo_, PrescribedMotionOperator* pmo_)
                  : comm(comm_), iod(iod_), spo(spo_), lso(lso_), mpo(mpo_), laser(laser_), embed(embed_),
                    heo(heo_), pmo(pmo_), IDn(comm_, &(dms_.ghosted1_1dof)), sso(NULL),
                    local_time_stepping(iod.ts.local_dt == TsData::YES)
{

  type = NONE;

  for(int i=0; i<(int)lso.size(); i++) {
    ls_mat_id.push_back(lso[i]->GetMaterialID());

    Phi_tmp.push_back(new SpaceVariable3D(comm_, &(dms_.ghosted1_1dof)));
  }

  if(local_time_stepping)
    assert(lso.size()==0);

  if(iod.ts.convergence_tolerance>0.0) //steady-state analysis
    sso = new SteadyStateOperator(comm_, iod.ts, spo.GetGlobalMeshInfo());
}

//----------------------------------------------------------------------------

TimeIntegratorBase::~TimeIntegratorBase()
{
  if(sso)
    delete sso;
} 

//----------------------------------------------------------------------------

void
TimeIntegratorBase::Destroy()
{
  IDn.Destroy();

  for(int i=0; i<(int)Phi_tmp.size(); i++) {
    Phi_tmp[i]->Destroy(); 
    delete Phi_tmp[i];
  }

  if(sso)
    sso->Destroy();
}

//----------------------------------------------------------------------------

void
TimeIntegratorBase::AddFluxWithLocalTimeStep(SpaceVariable3D &U, double alpha,
                                             SpaceVariable3D *Dt, SpaceVariable3D &R)
{
  int dof = U.NumDOF();
  assert(dof == R.NumDOF());
  assert(Dt);
  assert(Dt->NumDOF()==1);

  //add flux within domain interior
  int i0, j0, k0, imax, jmax, kmax;
  U.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);

  double*** u  = U.GetDataPointer();
  double*** dt = Dt->GetDataPointer();
  double*** r  = R.GetDataPointer();

  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++)
        for(int p=0; p<dof; p++)
          u[k][j][i*dof+p] += alpha*dt[k][j][i]*r[k][j][i*dof+p];        

  R.RestoreDataPointerToLocalVector();
  Dt->RestoreDataPointerToLocalVector();

  U.RestoreDataPointerAndInsert();
}

//----------------------------------------------------------------------------



//----------------------------------------------------------------------------
// FORWARD EULER
//----------------------------------------------------------------------------

TimeIntegratorFE::TimeIntegratorFE(MPI_Comm &comm_, IoData& iod_, DataManagers3D& dms_, 
                      SpaceOperator& spo_, vector<LevelSetOperator*>& lso_, MultiPhaseOperator& mpo_,
                      LaserAbsorptionSolver* laser_, EmbeddedBoundaryOperator* embed_,
                      HyperelasticityOperator* heo_, PrescribedMotionOperator* pmo_)
                : TimeIntegratorBase(comm_, iod_, dms_, spo_, lso_, mpo_, laser_, embed_, heo_, pmo_),
                  Un(comm_, &(dms_.ghosted1_5dof)),
                  Rn(comm_, &(dms_.ghosted1_5dof)), Rn_xi(NULL)
{
  type = FORWARD_EULER;

  for(int i=0; i<(int)lso.size(); i++) {
    Rn_ls.push_back(new SpaceVariable3D(comm_, &(dms_.ghosted1_1dof)));
  }

  if(heo_)
    Rn_xi = new SpaceVariable3D(comm_, &(dms_.ghosted1_3dof));
}

//----------------------------------------------------------------------------

TimeIntegratorFE::~TimeIntegratorFE()
{
  for(int i=0; i<(int)Rn_ls.size(); i++)
    delete Rn_ls[i];

  if(Rn_xi)
    delete Rn_xi; 

  //base class destructor is called automatically
}

//----------------------------------------------------------------------------

void TimeIntegratorFE::Destroy()
{
  Un.Destroy();
  Rn.Destroy();

  for(int i=0; i<(int)Rn_ls.size(); i++) {
    Rn_ls[i]->Destroy();
  }

  if(Rn_xi)
    Rn_xi->Destroy();

  TimeIntegratorBase::Destroy();
}

//----------------------------------------------------------------------------

void
TimeIntegratorFE::AdvanceOneTimeStep(SpaceVariable3D &V, SpaceVariable3D &ID, 
                                     vector<SpaceVariable3D*>& Phi, vector<SpaceVariable3D*>& NPhi,
                                     vector<SpaceVariable3D*>& KappaPhi,
                                     SpaceVariable3D *L, SpaceVariable3D *Xi,
                                     [[maybe_unused]] SpaceVariable3D *Vturb,
                                     SpaceVariable3D *Dt,
                                     double time, double dt, int time_step, int subcycle, double dts)
{

  // Store a copy of V at OVERSET ghost nodes (for boundary condition update)
  spo.UpdateOversetGhostNodes(V);

  // Make a copy of Phi for update of material ID. 
  if(time_step == 1) { // Copy entire domain, even in the case of narrow-band LS
    for(int i=0; i<(int)Phi.size(); i++)
      Phi_tmp[i]->AXPlusBY(0.0, 1.0, *Phi[i], true); // setting Phi_tmp[i] = Phi[i], including external ghosts
  } else {
    for(int i=0; i<(int)Phi.size(); i++)
      lso[i]->AXPlusBY(0.0, *Phi_tmp[i], 1.0, *Phi[i], true); //in case of narrow-band, go over only useful nodes
  }

  // Get embedded boundary data
  unique_ptr<vector<unique_ptr<EmbeddedBoundaryDataSet> > > EBDS 
    = embed ? embed->GetPointerToEmbeddedBoundaryData() : nullptr;

   //de_file << std::setw(8) << std::scientific << time << "    " << std::setw(8) << std::scientific << dt  << "       "  << std::setw(8);
   //residual_file << std::setw(8) << std::scientific << time << "    " << std::setw(8) << std::scientific << dt  << "       "  << std::setw(8);

  // -------------------------------------------------------------------------------
  // Forward Euler step for the N-S equations: U(n+1) = U(n) + dt*R(V(n))
  // -------------------------------------------------------------------------------
  spo.ComputeResidual(V, ID, Rn, time-dt, &riemann_solutions, &ls_mat_id, &Phi, &KappaPhi, EBDS.get(),
                      Xi); //compute Rn
  if(laser) laser->AddHeatToNavierStokesResidual(Rn, *L, ID);

  spo.PrimitiveToConservative(V, ID, Un); // get Un
  if(local_time_stepping)
    AddFluxWithLocalTimeStep(Un, 1.0, Dt, Rn);
  else
    Un.AXPlusBY(1.0, dt, Rn);
  spo.ConservativeToPrimitive(Un, ID, V); //updates V = V(n+1)
  spo.ClipDensityAndPressure(V, ID);
  spo.ApplyBoundaryConditions(V);

  // -------------------------------------------------------------------------------
  // Forward Euler step for the level set equation(s): Phi(n+1) = Phi(n) + dt*R(Phi(n))
  // -------------------------------------------------------------------------------
  for(int i=0; i<(int)Phi.size(); i++) {
    lso[i]->ComputeResidual(V, *Phi[i], *Rn_ls[i], time-dt); //compute Rn_ls (level set)
    lso[i]->AXPlusBY(1.0, *Phi[i], dt, *Rn_ls[i]); //in case of narrow-band, go over only useful nodes
    lso[i]->ApplyBoundaryConditions(*Phi[i]);

    if (iod.exact_riemann.surface_tension == ExactRiemannSolverData::YES) {
      lso[i]->ComputeNormal(*Phi[i], *NPhi[i]);
      lso[i]->ApplyBoundaryConditionsNPhi(*NPhi[i]);

      lso[i]->ComputeUnitNormalAndCurvature(*Phi[i], *NPhi[i], *KappaPhi[i]);
      lso[i]->ApplyBoundaryConditionsNPhi(*NPhi[i]);
      lso[i]->ApplyBoundaryConditionsKappaPhi(*KappaPhi[i]); 
    }
  }

  // -------------------------------------------------------------------------------
  // Forward Euler step for the reference map equation: Xi(n+1) = Xi(n) + dt*R(Xi(n))
  // -------------------------------------------------------------------------------
  if(Xi) {
    assert(heo);
    heo->ComputeReferenceMapResidual(V, *Xi, *Rn_xi, time-dt);
    if(local_time_stepping)
      AddFluxWithLocalTimeStep(*Xi, 1.0, Dt, *Rn_xi);
    else
      Xi->AXPlusBY(1.0, dt, *Rn_xi); 
    heo->ApplyBoundaryConditionsToReferenceMap(*Xi);
  }


  // Check of convergence (for steady-state computations)
  if(sso)
    sso->MonitorConvergence(Rn,ID); //Strictly speaking, should recompute R using updated V. But this is OK.
                                    //Rn has been divided by dxdydz


  // -------------------------------------------------------------------------------
  // End-of-step tasks
  // -------------------------------------------------------------------------------
  UpdateSolutionAfterTimeStepping(V, ID, Phi, EBDS.get(), L, time, time_step, subcycle, dts);
}

//----------------------------------------------------------------------------
// SECOND-ORDER RUNGE-KUTTA (HEUN'S METHOD, TVD)
//----------------------------------------------------------------------------

TimeIntegratorRK2::TimeIntegratorRK2(MPI_Comm &comm_, IoData& iod_, DataManagers3D& dms_, 
                                     SpaceOperator& spo_ ,vector<LevelSetOperator*>& lso_,
                                     MultiPhaseOperator& mpo_, LaserAbsorptionSolver* laser_,
                                     EmbeddedBoundaryOperator* embed_,
                                     HyperelasticityOperator* heo_,
                                     PrescribedMotionOperator* pmo_)
                 : TimeIntegratorBase(comm_, iod_, dms_, spo_, lso_, mpo_, laser_, embed_, heo_, pmo_),
                   Un(comm_, &(dms_.ghosted1_5dof)), 
                   U1(comm_, &(dms_.ghosted1_5dof)),
                   V1(comm_, &(dms_.ghosted1_5dof)), 
                   R(comm_, &(dms_.ghosted1_5dof)),
                   Xi1(NULL), Rxi(NULL)
{

  type = RUNGE_KUTTA_2;

  for(int i=0; i<(int)lso.size(); i++) {
    Phi1.push_back(new SpaceVariable3D(comm_, &(dms_.ghosted1_1dof)));
    Rls.push_back(new SpaceVariable3D(comm_, &(dms_.ghosted1_1dof)));
  }

  if(heo_) {
    Xi1 = new SpaceVariable3D(comm_, &(dms_.ghosted1_3dof));
    Rxi = new SpaceVariable3D(comm_, &(dms_.ghosted1_3dof));
  }
}

//----------------------------------------------------------------------------

TimeIntegratorRK2::~TimeIntegratorRK2()
{
  for(int i=0; i<(int)Rls.size(); i++) {
    delete Phi1[i];
    delete Rls[i];
  }

  if(Xi1) delete Xi1;
  if(Rxi) delete Rxi;
}

//----------------------------------------------------------------------------

void TimeIntegratorRK2::Destroy() 
{
  Un.Destroy(); 
  U1.Destroy(); 
  V1.Destroy(); 
  R.Destroy();

  for(int i=0; i<(int)Rls.size(); i++) {
    Phi1[i]->Destroy(); 
    Rls[i]->Destroy(); 
  }

  if(Xi1) Xi1->Destroy();
  if(Rxi) Rxi->Destroy();

  TimeIntegratorBase::Destroy();
}

//----------------------------------------------------------------------------

void
TimeIntegratorRK2::AdvanceOneTimeStep(SpaceVariable3D &V, SpaceVariable3D &ID, 
                                      vector<SpaceVariable3D*>& Phi,
                                      [[maybe_unused]] vector<SpaceVariable3D*>& NPhi,
                                      vector<SpaceVariable3D*>& KappaPhi,
                                      SpaceVariable3D* L, SpaceVariable3D *Xi,
                                      [[maybe_unused]] SpaceVariable3D* Vturb, SpaceVariable3D *Dt,
                                      double time, double dt, 
                                      int time_step, int subcycle, double dts)
{

  // Store a copy of V at OVERSET ghost nodes (for boundary condition update)
  spo.UpdateOversetGhostNodes(V);

  // Make a copy of Phi for update of material ID. 
  if(time_step == 1) { // Copy entire domain, even in the case of narrow-band LS
    for(int i=0; i<(int)Phi.size(); i++)
      Phi_tmp[i]->AXPlusBY(0.0, 1.0, *Phi[i], true); // setting Phi_tmp[i] = Phi[i], including external ghosts
  } else {
    for(int i=0; i<(int)Phi.size(); i++)
      lso[i]->AXPlusBY(0.0, *Phi_tmp[i], 1.0, *Phi[i], true); //in case of narrow-band, go over only useful nodes
  }

  // Get embedded boundary data
  unique_ptr<vector<unique_ptr<EmbeddedBoundaryDataSet> > > EBDS 
    = embed ? embed->GetPointerToEmbeddedBoundaryData() : nullptr;

  //Temporary feature for the stability of Heat Diffusion solver
  //bool run_heat = time < 41; //Xuning: Can be used to skip heat diffusion after a certain time

  //****************** STEP 1 FOR NS ******************
  // Forward Euler step for the N-S equations: U1 = U(n) + dt*R(V(n))
  spo.ComputeResidual(V, ID, R, time-dt, &riemann_solutions, &ls_mat_id, &Phi, &KappaPhi, EBDS.get(),
                      Xi/*, run_heat*/); //->R(V(n))

  if(laser) laser->AddHeatToNavierStokesResidual(R, *L, ID);

  spo.PrimitiveToConservative(V, ID, Un); // get U(n)
  U1.AXPlusBY(0.0, 1.0, Un); //U1 = U(n)
  if(local_time_stepping)
    AddFluxWithLocalTimeStep(U1, 1.0, Dt, R);
  else
    U1.AXPlusBY(1.0, dt, R); //U1 = U1 + dt*R(V(n))

  // Check & clip the intermediate state (U1/V1)
  spo.ConservativeToPrimitive(U1, ID, V1); //get V1
  int clipped = spo.ClipDensityAndPressure(V1, ID);
  if(clipped)
    spo.PrimitiveToConservative(V1, ID, U1); //update U1 after clipping

  // Apply B.C. to the intermediate state (fill ghost cells)
  spo.ApplyBoundaryConditions(V1); 
  //***************************************************


  //****************** STEP 1 FOR LS ****************** 
  // Forward Euler step for the level set equation(s): Phi1 = Phi(n) + dt*R(Phi(n))
  for(int i=0; i<(int)Phi.size(); i++) {
    lso[i]->ComputeResidual(V, *Phi[i], *Rls[i], time-dt); //compute R(Phi(n))
    lso[i]->AXPlusBY(0.0, *Phi1[i], 1.0, *Phi[i]); //in case of narrow-band, go over only useful nodes
    lso[i]->AXPlusBY(1.0, *Phi1[i], dt, *Rls[i]); //in case of narrow-band, go over only useful nodes
    lso[i]->ApplyBoundaryConditions(*Phi1[i]);
  }
  //***************************************************


  //****************** STEP 1 FOR Xi ****************** 
  // Forward Euler step for the reference map equation: Xi1 = Xi(n) + dt*R(Xi(n))
  if(Xi) {
    assert(heo);
    heo->ComputeReferenceMapResidual(V, *Xi, *Rxi, time-dt);
    Xi1->AXPlusBY(0.0, 1.0, *Xi); //Xi1 = Xi(n)
    if(local_time_stepping)
      AddFluxWithLocalTimeStep(*Xi1, 1.0, Dt, *Rxi);
    else
      Xi1->AXPlusBY(1.0, dt, *Rxi); //Xi1 = Xi(n) + dt*R(Xi(n))
    heo->ApplyBoundaryConditionsToReferenceMap(*Xi1);
  }
  //***************************************************



  //****************** STEP 2 FOR NS ******************
  // Step 2: U(n+1) = 0.5*U(n) + 0.5*U1 + 0.5*dt*R(V1)
  //compute R(V1) using prev.Phi, "loose coupling"
  spo.ComputeResidual(V1, ID, R, time, NULL, &ls_mat_id, &Phi, &KappaPhi, EBDS.get(), Xi1/*, run_heat*/);

  if(laser) {
    laser->ComputeLaserRadiance(V1,ID,*L,time,time_step);
    laser->AddHeatToNavierStokesResidual(R, *L, ID);
  }
  U1.AXPlusBY(0.5, 0.5, Un); //U(n+1) = 0.5*U(n) + 0.5*U1;
  if(local_time_stepping)
    AddFluxWithLocalTimeStep(U1, 0.5, Dt, R);
  else
    U1.AXPlusBY(1.0, 0.5*dt, R); //U(n+1) = U(n+1) + 0.5*dt*R(V1)
  
  spo.ConservativeToPrimitive(U1, ID, V); //updates V = V(n+1)
  spo.ClipDensityAndPressure(V, ID);
  spo.ApplyBoundaryConditions(V);
  //***************************************************


  //****************** STEP 2 FOR LS ******************
  // Step 2 for the level set equations: Phi(n+1) = 0.5*Phi(n) + 0.5*Phi1 + 0.5*dt*R(Phi1)
  for(int i=0; i<(int)Phi.size(); i++) {
    lso[i]->ComputeResidual(V1, *Phi1[i], *Rls[i], time);
    lso[i]->AXPlusBY(0.5, *Phi[i], 0.5, *Phi1[i]); //in case of narrow-band, go over only useful nodes
    lso[i]->AXPlusBY(1.0, *Phi[i], 0.5*dt, *Rls[i]); //in case of narrow-band, go over only useful nodes
    lso[i]->ApplyBoundaryConditions(*Phi[i]);
  }
  //***************************************************


  //****************** STEP 2 FOR Xi ****************** 
  // Step 2 for the reference map equation: Xi(n+1) = 0.5*Xi(n) + 0.5*Xi1 + 0.5*dt*R(Xi1)
  if(Xi) {
    assert(heo);
    heo->ComputeReferenceMapResidual(V1, *Xi1, *Rxi, time);
    Xi->AXPlusBY(0.5, 0.5, *Xi1); 
    if(local_time_stepping)
      AddFluxWithLocalTimeStep(*Xi, 0.5, Dt, *Rxi);
    else
      Xi->AXPlusBY(1.0, 0.5*dt, *Rxi); 
    heo->ApplyBoundaryConditionsToReferenceMap(*Xi); //pass t(n+1)
  }
  //***************************************************


  // Check of convergence (for steady-state computations)
  if(sso)
    sso->MonitorConvergence(R,ID); //Strictly speaking, should recompute R using updated V. But this is OK.
                                   //Rn has been divided by dxdydz


  // End-of-step tasks
  UpdateSolutionAfterTimeStepping(V, ID, Phi, EBDS.get(), L, time, time_step, subcycle, dts);
}


//----------------------------------------------------------------------------
// THIRD-ORDER RUNGE-KUTTA (GOTTLIEB & SHU, TVD)
//----------------------------------------------------------------------------
TimeIntegratorRK3::TimeIntegratorRK3(MPI_Comm &comm_, IoData& iod_, DataManagers3D& dms_, 
                                     SpaceOperator& spo_, vector<LevelSetOperator*>& lso_,
                                     MultiPhaseOperator& mpo_, LaserAbsorptionSolver* laser_,
                                     EmbeddedBoundaryOperator* embed_,
                                     HyperelasticityOperator* heo_,
                                     PrescribedMotionOperator* pmo_)
                 : TimeIntegratorBase(comm_, iod_, dms_, spo_, lso_, mpo_, laser_, embed_, heo_, pmo_),
                   Un(comm_, &(dms_.ghosted1_5dof)), 
                   U1(comm_, &(dms_.ghosted1_5dof)),
                   V1(comm_, &(dms_.ghosted1_5dof)), 
                   V2(comm_, &(dms_.ghosted1_5dof)), 
                   R(comm_, &(dms_.ghosted1_5dof)),
                   Xi1(NULL), Rxi(NULL)
                
{

  type = RUNGE_KUTTA_3;

  for(int i=0; i<(int)lso.size(); i++) {
    Phi1.push_back(new SpaceVariable3D(comm_, &(dms_.ghosted1_1dof)));
    Rls.push_back(new SpaceVariable3D(comm_, &(dms_.ghosted1_1dof)));
  }

  if(heo_) {
    Xi1 = new SpaceVariable3D(comm_, &(dms_.ghosted1_3dof));
    Rxi = new SpaceVariable3D(comm_, &(dms_.ghosted1_3dof));
  }
}

//----------------------------------------------------------------------------

TimeIntegratorRK3::~TimeIntegratorRK3()
{
  for(int i=0; i<(int)Rls.size(); i++) {
    delete Phi1[i];
    delete Rls[i];
  }

  if(Xi1) delete Xi1;
  if(Rxi) delete Rxi;
}

//----------------------------------------------------------------------------

void TimeIntegratorRK3::Destroy()
{
  Un.Destroy();
  U1.Destroy();
  V1.Destroy();
  V2.Destroy();
  R.Destroy();

  for(int i=0; i<(int)Rls.size(); i++) {
    Phi1[i]->Destroy();
    Rls[i]->Destroy();
  }

  if(Xi1) Xi1->Destroy();
  if(Rxi) Rxi->Destroy();

  TimeIntegratorBase::Destroy();
}

//----------------------------------------------------------------------------

void
TimeIntegratorRK3::AdvanceOneTimeStep(SpaceVariable3D &V, SpaceVariable3D &ID, 
                                      vector<SpaceVariable3D*>& Phi, 
                                      [[maybe_unused]] vector<SpaceVariable3D*>& NPhi,
                                      vector<SpaceVariable3D*>& KappaPhi, 
                                      SpaceVariable3D* L, SpaceVariable3D* Xi,
                                      [[maybe_unused]] SpaceVariable3D* Vturb, SpaceVariable3D *Dt,
                                      double time, double dt,
                                      int time_step, int subcycle, double dts)
{ 

  // Store a copy of V at OVERSET ghost nodes (for boundary condition update)
  spo.UpdateOversetGhostNodes(V);

  // Make a copy of Phi for update of material ID. 
  if(time_step == 1) { // Copy entire domain, even in the case of narrow-band LS
    for(int i=0; i<(int)Phi.size(); i++)
      Phi_tmp[i]->AXPlusBY(0.0, 1.0, *Phi[i], true); // setting Phi_tmp[i] = Phi[i], including external ghosts
  } else {
    for(int i=0; i<(int)Phi.size(); i++)
      lso[i]->AXPlusBY(0.0, *Phi_tmp[i], 1.0, *Phi[i], true); //in case of narrow-band, go over only useful nodes
  }

  // Get embedded boundary data
  unique_ptr<vector<unique_ptr<EmbeddedBoundaryDataSet> > > EBDS 
    = embed ? embed->GetPointerToEmbeddedBoundaryData() : nullptr;


  //****************** STEP 1 FOR NS ******************
  // Forward Euler step: U1 = U(n) + dt*R(V(n))
  spo.ComputeResidual(V, ID, R, time-dt, &riemann_solutions, &ls_mat_id, &Phi, &KappaPhi, EBDS.get(),
                      Xi); //->R(V(n))

  if(laser) laser->AddHeatToNavierStokesResidual(R, *L, ID);

  spo.PrimitiveToConservative(V, ID, Un); // get U(n)
  U1.AXPlusBY(0.0, 1.0, Un); //U1 = U(n)
  if(local_time_stepping)
    AddFluxWithLocalTimeStep(U1, 1.0, Dt, R);
  else
    U1.AXPlusBY(1.0, dt, R); //U1 = U1 + dt*R(V(n))

  // Check & clip the intermediate state (U1/V1)
  spo.ConservativeToPrimitive(U1, ID, V1); //get V1
  int clipped = spo.ClipDensityAndPressure(V1, ID);
  if(clipped)
    spo.PrimitiveToConservative(V1, ID, U1); //update U1 after clipping

  // Apply B.C. to the intermediate state (fill ghost cells)
  spo.ApplyBoundaryConditions(V1); 
  //***************************************************


  //****************** STEP 1 FOR LS ******************
  // Forward Euler step for the level set equation(s): Phi1 = Phi(n) + dt*R(Phi(n))
  for(int i=0; i<(int)Phi.size(); i++) {
    lso[i]->ComputeResidual(V, *Phi[i], *Rls[i], time-dt); //compute R(Phi(n))
    lso[i]->AXPlusBY(0.0, *Phi1[i], 1.0, *Phi[i]); //in case of narrow-band, go over only useful nodes
    lso[i]->AXPlusBY(1.0, *Phi1[i], dt, *Rls[i]); //in case of narrow-band, go over only useful nodes
    lso[i]->ApplyBoundaryConditions(*Phi1[i]);
  }
  //***************************************************


  //****************** STEP 1 FOR Xi ****************** 
  // Forward Euler step for the reference map equation: Xi1 = Xi(n) + dt*R(Xi(n))
  if(Xi) {
    assert(heo);
    heo->ComputeReferenceMapResidual(V, *Xi, *Rxi, time-dt);
    Xi1->AXPlusBY(0.0, 1.0, *Xi); //Xi1 = Xi(n)
    if(local_time_stepping)
      AddFluxWithLocalTimeStep(*Xi1, 1.0, Dt, *Rxi);
    else
      Xi1->AXPlusBY(1.0, dt, *Rxi); //Xi1 = Xi(n) + dt*R(Xi(n))
    heo->ApplyBoundaryConditionsToReferenceMap(*Xi1);
  }
  //***************************************************



  //****************** STEP 2 FOR NS ******************
  // Step 2: U2 = 0.75*U(n) + 0.25*U1 + 0.25*dt*R(V1))
  //compute R(V1) using prev.Phi, "loose coupling"
  spo.ComputeResidual(V1, ID, R, time, NULL, &ls_mat_id, &Phi, &KappaPhi, EBDS.get(), Xi1);

  if(laser) {
    laser->ComputeLaserRadiance(V1,ID,*L,time,time_step);
    laser->AddHeatToNavierStokesResidual(R, *L, ID);
  }

  U1.AXPlusBY(0.25, 0.75, Un); //U2 = 0.75*U(n) + 0.25*U1;
  if(local_time_stepping)
    AddFluxWithLocalTimeStep(U1, 0.25, Dt, R);
  else
    U1.AXPlusBY(1.0, 0.25*dt, R); //U2 = U2 + 0.25*dt*R(V1)
  
  // Check & clip the intermediate state (U2/V2)
  spo.ConservativeToPrimitive(U1, ID, V2); //get V2
  clipped = spo.ClipDensityAndPressure(V2, ID);
  if(clipped)
    spo.PrimitiveToConservative(V2, ID, U1); //update U2 after clipping

  // Apply B.C. to the intermediate state (fill ghost cells)
  spo.ApplyBoundaryConditions(V2); //apply B.C. by populating the ghost layer
  //***************************************************


  //****************** STEP 2 FOR LS ******************
  // Step 2: Phi2 = 0.75*Phi(n) + 0.25*Phi1 + 0.25*dt*R(Phi1)
  for(int i=0; i<(int)Phi.size(); i++) {
    lso[i]->ComputeResidual(V1, *Phi1[i], *Rls[i], time);
    lso[i]->AXPlusBY(0.25, *Phi1[i], 0.75, *Phi[i]); //in case of narrow-band, go over only useful nodes
    lso[i]->AXPlusBY(1.0, *Phi1[i], 0.25*dt, *Rls[i]); //in case of narrow-band, go over only useful nodes
    lso[i]->ApplyBoundaryConditions(*Phi1[i]);
  }
  //***************************************************


  //****************** STEP 2 FOR Xi ****************** 
  // Step 2: Xi2 = 0.75*Xi(n) + 0.25*Xi1 + 0.25*dt*R(Xi1)
  if(Xi) {
    assert(heo);
    heo->ComputeReferenceMapResidual(V1, *Xi1, *Rxi, time);
    Xi1->AXPlusBY(0.25, 0.75, *Xi);  //re-use Xi1 to store Xi2
    if(local_time_stepping)
      AddFluxWithLocalTimeStep(*Xi1, 0.25, Dt, *Rxi);
    else
      Xi1->AXPlusBY(1.0, 0.25*dt, *Rxi); 
    heo->ApplyBoundaryConditionsToReferenceMap(*Xi1); //t(n+1)
  }
  //***************************************************



  //****************** STEP 3 FOR NS ******************
  // Step 3: U(n+1) = 1/3*U(n) + 2/3*U2 + 2/3*dt*R(V2)
  //compute R(V2) using prev.Phi,"loose coupling"
  spo.ComputeResidual(V2, ID, R, time-0.5*dt, NULL, &ls_mat_id, &Phi, &KappaPhi, EBDS.get(), Xi1);

  if(laser) {
    laser->ComputeLaserRadiance(V2,ID,*L,time-0.5*dt,time_step);
    laser->AddHeatToNavierStokesResidual(R, *L, ID);
  }
  U1.AXPlusBY(2.0/3.0, 1.0/3.0, Un); //U2 = 1/3*U(n) + 2/3*U2;
  if(local_time_stepping)
    AddFluxWithLocalTimeStep(U1, 2.0/3.0, Dt, R);
  else
    U1.AXPlusBY(1.0, 2.0/3.0*dt, R); //U2 = U2 + 2/3*dt*R(V2)

  spo.ConservativeToPrimitive(U1, ID, V); //updates V = V(n+1)
  spo.ClipDensityAndPressure(V, ID);
  spo.ApplyBoundaryConditions(V);
  //***************************************************


  //****************** STEP 3 FOR LS ******************
  // Step 3: Phi(n+1) = 1/3*Phi(n) + 2/3*Phi2 + 2/3*dt*R(Phi2)
  for(int i=0; i<(int)Phi.size(); i++) {
    lso[i]->ComputeResidual(V2, *Phi1[i], *Rls[i], time-0.5*dt);
    lso[i]->AXPlusBY(1.0/3.0, *Phi[i], 2.0/3.0, *Phi1[i]); //in case of narrow-band, go over only useful nodes
    lso[i]->AXPlusBY(1.0, *Phi[i], 2.0/3.0*dt, *Rls[i]); //in case of narrow-band, go over only useful nodes
    lso[i]->ApplyBoundaryConditions(*Phi[i]);
  }
  //***************************************************


  //****************** STEP 3 FOR Xi ****************** 
  // Step 3: Xi(n+1) = 1/3*Xi(n) + 2/3*Xi2 + 2/3*dt*R(Xi2)
  if(Xi) {
    assert(heo);
    heo->ComputeReferenceMapResidual(V2, *Xi1, *Rxi, time-0.5*dt);
    Xi->AXPlusBY(1.0/3.0, 2.0/3.0, *Xi1); 
    if(local_time_stepping)
      AddFluxWithLocalTimeStep(*Xi, 2.0/3.0, Dt, *Rxi);
    else
      Xi->AXPlusBY(1.0, 2.0/3.0*dt, *Rxi); 
    heo->ApplyBoundaryConditionsToReferenceMap(*Xi); //t(n+1)
  }
  //***************************************************


  // Check of convergence (for steady-state computations)
  if(sso)
    sso->MonitorConvergence(R,ID); //Strictly speaking, should recompute R using updated V. But this is OK.
                                   //Rn has been divided by dxdydz


  // End-of-step tasks
  UpdateSolutionAfterTimeStepping(V, ID, Phi, EBDS.get(), L, time, time_step, subcycle, dts);
}

//----------------------------------------------------------------------------

void
TimeIntegratorBase::UpdateSolutionAfterTimeStepping(SpaceVariable3D &V, SpaceVariable3D &ID,
                                                    vector<SpaceVariable3D*> &Phi,
                                                    vector<unique_ptr<EmbeddedBoundaryDataSet> > *EBDS,
                                                    SpaceVariable3D *L,
                                                    double time, int time_step, int subcycle, double dts)
{
  // Assumes Phi_tmp carries the previous values of Phi (before time-stepping)
   
  if(lso.size()) {

    int resolved_conflicts = 0; //if non-zero, force reinitialization of (all) the level sets

    // Check & fix two things: (1) cells belonging to more than 1 subdomain; (2) cells isolated between
    // material boundaries. (2) is optional, and not done by default (frequency controlled by user).
    resolved_conflicts = mpo.ResolveConflictsInLevelSets(time_step, Phi);
    if(resolved_conflicts) { //enforce b.c. for Phi (only necessary when boundary is touched, can be improved later)
      for(int i=0; i<(int)Phi.size(); i++)
        lso[i]->ApplyBoundaryConditions(*Phi[i]);
    }

    // Update ID, V, and possibly also Phi (skipped for single-material simulations)
    IDn.AXPlusBY(0.0, 1.0, ID);  //IDn = ID

    mpo.UpdateMaterialIDByLevelSet(Phi_tmp, Phi, embed ? embed->GetPointerToIntersectors() : nullptr,
                                   ID); //update mat. id. (including the ghost layer outside the physical domain)

    // Correct ID and Phi to be consistent with embedded surfaces (to avoid "leaking")
    if(EBDS) {
      assert(embed);
      resolved_conflicts += mpo.ResolveConflictsWithEmbeddedSurfaces(Phi, IDn, ID, EBDS, embed->GetPointerToIntersectors());
    }

    // Reinitialize level set (also apply boundary conditions)
    if(resolved_conflicts>0) { //must reinitialize
      if(verbose>=1) {
        print("- Resolved %d conflicts between level sets and embedded boundaries. "
              "Reinitialize level sets...\n", resolved_conflicts);
      }
      for(int i=0; i<(int)Phi.size(); i++) 
        lso[i]->Reinitialize(time, dts, time_step, *Phi[i], 0/*no-special-maxIts*/, true/*"must_do"*/);
    } else {
      if(subcycle==0) { //frequency specified by user
        for(int i=0; i<(int)Phi.size(); i++) 
          lso[i]->Reinitialize(time, dts, time_step, *Phi[i], 0/*no-special-maxIts*/, false/*"must_do"*/);
      }
    }


    vector<Int3> unresolved;
    int nUnresolved = mpo.UpdateStateVariablesAfterInterfaceMotion(IDn, ID, V, riemann_solutions,
                              embed ? embed->GetPointerToIntersectors() : nullptr, unresolved); //update V
    if(nUnresolved) {//note that "unresolved" is not combined over all the subdomains, could be empty for some

      //update phi(s) at the unresolved cells
      if(!embed) {//if embed, directly jump to fixing resolved nodes with "apply_failsafe_density = true"!

        // Make a copy of Phi for update of material ID. 
        if(time_step == 1) { // Copy entire domain, even in the case of narrow-band LS
          for(int i=0; i<(int)Phi.size(); i++)
            Phi_tmp[i]->AXPlusBY(0.0, 1.0, *Phi[i], true); // setting Phi_tmp[i] = Phi[i], including external ghosts
        } else {
          for(int i=0; i<(int)Phi.size(); i++)
            lso[i]->AXPlusBY(0.0, *Phi_tmp[i], 1.0, *Phi[i], true); //in case of narrow-band, go over only useful nodes
        }

        mpo.UpdateLevelSetsInUnresolvedCells(Phi, unresolved);
        for(int i=0; i<(int)Phi.size(); i++)
          lso[i]->Reinitialize(time, dts, time_step, *Phi[i], 0, true); //will NOT change sign of phi
        mpo.UpdateMaterialIDByLevelSet(Phi_tmp, Phi, embed ? embed->GetPointerToIntersectors() : nullptr,
                                       ID); //should only update ID of unresolved cells      
      }

      mpo.FixUnresolvedNodes(unresolved, IDn, ID, V, embed ? embed->GetPointerToIntersectors() : nullptr,
                             unresolved/*not used*/, true); 
    }

    // add stored latent heat (Lambda) to cells that changed phase due to interface motion
    mpo.AddLambdaToEnthalpyAfterInterfaceMotion(IDn, ID, V);

    spo.ClipDensityAndPressure(V, ID);
    spo.ApplyBoundaryConditions(V);
  }


  // Check for phase transitions
  //lam_file << std::setw(8) << std::scientific << time << "    " << std::setw(8);
  if(lso.size()) {
    vector<int> phi_updated(lso.size(), 0); //0 or 1
    vector<Int3> new_useful_nodes[lso.size()];
    if(mpo.UpdatePhaseTransitions(Phi, ID, V, phi_updated, new_useful_nodes)>0) {//detected phase transitions
      for(int i=0; i<(int)Phi.size(); i++) {
        if(phi_updated[i] != 0)
          lso[i]->ReinitializeAfterPhaseTransition(*Phi[i], new_useful_nodes[i]);
      }
      spo.ClipDensityAndPressure(V, ID);
      spo.ApplyBoundaryConditions(V);
    }
  }


  // Apply smoothing to U (if specified by user)
  if(subcycle==0) //only consider doing this in the first subcycle
    spo.ApplySmoothingFilter(time, dts, time_step, V, ID);


  // Solve laser radiation equation
  if(laser)
    laser->ComputeLaserRadiance(V,ID,*L,time,time_step);

  // Prescribe velocity (if specified by user)
  if(pmo) {
    pmo->UpdateVelocity(V,ID,time);
    spo.ApplyBoundaryConditions(V); //to be safe...
  }

}

//----------------------------------------------------------------------------

