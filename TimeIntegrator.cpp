#include<TimeIntegrator.h>
using std::cout;
using std::endl;

//----------------------------------------------------------------------------
// FORWARD EULER
//----------------------------------------------------------------------------

TimeIntegratorFE::TimeIntegratorFE(MPI_Comm &comm_, IoData& iod_, DataManagers3D& dms_, 
                      SpaceOperator& spo_, vector<LevelSetOperator*>& lso_, MultiPhaseOperator& mpo_)
                : TimeIntegratorBase(comm_, iod_, spo_, lso_, mpo_),
                  Un(comm_, &(dms_.ghosted1_5dof)),
                  Rn(comm_, &(dms_.ghosted1_5dof)),
                  IDn(comm_, &(dms_.ghosted1_1dof))
{
  for(int i=0; i<lso.size(); i++) {
    Rn_ls.push_back(new SpaceVariable3D(comm_, &(dms_.ghosted1_1dof)));
  }
}

//----------------------------------------------------------------------------

void TimeIntegratorFE::Destroy()
{
  Un.Destroy();
  Rn.Destroy();
  IDn.Destroy();

  for(int i=0; i<Rn_ls.size(); i++) {
    Rn_ls[i]->Destroy(); delete Rn_ls[i];
  }
}

//----------------------------------------------------------------------------

void TimeIntegratorFE::AdvanceOneTimeStep(SpaceVariable3D &V, SpaceVariable3D &ID, 
                                          vector<SpaceVariable3D*>& Phi, double time,
                                          double dt, int time_step)
{
  // Forward Euler step for the N-S equations: U(n+1) = U(n) + dt*R(V(n))
  spo.ComputeResidual(V, ID, Rn, &riemann_solutions); // compute Rn
  spo.PrimitiveToConservative(V, ID, Un); // get Un
  Un.AXPlusBY(1.0, dt, Rn);
  spo.ConservativeToPrimitive(Un, ID, V); //updates V = V(n+1)
  spo.ClipDensityAndPressure(V, ID);
  spo.ApplyBoundaryConditions(V);

  // Forward Euler step for the level set equation(s): Phi(n+1) = Phi(n) + dt*R(Phi(n))
  for(int i=0; i<Phi.size(); i++) {
    lso[i]->ComputeResidual(V, *Phi[i], *Rn_ls[i], time, dt); //compute Rn_ls (level set)
    lso[i]->AXPlusBY(1.0, *Phi[i], dt, *Rn_ls[i]); //in case of narrow-band, go over only useful nodes
    lso[i]->ApplyBoundaryConditions(*Phi[i]);
  }


  // Reinitialize level set (frequency specified by user)
  for(int i=0; i<Phi.size(); i++) {
    lso[i]->Reinitialize(time, dt, time_step, *Phi[i]);
  }


  // Update ID and V (skipped for single-material simulations)
  if(lso.size()) {
    IDn.AXPlusBY(0.0, 1.0, ID);  //IDn = ID
    mpo.UpdateMaterialID(Phi, ID); //update mat. id. (including the ghost layer outside the physical domain)
    mpo.UpdateStateVariablesAfterInterfaceMotion(IDn, ID, V, riemann_solutions); //update V
    spo.ClipDensityAndPressure(V, ID);
    spo.ApplyBoundaryConditions(V);
  }

  // Check for phase transitions
  if(lso.size()) {
    vector<bool> phi_updated(lso.size(), false);
    if(mpo.UpdatePhaseTransitions(Phi, ID, V, phi_updated)>0) {//detected phase transitions
      for(int i=0; i<Phi.size(); i++) {
        if(phi_updated[i])
          lso[i]->Reinitialize(time, dt, time_step, *Phi[i], true/*must do*/);
      }
      spo.ClipDensityAndPressure(V, ID);
      spo.ApplyBoundaryConditions(V);
    }
  }

  // Apply smoothing to U (if specified by user)
  spo.ApplySmoothingFilter(time, dt, time_step, V, ID);

}

//----------------------------------------------------------------------------
// SECOND-ORDER RUNGE-KUTTA (HEUN'S METHOD, TVD)
//----------------------------------------------------------------------------

TimeIntegratorRK2::TimeIntegratorRK2(MPI_Comm &comm_, IoData& iod_, DataManagers3D& dms_, 
                                     SpaceOperator& spo_ ,vector<LevelSetOperator*>& lso_,
                                     MultiPhaseOperator& mpo_)
                 : TimeIntegratorBase(comm_, iod_, spo_, lso_, mpo_),
                   Un(comm_, &(dms_.ghosted1_5dof)), 
                   U1(comm_, &(dms_.ghosted1_5dof)),
                   V1(comm_, &(dms_.ghosted1_5dof)), 
                   R(comm_, &(dms_.ghosted1_5dof)),
                   IDn(comm_, &(dms_.ghosted1_1dof))
{
  for(int i=0; i<lso.size(); i++) {
    Phi1.push_back(new SpaceVariable3D(comm_, &(dms_.ghosted1_1dof)));
    Rls.push_back(new SpaceVariable3D(comm_, &(dms_.ghosted1_1dof)));
  }
}

//----------------------------------------------------------------------------

void TimeIntegratorRK2::Destroy() 
{
  Un.Destroy(); 
  U1.Destroy(); 
  V1.Destroy(); 
  R.Destroy();
  IDn.Destroy();

  for(int i=0; i<Rls.size(); i++) {
    Phi1[i]->Destroy(); delete Phi1[i];
    Rls[i]->Destroy(); delete Rls[i];
  }
}

//----------------------------------------------------------------------------

void TimeIntegratorRK2::AdvanceOneTimeStep(SpaceVariable3D &V, SpaceVariable3D &ID, 
                            vector<SpaceVariable3D*>& Phi, double time, double dt, int time_step)
{

  //****************** STEP 1 FOR NS ******************
  // Forward Euler step for the N-S equations: U1 = U(n) + dt*R(V(n))
  spo.ComputeResidual(V, ID, R, &riemann_solutions); // compute R = R(V(n))
  spo.PrimitiveToConservative(V, ID, Un); // get U(n)
  U1.AXPlusBY(0.0, 1.0, Un); //U1 = U(n)
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
  for(int i=0; i<Phi.size(); i++) {
    lso[i]->ComputeResidual(V, *Phi[i], *Rls[i], time, dt); //compute R(Phi(n))
    lso[i]->AXPlusBY(0.0, *Phi1[i], 1.0, *Phi[i]); //in case of narrow-band, go over only useful nodes
    lso[i]->AXPlusBY(1.0, *Phi1[i], dt, *Rls[i]); //in case of narrow-band, go over only useful nodes
    lso[i]->ApplyBoundaryConditions(*Phi1[i]);
  }
  //***************************************************


  //****************** STEP 2 FOR NS ******************
  // Step 2: U(n+1) = 0.5*U(n) + 0.5*U1 + 0.5*dt*R(V1)
  spo.ComputeResidual(V1, ID, R); //R = R(V1)
  U1.AXPlusBY(0.5, 0.5, Un); //U(n+1) = 0.5*U(n) + 0.5*U1;
  U1.AXPlusBY(1.0, 0.5*dt, R); //U(n+1) = U(n+1) + 0.5*dt*R(V1)
  
  spo.ConservativeToPrimitive(U1, ID, V); //updates V = V(n+1)
  spo.ClipDensityAndPressure(V, ID);
  spo.ApplyBoundaryConditions(V);
  //***************************************************


  //****************** STEP 2 FOR LS ******************
  // Step 2 for the level set equations: Phi(n+1) = 0.5*Phi(n) + 0.5*Phi1 + 0.5*dt*R(Phi1)
  for(int i=0; i<Phi.size(); i++) {
    lso[i]->ComputeResidual(V1, *Phi1[i], *Rls[i], time, dt);
    lso[i]->AXPlusBY(0.5, *Phi[i], 0.5, *Phi1[i]); //in case of narrow-band, go over only useful nodes
    lso[i]->AXPlusBY(1.0, *Phi[i], 0.5*dt, *Rls[i]); //in case of narrow-band, go over only useful nodes
    lso[i]->ApplyBoundaryConditions(*Phi[i]);
  }
  //***************************************************

  // Reinitialize level set (frequency specified by user)
  for(int i=0; i<Phi.size(); i++) {
    lso[i]->Reinitialize(time, dt, time_step, *Phi[i]);
  }


  // Update ID and V (skipped for single-material simulations)
  if(lso.size()) {
    IDn.AXPlusBY(0.0, 1.0, ID);  //IDn = ID
    mpo.UpdateMaterialID(Phi, ID); //update mat. id. (including the ghost layer outside the physical domain)
    mpo.UpdateStateVariablesAfterInterfaceMotion(IDn, ID, V, riemann_solutions); //update V
    spo.ClipDensityAndPressure(V, ID);
    spo.ApplyBoundaryConditions(V);
  }


  // Check for phase transitions
  if(lso.size()) {
    vector<bool> phi_updated(lso.size(), false);
    if(mpo.UpdatePhaseTransitions(Phi, ID, V, phi_updated)>0) {//detected phase transitions
      for(int i=0; i<Phi.size(); i++) {
        if(phi_updated[i])
          lso[i]->Reinitialize(time, dt, time_step, *Phi[i], true/*must do*/);
      }
      spo.ClipDensityAndPressure(V, ID);
      spo.ApplyBoundaryConditions(V);
    }
  }


  // Apply smoothing to U (if specified by user)
  spo.ApplySmoothingFilter(time, dt, time_step, V, ID);

}

//----------------------------------------------------------------------------
// THIRD-ORDER RUNGE-KUTTA (GOTTLIEB & SHU, TVD)
//----------------------------------------------------------------------------
TimeIntegratorRK3::TimeIntegratorRK3(MPI_Comm &comm_, IoData& iod_, DataManagers3D& dms_, 
                                     SpaceOperator& spo_, vector<LevelSetOperator*>& lso_,
                                     MultiPhaseOperator &mpo_)
                 : TimeIntegratorBase(comm_, iod_, spo_, lso_, mpo_),
                   Un(comm_, &(dms_.ghosted1_5dof)), 
                   U1(comm_, &(dms_.ghosted1_5dof)),
                   V1(comm_, &(dms_.ghosted1_5dof)), 
                   R(comm_, &(dms_.ghosted1_5dof)),
                   IDn(comm_, &(dms_.ghosted1_1dof))
                
{
  for(int i=0; i<lso.size(); i++) {
    Phi1.push_back(new SpaceVariable3D(comm_, &(dms_.ghosted1_1dof)));
    Rls.push_back(new SpaceVariable3D(comm_, &(dms_.ghosted1_1dof)));
  }
}

//----------------------------------------------------------------------------

void TimeIntegratorRK3::Destroy()
{
  Un.Destroy();
  U1.Destroy();
  V1.Destroy();
  R.Destroy();
  IDn.Destroy();

  for(int i=0; i<Rls.size(); i++) {
    Phi1[i]->Destroy(); delete Phi1[i];
    Rls[i]->Destroy(); delete Rls[i];
  }
}

//----------------------------------------------------------------------------

void TimeIntegratorRK3::AdvanceOneTimeStep(SpaceVariable3D &V, SpaceVariable3D &ID, 
                                           vector<SpaceVariable3D*>& Phi, double time, 
                                           double dt, int time_step)
{ 

  //****************** STEP 1 FOR NS ******************
  // Forward Euler step: U1 = U(n) + dt*R(V(n))
  spo.ComputeResidual(V, ID, R, &riemann_solutions); // get R = R(V(n))
  spo.PrimitiveToConservative(V, ID, Un); // get U(n)
  U1.AXPlusBY(0.0, 1.0, Un); //U1 = U(n)
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
  for(int i=0; i<Phi.size(); i++) {
    lso[i]->ComputeResidual(V, *Phi[i], *Rls[i], time, dt); //compute R(Phi(n))
    lso[i]->AXPlusBY(0.0, *Phi1[i], 1.0, *Phi[i]); //in case of narrow-band, go over only useful nodes
    lso[i]->AXPlusBY(1.0, *Phi1[i], dt, *Rls[i]); //in case of narrow-band, go over only useful nodes
    lso[i]->ApplyBoundaryConditions(*Phi1[i]);
  }
  //***************************************************


  //****************** STEP 2 FOR NS ******************
  // Step 2: U2 = 0.75*U(n) + 0.25*U1 + 0.25*dt*R(V1))
  spo.ComputeResidual(V1, ID, R); //R = R(V1)
  U1.AXPlusBY(0.25, 0.75, Un); //U2 = 0.75*U(n) + 0.25*U1;
  U1.AXPlusBY(1.0, 0.25*dt, R); //U2 = U2 + 0.25*dt*R(V1)
  
  // Check & clip the intermediate state (U2/V2)
  spo.ConservativeToPrimitive(U1, ID, V1); //get V2
  clipped = spo.ClipDensityAndPressure(V1, ID);
  if(clipped)
    spo.PrimitiveToConservative(V1, ID, U1); //update U2 after clipping

  // Apply B.C. to the intermediate state (fill ghost cells)
  spo.ApplyBoundaryConditions(V1); //apply B.C. by populating the ghost layer
  //***************************************************


  //****************** STEP 2 FOR LS ******************
  // Step 2: Phi2 = 0.75*Phi(n) + 0.2*Phi1 + 0.25*dt*R(Phi1)
  for(int i=0; i<Phi.size(); i++) {
    lso[i]->ComputeResidual(V1, *Phi1[i], *Rls[i], time, dt);
    lso[i]->AXPlusBY(0.25, *Phi1[i], 0.75, *Phi[i]); //in case of narrow-band, go over only useful nodes
    lso[i]->AXPlusBY(1.0, *Phi1[i], 0.25*dt, *Rls[i]); //in case of narrow-band, go over only useful nodes
    lso[i]->ApplyBoundaryConditions(*Phi1[i]);
  }
  //***************************************************


  //****************** STEP 3 FOR NS ******************
  // Step 3: U(n+1) = 1/3*U(n) + 2/3*U2 + 2/3*dt*R(V2)
  spo.ComputeResidual(V1, ID, R); //R = R(V2)
  U1.AXPlusBY(2.0/3.0, 1.0/3.0, Un); //U2 = 1/3*U(n) + 2/3*U2;
  U1.AXPlusBY(1.0, 2.0/3.0*dt, R); //U2 = U2 + 2/3*dt*R(V2)

  spo.ConservativeToPrimitive(U1, ID, V); //updates V = V(n+1)
  spo.ClipDensityAndPressure(V, ID);
  spo.ApplyBoundaryConditions(V);
  //***************************************************


  //****************** STEP 3 FOR LS ******************
  // Step 3: Phi(n+1) = 1/3*Phi(n) + 2/3*Phi2 + 2/3*dt*R(Phi2)
  for(int i=0; i<Phi.size(); i++) {
    lso[i]->ComputeResidual(V1, *Phi1[i], *Rls[i], time, dt);
    lso[i]->AXPlusBY(1.0/3.0, *Phi[i], 2.0/3.0, *Phi1[i]); //in case of narrow-band, go over only useful nodes
    lso[i]->AXPlusBY(1.0, *Phi[i], 2.0/3.0*dt, *Rls[i]); //in case of narrow-band, go over only useful nodes
    lso[i]->ApplyBoundaryConditions(*Phi[i]);
  }
  //***************************************************


  // Reinitialize level set (frequency specified by user)
  for(int i=0; i<Phi.size(); i++) {
    lso[i]->Reinitialize(time, dt, time_step, *Phi[i]);
  }


  // Update ID and V (skipped for single-material simulations)
  if(lso.size()) {
    IDn.AXPlusBY(0.0, 1.0, ID);  //IDn = ID
    mpo.UpdateMaterialID(Phi, ID); //update mat. id. (including the ghost layer outside the physical domain)
    mpo.UpdateStateVariablesAfterInterfaceMotion(IDn, ID, V, riemann_solutions); //update V
    spo.ClipDensityAndPressure(V, ID);
    spo.ApplyBoundaryConditions(V);
  }


  // Check for phase transitions
  if(lso.size()) {
    vector<bool> phi_updated(lso.size(), false);
    if(mpo.UpdatePhaseTransitions(Phi, ID, V, phi_updated)>0) {//detected phase transitions
      for(int i=0; i<Phi.size(); i++) {
        if(phi_updated[i])
          lso[i]->Reinitialize(time, dt, time_step, *Phi[i], true/*must do*/);
      }
      spo.ClipDensityAndPressure(V, ID);
      spo.ApplyBoundaryConditions(V);
    }
  }


  // Apply smoothing to U (if specified by user)
  spo.ApplySmoothingFilter(time, dt, time_step, V, ID);

}

//----------------------------------------------------------------------------

