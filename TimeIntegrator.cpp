#include<TimeIntegrator.h>
using std::cout;
using std::endl;

//----------------------------------------------------------------------------

void TimeIntegratorFE::AdvanceOneTimeStep(SpaceVariable2D &V, double dt)
{ 
  spo.ComputeResidual(V, Rn); // compute Rn

  spo.PrimitiveToConservative(V, Un); // get Un

  // Forward Euler step: U(n+1) = U(n) + dt*R(V(n))
  Un.AXPlusBY(1.0, dt, Rn);

  spo.ConservativeToPrimitive(Un,V); //updates V = V(n+1)
}

//----------------------------------------------------------------------------

void TimeIntegratorRK2::AdvanceOneTimeStep(SpaceVariable2D &V, double dt)
{
  // Forward Euler step: U1 = U(n) + dt*R(V(n))
  spo.ComputeResidual(V, R); // compute R = R(V(n))
  spo.PrimitiveToConservative(V, Un); // get U(n)
  U1.AXPlusBY(0.0, 1.0, Un); //U1 = U(n)
  U1.AXPlusBY(1.0, dt, R); //U1 = U1 + dt*R(V(n))

  // Check & clip the intermediate state (U1/V1)
  spo.ConservativeToPrimitive(U1, V1); //get V1
  int clipped = spo.ClipDensityAndPressure(V1);
  if(clipped)
    spo.PrimitiveToConservative(V1,U1); //update U1 after clipping

  // Apply B.C. to the intermediate state (fill ghost cells)
  spo.ApplyBoundaryConditions(V1); 

  // Step 2: U(n+1) = 0.5*U(n) + 0.5*U1 + 0.5*dt*R(V1))
  spo.ComputeResidual(V1, R); //R = R(V1)
  U1.AXPlusBY(0.5, 0.5, Un); //U(n+1) = 0.5*U(n) + 0.5*U1;
  U1.AXPlusBY(1.0, 0.5*dt, R); //U(n+1) = U(n+1) + 0.5*dt*R(V1)
  
  spo.ConservativeToPrimitive(U1,V); //updates V = V(n+1)
}

//----------------------------------------------------------------------------

void TimeIntegratorRK3::AdvanceOneTimeStep(SpaceVariable2D &V, double dt)
{

  // Forward Euler step: U1 = U(n) + dt*R(V(n))
  spo.ComputeResidual(V, R); // get R = R(V(n))
  spo.PrimitiveToConservative(V, Un); // get U(n)
  U1.AXPlusBY(0.0, 1.0, Un); //U1 = U(n)
  U1.AXPlusBY(1.0, dt, R); //U1 = U1 + dt*R(V(n))

  // Check & clip the intermediate state (U1/V1)
  spo.ConservativeToPrimitive(U1, V1); //get V1
  int clipped = spo.ClipDensityAndPressure(V1);
  if(clipped)
    spo.PrimitiveToConservative(V1,U1); //update U1 after clipping

  // Apply B.C. to the intermediate state (fill ghost cells)
  spo.ApplyBoundaryConditions(V1); 

  // Step 2: U2 = 0.75*U(n) + 0.25*U1 + 0.25*dt*R(V1))
  spo.ComputeResidual(V1, R); //R = R(V1)
  U1.AXPlusBY(0.25, 0.75, Un); //U2 = 0.75*U(n) + 0.25*U1;
  U1.AXPlusBY(1.0, 0.25*dt, R); //U2 = U2 + 0.25*dt*R(V1)
  
  // Check & clip the intermediate state (U2/V2)
  spo.ConservativeToPrimitive(U1, V1); //get V2
  clipped = spo.ClipDensityAndPressure(V1);
  if(clipped)
    spo.PrimitiveToConservative(V1,U1); //update U2 after clipping

  // Apply B.C. to the intermediate state (fill ghost cells)
  spo.ApplyBoundaryConditions(V1); //apply B.C. by populating the ghost layer

  // Step 3: U(n+1) = 1/3*U(n) + 2/3*U2 + 2/3*dt*R(V2)
  spo.ComputeResidual(V1, R); //R = R(V2)
  U1.AXPlusBY(2.0/3.0, 1.0/3.0, Un); //U2 = 1/3*U(n) + 2/3*U2;
  U1.AXPlusBY(1.0, 2.0/3.0*dt, R); //U2 = U2 + 2/3*dt*R(V2)

  spo.ConservativeToPrimitive(U1,V); //updates V = V(n+1)
}

//----------------------------------------------------------------------------

