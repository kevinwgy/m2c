/**********************************************************************************
 * Copyright © Kevn Wang, 2020
 * (1) Redistribution and use in source and binary forms, with or without modification,
 *     are permitted, provided that this copyright notice is retained.
 * (2) Use at your own risk.
 **********************************************************************************/
#include <time.h>
#include <petscdmda.h> //PETSc
#include <Utils.h>
#include <IoData.h>
#include <Output.h>
#include <SpaceVariable.h>
#include <VarFcnSGEuler.h>
#include <SpaceOperator.h>
using std::cout;
using std::endl;
/*************************************
 * Main Function
 ************************************/
int main(int argc, char* argv[])
{
  clock_t start_time = clock(); //for timing purpose only

  //! Initialize PETSc and MPI 
  PetscInitialize(&argc, &argv, argc>=3 ? argv[2] : (char*)0, (char*)0);
  MPI_Comm comm = PETSC_COMM_WORLD; //be default, this is MPI_COMM_WORLD
  printLogo();

  print("\033[0;32m==========================================\n\033[0m");
  print("\033[0;32m                 START \n\033[0m"); 
  print("\033[0;32m==========================================\n\033[0m");
  print("\n");

  //! Read user's input file
  IoData iod(argc, argv);

  //! Setup PETSc data array (da) structure for nodal variables
  DataManagers2D dms(comm, iod.mesh.Nx, iod.mesh.Ny);

  //! Initialize VarFcn (EOS, etc.)
  VarFcnSGEuler vf(iod.eqs.fluid1);

  //! Initialize space operator
  SpaceOperator spo(comm, dms, iod, vf);

  //! Initialize State Variable and fluxes 
  SpaceVariable2D U(comm, &(dms.ghosted1_5dof)); //!< conservative state variables
  SpaceVariable2D V(comm, &(dms.ghosted1_5dof)); //!< primitive state variables

  SpaceVariable2D F(comm, &(dms.ghosted1_5dof)); //!< advection fluxes

  //! Allocate memory for reconstructed cons. s.v. at l,r,t,b faces of each cell
  SpaceVariable2D Ul(comm, &(dms.ghosted1_5dof)); 
  SpaceVariable2D Ur(comm, &(dms.ghosted1_5dof));
  SpaceVariable2D Ut(comm, &(dms.ghosted1_5dof));
  SpaceVariable2D Ub(comm, &(dms.ghosted1_5dof));

  //! Impose initial condition
  spo.SetInitialCondition(V);
  spo.PrimitiveToConservative(V,U);

  //! Initialize output
  Output out(comm, dms, iod, vf);
  out.InitializeOutput(spo.GetMeshCoordinates());

  /*************************************
   * Main Loop 
   ************************************/
  print("\n");
  print("\033[0;32m----------------------------\n\033[0m");
  print("\033[0;32m--       Main Loop        --\n\033[0m");
  print("\033[0;32m----------------------------\n\033[0m");
  double t = 0.0; //!< simulation (i.e. physical) time
  double dt = 0.0;
  double cfl = 0.0;
  int time_step = 0;

  //! write initial condition to file
  out.WriteSolutionSnapshot(t, time_step, V);

  while(t<iod.ts.maxTime && time_step<iod.ts.maxIts) {

    time_step++;
    print("Time step %d: t = %e, dt = %e, cfl = %e.\n", time_step, t, dt, cfl);

    spo.ComputeTimeStepSize(V, dt, cfl); 
    //spo.ReconstructConservativeStateVariables(U, Ul, Ur, Ut, Ub);
    //spo.ComputeConvectionFluxes(Ul, Ur, Ut, Ub, F);
    //integreter.integrate
    //TODO

    t += dt;

    if(out.ToWriteSolutionSnapshot(t, dt, time_step))
      out.WriteSolutionSnapshot(t, time_step, V);
  }

  if(t > out.GetLastSnapshotTime()+0.1*dt) //!< write final solution to file (if it has not been written)
    out.WriteSolutionSnapshot(t, time_step, V);

  print("\n");
  print("\033[0;32m==========================================\n\033[0m");
  print("\033[0;32m   NORMAL TERMINATION (t = %e) \n\033[0m", t); 
  print("\033[0;32m==========================================\n\033[0m");
  print("Total Computation Time: %f sec.\n", ((double)(clock()-start_time))/CLOCKS_PER_SEC);
  print("\n");

  //! finalize 
  U.Destroy();
  Ul.Destroy();
  Ur.Destroy();
  Ut.Destroy();
  Ub.Destroy();
  V.Destroy();
  F.Destroy();

  out.FinalizeOutput();
  spo.Destroy();
  dms.DestroyAllDataManagers();
  PetscFinalize();

  return 0;
}

//--------------------------------------------------------------
