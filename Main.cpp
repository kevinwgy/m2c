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

  print("\033[0;32m==========================================\033[0m\n");
  print("\033[0;32m                 START                    \033[0m\n"); 
  print("\033[0;32m==========================================\033[0m\n");
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
  print("----------------------------\n");
  print("--       Main Loop        --\n");
  print("----------------------------\n");
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
    spo.ComputeAdvectionFluxes(U, V, F); //pass in both the conservative and the primitive s.v.
    //integreter.integrate
    //TODO

    t += dt;

    if(out.ToWriteSolutionSnapshot(t, dt, time_step))
      out.WriteSolutionSnapshot(t, time_step, V);
  }

  if(t > out.GetLastSnapshotTime()+0.1*dt) //!< write final solution to file (if it has not been written)
    out.WriteSolutionSnapshot(t, time_step, V);

  print("\n");
  print("\033[0;32m==========================================\033[0m\n");
  print("\033[0;32m   NORMAL TERMINATION (t = %e)  \033[0m\n", t); 
  print("\033[0;32m==========================================\033[0m\n");
  print("Total Computation Time: %f sec.\n", ((double)(clock()-start_time))/CLOCKS_PER_SEC);
  print("\n");

  //! finalize 
  U.Destroy();
  V.Destroy();
  F.Destroy();

  out.FinalizeOutput();
  spo.Destroy();
  dms.DestroyAllDataManagers();
  PetscFinalize();

  return 0;
}

//--------------------------------------------------------------
