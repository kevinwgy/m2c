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
#include <FluxFcnGenRoe.h>
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
  VarFcnBase *vf = NULL;
  if(iod.eqs.fluid1.fluid == FluidModelData::STIFFENED_GAS)
    vf = new VarFcnSGEuler(iod.eqs.fluid1, iod.output.verbose);
  else {
    print_error("Error: Unable to initialize variable functions (VarFcn) for the specified fluid model.\n");
    exit_mpi();
  }

  //! Initialize FluxFcn
  FluxFcnBase *ff = NULL;
  if(iod.schemes.ns.flux == SchemeData::ROE)
    ff = new FluxFcnGenRoe(*vf, iod);
  else {
    print_error("Error: Unable to initialize flux calculator (FluxFcn) for the specified numerical method.\n");
    exit_mpi();
  }

  //! Initialize space operator
  SpaceOperator spo(comm, dms, iod, *vf, *ff);

  //! Initialize State Variables
  SpaceVariable2D V(comm, &(dms.ghosted1_5dof)); //!< primitive state variables

  //! Impose initial condition
  spo.SetInitialCondition(V);

  //! Initialize output
  Output out(comm, dms, iod, *vf);
  out.InitializeOutput(spo.GetMeshCoordinates());

  //! Initialize time integrator
  TimeIntegratorBase *integrator = NULL;
  if(iod.ts.type == TsData::EXPLICIT) {
    if(iod.ts.expl.type == FORWARD_EULER)
      integrator = new TimeIntegratorFE(comm, iod, dms, spo);
    else if(iod.ts.expl.type == RUNGE_KUTTA_2)
      integrator = new TimeIntegratorRK2(comm, iod, dms, spo);
    else if(iod.ts.expl.type == RUNGE_KUTTA_3)
      integrator = new TimeIntegratorRK3(comm, iod, dms, spo);
    else {
      print_error("Error: Unable to initialize time integrator for the specified (explicit) method.\n");
      exit_mpi();
    }
  } else {
    print_error("Error: Unable to initialize time integrator for the specified method.\n");
    exit_mpi();
  }



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

    // Apply boundary conditions by filling the ghost layer (outside the physical domain) of V
    spo.ApplyBoundaryConditions(V);

    // Compute time step size
    spo.ComputeTimeStepSize(V, dt, cfl); 

    if(t+dt >= iod.ts.maxTime) { //update dt at the LAST time step so it terminates at maxTime
      cfl *= (iod.ts.maxTime - t)/dt;
      dt = iod.ts.maxTime - t;
    }
    print("Time step %d: t = %e, dt = %e, cfl = %e.\n", time_step, t, dt, cfl);

    //----------------------------------------------------
    // Move forward by one time-step: Update V
    //----------------------------------------------------
    integrator.AdvanceOneTimeStep(V,dt); 
    spo.ClipDensityAndPressure(V);
    //----------------------------------------------------

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
  //! In general, "Destroy" should be called for classes that store Petsc DMDA data (which need to be "destroyed").
  V.Destroy();

  out.FinalizeOutput();
  integrator.Destroy();
  spo.Destroy();
  dms.DestroyAllDataManagers();
  PetscFinalize();

  delete integrator;
  delete ff;
  delete vf;

  return 0;
}

//--------------------------------------------------------------
