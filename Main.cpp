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
#include <SpaceVariable.h>
#include <VarFcnSGEuler.h>
#include <SpaceOperator.h>

//--------------------------------------------------------------
// Main Function
//--------------------------------------------------------------
int main(int argc, char* argv[])
{
  clock_t start_time = clock(); //for timing purpose only

  //--------------------------------------------------------------
  // Initialize PETSc and MPI 
  //--------------------------------------------------------------
  PetscInitialize(&argc, &argv, argc>=3 ? argv[2] : (char*)0, (char*)0);
  MPI_Comm comm = PETSC_COMM_WORLD; //be default, this is MPI_COMM_WORLD
  printLogo();

  //--------------------------------------------------------------
  // Read user's input file
  //--------------------------------------------------------------
  IoData iod(argc, argv);

  //--------------------------------------------------------------
  // Setup PETSc data array (da) structure for nodal variables
  //--------------------------------------------------------------
  DataManagers2D dms(comm, iod.mesh.Nx, iod.mesh.Ny);

  //--------------------------------------------------------------
  // Initialize VarFcn (EOS, etc.)
  //--------------------------------------------------------------
  VarFcnSGEuler vf(iod.eqs.fluid1);

  //--------------------------------------------------------------
  // Initialize space operator
  //--------------------------------------------------------------
  SpaceOperator spo(comm, dms, iod, vf);

  //--------------------------------------------------------------
  // Initialize State Variable and Residual
  //--------------------------------------------------------------

  //--------------------------------------------------------------
  // Main Loop
  //--------------------------------------------------------------
  double t = 0.0; //simulation (i.e. physical) time
  print("\n");
  print("----------------------------\n");
  print("--       Main Loop        --\n");
  print("----------------------------\n");


  //TODO





  print("\n");
  print("======================================\n");
  print("   NORMAL TERMINATION (t = %f) \n", t); 
  print("======================================\n");
  print("Total Computation Time: %f sec.\n", ((double)(clock()-start_time))/CLOCKS_PER_SEC);

  // finalize 
  spo.Destroy();
  dms.DestroyAllDataManagers();
  PetscFinalize();

  return 0;
}

//--------------------------------------------------------------
