/**********************************************************************************
 * Copyright © Kevi Wang, 2020
 * (1) Redistribution and use in source and binary forms, with or without modification,
 *     are permitted, provided that this copyright notice is retained.
 * (2) Use at your own risk.
 **********************************************************************************/
#include <stdio.h>
#include <iostream>
#include <vector>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include <petscdmda.h> //PETSc
#include "hgversion.h" //Version control
#include "input.h"
using namespace std;

//--------------------------------------------------------------
void printLogo(MPI_Comm comm);
//--------------------------------------------------------------

typedef struct {
  PetscReal rho, rhou, rhov, rhoe;
} StateVariables;

//--------------------------------------------------------------
// Main Function
//--------------------------------------------------------------
int main(int argc, char* argv[])
{
  clock_t start_time = clock(); //for timing purpose only

  //--------------------------------------------------------------
  // Read user's input file
  //--------------------------------------------------------------
  IoData iod; 
  iod.readCmdLine(argc, argv);
  iod.readCmdFile();

  //--------------------------------------------------------------
  // Initialize PETSc and MPI 
  //--------------------------------------------------------------
  PetscInitialize(&argc, &argv, argc>=3 ? argv[2] : (char*)0, (char*)0);
  MPI_Comm comm = PETSC_COMM_WORLD; //be default, this is MPI_COMM_WORLD
  int MPI_rank, MPI_size;
  MPI_Comm_rank(comm, &MPI_rank);
  MPI_Comm_size(comm, &MPI_size);
  MPI_Barrier(comm);
  printLogo(comm);
  
     
  //--------------------------------------------------------------
  // Setup PETSc data array (da) structure for nodal variables
  //--------------------------------------------------------------
  DM da;
  int dim = 2; //TODO: use iodata
  auto ierr = DMDACreate2D(comm, DM_BOUNDARY_MIRROR, DM_BOUNDARY_MIRROR, DMDA_STENCIL_BOX, 
                           iod.mesh.Nx, input.mesh.Ny, PETSC_DECIDE, PETSC_DECIDE, 
                           dim+2, 1, 0, 0,
                           &da);
  CHKERRQ(ierr);
  DMSetFromOptions(da);
  DMSetUp(da);
  //DMSetApplicationContext(da, &iod.eqs);

  //--------------------------------------------------------------
  // Initialize space operator
  //--------------------------------------------------------------
  

  //--------------------------------------------------------------
  // Initialize State Variable and Residual
  //--------------------------------------------------------------
  Vec U, F;
  DMCreateGlobalVector(da, &U); U = 0.0;
  InitializeStateVaraibles(U);


  //--------------------------------------------------------------
  // Main Loop
  //--------------------------------------------------------------
  if(!MPI_rank) {
    cout << endl;
    cout << "----------------------------" << endl;
    cout << "--       Main Loop        --" << endl;
    cout << "----------------------------" << endl;
    cout.flush();
  } 
  // calculate time step size
  // KEVIN IS HERE
  while(!lastStep) {
    t += dt;
    iTimeStep++;
    if(!MPI_rank) {
      cout << endl;                                                                                                            cout << "* Time-Step " << iTimeStep << ". Time: " << t << endl; cout.flush();}

    // check if this is the last time step
    //  double dt = input.file.dt;
    if(t >= input.file.t_final - 0.001*dt) {//last time step
      lastStep = true;
      //dt = input.file.t_final - t;
    } 

  DMCreateGlobalVector(da, &F); F = 0.0; 
   

  if(!MPI_rank) {
    cout << "======================================" << endl;
    cout << "   NORMAL TERMINATION (t = " << t << ") " << endl;
    cout << "======================================" << endl;
  }
  MPI_Barrier(comm);
  if(!MPI_rank) cout << "Total Computation Time: " << ((double)(clock()-start_time))/CLOCKS_PER_SEC << " sec." << endl;
//  MPI_Finalize(); //will be called in the destructor of Tao

  return 0;
}

void printLogo(MPI_Comm comm)
{
  int MPI_rank = 0;
  MPI_Comm_rank(comm, &MPI_rank);
  if(!MPI_rank) {
    cout << endl;
    cout << " _____                      " << endl;
    cout << "/__  /  ____  ___  __  __   " << endl;
    cout << "  / /  / __ \\/ _ \\/ / / / " << endl;
    cout << " / /__/ /_/ /  __/ /_/ /    " << endl;
    cout << "/____/\\____/\\___/\\__, /  " << endl;
    cout << "                /____/      " << endl;
    cout << endl;
    cout << " Revision: " << hgRevisionNo << " | " << hgRevisionHash <<endl;
    cout << endl;
    cout.flush();
  }
}
