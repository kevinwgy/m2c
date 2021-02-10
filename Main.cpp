#include <time.h>
#include <petscdmda.h> //PETSc
#include <Utils.h>
#include <IoData.h>
#include <Output.h>
#include <SpaceVariable.h>
#include <VarFcnSG.h>
#include <VarFcnMG.h>
#include <FluxFcnGenRoe.h>
#include <FluxFcnLLF.h>
#include <SpaceOperator.h>
#include <TimeIntegrator.h>
#include <set>
#include <vector>
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
  printHeader(argc, argv);

  print("\033[0;32m==========================================\033[0m\n");
  print("\033[0;32m                 START                    \033[0m\n"); 
  print("\033[0;32m==========================================\033[0m\n");
  print("\n");

  //! Read user's input file
  IoData iod(argc, argv);

  //! Setup PETSc data array (da) structure for nodal variables
  DataManagers3D dms(comm, iod.mesh.Nx, iod.mesh.Ny, iod.mesh.Nz);

  //! Initialize VarFcn (EOS, etc.)
  VarFcnBase *vf = NULL;
  if(iod.eqs.material1.eos == MaterialModelData::STIFFENED_GAS)
    vf = new VarFcnSG(iod.eqs.material1, iod.output.verbose);
  else if(iod.eqs.material1.eos == MaterialModelData::MIE_GRUNEISEN)
    vf = new VarFcnMG(iod.eqs.material1, iod.output.verbose);
  else {
    print_error("Error: Unable to initialize variable functions (VarFcn) for the specified material model.\n");
    exit_mpi();
  }

  //! Initialize FluxFcn
  FluxFcnBase *ff = NULL;
  if(iod.schemes.ns.flux == SchemeData::ROE)
    ff = new FluxFcnGenRoe(vf, iod);
  else if(iod.schemes.ns.flux == SchemeData::LOCAL_LAX_FRIEDRICHS)
    ff = new FluxFcnLLF(vf, iod);
  else {
    print_error("Error: Unable to initialize flux calculator (FluxFcn) for the specified numerical method.\n");
    exit_mpi();
  }

  //! Initialize space operator
  SpaceOperator spo(comm, dms, iod, *vf, *ff);

  //! Initialize State Variables
  SpaceVariable3D V(comm, &(dms.ghosted1_5dof)); //!< primitive state variables

  //! Impose initial condition
  spo.SetInitialCondition(V);

  //! Initialize Levelset(s)
  std::vector<LevelSetOperator*> lso;
  std::vector<SpaceVariable3D*>  Phi;
  std::set<int> ls_tracker;
  for(int i=0; i<SchemesData::MAXLS; i++) {
    int matid = iod.schemes.ls[i].materialid;
    if(matid<=0)
      continue; //levelsets are for tracking subdomains w/ materialid>0
    if(ls_tracker.find(matid) != ls_tracker.end()) {
      print_error("Error: Cannot initialize multiple level sets for the same material (id=%d).\n", matid);
      exit_mpi();
    }
    ls_tracker.insert(matid);    
    lso.push_back(new LevelSetOperator(comm, dms, iod, iod.schemes.ls[i], spo));
    Phi.push_back(new SpaceVariable3D(comm, &(dms.ghosted1_1dof)));

    lso.back()->SetInitialCondition(*Phi.back());

    print("- Initialized level set function (%d) for tracking the boundary of material %d.\n", 
          lso.size()-1, matid);
  }  
  
  //! Initialize output
  Output out(comm, dms, iod, *vf);
  out.InitializeOutput(spo.GetMeshCoordinates());

  //! Initialize time integrator
  TimeIntegratorBase *integrator = NULL;
  if(iod.ts.type == TsData::EXPLICIT) {
    if(iod.ts.expl.type == ExplicitData::FORWARD_EULER)
      integrator = new TimeIntegratorFE(comm, iod, dms, spo, lso);
    else if(iod.ts.expl.type == ExplicitData::RUNGE_KUTTA_2)
      integrator = new TimeIntegratorRK2(comm, iod, dms, spo, lso);
    else if(iod.ts.expl.type == ExplicitData::RUNGE_KUTTA_3)
      integrator = new TimeIntegratorRK3(comm, iod, dms, spo, lso);
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
  out.WriteSolutionSnapshot(t, time_step, V, Phi);

  while(t<iod.ts.maxTime && time_step<iod.ts.maxIts) {

    time_step++;

    // Apply boundary conditions by filling the ghost layer (outside the physical domain).
    spo.ApplyBoundaryConditions(V);
    for(int ls=0; ls<lso.size(); ls++) 
      lso[ls]->ApplyBoundaryConditions(*Phi[ls]);
    
    // Compute time step size
    spo.ComputeTimeStepSize(V, dt, cfl); 

    if(t+dt >= iod.ts.maxTime) { //update dt at the LAST time step so it terminates at maxTime
      cfl *= (iod.ts.maxTime - t)/dt;
      dt = iod.ts.maxTime - t;
    }
    
    print("Step %d: t = %e, dt = %e, cfl = %.4e. Computation time: %.4e s.\n", time_step, t, dt, cfl, 
          ((double)(clock()-start_time))/CLOCKS_PER_SEC);

    //----------------------------------------------------
    // Move forward by one time-step: Update V and Phi
    //----------------------------------------------------
    integrator->AdvanceOneTimeStep(V, Phi, dt); 
    spo.ClipDensityAndPressure(V);
    //----------------------------------------------------

    t += dt;

    if(out.ToWriteSolutionSnapshot(t, dt, time_step))
      out.WriteSolutionSnapshot(t, time_step, V, Phi);

  }

  if(t > out.GetLastSnapshotTime()+0.1*dt) //!< write final solution to file (if it has not been written)
    out.WriteSolutionSnapshot(t, time_step, V, Phi);

  print("\n");
  print("\033[0;32m==========================================\033[0m\n");
  print("\033[0;32m   NORMAL TERMINATION (t = %e)  \033[0m\n", t); 
  print("\033[0;32m==========================================\033[0m\n");
  print("Total Computation Time: %f sec.\n", ((double)(clock()-start_time))/CLOCKS_PER_SEC);
  print("\n");

  //! finalize 
  //! In general, "Destroy" should be called for classes that store Petsc DMDA data (which need to be "destroyed").
  V.Destroy();


  //! Detroy the levelsets
  for(int ls = 0; ls<lso.size(); ls++) {
    Phi[ls]->Destroy(); delete Phi[ls];
    lso[ls]->Destroy(); delete lso[ls];
  }

  out.FinalizeOutput();
  integrator->Destroy();
  spo.Destroy();
  dms.DestroyAllDataManagers();
  PetscFinalize();

  delete integrator;
  delete ff;
  delete vf;

  return 0;
}

//--------------------------------------------------------------
