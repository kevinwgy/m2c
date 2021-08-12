#include <time.h>
#include <petscdmda.h> //PETSc
#include <MeshGenerator.h>
#include <Output.h>
#include <VarFcnSG.h>
#include <VarFcnMG.h>
#include <VarFcnJWL.h>
#include <FluxFcnGenRoe.h>
#include <FluxFcnLLF.h>
#include <FluxFcnHLLC.h>
#include <FluxFcnGodunov.h>
#include <SpaceOperator.h>
#include <TimeIntegrator.h>
#include <MultiPhaseOperator.h>
#include <GradientCalculatorCentral.h>
#include <set>
using std::cout;
using std::endl;
int verbose;

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
  verbose = iod.output.verbose;

  //! Calculate mesh coordinates
  vector<double> xcoords, dx, ycoords, dy, zcoords, dz;
  MeshGenerator meshgen;
  meshgen.ComputeMeshCoordinatesAndDeltas(iod.mesh, xcoords, ycoords, zcoords, dx, dy, dz);
  
  //! Setup PETSc data array (da) structure for nodal variables
  DataManagers3D dms(comm, xcoords.size(), ycoords.size(), zcoords.size());

  //! Initialize VarFcn (EOS, etc.) 

  std::vector<VarFcnBase *> vf;
  for(int i=0; i<iod.eqs.materials.dataMap.size(); i++)
    vf.push_back(NULL); //allocate space for the VarFcn pointers

  for(auto it = iod.eqs.materials.dataMap.begin(); it != iod.eqs.materials.dataMap.end(); it++) {
    int matid = it->first;
    if(matid < 0 || matid >= vf.size()) {
      print_error("*** Error: Detected error in the specification of material indices (id = %d).\n", matid);
      exit_mpi();
    }
    if(it->second->eos == MaterialModelData::STIFFENED_GAS)
      vf[matid] = new VarFcnSG(*it->second);
    else if(it->second->eos == MaterialModelData::MIE_GRUNEISEN)
      vf[matid] = new VarFcnMG(*it->second);
    else if(it->second->eos == MaterialModelData::JWL)
      vf[matid] = new VarFcnJWL(*it->second);
    else {
      print_error("*** Error: Unable to initialize variable functions (VarFcn) for the specified material model.\n");
      exit_mpi();
    }
  }

  //! Initialize the exact Riemann problem solver.
  ExactRiemannSolverBase riemann(vf, iod.exact_riemann);


  //! Initialize FluxFcn for the advector flux of the N-S equations
  FluxFcnBase *ff = NULL;
  if(iod.schemes.ns.flux == SchemeData::ROE)
    ff = new FluxFcnGenRoe(vf, iod);
  else if(iod.schemes.ns.flux == SchemeData::LOCAL_LAX_FRIEDRICHS)
    ff = new FluxFcnLLF(vf, iod);
  else if(iod.schemes.ns.flux == SchemeData::HLLC)
    ff = new FluxFcnHLLC(vf, iod);
  else if(iod.schemes.ns.flux == SchemeData::GODUNOV)
    ff = new FluxFcnGodunov(vf, iod);
  else {
    print_error("*** Error: Unable to initialize flux calculator (FluxFcn) for the specified numerical method.\n");
    exit_mpi();
  }

  //! Initialize space operator
  SpaceOperator spo(comm, dms, iod, vf, *ff, riemann, xcoords, ycoords, zcoords, dx, dy, dz);

  //! Initialize interpolator
  InterpolatorBase *interp = NULL;
  if(true) //may add more choices later
    interp = new InterpolatorLinear(comm, dms, spo.GetMeshCoordinates(), spo.GetMeshDeltaXYZ());

  //! Initialize (sptial) gradient calculator
  GradientCalculatorBase *grad = NULL;
  if(true) //may add more choices later
    grad = new GradientCalculatorCentral(comm, dms, spo.GetMeshCoordinates(), spo.GetMeshDeltaXYZ(), *interp);
  
  //! Setup viscosity operator in spo (if viscosity model is not NONE)
  spo.SetupViscosityOperator(interp, grad);

  //! Initialize State Variables
  SpaceVariable3D V(comm, &(dms.ghosted1_5dof)); //!< primitive state variables
  SpaceVariable3D ID(comm, &(dms.ghosted1_1dof)); //!< material id

  //! Impose initial condition
  spo.SetInitialCondition(V, ID);

  //! Initialize Levelset(s)
  std::vector<LevelSetOperator*> lso;
  std::vector<SpaceVariable3D*>  Phi;
  std::set<int> ls_tracker;
  for(auto it = iod.schemes.ls.dataMap.begin(); it != iod.schemes.ls.dataMap.end(); it++) {
    int matid = it->second->materialid;
    if(matid<=0 || matid>=vf.size()) { //cannot use ls to track material 0
      print_error("*** Error: Cannot initialize a level set for tracking material %d.\n", matid);
      exit_mpi();
    }
    if(ls_tracker.find(matid) != ls_tracker.end()) {
      print_error("*** Error: Cannot initialize multiple level sets for the same material (id=%d).\n", matid);
      exit_mpi();
    }
    ls_tracker.insert(matid);    
    lso.push_back(new LevelSetOperator(comm, dms, iod, *it->second, spo));
    Phi.push_back(new SpaceVariable3D(comm, &(dms.ghosted1_1dof)));

    lso.back()->SetInitialCondition(*Phi.back());

    print("- Initialized level set function (%d) for tracking the boundary of material %d.\n", 
          lso.size()-1, matid);
  }  

#ifdef LEVELSET_TEST
  if(!lso.empty()) {
    print("\n");
    print("\033[0;32m- Testing the Level Set Solver using a prescribed velocity field (%d). "
          "N-S solver not activated.\033[0m\n", (int)LEVELSET_TEST);
  }
#endif
  
  //! Initialize multiphase operator (for updating "phase change")
  MultiPhaseOperator mpo(comm, dms, iod, vf, spo, lso);
  mpo.UpdateMaterialID(Phi,ID); //populate the ghost layer of ID (outside domain boundary)

  //! Initialize output
  Output out(comm, dms, iod, vf, spo.GetMeshCellVolumes()); 
  out.InitializeOutput(spo.GetMeshCoordinates());

  //! Initialize time integrator
  TimeIntegratorBase *integrator = NULL;
  if(iod.ts.type == TsData::EXPLICIT) {
    if(iod.ts.expl.type == ExplicitData::FORWARD_EULER)
      integrator = new TimeIntegratorFE(comm, iod, dms, spo, lso, mpo);
    else if(iod.ts.expl.type == ExplicitData::RUNGE_KUTTA_2)
      integrator = new TimeIntegratorRK2(comm, iod, dms, spo, lso, mpo);
    else if(iod.ts.expl.type == ExplicitData::RUNGE_KUTTA_3)
      integrator = new TimeIntegratorRK3(comm, iod, dms, spo, lso, mpo);
    else {
      print_error("*** Error: Unable to initialize time integrator for the specified (explicit) method.\n");
      exit_mpi();
    }
  } else {
    print_error("*** Error: Unable to initialize time integrator for the specified method.\n");
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
  out.OutputSolutions(t, dt, time_step, V, ID, Phi, true/*force_write*/);

  while(t<iod.ts.maxTime && time_step<iod.ts.maxIts) {

    time_step++;

    // Compute time step size
    spo.ComputeTimeStepSize(V, ID, dt, cfl); 

    if(t+dt >= iod.ts.maxTime) { //update dt at the LAST time step so it terminates at maxTime
      cfl *= (iod.ts.maxTime - t)/dt;
      dt = iod.ts.maxTime - t;
    }
    
    print("Step %d: t = %e, dt = %e, cfl = %.4e. Computation time: %.4e s.\n", time_step, t, dt, cfl, 
          ((double)(clock()-start_time))/CLOCKS_PER_SEC);

    //----------------------------------------------------
    // Move forward by one time-step: Update V, Phi, and ID
    //----------------------------------------------------
    t += dt;
    integrator->AdvanceOneTimeStep(V, ID, Phi, t, dt, time_step); 
    //----------------------------------------------------


    out.OutputSolutions(t, dt, time_step, V, ID, Phi, false/*force_write*/);

  }

  out.OutputSolutions(t, dt, time_step, V, ID, Phi, true/*force_write*/);

  print("\n");
  print("\033[0;32m==========================================\033[0m\n");
  print("\033[0;32m   NORMAL TERMINATION (t = %e)  \033[0m\n", t); 
  print("\033[0;32m==========================================\033[0m\n");
  print("Total Computation Time: %f sec.\n", ((double)(clock()-start_time))/CLOCKS_PER_SEC);
  print("\n");

  //! finalize 
  //! In general, "Destroy" should be called for classes that store Petsc DMDA data (which need to be "destroyed").
  V.Destroy();
  ID.Destroy();


  //! Detroy the levelsets
  for(int ls = 0; ls<lso.size(); ls++) {
    Phi[ls]->Destroy(); delete Phi[ls];
    lso[ls]->Destroy(); delete lso[ls];
  }

  out.FinalizeOutput();
  integrator->Destroy();
  mpo.Destroy();
  spo.Destroy();
  if(grad) {
    grad->Destroy();
    delete grad;
  }
  if(interp) {
    interp->Destroy();
    delete interp;
  }

  dms.DestroyAllDataManagers();
  PetscFinalize();

  delete integrator;
  delete ff;

  for(int i=0; i<vf.size(); i++)
    delete vf[i];

  return 0;
}

//--------------------------------------------------------------
