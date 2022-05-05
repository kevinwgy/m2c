#include <time.h>
#include <petscdmda.h> //PETSc
#include <ConcurrentProgramsHandler.h>
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
#include <LaserAbsorptionSolver.h>
#include <IonizationOperator.h>
#include <TriangulatedSurface.h>
#include <EmbeddedBoundaryOperator.h>
#include <SpecialToolsDriver.h>
#include <set>
using std::cout;
using std::endl;
int verbose;
double domain_diagonal;
clock_t start_time;
MPI_Comm m2c_comm;

int INACTIVE_MATERIAL_ID;

/*************************************
 * Main Function
 ************************************/
int main(int argc, char* argv[])
{
  start_time = clock(); //for timing purpose only

  //! Initialize MPI 
  MPI_Init(NULL,NULL); //called together with all concurrent programs -> MPI_COMM_WORLD

  //! Print header (global proc #1, assumed to be a M2C proc)
  m2c_comm = MPI_COMM_WORLD; //temporary, just for the next few lines of code
  printHeader(argc, argv);

  //! Read user's input file (read the parameters)
  IoData iod(argc, argv);
  verbose = iod.output.verbose;

  //! Partition MPI, if there are concurrent programs
  MPI_Comm comm; //this is going to be the M2C communicator
  ConcurrentProgramsHandler concurrent(iod, MPI_COMM_WORLD, comm);
  m2c_comm = comm; //correct it
 
  //! Finalize IoData (read additional files and check for errors)
  iod.finalize();


  /********************************************************
   *                   Special Tools                      *
   *******************************************************/
  if(iod.special_tools.type != SpecialToolsData::NONE) {
    SpecialToolsDriver special_tools_driver(iod, comm, concurrent);
    special_tools_driver.Run();
    return 0;
  }
  /*******************************************************/


  //! Initialize Embedded Boundary Operator, if needed
  EmbeddedBoundaryOperator *embed = NULL;
  if(iod.concurrent.aeros.fsi_algo != AerosCouplingData::NONE) {
    embed = new EmbeddedBoundaryOperator(comm, iod, true); 
    concurrent.InitializeMessengers(embed->GetPointerToSurface(0),
                                    embed->GetPointerToForcesOnSurface(0));
  } else if(iod.ebm.embed_surfaces.surfaces.dataMap.size() != 0)
    embed = new EmbeddedBoundaryOperator(comm, iod, false);


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

  //! Calculate mesh coordinates
  vector<double> xcoords, dx, ycoords, dy, zcoords, dz;
  MeshGenerator meshgen;
  meshgen.ComputeMeshCoordinatesAndDeltas(iod.mesh, xcoords, ycoords, zcoords, dx, dy, dz);
  domain_diagonal = sqrt(pow(iod.mesh.xmax - iod.mesh.x0, 2) +
                         pow(iod.mesh.ymax - iod.mesh.y0, 2) +
                         pow(iod.mesh.zmax - iod.mesh.z0, 2));
  
  //! Initialize PETSc
  PETSC_COMM_WORLD = comm;
  PetscInitialize(&argc, &argv, argc>=3 ? argv[2] : (char*)0, (char*)0);

  //! Setup PETSc data array (da) structure for nodal variables
  DataManagers3D dms(comm, xcoords.size(), ycoords.size(), zcoords.size());

  //! Initialize space operator
  SpaceOperator spo(comm, dms, iod, vf, *ff, riemann, xcoords, ycoords, zcoords, dx, dy, dz);

  //! Track the embedded boundaries
  if(embed) {
    embed->SetCommAndMeshInfo(dms, spo.GetMeshCoordinates(), 
                              *(spo.GetPointerToInnerGhostNodes()), *(spo.GetPointerToOuterGhostNodes()),
                              xcoords, ycoords, zcoords, dx, dy, dz);
    embed->SetupIntersectors();
    embed->TrackSurfaces();
    return 0;
  }

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

  //! Setup heat diffusion operator in spo (if heat diffusion model is not NONE)
  spo.SetupHeatDiffusionOperator(interp, grad);
 
  //! Initialize State Variables
  SpaceVariable3D V(comm, &(dms.ghosted1_5dof)); //!< primitive state variables
  SpaceVariable3D ID(comm, &(dms.ghosted1_1dof)); //!< material id
  INACTIVE_MATERIAL_ID = -1;

  //! Impose initial condition
  std::map<int, std::pair<int,int> >
  id2closure = spo.SetInitialCondition(V, ID, embed ? embed->GetPointerToIntersectorResults() : nullptr);

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

    if(id2closure.find(matid) != id2closure.end()) {
      lso.back()->SetInitialCondition(*Phi.back(), 
                                      embed->GetPointerToIntersectoResultsOnSurface(id2closure[matid].first),
                                      id2closure[matid].second);
    } else
      lso.back()->SetInitialCondition(*Phi.back());

    print("- Initialized level set function (%d) for tracking the boundary of material %d.\n\n", 
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


  //! Initialize laser radiation solver (if needed)
  LaserAbsorptionSolver* laser = NULL;
  SpaceVariable3D* L = NULL;  //laser radiance
  if(iod.laser.source_power>0.0 || iod.laser.source_intensity>0.0 ||
     strcmp(iod.laser.source_power_timehistory_file, "") != 0) {//laser source is specified

    if(iod.laser.parallel == LaserData::BALANCED) {//re-balance the load
      print("- Initializing the laser radiation solver on a re-partitioned sub-mesh.\n");
      laser = new LaserAbsorptionSolver(comm, dms, iod, vf, spo.GetMeshCoordinates(), 
                                        //the following inputs are used for creating a new dms/spo
                                        *ff, riemann, xcoords, ycoords, zcoords, dx, dy, dz);
    } else
      laser = new LaserAbsorptionSolver(comm, dms, iod, vf, spo.GetMeshCoordinates(), 
                                        spo.GetMeshDeltaXYZ(), spo.GetMeshCellVolumes(),
                                        *(spo.GetPointerToInnerGhostNodes()),
                                        *(spo.GetPointerToOuterGhostNodes()));

    L = new SpaceVariable3D(comm, &(dms.ghosted1_1dof)); 
  }
 

  //! Initialize ionization solver (if needed)
  IonizationOperator* ion = NULL;
  if(iod.ion.materialMap.dataMap.size() != 0)
    ion = new IonizationOperator(comm, dms, iod, vf);


  //! Initialize output
  Output out(comm, dms, iod, vf, spo.GetMeshCellVolumes(), ion); 
  out.InitializeOutput(spo.GetMeshCoordinates());

  //! Initialize time integrator
  TimeIntegratorBase *integrator = NULL;
  if(iod.ts.type == TsData::EXPLICIT) {
    if(iod.ts.expl.type == ExplicitData::FORWARD_EULER)
      integrator = new TimeIntegratorFE(comm, iod, dms, spo, lso, mpo, laser);
    else if(iod.ts.expl.type == ExplicitData::RUNGE_KUTTA_2)
      integrator = new TimeIntegratorRK2(comm, iod, dms, spo, lso, mpo, laser);
    else if(iod.ts.expl.type == ExplicitData::RUNGE_KUTTA_3)
      integrator = new TimeIntegratorRK3(comm, iod, dms, spo, lso, mpo, laser);
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

  if(laser) //initialize L (otherwise the initial output will only have 0s)
    laser->ComputeLaserRadiance(V, ID, *L, t);

  //! Compute force on embedded surfaces (if any) using initial state
  if(embed) 
    embed->ComputeForces(V, ID);

  //! write initial condition to file
  out.OutputSolutions(t, dt, time_step, V, ID, Phi, L, true/*force_write*/);

  if(concurrent.Coupled()) {
    concurrent.CommunicateBeforeTimeStepping(); 
    if(embed)
      embed->TrackUpdatedSurfaceFromOtherSolver();
  }

  // find maxTime, and dts (meaningful only with concurrent programs)
  double tmax = iod.ts.maxTime;
  double dts = 0.0;
  if(concurrent.Coupled()) {
    dts =  concurrent.GetTimeStepSize();
    tmax = concurrent.GetMaxTime(); //std::max(iod.ts.maxTime, concurrent.GetMaxTime());
  }

  while(t<tmax && time_step<iod.ts.maxIts) {

    double dtleft = dts;

    time_step++;
    int subcycle = 0;

    do { //subcycling w.r.t. concurrent programs

      // Compute time step size
      spo.ComputeTimeStepSize(V, ID, dt, cfl); 

      if(concurrent.Coupled() && dt>dtleft) { //reduce dt to match the structure
        cfl *= dtleft/dt;
        dt = dtleft;
      }

      if(t+dt >= tmax) { //update dt at the LAST time step so it terminates at tmax
        cfl *= (tmax - t)/dt;
        dt = tmax - t;
      }

      if(!concurrent.Coupled())
        dts = dt;
 
      if(dts<=dt)
        print("Step %d: t = %e, dt = %e, cfl = %.4e. Computation time: %.4e s.\n", time_step, t, dt, cfl, 
              ((double)(clock()-start_time))/CLOCKS_PER_SEC);
      else
        print("Step %d(%d): t = %e, dt = %e, cfl = %.4e. Computation time: %.4e s.\n", time_step, subcycle+1, t, dt, cfl, 
              ((double)(clock()-start_time))/CLOCKS_PER_SEC);

      //----------------------------------------------------
      // Move forward by one time-step: Update V, Phi, and ID
      //----------------------------------------------------
      t      += dt;
      dtleft -= dt;
      integrator->AdvanceOneTimeStep(V, ID, Phi, L, t, dt, time_step, subcycle, dts); 
      subcycle++; //do this *after* AdvanceOneTimeStep.
      //----------------------------------------------------

    } while (concurrent.Coupled() && dtleft != 0.0);

    if(embed) 
      embed->ComputeForces(V, ID);

    //Exchange data with concurrent programs (Note: This chunk should be at the end of each time-step.)
    if(concurrent.Coupled()) {
      if(t<tmax && time_step<iod.ts.maxIts) {//not the last time-step
        if(time_step==1)
          concurrent.FirstExchange();
        else
          concurrent.Exchange();
      } 

      dts =  concurrent.GetTimeStepSize();
      tmax = concurrent.GetMaxTime(); //at final time-step, tmax is set to a very small number

      if(embed)
        embed->TrackUpdatedSurfaceFromOtherSolver();
    }

    out.OutputSolutions(t, dts, time_step, V, ID, Phi, L, false/*force_write*/);

  }

  if(concurrent.Coupled())
    concurrent.FinalExchange();

  out.OutputSolutions(t, dts, time_step, V, ID, Phi, L, true/*force_write*/);

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

  if(laser) {laser->Destroy(); delete laser;}
  if(L)     {L->Destroy(); delete L;}

  if(ion) {ion->Destroy(); delete ion;}

  if(embed) {embed->Destroy(); delete embed;}

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

  delete integrator;
  delete ff;

  for(int i=0; i<vf.size(); i++)
    delete vf[i];

  PetscFinalize();
  MPI_Finalize();

  return 0;
}

//--------------------------------------------------------------

