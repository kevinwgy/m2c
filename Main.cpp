/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <time.h>
#include <petscdmda.h> //PETSc
#include <ConcurrentProgramsHandler.h>
#include <MeshGenerator.h>
#include <Output.h>
#include <VarFcnSG.h>
#include <VarFcnNASG.h>
#include <VarFcnMG.h>
#include <VarFcnMGExt.h>
#include <VarFcnTillot.h>
#include <VarFcnJWL.h>
#include <VarFcnANEOSEx1.h>
#include <VarFcnHomoIncomp.h>
#include <VarFcnDummy.h>
#include <FluxFcnGenRoe.h>
#include <FluxFcnLLF.h>
#include <FluxFcnHLLC.h>
#include <FluxFcnGodunov.h>
#include <GravityHandler.h>
#include <TimeIntegratorSemiImp.h>
#include <SpaceInitializer.h>
#include <GradientCalculatorCentral.h>
#include <IonizationOperator.h>
#include <PrescribedMotionOperator.h>
#include <SpecialToolsDriver.h>
#include <ExactRiemannSolverInterfaceJump.h>
#include <set>
#include <string>
using std::to_string;
#include <limits>

// for timing
//using std::chrono::high_resolution_clock;
//using std::chrono::duration_cast;
//using std::chrono::duration;
//using std::chrono::milliseconds;



int verbose;
double domain_diagonal;
double start_time; //start time in seconds
MPI_Comm m2c_comm;

int INACTIVE_MATERIAL_ID;

/*************************************
 * Main Function
 ************************************/
int main(int argc, char* argv[])
{
  
  //! Initialize MPI 
  MPI_Init(NULL,NULL); //called together with all concurrent programs -> MPI_COMM_WORLD
  start_time = walltime(); //for timing purpose only (calls MPI_Wtime)

  //! Print header (global proc #0, assumed to be a M2C proc)
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



  //! Initialize VarFcn (EOS, etc.) 
  std::set<int> vf_tracker;
  std::vector<VarFcnBase *> vf;
  bool incompressible = false; // whether the flow is incompressible
  int num_incompressible = 0;
  for(int i=0; i<(int)iod.eqs.materials.dataMap.size(); i++)
    vf.push_back(NULL); //allocate space for the VarFcn pointers

  for(auto it = iod.eqs.materials.dataMap.begin(); it != iod.eqs.materials.dataMap.end(); it++) {
    int matid = it->first;
    vf_tracker.insert(matid);
    if(matid < 0 || matid >= (int)vf.size()) {
      print_error("*** Error: Detected error in the specification of material indices (id = %d).\n", matid);
      exit_mpi();
    }
    if(it->second->eos == MaterialModelData::STIFFENED_GAS)
      vf[matid] = new VarFcnSG(*it->second);
    else if(it->second->eos == MaterialModelData::NOBLE_ABEL_STIFFENED_GAS)
      vf[matid] = new VarFcnNASG(*it->second);
    else if(it->second->eos == MaterialModelData::MIE_GRUNEISEN)
      vf[matid] = new VarFcnMG(*it->second);
    else if(it->second->eos == MaterialModelData::EXTENDED_MIE_GRUNEISEN)
      vf[matid] = new VarFcnMGExt(*it->second);
    else if(it->second->eos == MaterialModelData::TILLOTSON)
      vf[matid] = new VarFcnTillot(*it->second);
    else if(it->second->eos == MaterialModelData::JWL)
      vf[matid] = new VarFcnJWL(*it->second);
    else if(it->second->eos == MaterialModelData::ANEOS_BIRCH_MURNAGHAN_DEBYE)
      vf[matid] = new VarFcnANEOSEx1(*it->second);
    else if(it->second->eos == MaterialModelData::HOMOGENEOUS_INCOMPRESSIBLE) {
      vf[matid] = new VarFcnHomoIncomp(*it->second);
      num_incompressible++;
    } else {
      print_error("*** Error: Unable to initialize variable functions (VarFcn) for the "
                  "specified material model.\n");
      exit_mpi();
    }
  }
  if(vf_tracker.size() != vf.size()) {
    print_error("*** Error: Detected error in the specification of material IDs.\n");
    exit_mpi();
  }
  vf_tracker.clear();   

  if(num_incompressible>0 && num_incompressible<(int)iod.eqs.materials.dataMap.size()) {
    print_error("*** Error: Currently, M2C does not support compressible and incompressible materials "
                "in a single simulation.");
    exit_mpi();
  }
  else if (num_incompressible == (int)iod.eqs.materials.dataMap.size())
    incompressible = true;
    
  vf.push_back(new VarFcnDummy(iod.eqs.dummy_state)); //for "inactive" nodes, e.g., occluded or inside a solid body
  INACTIVE_MATERIAL_ID = vf.size() - 1;



  /********************************************************
   *                   Special Tools                      *
   *******************************************************/
  if(iod.special_tools.type != SpecialToolsData::NONE) {
    SpecialToolsDriver special_tools_driver(iod, vf, comm, concurrent);
    special_tools_driver.Run();
    concurrent.Destroy();
    MPI_Finalize();
    return 0;
  }
  /*******************************************************/



  //! Initialize Embedded Boundary Operator, if needed
  EmbeddedBoundaryOperator *embed = NULL;
  if(iod.concurrent.aeros.fsi_algo != AerosCouplingData::NONE)
    embed = new EmbeddedBoundaryOperator(comm, iod, true); 
  else if(iod.ebm.embed_surfaces.surfaces.dataMap.size() != 0)
    embed = new EmbeddedBoundaryOperator(comm, iod, false);

  //! Call Messengers' initializers (which may require different inputs)
  if(concurrent.Coupled()) {
    if(iod.concurrent.aeros.fsi_algo != AerosCouplingData::NONE)
      concurrent.InitializeMessengers(embed->GetPointerToSurface(0),
                                      embed->GetPointerToForcesOnSurface(0));
    else
      concurrent.InitializeMessengers(NULL, NULL);
  }


  // -------------------------------------------------------------------------------------------------
  //! Note: If incompresible == true, riemann and some other classes are not needed. But to avoid too
  //!       many if-statements which complicates the code, we still construct these classes.
  // -------------------------------------------------------------------------------------------------


  //! Initialize the exact Riemann problem solver.
  ExactRiemannSolverBase *riemann = NULL;
  if(iod.exact_riemann.surface_tension == ExactRiemannSolverData::NO) {
    riemann = new ExactRiemannSolverBase(vf, iod.exact_riemann);
  }
  else {// with surface tension (an experimental feature)
    riemann = new ExactRiemannSolverInterfaceJump(vf, iod.exact_riemann);
    if(incompressible || iod.ts.expl.type != ExplicitTsData::FORWARD_EULER) {
      print_error("*** Error: Currently, only the Forward Euler time integrator for compressible flows "
                  "supports surface tension.\n");
      exit_mpi();
    }
  }

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
  
  //! Setup global mesh info
  GlobalMeshInfo global_mesh(xcoords, ycoords, zcoords, dx, dy, dz, incompressible);

  //! Initialize PETSc
  PETSC_COMM_WORLD = comm;
  PetscInitialize(&argc, &argv, argc>=3 ? argv[2] : (char*)0, (char*)0);

  //! Setup PETSc data array (da) structure for nodal variables
  DataManagers3D dms(comm, xcoords.size(), ycoords.size(), zcoords.size());

  //! Let global_mesh find subdomain boundaries and neighbors
  global_mesh.FindSubdomainInfo(comm, dms);

  //! Initialize space operator
  SpaceOperator spo(comm, dms, iod, vf, *ff, *riemann, global_mesh);

  //! Track the embedded boundaries
  if(embed) {
    // determine whether force should be "spread out" to a 3D structure
    embed->SetCommAndMeshInfo(dms, spo.GetMeshCoordinates(), 
                              *(spo.GetPointerToInnerGhostNodes()), *(spo.GetPointerToOuterGhostNodes()),
                              global_mesh);
    embed->SetupIntersectors();
    embed->TrackSurfaces();
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
  spo.SetupViscosityOperator(interp, grad, embed!=NULL);

  //! Setup heat diffusion operator in spo (if heat diffusion model is not NONE)
  spo.SetupHeatDiffusionOperator(interp, grad);
 
  /** Create a set that stores additional nodes/cells where solutions should *not* be updated
    * In the vast majority of cases, this set should be empty. Use it only when other options
    * are impossible or a lot more intrusive.*/
  std::set<Int3> spo_frozen_nodes;
  spo.SetPointerToFrozenNodes(&spo_frozen_nodes); 


  //! Create turbulence working variable(s)
  SpaceVariable3D* Vturb = NULL;
  if(incompressible) {
    if(iod.rans.model == RANSTurbulenceModelData::SPALART_ALLMARAS)
      Vturb = new SpaceVariable3D(comm, &(dms.ghosted1_1dof));
  } else {
    if(iod.rans.model != RANSTurbulenceModelData::NONE) {
      print_error("*** Error: Currently, only the incompressible flow solver supports turbulence models.\n"); 
      exit_mpi();
    }
  }


  //! Initialize incompressible flow space operator
  IncompressibleOperator* inco = NULL;
  if(incompressible) {
    assert(interp);
    inco = new IncompressibleOperator(comm, dms, iod, vf, spo, *interp);
    if(Vturb)
      inco->InitializeTurbulenceVariables(*Vturb);
  }

  //! Allocate memory for V and ID 
  SpaceVariable3D V(comm, &(dms.ghosted1_5dof)); //!< primitive state variables
  SpaceVariable3D ID(comm, &(dms.ghosted1_1dof)); //!< material id

  //! Allocate memory for Phi. Initialize LevelSetOperators
  std::vector<LevelSetOperator*> lso;
  std::vector<SpaceVariable3D*>  Phi;
  std::vector<SpaceVariable3D*> NPhi; // unit normal, first derivative of Phi
  std::vector<SpaceVariable3D*> KappaPhi; // curvature, based on the second derivative of Phi

  std::set<int> ls_tracker;
  int ls_input_id_min = 9999, ls_input_id_max = -9999;
  for(auto it = iod.schemes.ls.dataMap.begin(); it != iod.schemes.ls.dataMap.end(); it++) {
    if(ls_tracker.find(it->first) != ls_tracker.end()){
      print_error("*** Error: Detected two level sets with the same id (%d).\n", it->first);
      exit_mpi();
    }
    ls_tracker.insert(it->first);
    ls_input_id_min = std::min(ls_input_id_min, it->first);
    ls_input_id_max = std::max(ls_input_id_max, it->first);
  } 
  if(ls_input_id_min<0 || ls_input_id_max>=(int)ls_tracker.size()){
    print_error("*** Error: Level set ids should start from 0 and have no gaps.\n"); 
    exit_mpi();
  }
  lso.resize(ls_tracker.size(),NULL);
  Phi.resize(ls_tracker.size(),NULL);
  NPhi.resize(ls_tracker.size(),NULL);
  KappaPhi.resize(ls_tracker.size(),NULL);
  ls_tracker.clear(); //for re-use
  for(auto it = iod.schemes.ls.dataMap.begin(); it != iod.schemes.ls.dataMap.end(); it++) {
    int matid = it->second->materialid;
    if(matid<=0 || matid>=(int)vf.size()) { //cannot use ls to track material 0
      print_error("*** Error: Cannot initialize a level set for tracking material %d.\n", matid);
      exit_mpi();
    }
    if(ls_tracker.find(matid) != ls_tracker.end()) {
      print_error("*** Error: Cannot initialize multiple level sets for the same material (id=%d).\n", matid);
      exit_mpi();
    }
    ls_tracker.insert(matid);    
    lso[it->first] = new LevelSetOperator(comm, dms, iod, *it->second, spo);
    Phi[it->first] = new SpaceVariable3D(comm, &(dms.ghosted1_1dof));
    NPhi[it->first] = new SpaceVariable3D(comm, &(dms.ghosted1_3dof));
    KappaPhi[it->first] = new SpaceVariable3D(comm, &(dms.ghosted1_1dof));
  }
  // check for user error
  for(int ls=0; ls<OutputData::MAXLS; ls++)
    if(iod.output.levelset[ls] == OutputData::ON && ls>=(int)Phi.size()) {
      print_error("*** Error: Cannot output level set %d, which is undefined.\n"); exit_mpi();}


  // Create gravity handler if specified by user
  GravityHandler* ghand = NULL;
  if(iod.ic.floodIc.source_x<DBL_MAX && iod.ic.floodIc.source_y<DBL_MAX &&
     iod.ic.floodIc.source_z<DBL_MAX)
    ghand = new GravityHandler(comm, dms, iod, spo.GetMeshCoordinates(),
                               *(spo.GetPointerToInnerGhostNodes()), *(spo.GetPointerToOuterGhostNodes()),
                               global_mesh);


#ifdef LEVELSET_TEST
  print("\n");
  if(lso.empty()) {
    print_error("*** Error: Activated LEVELSET_TEST (%d), but LevelSetOperator is not specified.\n",
                (int)LEVELSET_TEST);
    exit_mpi();
  }
  print("\033[0;32m- Testing the Level Set Solver using a prescribed velocity field (%d). "
        "N-S solver not activated.\033[0m\n", (int)LEVELSET_TEST);
#endif
  

  // ------------------------------------------------------------------------
  //! Initialize V, ID, Phi. 
  SpaceInitializer spinit(comm, dms, iod, global_mesh, spo.GetMeshCoordinates());
  std::multimap<int, std::pair<int,int> >
  id2closure = spinit.SetInitialCondition(V, ID, Phi, NPhi, KappaPhi, spo, lso, ghand,
                                          embed ? embed->GetPointerToEmbeddedBoundaryData() : nullptr,
                                          embed ? embed->GetPointerToMultiSurfaceIntersectors() : nullptr);
  
  // Boundary conditions are applied to V and Phi. But the ghost nodes of ID have not been populated.
  // ------------------------------------------------------------------------

  if(embed) //even if id2closure is empty, we must still call this function to set "inactive_elem_status"
    embed->FindSolidBodies(id2closure);  //tracks the colors of solid bodies


  //! Initialize multiphase operator (for updating "phase change")
  MultiPhaseOperator mpo(comm, dms, iod, vf, global_mesh, spo, lso);
  if((int)lso.size()>1) { //at each node, at most one "phi" can be negative
    int overlap = mpo.CheckLevelSetOverlapping(Phi);
    if(overlap>0) {
      print_error("*** Error: Found overlapping material subdomains. Number of overlapped cells "
                  "(including duplications): %d.\n", overlap);
      exit_mpi();
    }
  }
  mpo.UpdateMaterialIDAtGhostNodes(ID); //ghost nodes (outside domain) get the ID of their image nodes

  if(incompressible) {
    inco->FinalizeInitialCondition(V, ID); //Shift vel to cell faces, set rho=rho0, p=0 (must be followed by ApplyBC)
    inco->ApplyBoundaryConditions(V);
  }

  //! Initialize laser radiation solver (if needed)
  LaserAbsorptionSolver* laser = NULL;
  SpaceVariable3D* L = NULL;  //laser radiance
  if(iod.laser.source_power>0.0 || iod.laser.source_intensity>0.0 ||
     strcmp(iod.laser.source_power_timehistory_file, "") != 0) {//laser source is specified

    if(iod.laser.parallel == LaserData::BALANCED) {//re-balance the load
      print("- Initializing the laser radiation solver on a re-partitioned sub-mesh.\n");
      laser = new LaserAbsorptionSolver(comm, dms, iod, vf, spo.GetMeshCoordinates(), 
                                        //the following inputs are used for creating a new dms/spo
                                        *ff, *riemann, xcoords, ycoords, zcoords, dx, dy, dz);
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


  //! Initialize hyperelasticity solver and reference map (if needed)
  HyperelasticityOperator* heo = NULL;
  SpaceVariable3D* Xi = NULL; //reference map
  bool activate_heo = false;
  for(auto&& material : iod.eqs.materials.dataMap)
    if(material.second->hyperelasticity.type != HyperelasticityModelData::NONE) {
      activate_heo = true;
      break;
    }
  if(activate_heo) {
    heo = new HyperelasticityOperator(comm, dms, iod, vf, spo.GetMeshCoordinates(),
                                      spo.GetMeshDeltaXYZ(), spo.GetMeshCellVolumes(),
                                      global_mesh, *interp, *grad,
                                      *(spo.GetPointerToInnerGhostNodes()),
                                      *(spo.GetPointerToOuterGhostNodes()));
    spo.SetHyperelasticityOperatorPointer(heo);
    Xi = new SpaceVariable3D(comm, &(dms.ghosted1_3dof));
    heo->InitializeReferenceMap(*Xi);
  }


#ifdef HYPERELASTICITY_TEST 
  print("\n");
  if(!heo || !Xi) {
    print_error("*** Error: Activated HYPERELASTICITY_TEST (%d), but HyperelasticityOperator and/or "
                "ReferenceMapOperator are not defined.\n", (int)HYPERELASTICITY_TEST);
    exit_mpi();
  }
  print("\033[0;32m- Testing the Hyperelasticity Solver using a prescribed velocity field (%d).\033[0m\n",
        (int)HYPERELASTICITY_TEST);
#endif


  //! Create prescribed motion operator (if needed)
  PrescribedMotionOperator* pmo = NULL;
  if(!iod.schemes.pm.dataMap.empty())
    pmo = new PrescribedMotionOperator(iod.schemes.pm, vf.size()); //vf.size includes INACTIVE

/*
  ID.StoreMeshCoordinates(spo.GetMeshCoordinates());
  V.StoreMeshCoordinates(spo.GetMeshCoordinates());
  ID.WriteToVTRFile("ID.vtr","id");
  V.WriteToVTRFile("V.vtr", "sol");
  if(Phi.size()>0) {
    Phi[0]->StoreMeshCoordinates(spo.GetMeshCoordinates());
    Phi[0]->WriteToVTRFile("Phi0.vtr", "phi0");
  }
  print("I am here!\n");
  exit_mpi();
*/

  //! Initialize output
  Output out(comm, dms, iod, global_mesh, spo.GetPointerToOuterGhostNodes(), vf, laser, spo.GetMeshCoordinates(),
             spo.GetMeshDeltaXYZ(), spo.GetMeshCellVolumes(), ion, heo, inco); 
  out.InitializeOutput(spo.GetMeshCoordinates());
  mpi_barrier();


  //! Initialize time integrator
  TimeIntegratorBase *integrator = NULL;
  if(!incompressible) { //compressible
    if(iod.ts.type != TsData::EXPLICIT) {
      print_error("*** Error: Compressible flows require an explicit time-integrator.\n");
      exit_mpi();
    }
    if(iod.ts.expl.type == ExplicitTsData::FORWARD_EULER)
      integrator = new TimeIntegratorFE(comm, iod, dms, spo, lso, mpo, laser, embed, heo, pmo);
    else if(iod.ts.expl.type == ExplicitTsData::RUNGE_KUTTA_2)
      integrator = new TimeIntegratorRK2(comm, iod, dms, spo, lso, mpo, laser, embed, heo, pmo);
    else if(iod.ts.expl.type == ExplicitTsData::RUNGE_KUTTA_3)
      integrator = new TimeIntegratorRK3(comm, iod, dms, spo, lso, mpo, laser, embed, heo, pmo);
    else {
      print_error("*** Error: Unable to initialize time integrator for the specified (explicit) method.\n");
      exit_mpi();
    }
  } 
  else { //incompressible
    if(iod.ts.type != TsData::SEMI_IMPLICIT) {
      print_error("*** Error: Incompressible flows require a semi-implicit time-integrator.\n");
      exit_mpi();
    }
    if(iod.ts.semi_impl.type == SemiImplicitTsData::SIMPLE)
      integrator = new TimeIntegratorSIMPLE(comm, iod, dms, spo, *inco, lso, mpo, laser, embed, heo, pmo);
    else if(iod.ts.semi_impl.type == SemiImplicitTsData::SIMPLER)
      integrator = new TimeIntegratorSIMPLER(comm, iod, dms, spo, *inco, lso, mpo, laser, embed, heo, pmo);
    else if(iod.ts.semi_impl.type == SemiImplicitTsData::SIMPLEC)
      integrator = new TimeIntegratorSIMPLEC(comm, iod, dms, spo, *inco, lso, mpo, laser, embed, heo, pmo);
    else if(iod.ts.semi_impl.type == SemiImplicitTsData::SIMPLEC)
      integrator = new TimeIntegratorPISO(comm, iod, dms, spo, *inco, lso, mpo, laser, embed, heo, pmo);
    else {
      print_error("*** Error: Unable to initialize time integrator for the specified (semi-implicit)"
                  " method.\n");
      exit_mpi();
    }
  }


  //! Setup for steady-state computations
  SpaceVariable3D *LocalDt = NULL;
  bool steady_state = iod.ts.convergence_tolerance>0.0;
  if(iod.ts.local_dt == TsData::YES) {//local time-stepping
    if(!steady_state) {
      print_error("*** Error: Local time-stepping can be used only for steady-state computations.\n");
      exit_mpi();
    }
    if(iod.ts.timestep > 0.0) {
      print_error("*** Error: Local time-stepping conflicts with a constant user-specified time step.\n");
      exit_mpi();
    }
    LocalDt = new SpaceVariable3D(comm, &(dms.ghosted1_1dof)); //!< local dt for steady-state
  }
  if(steady_state) { //steady-state analysis
    if(concurrent.Coupled()) {
      print_error("*** Error: Unable to perform steady-state analysis with concurrent programs.\n");
      exit_mpi();
    }
    if(lso.size()>0) {
      print_error("*** Error: Unable to perform steady-state analysis with level set(s).\n");
      exit_mpi();
    }

    print("- Performing a steady-state analysis. Tolerance (rel. residual): %e.\n\n",
          iod.ts.convergence_tolerance);
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

  // In the case of steady-state simulation with local time-stepping, the constant "dt" is not
  // actually used. It represents the smallest dt in all the cells.

  if(laser) //initialize L (otherwise the initial output will only have 0s)
    laser->ComputeLaserRadiance(V, ID, *L, t, time_step);

  print("I got here.\n");
  //! Compute force on embedded surfaces (if any) using initial state
  if(embed) {
    embed->ComputeForces(V, ID);
    embed->UpdateSurfacesPrevAndFPrev();

    embed->OutputSurfaces(); //!< write the mesh(es) to file
    embed->OutputResults(t, dt, time_step, true/*force_write*/); //!< write displacement and nodal loads to file
  }

  print("I got here 2.\n");
  //! write initial condition to file
  out.OutputSolutions(t, dt, time_step, V, ID, Phi, NPhi, KappaPhi, L, Xi, Vturb, true/*force_write*/);

  print("I got here 3.\n");
  if(concurrent.Coupled()) {
    concurrent.CommunicateBeforeTimeStepping(&spo.GetMeshCoordinates(), &dms,
                                             spo.GetPointerToInnerGhostNodes(),
                                             spo.GetPointerToOuterGhostNodes(),
                                             &global_mesh, &V,
                                             &ID, //updated in inactive regions (overset grids)
                                             &spo_frozen_nodes); //updated w/ "green boxes" (overset grids)
    spo.UpdateOversetGhostNodes(V);
    
/*
    SpaceVariable3D Test(comm, &(dms.ghosted1_1dof));
    double*** test = Test.GetDataPointer();
    SpaceVariable3D& coordinates(spo.GetMeshCoordinates());
    Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
    int ii0,jj0,kk0,iimax,jjmax,kkmax;
    coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax); 
    double PI = 2.0*acos(0.0);
    for(int k=kk0; k<kkmax; k++)
      for(int j=jj0; j<jjmax; j++)
        for(int i=ii0; i<iimax; i++) {
          if(concurrent.GetTwinningStatus() == ConcurrentProgramsHandler::LEADER)
            test[k][j][i] = cos(coords[k][j][i][0]*10.0*PI)*sin(coords[k][j][i][1]*5.0*PI);
          else if(fabs(coords[k][j][i][0])>=2.0 || coords[k][j][i][1]>=2.0)
            test[k][j][i] = cos(coords[k][j][i][0]*10.0*PI)*sin(coords[k][j][i][1]*5.0*PI);
        }
    coordinates.RestoreDataPointerToLocalVector();
    Test.RestoreDataPointerAndInsert();
    concurrent.FirstExchange(&Test, 1.0e-8, 1.0, false);
    Test.StoreMeshCoordinates(coordinates);
    if(concurrent.GetTwinningStatus() == ConcurrentProgramsHandler::LEADER)
      Test.WriteToVTRFile("Solution1.vtr","Test");
    else
      Test.WriteToVTRFile("Solution2.vtr","Test");


    MPI_Barrier(MPI_COMM_WORLD); exit(-1);
*/
  }

  if(embed) {
    print("I got here 4.\n");
    embed->ApplyUserDefinedSurfaceDynamics(t, dt); //update surfaces provided through input (not concurrent solver)
    print("I got here 5.\n");
    embed->TrackUpdatedSurfaces();
    print("I got here 6.\n");
    mpi_barrier();
    exit_mpi(); 
    int boundary_swept = mpo.UpdateCellsSweptByEmbeddedSurfaces(V, ID, Phi,
                                 embed->GetPointerToEmbeddedBoundaryData(),
                                 embed->GetPointerToIntersectors()); //update V, ID, Phi
    spo.ClipDensityAndPressure(V,ID);
    if(boundary_swept) {
      if(incompressible)
        inco->ApplyBoundaryConditions(V);
      else
        spo.ApplyBoundaryConditions(V);

      for(int i=0; i<(int)Phi.size(); i++) 
        lso[i]->ApplyBoundaryConditions(*Phi[i]);
    }
  }


  // Enforce user-prescribed velocity
  if(pmo) {
    pmo->UpdateVelocity(V,ID,t);
    spo.ApplyBoundaryConditions(V);
  }


  // find maxTime, and dts (meaningful only with concurrent programs)
  double tmax = iod.ts.maxTime;
  double dts = 0.0;
  if(concurrent.Coupled()) {
    dts =  concurrent.GetTimeStepSize(); //should return negative if not specified
    double tmaxs = concurrent.GetMaxTime(); //should return negative if not specified
    if(tmaxs>0) //if "concurrent" returns a tmax, use it.
      tmax = tmaxs;
  }

  // set max time-step number to user input, or INF if it is a follower in Chimera
  int maxIts = concurrent.GetTwinningStatus() == ConcurrentProgramsHandler::FOLLOWER ? INT_MAX 
             : iod.ts.maxIts;

  // Time-Stepping
  while(t<tmax-1.0e-6*dts && time_step<maxIts && !integrator->Converged()) {// the last one is for steady-state

    double dtleft = dts;

    time_step++;
    int subcycle = 0;

    do { //subcycling w.r.t. concurrent programs

      // Compute time step size
      if(!incompressible)
        spo.ComputeTimeStepSize(V, ID, dt, cfl, LocalDt); 
      else 
        inco->ComputeTimeStepSize(V, ID, dt, cfl, LocalDt);

      // Modify dt and cfl, dtleft if needed
      if(concurrent.GetTimeStepSize()>0) {//concurrent solver provides a "dts"
        if(dt>dtleft) { //dt should be reduced to dtleft
          cfl *= dtleft/dt;
          dt = dtleft;
        }
      } else
        dtleft = dt;

      if(t+dt > tmax) { //update dt at the LAST time step so it terminates at tmax
        cfl *= (tmax - t)/dt;
        dt = tmax - t;
        dtleft = dt;
      }

      // If concurrent programs do not specify dt, set dts = dt
      if(concurrent.GetTimeStepSize()<=0.0)
        dts = dt;

/*
      if(concurrent.GetTwinningStatus() == ConcurrentProgramsHandler::LEADER)
        fprintf(stdout,"[Leader] dts = %e, dt = %e, dtleft = %e.\n", dts, dt, dtleft);
      else
        fprintf(stdout,"[Follower] dts = %e, dt = %e, dtleft = %e.\n", dts, dt, dtleft);
*/
 
      if(steady_state)  {
        print("Step %d: t = %e, dt = %e, cfl = %.4e. Computation time: %.4e s.\n", 
              time_step, t, dt, cfl, walltime()-start_time);
      }
      else { //unsteady
        if(dts<=dt)
          print("Step %d: t = %e, dt = %e, cfl = %.4e. Computation time: %.4e s.\n", 
                time_step, t, dt, cfl, walltime()-start_time);
        else
          print("Step %d(%d): t = %e, dt = %e, cfl = %.4e. Computation time: %.4e s.\n", 
                time_step, subcycle+1, t, dt, cfl, walltime()-start_time);
      }

      //----------------------------------------------------
      // Move forward by one time-step: Update V, Phi, and ID
      //----------------------------------------------------
      t      += dt;
      dtleft -= dt;
      integrator->AdvanceOneTimeStep(V, ID, Phi, NPhi, KappaPhi, L, Xi, Vturb, LocalDt, t, dt, time_step,
                                     subcycle, dts); 
      subcycle++; //do this *after* AdvanceOneTimeStep.
      //----------------------------------------------------

      if(steady_state)
        print("  => Residual: %.4e (func 1 norm) |  %.4e (2 norm)  |  %.4e (max).\n",
              integrator->GetRelativeResidual1Norm(), integrator->GetRelativeResidual2Norm(),
              integrator->GetRelativeResidualInfNorm());

    } while (concurrent.Coupled() && dtleft != 0.0);


    if(embed) {
      embed->ComputeForces(V, ID);
      embed->UpdateSurfacesPrevAndFPrev();

      embed->OutputResults(t, dts, time_step, false/*force_write*/); //!< write displacement and nodal loads to file
    }


    //Exchange data with concurrent programs (Note: This chunk should be at the end of each time-step.)
    double dts0 = dts; //the previous time step size
    if(concurrent.Coupled()) {

      if(t<tmax && time_step<maxIts) {//not the last time-step
        if(time_step==1)
          concurrent.FirstExchange(&V, dts, tmax);
        else
          concurrent.Exchange(&V, dts, tmax);
      } 

      double dts_tmp = concurrent.GetTimeStepSize();
      if(dts_tmp>0)
        dts = dts_tmp;

      double tmax_tmp = concurrent.GetMaxTime(); //at final time-step, tmax is set to a very small number
      if(tmax_tmp>0)
        tmax = tmax_tmp;
    }

    if(embed) {
      embed->ApplyUserDefinedSurfaceDynamics(t, dts0); //update surfaces provided through input (not concurrent solver)
      embed->TrackUpdatedSurfaces();
      int boundary_swept = mpo.UpdateCellsSweptByEmbeddedSurfaces(V, ID, Phi,
                                   embed->GetPointerToEmbeddedBoundaryData(),
                                   embed->GetPointerToIntersectors()); //update V, ID, Phi
      spo.ClipDensityAndPressure(V,ID);
      if(boundary_swept) {
        spo.ApplyBoundaryConditions(V);
        for(int i=0; i<(int)Phi.size(); i++) 
          lso[i]->ApplyBoundaryConditions(*Phi[i]);
      }
    }

    

    out.OutputSolutions(t, dts0, time_step, V, ID, Phi, NPhi, KappaPhi, L, Xi, Vturb, false/*force_write*/);

  }

  if(concurrent.Coupled())
    concurrent.FinalExchange(&V);


  // Final outputs

  if(embed)
    embed->OutputResults(t, dt, time_step, true/*force_write*/);

  out.OutputSolutions(t, dts, time_step, V, ID, Phi, NPhi, KappaPhi, L, Xi, Vturb, true/*force_write*/);

  print("\n");
  print("\033[0;32m==========================================\033[0m\n");
  print("\033[0;32m   NORMAL TERMINATION (t = %e)  \033[0m\n", t); 
  print("\033[0;32m==========================================\033[0m\n");
  print("Total Computation Time: %f sec.\n", walltime()-start_time);
  print("\n");



  //! finalize 
  //! In general, "Destroy" should be called for classes that store Petsc DMDA data (which need to be "destroyed").
  
  concurrent.Destroy();

  V.Destroy();
  ID.Destroy();

  spinit.Destroy();

  if(LocalDt) {LocalDt->Destroy(); delete LocalDt;}

  if(laser) {laser->Destroy(); delete laser;}
  if(L)     {L->Destroy(); delete L;}

  if(ion) {ion->Destroy(); delete ion;}

  if(embed) {embed->Destroy(); delete embed;}

  //! Detroy the levelsets
  for(int ls = 0; ls<(int)lso.size(); ls++) {
    Phi[ls]->Destroy(); delete Phi[ls];
    lso[ls]->Destroy(); delete lso[ls];
    NPhi[ls]->Destroy(); delete NPhi[ls];
    KappaPhi[ls]->Destroy(); delete KappaPhi[ls];
  }

  if(heo) {heo->Destroy(); delete heo;}
  if(Xi) {Xi->Destroy(); delete Xi;}

  if(pmo) {pmo->Destroy(); delete pmo;}

  if(ghand) {ghand->Destroy(); delete ghand;}

  if(Vturb) {Vturb->Destroy(); delete Vturb;}

  out.FinalizeOutput();
  integrator->Destroy();
  mpo.Destroy();
  spo.Destroy();

  if(inco) {inco->Destroy(); delete inco;}

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
  delete riemann;

  for(int i=0; i<(int)vf.size(); i++)
    delete vf[i];

  PetscFinalize();
  MPI_Finalize();

  return 0;
}

//--------------------------------------------------------------

