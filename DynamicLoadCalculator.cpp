/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<DynamicLoadCalculator.h>
#include<EmbeddedBoundaryOperator.h>
#include<rbf_interp.hpp>
#include<cfloat> //DBL_MAX

using std::string;
using std::vector;
using std::shared_ptr;

extern double start_time; //time computation started (in Main.cpp)

//-----------------------------------------------------------------

DynamicLoadCalculator::DynamicLoadCalculator(IoData &iod_, MPI_Comm &comm_, 
                                             ConcurrentProgramsHandler &concurrent_)
                     : comm(comm_), iod(iod_), concurrent(concurrent_), 
                       lagout(comm_, iod_.special_tools.transient_input.output),
                       S0(NULL), S1(NULL), tree0(NULL), tree1(NULL), id0(-INT_MAX), id1(-INT_MAX),
                       force_calculator(std::make_tuple(nullptr,nullptr,nullptr))

{ }

//-----------------------------------------------------------------

DynamicLoadCalculator::~DynamicLoadCalculator()
{ 
  //the smart pointers ("shared_ptr") are deleted automatically

  //delete user-specified force calculator
  if(std::get<0>(force_calculator)) {
    assert(std::get<2>(force_calculator)); //this is the destruction function
    (std::get<2>(force_calculator))(std::get<0>(force_calculator)); //use the destruction function to destroy the calculator
    assert(std::get<1>(force_calculator));
    dlclose(std::get<1>(force_calculator));
  }
}

//-----------------------------------------------------------------

void
DynamicLoadCalculator::Run()
{

  prefix = string(iod.special_tools.transient_input.snapshot_file_prefix);
  suffix = string(iod.special_tools.transient_input.snapshot_file_suffix);

  if(strcmp(iod.special_tools.transient_input.metafile, "") != 0)
    ReadMetaFile(string(iod.special_tools.transient_input.metafile));

  if(iod.concurrent.aeros.fsi_algo != AerosCouplingData::NONE)
    RunForAeroS();
  else {
    print_error("*** Error: No concurrent program specified.\n");
    exit_mpi();
  }

  print("\n");
  print("\033[0;32m==========================================\033[0m\n");
  print("\033[0;32m            NORMAL TERMINATION            \033[0m\n");
  print("\033[0;32m==========================================\033[0m\n");
  print("Total Computation Time: %f sec.\n", walltime() - start_time);
  print("\n");
}

//-----------------------------------------------------------------

void
DynamicLoadCalculator::RunForAeroS()
{

  // Initialize Embedded Boundary Operator
  EmbeddedBoundaryOperator *embed = NULL;
  embed = new EmbeddedBoundaryOperator(comm, iod, true);
  concurrent.InitializeMessengers(embed->GetPointerToSurface(0),
                                  embed->GetPointerToForcesOnSurface(0));
  
  // Get pointer to embedded surface and force vector
  TriangulatedSurface *surface = embed->GetPointerToSurface(0);
  vector<Vec3D>       *force   = embed->GetPointerToForcesOnSurface(0);

  // Allocate memory for the (internal) force vector
  F0.assign(surface->X.size(), Vec3D(0.0));
  F1.assign(surface->X.size(), Vec3D(0.0)); //note that X may contain unused nodes (e.g., in fracture)

  // Setup user-defined force calculator (if provided)
  SetupUserDefinedForces();

  // ---------------------------------------------------
  // Main time loop (mimics the time loop in Main.cpp
  // ---------------------------------------------------
  double t = 0.0, dt = 0.0, tmax = 0.0;
  int time_step = 0; 
  ComputeForces(surface, force, t);

  lagout.OutputTriangulatedMesh(surface->X0, surface->elems);
  lagout.OutputResults(t, dt, time_step, surface->X0, surface->X, *force, NULL, true);

  concurrent.CommunicateBeforeTimeStepping(); 
  dt   = concurrent.GetTimeStepSize();
  tmax = concurrent.GetMaxTime();

  while(t<tmax) {

    time_step++;

    if(t+dt >= tmax)
      dt = tmax - t;

    //---------------------------------
    // Move forward by one time-step
    //---------------------------------
    t += dt;
    print("Step %d: t = %e, dt = %e. Computation time: %.4e s.\n", time_step, t, dt, walltime()-start_time);
 
    ComputeForces(surface, force, t); 

    if(t<tmax) {
      if(time_step==1)
        concurrent.FirstExchange();
      else
        concurrent.Exchange();
    } 

    dt   = concurrent.GetTimeStepSize();
    tmax = concurrent.GetMaxTime(); //set to a small number at final time-step

    lagout.OutputResults(t, dt, time_step, surface->X0, surface->X, *force, NULL, false);
  }

  concurrent.FinalExchange();

  lagout.OutputResults(t, dt, time_step, surface->X0, surface->X, *force, NULL, true);

  if(embed)
    delete embed;
}

//-----------------------------------------------------------------

void
DynamicLoadCalculator::SetupUserDefinedForces()
{
  
  for(auto it = iod.ebm.embed_surfaces.surfaces.dataMap.begin();
           it != iod.ebm.embed_surfaces.surfaces.dataMap.end(); it++) {
    int index = it->first;
    if(index != 0)
      continue; // Externally provided embedded surface must have id = 0

    if(strcmp(it->second->force_calculator, "") == 0)
      return; // user did not specify a force calculator

    // setup the calculator: dynamically load the object
    std::get<1>(force_calculator) = dlopen(it->second->force_calculator, RTLD_NOW);
    if(!std::get<1>(force_calculator)) {
      print_error("*** Error: Unable to load object %s.\n", it->second->force_calculator);
      exit_mpi();
    }
    dlerror();
    CreateUDF* create = (CreateUDF*) dlsym(std::get<1>(force_calculator), "Create");
    const char* dlsym_error = dlerror();
    if(dlsym_error) {
      print_error("*** Error: Unable to find function Create in %s.\n",
                  it->second->force_calculator);
      exit_mpi();
    }
    std::get<2>(force_calculator) = (DestroyUDF*) dlsym(std::get<1>(force_calculator), "Destroy");
    dlsym_error = dlerror();
    if(dlsym_error) {
      print_error("*** Error: Unable to find function Destroy in %s.\n",
                  it->second->force_calculator);
      exit_mpi();
    }

    //This is the actual calculator
    std::get<0>(force_calculator) = create();

    print("- Loaded user-defined force calculator for surface %d from %s.\n",
          index, it->second->force_calculator);

    return; //done!
  } 

}

//-----------------------------------------------------------------

void
DynamicLoadCalculator::ReadMetaFile(string filename)
{

  if(!strcmp(filename.c_str(), "")) {//file not specified
    print_error("*** Error: Metafile is not specified.\n");
    exit_mpi();
  }

  std::fstream input;
  input.open(filename.c_str(), std::fstream::in);
  if (!input.is_open()) {
    print_error("*** Error: Cannot open file %s.\n", filename.c_str());
    exit_mpi();
  } else
    print("- Reading user-specified metafile: %s.\n", filename.c_str());

  std::string word, line;

  //Line #1 --- user's comment (skip)
  input.ignore(512,'\n'); //done with line 1 (skip, user's comment)

  //Line #2 --- fields in data (e.g., ## Coordinates Density Velocity Pressure)
  input.ignore(2,' '); //This line should start with ##
  getline(input, line);

  std::istringstream iss(line);
  int column = 0;
  for(;iss>>word;) {
    // just check the first four letters
    if(!(word.compare(0,4,"Node",0,4) &&
         word.compare(0,4,"NODE",0,4) &&
         word.compare(0,4,"node",0,4) &&
         word.compare(0,4,"Index",0,4) &&
         word.compare(0,4,"INDEX",0,4) &&
         word.compare(0,4,"index",0,4))) {
      column2var[column]   = NODE_NUMBER;
      var2column[NODE_NUMBER] = column; 
      column += 1;
    } 
    else if(!(word.compare(0,4,"Coordinates",0,4) &&
              word.compare(0,4,"COORDINATES",0,4) &&
              word.compare(0,4,"coordinates",0,4))) {
      column2var[column]   = COORDINATES;
      column2var[column+1] = COORDINATES;
      column2var[column+2] = COORDINATES;
      var2column[COORDINATES] = column; 
      column += 3;
    } 
    else if(!(word.compare(0,4,"Density",0,4) &&
              word.compare(0,4,"DENSITY",0,4) &&
              word.compare(0,4,"density",0,4))) {
      column2var[column]   = DENSITY;
      var2column[DENSITY] = column; 
      column += 1;
    } 
    else if(!(word.compare(0,4,"Velocity",0,4) &&
              word.compare(0,4,"VELOCITY",0,4) &&
              word.compare(0,4,"velocity",0,4))) {
      column2var[column]   = VELOCITY;
      column2var[column+1] = VELOCITY;
      column2var[column+2] = VELOCITY;
      var2column[VELOCITY] = column; 
      column += 3;
    }
    else if(!(word.compare(0,4,"Pressure",0,4) &&
              word.compare(0,4,"PRESSURE",0,4) &&
              word.compare(0,4,"pressure",0,4))) {
      column2var[column]   = PRESSURE;
      var2column[PRESSURE] = column; 
      column += 1;
    } 
    else if(!(word.compare(0,4,"Force",0,4) &&
              word.compare(0,4,"FORCE",0,4) &&
              word.compare(0,4,"force",0,4))) {
      column2var[column]   = FORCE;
      column2var[column+1] = FORCE;
      column2var[column+2] = FORCE;
      var2column[FORCE]    = column; 
      column += 3;
    } else {
      print_error("*** Error: I do not understand the word '%s' in the metafile.\n", word.c_str());
      exit_mpi();
    }
  }

  //Line #3 (and onwards) - time stamp and file label (i.e. the string between preflix and suffix)
  stamp.clear();
  double time_stamp;
  string file_label;
  for(int i=0; i<INT_MAX; i++) {
    getline(input, line);
    std::istringstream is(line);
    is >> time_stamp;
    if(is.fail())
      break;
    is >> file_label;
    if(is.fail())
      break;
    stamp.push_back(std::make_pair(time_stamp, file_label));
  }    

  //Sort by first component (i.e. time)
  std::sort(stamp.begin(), stamp.end());
  tmin = stamp.front().first;
  tmax = stamp.back().first;

  input.close();

}

//-----------------------------------------------------------------

void
DynamicLoadCalculator::ReadSnapshot(string filename, vector<vector<double> >& S)
{

  assert(strcmp(filename.c_str(), ""));

  // open file
  std::fstream input;
  input.open(filename.c_str(), std::fstream::in);
  if (!input.is_open()) {
    print_error("*** Error: Cannot open file %s.\n", filename.c_str());
    exit_mpi();
  } 

  int nCol = column2var.size();

  // clear the vector
  int initial_size = 1000;
  if(S.empty())
    S.resize(initial_size, vector<double>(nCol,0.0));

  // Now, start reading the file
  std::string word, line;
  //Line #1 --- user's comment (skip)
  input.ignore(512,'\n'); 

  int r;
  double data[nCol];
  for(r=0; r<INT_MAX; r++) {
  
    getline(input, line);

    std::istringstream is(line);
    bool done = false;
    for(int i=0; i<nCol; i++) {
      is >> data[i];
      if(is.fail()) {
        done = true;
        break;
      }
    }
    if(done)
      break;
    if(r>=(int)S.size())
      S.resize(S.size() + initial_size, vector<double>(nCol,0.0));
    for(int i=0; i<nCol; i++)
      S[r][i] = data[i];
  }

  S.resize(r);
} 

//-----------------------------------------------------------------

void
DynamicLoadCalculator::BuildKDTree(vector<vector<double> >& S, KDTree<PointIn3D,3>* tree, vector<PointIn3D>& p)
{
  if(tree)
    delete tree; //build a new tree

  int N = S.size();

  p.resize(N);

  int ix = var2column[COORDINATES];
  for(int i=0; i<N; i++) {
    Vec3D xyz(S[i][ix],S[i][ix+1],S[i][ix+2]);
    p[i] = PointIn3D(i, xyz);
  }

  tree = new KDTree<PointIn3D, 3>(N, p.data());
} 

//-----------------------------------------------------------------

void
DynamicLoadCalculator::InterpolateInSpace(vector<vector<double> >& S, KDTree<PointIn3D,3>* tree,
                                          vector<Vec3D>& X, int active_nodes, Var var, int var_dim,
                                          double* output)
{
  //active_nodes is the number of nodes (in X) that need force. May be different from X.size() in case of fracture
  //do NOT use "X.size()"!
  //assuming "output" has the right size (active_nodes x var_dim, or bigger)
  
  assert(tree);
  assert(var2column.find(var) != var2column.end());
  int col = var2column[var];

  int numPoints = iod.special_tools.transient_input.numPoints; //this is the number of points for interpolation
  int maxCand = numPoints*10;
  PointIn3D candidates[10*maxCand]; //a lot more than enough. (maxCand may be locally increased for failsafe)

  // find an initial/tentative cut-off distance
  assert(var2column.find(COORDINATES) != var2column.end());
  int ix = var2column[COORDINATES];
  assert(S.size()>=2);
  Vec3D tmpx0(S[0][ix],S[0][ix+1],S[0][ix+2]); 
  Vec3D tmpx1(S[1][ix],S[1][ix+1],S[1][ix+2]); 
  double cutoff = 0.1*(tmpx1-tmpx0).norm(); //will be adjusted


  // split the job among the processors (for parallel interpolation)
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  int block_size = active_nodes/mpi_size;
  int remainder = active_nodes - block_size*mpi_size;
  assert(remainder>=0 && remainder<mpi_size);

  int* counts = new int[mpi_size];
  int* displacements = new int[mpi_size];
  for(int i=0; i<mpi_size; i++) {
    counts[i] = (i<remainder) ? block_size + 1 : block_size;
    displacements[i] = (i<remainder) ? (block_size+1)*i : block_size*i + remainder;
  }
  assert(displacements[mpi_size-1]+counts[mpi_size-1] == active_nodes);

  int my_start_id = displacements[mpi_rank];
  int my_block_size = counts[mpi_rank];  


  //fprintf(stdout,"[%d] my_start_id = %d, my_block_size = %d.\n", mpi_rank, my_start_id, my_block_size);
  //MPI_Barrier(comm);
 

  //choose a radial basis function for interpolation
  void (*phi)(int, double[], double, double[]); //a function pointer
  switch (iod.special_tools.transient_input.basis) {
    case TransientInputData::MULTIQUADRIC :
      phi = MathTools::phi1;  break;
    case TransientInputData::INVERSE_MULTIQUADRIC :
      phi = MathTools::phi2;  break;
    case TransientInputData::THIN_PLATE_SPLINE :
      phi = MathTools::phi3;  break;
    case TransientInputData::GAUSSIAN :
      phi = MathTools::phi4;  break;
    default :
      phi = MathTools::phi2;  break; //Inverse multi-quadric
  }

/*
  double cutoffs[2000];
  int founds[2000];
*/
  // interpolation
  int maxIt = 1000;
  int index = my_start_id;
  for(int i=0; i<my_block_size; i++) {

    Vec3D& pnode(X[index]);

    // find sample points using KDTree
    int nFound = 0, counter = 0;
    double low_cut = 0.0, high_cut = DBL_MAX;
    while(nFound<numPoints || nFound>maxCand) {
      if(++counter>maxIt) {
        fprintf(stdout,"\033[0;31m*** Error: Cannot find required number of sample points for "
                       "interpolation (by RBF) after %d iterations. "
                       "Coord(3D):%e %e %e, Candidates: %d, cutoff = %e.\n\033[0m",
                        counter, pnode[0], pnode[1], pnode[2], nFound, cutoff);
        for(int i=0; i<std::min(nFound,maxCand); i++)
          fprintf(stdout,"%d  %d  %e  %e  %e  d = %e.\n", i, candidates[i].id, candidates[i].x[0],
                  candidates[i].x[1], candidates[i].x[2], (candidates[i].x-pnode).norm());

        fprintf(stdout,"low_cut = %e, high_cut = %e.\n", low_cut, high_cut);
//        for(int i=0; i<counter-1; i++)
//          fprintf(stdout,"%d  %d  %e.\n", i, founds[i], cutoffs[i]);

        exit(-1);
      }
      nFound = tree->findCandidatesWithin(pnode, candidates, maxCand, cutoff);
/*
      cutoffs[counter-1] = cutoff;
      founds[counter-1] = nFound;
*/
      if(nFound<numPoints) {
        low_cut = std::max(low_cut, cutoff);
        if(high_cut>0.5*DBL_MAX)
          cutoff *= 4.0;
        else
          cutoff = 0.5*(low_cut + high_cut); 
      }
      else if(nFound>maxCand) {
        high_cut = std::min(high_cut, cutoff);
        cutoff = 0.5*(low_cut + high_cut); 
      }
    

      if((high_cut - low_cut)/high_cut<1e-6) { //fail-safe
        nFound = tree->findCandidatesWithin(pnode, candidates, 10*maxCand, high_cut);
        if(nFound>10*maxCand) {
          fprintf(stdout,"\033[0;31m*** Error: Cannot find required number of sample points at any"
                         " cutoff distance.\n\033[0m");
          exit(-1); 
        } 

        assert(nFound>=numPoints);
        fprintf(stdout,"\033[0;35mWarning: Unusual behavior. Found %d candidates with cutoff = %e "
                       "(node: %e %e %e).\n\033[0m", nFound, high_cut, pnode[0], pnode[1], pnode[2]);
        break; 
      }

    }

    //figure out the actual points for interpolation (numPoints)
    vector<pair<double,int> > dist2node;
    for(int i=0; i<nFound; i++)
      dist2node.push_back(std::make_pair((candidates[i].x-pnode).norm(), candidates[i].id));

    if(nFound>numPoints)
      sort(dist2node.begin(), dist2node.end());
    dist2node.resize(numPoints);

    // prepare to interpolate 
    double xd[3*numPoints];
    for(int i=0; i<numPoints; i++) {
      xd[3*i]   = S[dist2node[i].second][ix];
      xd[3*i+1] = S[dist2node[i].second][ix+1];
      xd[3*i+2] = S[dist2node[i].second][ix+2];
    }
    double r0; //smaller than maximum separation, larger than typical separation
    r0 = dist2node.front().first + dist2node.back().first;
    double fd[numPoints];

    // interpolate
    for(int dim=0; dim<var_dim; dim++) {
      for(int i=0; i<numPoints; i++) 
        fd[i] = S[dist2node[i].second][col+dim];
 
      double *rbf_weight = new double[numPoints];
      double *interp = new double[numPoints];
      MathTools::rbf_weight(3, numPoints, xd, r0, phi, fd, rbf_weight);
      MathTools::rbf_interp(3, numPoints, xd, r0, phi, rbf_weight, 1, pnode, interp);

      output[var_dim*index + dim] = interp[0]; 

      delete [] rbf_weight;
      delete [] interp;
    }

    index++;

  }


  //communication
  for(int i=0; i<mpi_size; i++) {
    counts[i] *= var_dim;
    displacements[i] *= var_dim;
  }
  MPI_Allgatherv(MPI_IN_PLACE, var_dim*my_block_size, MPI_DOUBLE, output, counts, 
                 displacements, MPI_DOUBLE, comm);


  delete [] counts;
  delete [] displacements; 
}

//-----------------------------------------------------------------

void
DynamicLoadCalculator::InterpolateInTime(double t1, double* input1, double t2, double* input2, 
                                         double t, double* output, int size)
{
  assert(t2>t1);
  double c1 = (t2-t)/(t2-t1);
  double c2 = 1.0 - c1;
  for(int i=0; i<size; i++)
    output[i] = c1*input1[i] + c2*input2[i];
}

//-----------------------------------------------------------------

void
DynamicLoadCalculator::ComputeForces(TriangulatedSurface *surface, vector<Vec3D> *force, double t)
{

  if(std::get<0>(force_calculator)) {//User has specified a force calculator
    ApplyUserDefinedForces(surface, force, t);
    return;
  }
  

  // find the time interval
  int k0(-1), k1(-1);
  if(t<tmin) {
    k0 = k1 = 0;
    print_warning("- Warning: Calculating forces at %e by const extrapolation "
                  "(outside data interval [%e,%e]).\n", t, tmin, tmax);
  }
  else if(t>tmax) {
    k0 = k1 = stamp.size() - 1;
    print_warning("- Warning: Calculating forces at %e by const extrapolation "
                  "(outside data interval [%e,%e]).\n", t, tmin, tmax);
  }
  else {
    for(int k=1; k<(int)stamp.size(); k++) {
      if(t<stamp[k].first) {
        k0 = k-1;
        k1 = k;
        break;
      }
    }
  }
  assert(k0>=0 && k1>=0);

  // check if need to be build new tree(s)  
  if(k0 != id0) {
    if(k0 == id1) {
      id0   = id1;
      S0    = S1;
      tree0 = tree1;
    }
    else { //read data from file

      S0.reset(new vector<vector<double> >);
      string file_to_read = prefix + stamp[k0].second + suffix;
      ReadSnapshot(file_to_read, *S0);

      KDTree<PointIn3D,3> *tr(NULL);
      BuildKDTree(*S0, tr, tree0_data);
      tree0.reset(tr);

      id0 = k0;
    }
  } 
  if(k1 != id1) { //read data from file

      S1.reset(new vector<vector<double> >);
      string file_to_read = prefix + stamp[k1].second + suffix;
      ReadSnapshot(file_to_read, *S1);

      KDTree<PointIn3D,3> *tr(NULL);
      BuildKDTree(*S1, tr, tree1_data);
      tree1.reset(tr);

      id1 = k1;
  }

  // at this point, id0 = k0, id1 = k1

  // interpolate, first in space, then in time
  if(var2column.find(FORCE) != var2column.end()) {

    InterpolateInSpace(*S0, tree0.get(), surface->X, surface->active_nodes, FORCE, 3, (double*)F0.data());
    InterpolateInSpace(*S1, tree1.get(), surface->X, surface->active_nodes, FORCE, 3, (double*)F1.data());
    InterpolateInTime(stamp[id0].first, (double*)F0.data(), stamp[id1].first, (double*)F1.data(),
                      t, (double*)(force->data()), 3*surface->active_nodes);

  }
  else {
    print_error("*** Error: Currently, can only interpolate force by user-specified force.\n");
    exit_mpi();
  }

}

//-----------------------------------------------------------------

void
DynamicLoadCalculator::ApplyUserDefinedForces(TriangulatedSurface *surface, 
                                              std::vector<Vec3D> *force, double t)
{
  UserDefinedForces *calculator(std::get<0>(force_calculator));
  if(!calculator)
    return;

  vector<Vec3D> &Xs(surface->X);
  vector<Vec3D> &X0(surface->X0);
  vector<Int3> &elems(surface->elems);

  calculator->GetUserDefinedForces(t, X0.size(), (double*)X0.data(), (double*)Xs.data(),
                                   elems.size(), (int*)elems.data(),
                                   (double*)force->data()); //output
}



//-----------------------------------------------------------------




//-----------------------------------------------------------------







