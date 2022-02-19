#include<DynamicLoadCalculator.h>
#include<KDTree.h>
#include<rbf_interp_2d.hpp>
//#include<bits/stdc++.h> //std::sort
//#include<time.h>

using std::string;
using std::vector;

extern clock_t start_time; //time computation started (in Main.cpp)

//-----------------------------------------------------------------

DynamicLoadCalculator::DynamicLoadCalculator(IoData &iod_, MPI_Comm &comm_, 
                                             ConcurrentProgramsHandler &concurrent_)
                     : comm(comm_), iod(iod_), concurrent(concurrent_)
{ }

//-----------------------------------------------------------------

DynamicLoadCalculator::~DynamicLoadCalculator()
{ }

//-----------------------------------------------------------------

void
DynamicLoadCalculator::Run()
{
  print("\n");
  print("----------------------------------------------------\n");
  print("- Activated special tool: Dynamic load calculator. -\n");
  print("----------------------------------------------------\n");
  print("\n");

  prefix = string(iod.special_tools.transient_input.snapshot_file_prefix);
  suffix = string(iod.special_tools.transient_input.snapshot_file_suffix);
  ReadMetaFile(string(iod.special_tools.transient_input.metafile));

  if(concurrent.aeros != NULL)
    RunForAeroS();
  else {
    print_error("*** Error: No concurrent program specified.\n");
    exit_mpi();
  }

  print("\n");
  print("\033[0;32m==========================================\033[0m\n");
  print("\033[0;32m            NORMAL TERMINATION            \033[0m\n", t);
  print("\033[0;32m==========================================\033[0m\n");
  print("Total Computation Time: %f sec.\n", ((double)(clock()-start_time))/CLOCKS_PER_SEC);
  print("\n");
}

//-----------------------------------------------------------------

void
DynamicLoadCalculator::RunForAeroS()
{

  // Initialize Embedded Boundary Operator
  EmbeddedBoundaryOperator *embed = NULL;
  if(iod.concurrent.aeros.fsi_algo != AerosCouplingData::NONE) {
    embed = new EmbeddedBoundaryOperator(iod, true);
    concurrent.InitializeMessengers(embed->GetPointerToSurface(0),
                                    embed->GetPointerToForcesOnSurface(0));
  }
  assert(embed); //shouldn't be NULL...
  



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
    else if(!(word.compare(0,4,"Coordinate",0,4) &&
              word.compare(0,4,"COORDINATE",0,4) &&
              word.compare(0,4,"coordinate",0,4))) {
      column2var[column]   = COORDINATE;
      column2var[column+1] = COORDINATE;
      column2var[column+2] = COORDINATE;
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
  tmin = stamp.front.first();
  tmax = stamp.back.first();

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
    if(r>=S.size())
      S.resize(S.size() + initial_size, vector<double>(nCol,0.0));
    for(int i=0; i<nCol; i++)
      S[r][i] = data[i];
  }

  S.resize(r);
} 
//-----------------------------------------------------------------

void
DynamicLoadCalculator::BuildKDTree(vector<vector<double> >& S, KDTree<PointIn3D,3> *tree)
{
  if(tree)
    delete tree; //build a new tree

  int N = S.size();
  vector<PointIn3D> p;
  p.resize(N);
  int ix = var2column[COORDINATES];
  for(int i=0; i<N; i++)
    p[i] = Vec3D(S[i][ix],S[i][ix+1],S[i][ix+2]);

  tree = new KDTree<PointIn3D, 3>(N, p);
} 

//-----------------------------------------------------------------

void
DynamicLoadCalculator::InterpolateInSpace(vector<vector<double> >& S, KDTree<PointIn3D,3>* tree,
                                          vector<Vec3D>& X, Var var, int var_dim,
                                          double* output)
{
  //assuming "output" has the right size
  
  assert(tree);
  assert(var2column.find(var) != var2column.end());
  int col = var2column[var];

  int numPoints = iod.special_tools.transient_input.numPoints; //this is the number of points for interpolation
  int maxCand = numPoints*10;
  PointIn3D candidates[maxCand];

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

  int block_size = floor((double)X.size()/(double)mpi_size); 
  int remainder = X.size() - block_size*mpi_size;

  int* counts = new int[mpi_size];
  int* displacements = new int[mpi_size];
  for(int i=0; i<mpi_size; i++) {
    counts[i] = (i<remainder) ? block_size + 1 : block_size;
    displacements[i] = (i<remainder) ? (block_size+1)*i : block_size*i + remainder;
  }

  int my_start_id = displacements[mpi_rank];
  int my_block_size = counts[mpi_rank];  


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


  // interpolation
  int index = my_start_id;
  for(int i=0; i<my_block_size; i++) {

    Vec3D& pnode(X[index]);

    // find sample points using KDTree
    int nFound = 0, counter = 0;
    while(nFound<numPoints || nFound>maxCand) {
      if(++counter>1000) {
        fprintf(stderr,"*** Error: Cannot find required number of sample points for "
                       "interpolation (by RBF) after 2000 iterations. "
                       "Coord(3D):%e %e %e, Candidates: %d, cutoff = %e.\n",
                        pnode[0], pnode[1], pnode[2], nFound, cutoff);
        exit(-1);
      }
      nFound = tree.findCandidatesWithin(pnode, candidates, maxCand, cutoff);
      if(nFound==0)              cutoff *= 4.0;
      else if(nFound<numPoints)  cutoff *= 1.5*sqrt((double)numPoints/(double)nFound);
      else if(nFound>maxCand)    cutoff /= 1.5*sqrt((double)nFound/(double)maxCand);
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
 
      double *rbf_weight, *interp;
      rbf_weight = MathTools::rbf_weight(3, numPoints, xd, r0, phi, fd);
      interp = MathTools::rbf_interp(3, numPoints, xd, r0, phi, rbf_weight, 1, pnode);

      output[var_dim*index + dim] = interp[0]; 
    }
    index++;

    delete [] rbf_weight;
    delete [] interp;
  }


  //communication
  MPI_Allgatherv(MPI_IN_PLACE, var_dim*my_block_size, MPI_DOUBLE, output, var_dim*counts, 
                 var_dim*displacements, MPI_DOUBLE, comm);


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

