#include<LaserAbsorptionSolver.h>
#include<GeoTools.h>
#include<algorithm> //std::sort
#include <numeric> //std::iota

//#include<chrono> //for timing only
//using namespace std::chrono;

using std::vector;
using std::pair;

extern int verbose;

//--------------------------------------------------------------------------

LaserAbsorptionSolver::LaserAbsorptionSolver(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_, 
                         std::vector<VarFcnBase*> &varFcn_, SpaceVariable3D &coordinates_, 
                         SpaceVariable3D &delta_xyz_, SpaceVariable3D &volume_,
                         std::vector<GhostPoint> &ghost_nodes_inner_, 
                         std::vector<GhostPoint> &ghost_nodes_outer_)
                     : comm(comm_), iod(iod_), varFcn(varFcn_), coordinates(coordinates_), delta_xyz(delta_xyz_),
                       volume(volume_), ghost_nodes_inner(ghost_nodes_inner_),
                       ghost_nodes_outer(ghost_nodes_outer_),
                       L0(comm_, &(dm_all_.ghosted1_1dof)),
                       FluxIn(comm_, &(dm_all_.ghosted1_1dof)),
                       FluxOut(comm_, &(dm_all_.ghosted1_1dof)),
                       Phi(comm_, &(dm_all_.ghosted1_1dof)),
                       Level(comm_, &(dm_all_.ghosted1_1dof)),
                       Tag(comm_, &(dm_all_.ghosted1_1dof))
{
  // Check input parameters
  CheckForInputErrors();

  // Mesh info
  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);
  CalculateGlobalMeshInfo();

  // Get absorption coefficient for each material
  int numMaterials = iod.eqs.materials.dataMap.size();
  abs.resize(numMaterials, std::make_pair(0.0, 0.0)); //by default, set coeff = 0
  for (auto it = iod.laser.abs.dataMap.begin(); it != iod.laser.abs.dataMap.end(); it++) {
    if(it->second->materialid < 0 || it->second->materialid >= numMaterials) {
      fprintf(stderr,"ERROR: Found laser absorption coefficients for an unknown material (id = %d).\n", it->first);
      exit_mpi();
    }
    abs[it->first].first = it->second->slope;
    abs[it->first].second = it->second->alpha0;
  }

  // Store time history of source power (if provided by user)
  if(strcmp(iod.laser.source_power_timehistory_file, ""))
    ReadUserSpecifiedPowerFile(iod.laser.source_power_timehistory_file);

  // Calculate laser info ("source" and "laser_range")
  CalculateLaserInfo();

  // Calculate distance from each node to source
  CalculateDistanceToSource(Phi); 

  // Sort nodes in scope (Not all the nodes) by queue-level (primary) and dist-to-source (secondary)
  BuildSortedNodeList();

  // Create a custom communicator for each level
  BuildCustomizedCommunicators();
  //VerifySortedNodesAndCommunicators(); //for debug purpose

  // Find ghost nodes outside laser boundary (These include nodes inside and outside the physical domain!)
  SetupLaserGhostNodes();
  
  exit_mpi();
}

//--------------------------------------------------------------------------

LaserAbsorptionSolver::~LaserAbsorptionSolver()
{
  for(int i=0; i<levelcomm.size(); i++)
    if(levelcomm[i])
      delete levelcomm[i];
}

//--------------------------------------------------------------------------

void
LaserAbsorptionSolver::Destroy()
{
  L0.Destroy();
  FluxIn.Destroy();
  FluxOut.Destroy();
  Phi.Destroy();
  Level.Destroy();
  Tag.Destroy();
}

//--------------------------------------------------------------------------

void
LaserAbsorptionSolver::CheckForInputErrors()
{

  LaserData &laser(iod.laser);

  int error = 0;
  if(laser.source_distribution == LaserData::GAUSSIAN &&
     laser.source_power<=0.0 && strcmp(laser.source_power_timehistory_file,"")==0 ) {
    print_error("*** Error: Laser power must be specified for Gaussian distribution.\n");
    error ++;
  }
  if(laser.alpha<0.5 || laser.alpha>1.0+1e-14) {
    print_error("*** Error: The numerical parameter 'alpha' in Laser Absorption should be between 0.5 and 1.0.\n");
    error ++;
  }
  if(laser.relax_coeff<=0) {
    print_error("*** Error: The numerical parameter 'relax_coeff' in Laser Absorption should be positive.\n");
    error ++;
  }
  if(laser.source_dir_x*laser.source_dir_x +
     laser.source_dir_y*laser.source_dir_y +
     laser.source_dir_z*laser.source_dir_z == 0.0) {
    print_error("*** Error: Laser direction is not specified.\n");
    error ++;
  }
  if(laser.source_depth<=0.0) {//the user did not specify a (valid) depth
    print_error("*** Error: 'source_depth' in Laser Absorption is not specified.\n");
    error ++;
  }
  if(laser.focusing_angle_degrees<=-90 || laser.focusing_angle_degrees>=90) {
    print_error("*** Error: 'focusing_angle_degrees' in Laser Absorption must be between -90 and 90.\n");
    error ++;
  }
  if(error>0)
    exit_mpi();
}

//--------------------------------------------------------------------------

void
LaserAbsorptionSolver::CalculateGlobalMeshInfo()
{
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
  Vec3D*** dxyz = (Vec3D***)delta_xyz.GetDataPointer();

  dxmin_glob = DBL_MAX; = dymin_glob = dzmin_glob = DBL_MAX;

  for(int i=ii0; i<iimax; i++)
    dxmin_glob[0] = std::min(dxmin_glob[0], dxyz[kk0][jj0][i][0]);
  for(int j=jj0; j<jjmax; j++)
    dxmin_glob[1] = std::min(dxmin_glob[1], dxyz[kk0][j][ii0][1]);
  for(int k=kk0; k<kkmax; k++)
    dxmin_glob[2] = std::min(dxmin_glob[2], dxyz[k][jj0][ii0][2]);

  MPI_Allreduce(MPI_IN_PLACE, dxmin_glob, 3, MPI_DOUBLE, MPI_MIN, comm);

  for(int i=0; i<3; i++) {
    xmin_glob[i]       = coords[k0][j0][i0][i];
    xmin_glob_ghost[i] = coords[kk0][jj0][ii0][i];
    xmax_glob[i]       = coords[kmax][jmax][imax][i];
    xmax_glob_ghost[i] = coords[kkmax][jjmax][iimax][i];
  }
 
  MPI_Allreduce(MPI_IN_PLACE, xmin_glob, 3, MPI_DOUBLE, MPI_MIN, comm);
  MPI_Allreduce(MPI_IN_PLACE, xmax_glob, 3, MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(MPI_IN_PLACE, xmin_glob_ghost, 3, MPI_DOUBLE, MPI_MIN, comm);
  MPI_Allreduce(MPI_IN_PLACE, xmax_glob_ghost, 3, MPI_DOUBLE, MPI_MAX, comm);

  delta_xyz.RestoreDataPointerToLocalVector();
  coordinates.RestoreDataPointerToLocalVector();
}

//--------------------------------------------------------------------------

void
LaserAbsorptionSolver::ReadUserSpecifiedPowerFile(const char *source_power_timehistory_file)
{
  std::fstream input;
  char *filename = new char[strlen(source_power_timehistory_file) + 5];
  sprintf(filename, "%s", iod.laser.source_power_timehistory_file);
  input.open(filename, std::fstream::in);
  if(!input.is_open()) {
    print_error("*** Error: Could not open user-specified laser power file %s.\n", filename);
    exit_mpi();
  }

  double time0(-DBL_MAX), time, power;
  while(input >> time >> power) {
    assert(time>time0);
    source_power_timehistory.push_back(std::make_pair(time, power));
    time0 = time;
  }
  print("- Loaded user-specified laser power file %s, final time: %e.\n", filename, time);
  if(source_power_timehistory[0].first>0) //if the first time stamp is greater than 0, insert one at 0.
    source_power_timehistory.insert(source_power_timehistory.begin(), 
                                    std::make_pair(0.0, source_power_timehistory[0].second));
  else if(source_power_timehistory[0].first<0) {
    print_error("*** Error: Detected negative time stamp (%e) in laser power file.\n",
                source_power_timehistory[0].first);
    exit_mpi();
  }
}

//--------------------------------------------------------------------------

void
LaserAbsorptionSolver::CalculateLaserInfo()
{

  LaserData &laser(iod.laser);

  source.radius = laser.source_radius;
  source.x0     = Vec3D(laser.source_center_x, laser.source_center_y, laser.source_center_z);
  source.dir    = Vec3D(laser.source_dir_x, laser.source_dir_y, laser.source_dir_z);
  source.dir   /= source.dir.norm(); //normalize direction

  double PI = 2.0*acos(0.0);
  source.angle = laser.focusing_angle_degrees/180.0*PI; //in radians

  if(angle>0.0){
    double d2f = laser.source_radius / sin(angle);
    double eps = std::max(dxmin_glob[0], std::max(dxmin_glob[1], dxmin_glob[2]));
    if(laser.range > d2f - 2.0*eps) {
      print_error("*** Error: Laser range must be less than the radius of curvature (by at least two cells).\n");
      exit_mpi();
    }
  }
  laser_range = laser.range;

  if(source.angle == 0) {//parallel laser
    R = DBL_MAX;
    O = Vec3D(DBL_MAX,DBL_MAX,DBL_MAX);
  }
  else if(source.angle>0) {//focusing laser 
    source.R = source.radius/sin(source.angle); //radius of curvature
    source.O = source.x0 + source.R*cos(source.angle)*source.dir; //origin of the circle/curve 
  }
  else if(source.angle<0) {//diverging laser
    double theta = -source.angle; //a positive angle
    source.R = source.radius/sin(theta); //radius of curvature
    source.O = source.x0 - source.R*cos(theta)*source.dir; //origin
  }

  print("- Laser absorption solver activated. Focusing angle: %f degrees; Source radius: %f; Range: %f\n",
        laser.focusing_angle_degrees, laser.source_radius, laser_range);

}

//--------------------------------------------------------------------------

void
LaserAbsorptionSolver::CalculateDistanceToSource(SpaceVariable3D &Phi)
{

  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
  double*** phi   = Phi.GetDataPointer();

  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        phi[k][j][i] = -1.0; //default (out of scope)

        int loc_info;
        double dist = DistanceToSource(coords[k][j][i], loc_info);
        if(loc_info == 0) //inside laser domain
          phi[k][j][i] = dist; 
      }

  coordinates.RestoreDataPointerToLocalVector();
  Phi.RestoreDataPointerAndInsert(); 

}

//--------------------------------------------------------------------------

double
LaserAbsorptionSolver::DistanceToSource(Vec3D &x, int &loc, double *image)
{
  //each digit of loc corresponds to a check
  //if "image" is not NULL, and if x is out of scope, "image" will get the coordinates
  //of the image point. (This function assumes that x is at least somewhat close to
  // the laser scope --- not an arbitrary point in R3.)
  
  double dist = DBL_MAX;
  loc = 0; //default: inside laser scope

  if(source.angle == 0) { //parallel laser

    double d2p  = (x - source.x0)*source.dir; //distance to source plane
    if(d2p<0)            loc += 1;
    if(d2p>laser_range)  loc += 10;
        
    Vec3D xp = x - source.x0 - d2p*source.dir; //projection of x on the plane of the source
    if(xp.norm()>source.radius) loc += 100;

    dist = d2p;

    if(image && loc) { //calculate image coordinates
      Vec3D xim = x;
      if(loc == 1 || loc == 101)
        xim -= 2.0*d2p*source.dir; //d2p<0
      if(loc == 10 || loc == 110)
        xim -= 2.0*(d2p-laser_range)*source.dir;
      if(loc >= 100) {//100, 101, 110
        Vec3D xp_dir = xp/xp.norm();  
        xim -= 2.0*(xp.norm()-source.radius)*xp_dir;
      } 
      for(int i=0; i<3; i++)
        image[i] = xim[i];
    }

  }
  else if(source.angle > 0) { //focusing laser

    Vec3D Ox = x - source.O;
    double d2O = Ox.norm();
    if(d2O>source.R)               loc += 1;
    if(d2O<(source.R-laser_range)) loc += 10;

    double r = Ox*(-source.dir);
    if(r<0) {
      fprintf(stderr,"*** Error: Detected ghost node on the wrong side of the origin. Reduce laser range.\n");
      exit(-1);
    }

    Vec3D xp = source.O - r*source.dir;
    if((x - xp).norm() > r*tan(source.angle))  loc += 100;

    dist = source.R - d2O;

    if(image && loc) { //calculate image coordinates
      Vec3D xim = x;
      if(loc == 1 || loc == 101)
        xim -= 2.0*(d2O-source.R)*Ox/Ox.norm();
      if(loc == 10 || loc == 110)
        xim += 2.0*(source.R-laser_range-d2O)*Ox/Ox.norm();
      if(loc >= 100) { //100, 101, 110
        Vec3D Oxim = xim - source.O;
        double d2Oxim = Oxim.norm();
        double theta = acos(fabs(Oxim*source.dir)/d2Oxim);
        assert(theta>source.angle);
        theta = 2.0*source.angle - theta;
        Vec3D ximp = source.O + d2Oxim*cos(theta)*source.dir;
        Vec3D vert_dir = x - xp;
        vert_dir /= vert_dir.norm();
        xim = ximp + d2Oxim*sin(theta)*vert_dir;
      }
      for(int i=0; i<3; i++)
        image[i] = xim[i];
    }
  }
  else if(source.angle < 0) { //focusing laser

    double theta = -source.angle; //a positive angle
    Vec3D Ox = x - source.O;
    double d2O = Ox.norm();

    if(d2O<source.R)               loc += 1;
    if(d2O>(source.R+laser_range)) loc += 10;

    double r = Ox*source.dir;
    Vec3D xp = source.O + r*source.dir;

    if((x - xp).norm() > r*tan(theta))  loc += 100;

    if(r<0) {
      fprintf(stderr,"*** Error: Detected ghost node on the wrong side of the origin (diverging laser).\n");
      exit(-1);
    }

    dist = d2O - source.R;

    if(image && loc) { //calculate image coordinates
      Vec3D xim = x;
      if(loc == 1 || loc == 101)
        xim += 2.0*(source.R-d2O)*Ox/Ox.norm();
      if(loc == 10 || loc == 110)
        xim -= 2.0*(d2O-source.R-laser_range)*Ox/Ox.norm();
      if(loc >= 100) { //100, 101, 110
        Vec3D Oxim = xim - source.O;
        double d2Oxim = Oxim.norm();
        theta = acos(fabs(Oxim*source.dir)/d2Oxim); //re-defining theta (the previous one is no longer needed)
        assert(theta>(-source.angle));
        theta = 2.0*(-source.angle) - theta;
        Vec3D ximp = source.O + d2Oxim*cos(theta)*source.dir;
        Vec3D vert_dir = x - xp;
        vert_dir /= vert_dir.norm();
        xim = ximp + d2Oxim*sin(theta)*vert_dir;
      }
      for(int i=0; i<3; i++)
        image[i] = xim[i];
    }
  }

  return dist;
}

//--------------------------------------------------------------------------

int
LaserAbsorptionSolver::FindImagePoint(Vec3D &x, Vec3D &image)
{
  int loc;
  double dist = DistanceToSource(x, loc, image);
  return loc;
}

//--------------------------------------------------------------------------

bool
LaserAbsorptionSolver::LaserDistCompare(NodalLaserInfo& node1, NodalLaserInfo& node2) 
{
  //compare queue-level, then dist-to-source  

  if(node1.level <  node2.level) return true;
  else if(node1.level > node2.level) return false;

  if(node1.phi < node2.phi) return true;
  else if (node1.phi > node2.phi) return false;

  if(node1.k < node2.k) return true;
  else if(node1.k > node2.k) return false;

  if(node1.j < node2.j) return true;
  else if(node1.j > node2.j) return false;

  if(node1.i < node2.i) return true;
  else if(node1.i > node2.i) return false;

  return false;
}

//--------------------------------------------------------------------------

void
LaserAbsorptionSolver::BuildSortedNodeList()
{
  ResetTag();

  int NX, NY, NZ;
  coordinates.GetGlobalSize(&NX, &NY, &NZ);

  double*** phi = Phi.GetDataPointer();

  //----------------------------------------------------------------------------------------
  //Step 1: Initialize the set of "sortedNodes" to be the nodes that NEED TO BE SORTED
  //----------------------------------------------------------------------------------------
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++)
        if(phi[k][j][i] >= 0.0)
          sortedNodes.push_back(NodalLaserInfo(i,j,k, phi[k][j][i], INT_MAX));
  numNodesInScope = sortedNodes.size();

  //----------------------------------------------------------------------------------------
  //Step 2: pre-sort "sortedNodes" based on distance to source (i.e. phi)
  //----------------------------------------------------------------------------------------
  std::sort(sortedNodes.begin(), sortedNodes.end(), LaserAbsorptionSolver::LaserDistCompare); 


  //----------------------------------------------------------------------------------------
  //Step 3: Get complete info for sortedNodes, and order them
  //----------------------------------------------------------------------------------------
  int tip, level, nNodesOnLevel;
  double*** tag = Tag.GetDataPointer();

  //Step 3.1. Create Level 0, i.e., nodes within the ``depth'' of the laser source
  tip = 0;
  level = 0;
  nNodesOnLevel = 0;
  for(auto it = sortedNodes.begin(); it != sortedNodes.end(); it++) {
    if(it->phi <= iod.laser.source_depth) {
      int i(it->i), j(it->j), k(it->k);
      it->level = level;
      tag[k][j][i] = 1; //computed 
      nNodesOnLevel++;
    } else
      break; //because sortedNodes are already sorted by phi
  }
  queueCounter.push_back(nNodesOnLevel);
  std::sort(sortedNodes.begin(), sortedNodes.end(), LaserDistCompare); 
  Tag.RestoreDataPointerAndInsert();

  //Step 3.2. Create the other levels one-by-one
  tip += nNodesOnLevel;
  level++;
  int done = 0;
  while(!done) {

    nNodesOnLevel = 0;
    tag = Tag.GetDataPointer();
    for(auto it = sortedNodes.begin() + tip;  it != sortedNodes.end();  it++) {
      int i(it->i), j(it->j), k(it->k);

      // if a neighbor is tagged (meaning on a previous level), set this node
      if(!( (i-1>=0 && phi[k][j][i-1]<phi[k][j][i] && tag[k][j][i-1]==0) || 
            (i+1<NX && phi[k][j][i+1]<phi[k][j][i] && tag[k][j][i+1]==0) ||
            (j-1>=0 && phi[k][j-1][i]<phi[k][j][i] && tag[k][j-1][i]==0) || 
            (j+1<NY && phi[k][j+1][i]<phi[k][j][i] && tag[k][j+1][i]==0) || 
            (k-1>=0 && phi[k-1][j][i]<phi[k][j][i] && tag[k-1][j][i]==0) || 
            (k+1<NZ && phi[k+1][j][i]<phi[k][j][i] && tag[k+1][j][i]==0) ) ) {
        it->level = level;
        nNodesOnLevel++;
      }

    } 
    queueCounter.push_back(nNodesOnLevel);
    std::sort(sortedNodes.begin() + tip, sortedNodes.end(), LaserDistCompare); 

    //Update tag for nodes on the current level
    for(auto it = sortedNodes.begin() + tip; it != sortedNodes.begin() + tip + nNodesOnLevel; it++)
      tag[it->k][it->j][it->i] = 1;
    Tag.RestoreDataPointerAndInsert();

    tip += nNodesOnLevel;
    assert(tip<=sortedNodes.size());

    done = (tip == sortedNodes.size()) ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &done, 1, MPI_INT, MPI_MIN, comm);
    if(done)
      break;

    level++;
  }

  Phi.RestoreDataPointerToLocalVector();

  //Populate Level
  Level.SetConstantValue(-1, true);
  double*** lev = Level.GetDataPointer();
  int previous_node_level = 0; //for code verification
  for(auto it = sortedNodes.begin();  it != sortedNodes.end();  it++) {
    int i(it->i), j(it->j), k(it->k);
    assert(it->level >= previous_node_level); //for verification
    previous_node_level = it->level; 

    lev[k][j][i] = it->level;
  }
  Level.RestoreDataPointerAndInsert();
}

//--------------------------------------------------------------------------

void
LaserAbsorptionSolver::ResetTag()
{
  //-1: out of scope (no need to compute intensity); 
  // 0: intensity not computed yet, but need to be computed;  
  // 1: intensity has been computed

  double*** tag = Tag.GetDataPointer();
  double*** phi = Phi.GetDataPointer();

  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        if(phi[k][j][i] >= 0.0)
          tag[k][j][i] = 0;
        else
          tag[k][j][i] = -1;

      }
  
  Tag.RestoreDataPointerAndInsert();
  Phi.RestoreDataPointerToLocalVector();
}

//--------------------------------------------------------------------------

void
LaserAbsorptionSolver::BuildCustomizedCommunicators()
{

  // Build one custom communicator per level  
  vector<vector<Int3> > nodes_on_level(queueCounter.size());

  double*** level = Level.GetDataPointer();
  for(auto it = ghost_nodes_inner.begin(); it != ghost_nodes_inner.end(); it++) { 
    int i(it->ijk[0]), j(it->ijk[1]), k(it->ijk[2]); 
    if(level[k][j][i]<0)
      continue;
    assert(level[k][j][i]<queueCounter.size());
    nodes_on_level[level[k][j][i]].push_back(Int3(i,j,k));
  }
  Level.RestoreDataPointerToLocalVector();

  levelcomm.clear();
  for(int lvl = 0; lvl<queueCounter.size(); lvl++) {
    //fprintf(stderr,"Level[%d]: number of ghost nodes = %d.\n", lvl, (int)nodes_on_level[lvl].size());
    levelcomm.push_back(new CustomCommunicator(comm, L0, nodes_on_level[lvl]));
    //print("Good! Created comm level %d.\n", lvl);
  }

}

void
LaserAbsorptionSolver::VerifySortedNodesAndCommunicators()
{
  //---------------------------------------------------
  // Verify communicators
  double*** l0 = L0.GetDataPointer();
  double*** level = Level.GetDataPointer();
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++)
        l0[k][j][i] = level[k][j][i];
  L0.RestoreDataPointerToLocalVector();

  l0 = L0.GetDataPointer();
  for(int lvl = 0; lvl<queueCounter.size(); lvl++) {
    levelcomm[lvl]->ExchangeAndInsert(l0);
    for(auto it = ghost_nodes_inner.begin(); it != ghost_nodes_inner.end(); it++) {
      int i(it->ijk[0]), j(it->ijk[1]), k(it->ijk[2]);
      if(level[k][j][i]==lvl) 
        assert(l0[k][j][i] == lvl);
    }
  }
  L0.RestoreDataPointerToLocalVector();
  Level.RestoreDataPointerToLocalVector();


  //---------------------------------------------------
  // Verify sorted nodes
  double*** phi = Phi.GetDataPointer();
  level = Level.GetDataPointer();
  auto it = sortedNodes.begin();
  int level_previous = 0;
  double phi_previous = -1;
  int counter = 0;
  for(int lvl = 0; lvl<queueCounter.size(); lvl++) {
    for(int n = 0; n < queueCounter[lvl]; n++) {
      int i(it->i),j(it->j),k(it->k);
      assert(it->level == level[k][j][i]);
      assert(it->level == lvl);
      assert(it->phi == phi[k][j][i]);
      assert(it->level >= level_previous);
      if(it->level == level_previous)
        assert(it->phi >= phi_previous);
      counter++;
      it++;
      level_previous = it->level;
      phi_previous = it->phi;
    }
  }
  assert(counter == numNodesInScope);

  L0.SetConstantValue(0.0,true);
  Level.RestoreDataPointerToLocalVector();
  Phi.RestoreDataPointerToLocalVector();

  //---------------------------------------------------


/*
  //---------------------------------------------------
  // Timing 
  MPI_Barrier(comm);
  auto start = high_resolution_clock::now();
  l0 = L0.GetDataPointer();
  for(int counter = 0; counter < 100; counter++) {
    for(int lvl = 0; lvl<queueCounter.size(); lvl++) {
      levelcomm[lvl]->ExchangeAndInsert(l0);
    }
  }
  L0.RestoreDataPointerToLocalVector();
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);
  double cost = duration.count();
  MPI_Allreduce(MPI_IN_PLACE, &cost, 1, MPI_DOUBLE, MPI_MAX, comm);
  print(stderr,"cost is %e.\n", cost); 


  MPI_Barrier(comm);
  start = high_resolution_clock::now();
  for(int counter = 0; counter < 100; counter++) {
    for(int lvl = 0; lvl<queueCounter.size(); lvl++) {
      l0 = L0.GetDataPointer();
      L0.RestoreDataPointerAndInsert();
    }
  }
  stop = high_resolution_clock::now();
  duration = duration_cast<microseconds>(stop - start);
  cost = duration.count();
  MPI_Allreduce(MPI_IN_PLACE, &cost, 1, MPI_DOUBLE, MPI_MAX, comm);
  print(stderr,"cost is %e.\n", cost); 
  //---------------------------------------------------
*/

}

//--------------------------------------------------------------------------

void
LaserAbsorptionSolver::SetupLaserGhostNodes()
{

  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();

  //----------------------------------------------------------------
  // Step 1: Find ghost nodes (inside and outside physical domain)
  //----------------------------------------------------------------
  Tag.SetConstantValue(0, true);
  double*** level = Level.GetDataPointer(); //level only tags nodes inside subdomain
  double*** tag   = Tag.GetDataPointer(); 
  //loop through sorted nodes on level 1, 2, 3, ...
  for(int n = queueCounter[0]; n < sortedNodes.size(); n++) {
    int i(sortedNodes[n].i), j(sortedNodes[n].j), k(sortedNodes[n].k); //all sorted nodes are inside subdomain
    // Note: because of the Cartesian grid and the assuption that layer boundaries are straight lines, we only
    //       need to loop through the following 6 direct neighbors, instead of 26 neighbors based on element 
    //       connectivity.
    if(level[k][j][i-1]<0 && tag[k][j][i-1]==0)  tag[k][j][i-1] = 1;
    if(level[k][j][i+1]<0 && tag[k][j][i+1]==0)  tag[k][j][i+1] = 1;
    if(level[k][j-1][i]<0 && tag[k][j-1][i]==0)  tag[k][j-1][i] = 1;
    if(level[k][j+1][i]<0 && tag[k][j+1][i]==0)  tag[k][j+1][i] = 1;
    if(level[k-1][j][i]<0 && tag[k-1][j][i]==0)  tag[k-1][j][i] = 1;
    if(level[k+1][j][i]<0 && tag[k+1][j][i]==0)  tag[k+1][j][i] = 1;
  }  
  Tag.RestoreDataPointerAndInsert();
  Level.RestoreDataPointerToLocalVector();

  vector<Int3> ghost_tmp;

  tag = Tag.GetDataPointer();
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++)
        if(tag[k][j][i] == 1)
          ghost_tmp.push_back(Int3(i,j,k));
  Tag.RestoreDataPointerToLocalVector();

  //----------------------------------------------------------------
  // Step 2: Find the location (coordinates) of the image of each ghost node
  //----------------------------------------------------------------
  vector<Vec3D> image;
  vector<int>   status;
  vector<double> image_phi;
  // status: 
  // (1) If the image is outside laser domain,
  //     let Q be the closest point on the boundary of the laser domain
  // (1.1) if Q is on the source disk/sector                ---> status = 1
  // (1.2)       ...  the disk/sector of max. range         ---> status = 2
  // (1.3)       ...  the side boundary of the laser domain ---> status = 3
  // (1.4)       ...  an edge of the laser domain boundary  ---> status = 4
  // (2) If the image is inside laser domain, (the ghost node must be
  //                                           outside of the "physical domain")
  // (2.1) if it is next to a symmetry boundary (i.e. should be populated by mirroring) ---> status  = -1
  //       if it is "downstream" (i.e. should be populated by linear extrap.) --> report error!
  //
  // (status = -99) reserved for problematic nodes (see below).
  //
  // Note: Althought status 1 and 2 are tracked here. In the embedded boundary method for laser absorption,
  //       we don't really need to deal with these two cases, as long as (1) the scheme is strictly upwinding
  //       and (2) the source depth covers at least one (or 2?) layers of nodes.
  //
  for(auto it = ghost_tmp.begin(); it != ghost_tmp.end(); it++) {
    int i((*it)[0]), j((*it)[1]), k((*it)[2]);
    Vec3D& x(coords[k][j][i]);

    int location_info;
    Vec3D image;
    double dist = DistanceToSource(x, location_info, image);
    image_phi.push_back(dist);

    if(location_info) {//Category (1)
      image.push_back(image); //inserting the coordinates
      if(loation_info == 1)        status.push_back(1);
      else if(loation_info == 10)  status.push_back(2);
      else if(loation_info == 100) status.push_back(3);
      else                         status.push_back(4);
    }
    else { //Category (2)
      assert(coordinates.BoundaryType(i,j,k) == 1); //boundary-face
      bool found = false;
      for(auto it = ghost_nodes_outer.begin(); it != ghost_nodes_outer.end(); it++) {
        if(it->ijk[0] == i && it->ijk[1] == j && it->ijk[2] == k) {
        
          // verify that it is next to a symmetry boundary (or the laser)
          if(fabs(it->outward_normal*source.dir)>1.0e-12) {
            fprintf(stderr,"*** Error: Detected ghost node (%d,%d,%d)(%e,%e,%e) within the laser domain. Reduce laser range.\n",
                    i,j,k, x[0],x[1],x[2]);   
            exit(-1);
          }

          image.push_back(Vec3D(it->image_ijk[0], it->image_ijk[1], it->image_ijk[2])); //NOTE: inserting the indices here.
          status.push_back(-1);
          found = true;
          break;
        }
      }
      assert(found);
    }
  }

  //----------------------------------------------------------------
  // Step 3: Setup trilinear interpolation coeffs for each mirror image
  //----------------------------------------------------------------
  vector<Int3> image_ijk(image.size(), Int3(INT_MAX)); //lower-left corner of the interpolation box
  vector<Vec3D> image_xi(image.size(), Vec3D(0.0));
  vector<int> prob_nodes;
  for(int n=0; n<image.size(); n++) {
    if(status[n] == -1) {
      image_ijk[n] = Int3((int)image[n][0], (int)image[n][1], (int)image[n][2]);
    }
    else {
      int ihere(INT_MAX), jhere(INT_MAX), khere(INT_MAX);
      for(int i=ii0; i<iimax-1; i++) 
        if(image[n][0]>=coords[kk0][jj0][i][0] && image[n][0]<coords[kk0][jj0][i+1][0]) {ihere = i; break;} 
      for(int j=jj0; j<jjmax-1; j++) 
        if(image[n][1]>=coords[kk0][j][ii0][1] && image[n][1]<coords[kk0][j+1][ii0][1]) {jhere = j; break;} 
      for(int k=kk0; k<kkmax-1; k++) 
        if(image[n][2]>=coords[k][jj0][ii0][2] && image[n][2]<coords[k+1][jj0][ii0][2]) {jhere = j; break;} 

      if(ihere==INT_MAX || jhere==INT_MAX || khere==INT_MAX) {//not found
        //I think if this happens, the image must be outside of the expanded mesh that includes
        //the ghost boundary layer. If not, report error --- and we need to fix it!
        fprintf(stderr,"- The image of node (%d %d %d)(%e,%e,%e) is outside the computational domain.\n",
                ghost_tmp[n][0], ghost_tmp[n][1], ghost_tmp[n][2], 
                coords[ghost_tmp[n][2]][ghost_tmp[n][1]][ghost_tmp[n][0]][0],
                coords[ghost_tmp[n][2]][ghost_tmp[n][1]][ghost_tmp[n][0]][1],
                coords[ghost_tmp[n][2]][ghost_tmp[n][1]][ghost_tmp[n][0]][2]);
        assert(image[n][0]<xmin_glob_ghost[0] || image[n][0]>=xmax_glob_ghost[0] ||
               image[n][1]<xmin_glob_ghost[1] || image[n][1]>=xmax_glob_ghost[1] ||
               image[n][2]<xmin_glob_ghost[2] || image[n][2]>=xmax_glob_ghost[2]);
        status[n] = -99; //change its status to -99
        prob_nodes.push_back(n);
        continue;  //deal with all the problematic images later
      }

      image_ijk[n] = Int3(ihere,jhere,khere);
      for(int i=0; i<3; i++)
        image_xi[n][i] = (image[n][i] - coords[khere][jhere][ihere][i])
                       / (coords[khere+1][jhere+1][ihere+1][i] - coords[khere][jhere][ihere][i]);
    }
  }

  // Now, work on the problematic nodes
  if(prob_nodes.size()>0) {

    vector<pair<Int3,double> > closest(prob_nodes.size(), std::make_pair(Int3(0), DBL_MAX));
    for(auto it0 = sortedNodes.begin(); it0 != sortedNodes.end(); it0++) {
      int i(it0->i), j(it0->j), k(it0->k);
      for(int n=0; n<prob_nodes.size(); n++) {
        double dist = (image[prob_nodes[n]] - coords[k][j][i]).norm();
        if(dist < closest[n].second) {
          closest[n].second = dist;
          closest[n].first = Int3(i,j,k);
        } 
      }
    }
    double*** phi = Phi.GetDataPointer();
    for(int p=0; p<prob_nodes.size(); p++) {
      int n = prob_nodes[p];
      if(closest[p].second == DBL_MAX) {
        fprintf(stderr,"*** Error: Cannot find a valid image for ghost node (%d,%d,%d)(%e,%e,%e). "
                       "Original image: (%e,%e,%e).\n", ghost_tmp[n][0], ghost_tmp[n][1], ghost_tmp[n][2],
                       coords[ghost_tmp[n][2]][ghost_tmp[n][1]][ghost_tmp[n][0]][0],
                       coords[ghost_tmp[n][2]][ghost_tmp[n][1]][ghost_tmp[n][0]][1],
                       coords[ghost_tmp[n][2]][ghost_tmp[n][1]][ghost_tmp[n][0]][2],
                       image[n][0], image[n][1], image[n][2]);
        exit(-1);
      }
      image_ijk[n] = closest[n].first;
      image_xi[n]  = Vec3D(0,0,0);
      image_phi[n] = phi[(closest[p].first)[2]][(closest[p].first)[1]][(closest[p].first)[0]];
    }
    Phi.RestoreDataPointerToLocalVector();

  }

  coordinates.RestoreDataPointerToLocalVector();
  
  //----------------------------------------------------------------
  // Step 4: Sort the list of ghost nodes based on distance to source
  //         (image_phi).
  //----------------------------------------------------------------
  vector<int> ordered(image_phi.size());
  std::iota(std::begin(ordered), std::end(ordered), 0); //ordered = {0, 1, 2, 3, ...}
  std::sort(std::begin(ordered), std::end(ordered), 
            [&image_phi](int lhs, int rhs) {return image_phi[lhs]<image_phi[rhs];} );  

  //----------------------------------------------------------------
  // Step 5: Figure out formula for ghost node update
  //----------------------------------------------------------------
  double eps = std::min(dxmin_glob[0], std::min(dxmin_glob[1], dxmin_glob[2])); //geometry tolerance
  eps /= 1.0e-10;
  sortedGhostNodes.resize(ordered.size(), std::make_pair(Int3(0), EmbeddedBoundaryFormula(eps)));
  for(int o = 0; o < ordered.size(); o++) {
    int n = ordered[o]; 
    sortedGhostNodes[o].first = ghost_tmp[n];
    if(status[n] == -1 || status[n] == -99)
      sortedGhostNodes[o].second.SpecifyFormula(EmbeddedBoundaryFormula::MIRRORING, 
                                                EmbeddedBoundaryFormula::NODE, 
                                                vector<Int3>(1,image_ijk[n]), vector<double>(1,1.0));
    else
      sortedGhostNodes[o].second.BuildMirroringFormula(ghost_tmp[n], image_ijk[n], image_xi[n]);
  }
  
}

//--------------------------------------------------------------------------

void
LaserAbsorptionSolver::SetSourceRadiance(SpaceVariable3D &L, double t)
{
  // User-specified power time-history overrides (fixed) power, which then overrides (fixed) intensity.
  // If user-specified power time-history is available, find the current power by linear interpolation

  double power = -1.0, dt;
  if(!source_power_timehistory.empty()) {//user-specified power time-history is available
    for(auto it=source_power_timehistory.begin(); it!=source_power_timehistory.end(); it++) {
      auto next = it + 1;
      if(next==source_power_timehistory.end()) //reached the end
        power = it->second;
      else if(t>=it->first && t<next->first) {
        dt = next->first - it->first;
        power = (next->first - t)/dt*it->second + (t - it->first)/dt*next->second;
        break;
      }
    }
    if(power <0.0) {
      print_error("*** Error: Found negative laser source power %e at t = %e.\n", power, t);
      exit_mpi();
    }
  } else
    power = iod.laser.source_power;

  if(iod.laser.source_distribution == LaserData::CONSTANT){ //The Laser source is uniform

    double l0;
    if(power > 0.0) {//power (if specified) overrides source_intensity
      double PI = 2.0*acos(0.0);
      //area of sector surface, also applicable to parallel laser (angle = 0)
      double source_area = 2.0*PI*source.radius*source.radius/(1.0+cos(source.angle)); 
      l0 = power/source_area;
    } else
      l0 = iod.laser.source_intensity;

    double*** l = L.GetDataPointer();
    for(int n=0; n<queueCounter[0]; n++) //level 0
      l[sortedNodes[n].k][sortedNodes[n].j][sortedNodes[n].i] = l0;
    levelcomm[0]->ExchangeAndInsert(l);
    L.RestoreDataPointerToLocalVector();

  }
  else if(iod.laser.source_distribution == LaserData::GAUSSIAN){ //The Laser source is Gaussia

    //give the direction of X
    double beamwaist = iod.laser.source_beam_waist;
    double PI = 2.0*acos(0.0);

    double*** l = L.GetDataPointer();
    Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
    for(int n=0; n<queueCounter[0]; n++) {//level 0
      int i(sortedNodes[n].i), j(sortedNodes[n].j), k(sortedNodes[n].k);

      double ratio = 0.0;

      if(source.angle == 0.0) {
        double d2p = (coords[k][j][i]-source.x0)*source.dir;
        Vec3D xp = (coords[k][j][i]-source.x0) - d2p*source.dir;
        ratio = xp.norm()/beamwaist;
      }
      else {//focusing or diverging
        Vec3D O2p = coords[k][j][i] - source.O;
        double theta = acos(fabs(O2p*source.dir)/O2p.norm());
        ratio = theta/fabs(source.angle);
      }

      l[k][j][i] = (2.0*power/(PI*beamwaist*beamwaist))*exp(-2.0*ratio*ratio);

    }
    levelcomm[0]->ExchangeAndInsert(l);
    L.RestoreDataPointerToLocalVector();
    coordinates.RestoreDataPointerToLocalVector();

  }

}

//--------------------------------------------------------------------------





