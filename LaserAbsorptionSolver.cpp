#include<LaserAbsorptionSolver.h>
#include<GeoTools.h>
#include<algorithm> //std::sort
#include<numeric> //std::iota
#include<map>

//#include<chrono> //for timing only
//using namespace std::chrono;

using std::vector;
using std::pair;

extern int verbose;

//--------------------------------------------------------------------------

LaserAbsorptionSolver::LaserAbsorptionSolver(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_, 
                         vector<VarFcnBase*> &varFcn_, SpaceVariable3D &coordinates_, 
                         SpaceVariable3D &delta_xyz_, SpaceVariable3D &volume_,
                         vector<GhostPoint> &ghost_nodes_inner_, 
                         vector<GhostPoint> &ghost_nodes_outer_)
                     : comm(comm_), iod(iod_), varFcn(varFcn_), coordinates(coordinates_), delta_xyz(delta_xyz_),
                       volume(volume_), ghost_nodes_inner(ghost_nodes_inner_),
                       ghost_nodes_outer(ghost_nodes_outer_),
                       Temperature(comm_, &(dm_all_.ghosted1_1dof)),
                       L0(comm_, &(dm_all_.ghosted1_1dof)),
                       Lbk(comm_, &(dm_all_.ghosted1_1dof)),
                       Phi(comm_, &(dm_all_.ghosted1_1dof)),
                       Level(comm_, &(dm_all_.ghosted1_1dof)),
                       Tag(comm_, &(dm_all_.ghosted1_1dof))
{
  // Check input parameters
  CheckForInputErrors();

  //MPI info
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  // Mesh info
  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);
  CalculateGlobalMeshInfo();

  // cut-off radiance (L must be positive inside the laser domain)
  lmin = iod.laser.lmin;

  // parameters in mean flux method & SOR
  mfm_alpha = iod.laser.alpha;
  sor_relax = iod.laser.relax_coeff;

  // Get absorption coefficient for each material
  int numMaterials = iod.eqs.materials.dataMap.size();
  absorption.resize(numMaterials, std::make_tuple(0,0,0)); //by default, set coeff = 0
  for (auto it = iod.laser.abs.dataMap.begin(); it != iod.laser.abs.dataMap.end(); it++) {
    if(it->second->materialid < 0 || it->second->materialid >= numMaterials) {
      fprintf(stderr,"ERROR: Found laser absorption coefficients for an unknown material (id = %d).\n", it->first);
      exit_mpi();
    }
    std::get<0>(absorption[it->first]) = it->second->slope;
    std::get<1>(absorption[it->first]) = it->second->T0;
    std::get<2>(absorption[it->first]) = it->second->alpha0;
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

  // Sort nodes in scope (Not all the nodes) by queue-level (primary) and dist-to-source (secondary)
  // Create a custom communicator for each level
  BuildCustomizedCommunicators();
  //VerifySortedNodesAndCommunicators(); //for debug purpose

  // Find ghost nodes outside laser boundary (These include nodes inside and outside the physical domain!)
  SetupLaserGhostNodes();
  
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
  Temperature.Destroy();
  L0.Destroy();
  Lbk.Destroy();
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

  dxmin_glob = DBL_MAX;
  for(int i=ii0; i<iimax; i++)
    dxmin_glob[0] = std::min(dxmin_glob[0], dxyz[kk0][jj0][i][0]);
  for(int j=jj0; j<jjmax; j++)
    dxmin_glob[1] = std::min(dxmin_glob[1], dxyz[kk0][j][ii0][1]);
  for(int k=kk0; k<kkmax; k++)
    dxmin_glob[2] = std::min(dxmin_glob[2], dxyz[k][jj0][ii0][2]);

  MPI_Allreduce(MPI_IN_PLACE, dxmin_glob, 3, MPI_DOUBLE, MPI_MIN, comm);

  xmin_glob       = coords[k0][j0][i0];
  xmax_glob       = coords[kmax-1][jmax-1][imax-1];
  xmin_glob_ghost = coords[kk0][jj0][ii0];
  xmax_glob_ghost = coords[kkmax-1][jjmax-1][iimax-1];
 
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

  delete [] filename;
}

//--------------------------------------------------------------------------

void
LaserAbsorptionSolver::CalculateLaserInfo()
{

  LaserData &laser(iod.laser);

  laser_range = laser.range;

  source.radius = laser.source_radius;
  source.x0     = Vec3D(laser.source_center_x, laser.source_center_y, laser.source_center_z);
  source.dir    = Vec3D(laser.source_dir_x, laser.source_dir_y, laser.source_dir_z);
  source.dir   /= source.dir.norm(); //normalize

  double PI = 2.0*acos(0.0);
  source.angle = laser.focusing_angle_degrees/180.0*PI; //in radians

  if(source.angle == 0) {//parallel laser
    source.R = DBL_MAX;
    source.O = Vec3D(DBL_MAX,DBL_MAX,DBL_MAX);
  }
  else if(source.angle>0) {//focusing laser 
    source.R = source.radius/sin(source.angle); //radius of curvature
    double eps = std::max(dxmin_glob[0], std::max(dxmin_glob[1], dxmin_glob[2]));
    if(laser_range > source.R - 2.0*eps) {
      print_error("*** Error: Laser range must not include the focal point. range = %e, R = %e, hmin = %e. (Stay away by at least 2 cells.)\n",
                   laser_range, source.R, eps);
      exit_mpi();
    }
    source.O = source.x0 + source.R*cos(source.angle)*source.dir; //origin of the circle/curve 
  }
  else if(source.angle<0) {//diverging laser
    double theta = -source.angle; //a positive angle
    source.R = source.radius/sin(theta); //radius of curvature
    source.O = source.x0 - source.R*cos(theta)*source.dir; //origin
  }

  print("- Laser radiation solver activated. Focusing angle: %f degrees; Source radius: %f; Range: %f\n",
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

/*
  Phi.WriteToVTRFile("phi.vtr");
*/
}

//--------------------------------------------------------------------------

double
LaserAbsorptionSolver::DistanceToSource(Vec3D &x, int &loc, double *image)
{
  //each digit of loc corresponds to a check
  //if "image" is not NULL, and if x is out of scope, "image" will get the coordinates
  //of the image point. (This function assumes that x is somewhat close to
  //the laser scope --- not an arbitrary point in R3.)
  
  double dist = DBL_MAX;
  loc = 0; //default: inside laser scope

  if(source.angle == 0) { //parallel laser

    double d2p  = (x - source.x0)*source.dir; //distance to source plane
    if(d2p<0)            loc += 1;
    if(d2p>laser_range)  loc += 10;
        
    Vec3D xp = x - source.x0 - d2p*source.dir; 
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
    Vec3D xp = source.O - r*source.dir;
    if(r<0)
      loc += 1000; //wrong side! this function is not capable of calculating its image
    else
      if((x - xp).norm() > r*tan(source.angle))
        loc += 100;
    

    dist = source.R - d2O;

    if(image && loc) { //calculate image coordinates

      if(loc>=1000) {
        fprintf(stderr,"*** Error: Not able to find the image of a node on the wrong side of the focus (shouldn't need to).\n");
        exit(-1);
      }

      Vec3D xim = x;
      if(loc == 1 || loc == 101)
        xim -= 2.0*(d2O-source.R)/d2O*Ox;
      if(loc == 10 || loc == 110)
        xim += 2.0*(source.R-laser_range-d2O)/d2O*Ox;
      if(loc >= 100) { //100, 101, 110
        Vec3D Oxim = xim - source.O;
        double d2Oxim = Oxim.norm();
        double theta = acos(fabs(Oxim*source.dir)/d2Oxim);
        assert(theta>source.angle);
        theta = 2.0*source.angle - theta;
        Vec3D ximp = source.O - d2Oxim*cos(theta)*source.dir;
        Vec3D vert_dir = (x-xp)/((x-xp).norm());
        xim = ximp + d2Oxim*sin(theta)*vert_dir;
      }
      for(int i=0; i<3; i++)
        image[i] = xim[i];
    }
  }
  else if(source.angle < 0) { //diverging laser

    Vec3D Ox = x - source.O;
    double d2O = Ox.norm();

    if(d2O<source.R)               loc += 1;
    if(d2O>(source.R+laser_range)) loc += 10;

    double r = Ox*source.dir;
    Vec3D xp = source.O + r*source.dir;

    if(r<0) 
      loc += 1000; //wrong side! this function is not capable of calculating its image
    else if((x - xp).norm() > r*tan(-source.angle))  
      loc += 100;

    dist = d2O - source.R;

    if(image && loc) { //calculate image coordinates

      if(loc>=1000) {
        fprintf(stderr,"*** Error: Not able to find the image of a node on the wrong side of the origin (shouldn't need to).\n");
        exit(-1);
      }

      Vec3D xim = x;
      if(loc == 1 || loc == 101)
        xim += 2.0*(source.R-d2O)/d2O*Ox;
      if(loc == 10 || loc == 110)
        xim -= 2.0*(d2O-source.R-laser_range)/d2O*Ox;
      if(loc >= 100) { //100, 101, 110
        Vec3D Oxim = xim - source.O;
        double d2Oxim = Oxim.norm();
        double theta = acos(fabs(Oxim*source.dir)/d2Oxim); 
        assert(theta>(-source.angle));
        theta = 2.0*(-source.angle) - theta;
        Vec3D ximp = source.O + d2Oxim*cos(theta)*source.dir;
        Vec3D vert_dir = (x-xp)/((x-xp).norm());
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

      // If a neighbor is tagged (meaning on a previous level), set this node
      // (Because sortedNodes are all within subdomain interior, these neighbors
      //  are accessible from this subdomain.)
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


/*
  //DEBUG
  Level.WriteToVTRFile("Level.vtr");
  tag = Tag.GetDataPointer();
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++)
        tag[k][j][i] = mpi_rank;
  Tag.RestoreDataPointerAndInsert();
  Tag.WriteToVTRFile("Partition.vtr");
*/ 
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

  if(mpi_size>1)
    print("- Creating a customized communicator for the laser radiation solver.\n");

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
  for(int lvl = 0; lvl<queueCounter.size(); lvl++) 
    levelcomm.push_back(new CustomCommunicator(comm, L0, nodes_on_level[lvl]));

}

//--------------------------------------------------------------------------

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
  //geometry tolerance
  double eps = std::min(dxmin_glob[0], std::min(dxmin_glob[1], dxmin_glob[2]));
  eps /= 1.0e10;

  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();

  //----------------------------------------------------------------
  // Step 1: Find ghost nodes (inside and outside physical domain)
  //----------------------------------------------------------------
  Tag.SetConstantValue(0, true);
  double*** level = Level.GetDataPointer(); //level only tags nodes inside subdomain
  double*** tag   = Tag.GetDataPointer(); 
  //loop through sorted nodes on level 1, 2, 3, ... (Skipping level 0 as it does not need to be calculated)
  for(int n = queueCounter[0]; n < sortedNodes.size(); n++) {
    int i(sortedNodes[n].i), j(sortedNodes[n].j), k(sortedNodes[n].k); //all sorted nodes are inside subdomain
    // Note: because of the Cartesian grid and the assumption that laser boundaries are straight lines, we only
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
  vector<double> ghost_phi;
  // status: 
  // (1) If the ghost node is outside laser domain,
  //     let Q be the closest point on the boundary of the laser domain
  // (1.1) if Q is on the source disk/sector                ---> status = 1
  // (1.2)       ...  the disk/sector of max. range         ---> status = 2
  // (1.3)       ...  the side boundary of the laser domain ---> status = 3
  // (1.4)       ...  an edge of the laser domain boundary  ---> status = 4
  // (2) If the ghost node is inside laser domain, (the ghost node must be
  //                                           outside of the "physical domain")
  // (2.1) if it is next to a symmetry boundary (i.e. should be populated by mirroring) ---> status  = -1
  //       if it is "downstream" (i.e. should be populated by linear extrap.) --> report error!
  //
  // (status = -99 and -98) reserved for special cases (see below).
  //
  // Note: Althought status 1 and 2 are tracked here. In the embedded boundary method for laser absorption,
  //       we don't really need to deal with these two cases, as long as (1) the scheme is strictly upwinding
  //       and (2) the source depth covers at least one (or 2?) layers of nodes.
  //
  for(auto it = ghost_tmp.begin(); it != ghost_tmp.end(); it++) {
    int i((*it)[0]), j((*it)[1]), k((*it)[2]);
    Vec3D& x(coords[k][j][i]);

    int location_info;
    Vec3D im;
    double dist = DistanceToSource(x, location_info, im);
    ghost_phi.push_back(dist);

    if(location_info) {//Category (1)
      image.push_back(im); //inserting the coordinates
      if(location_info == 1)        status.push_back(1);
      else if(location_info == 10)  status.push_back(2);
      else if(location_info == 100) status.push_back(3);
      else                          status.push_back(4);
    }
    else { //Category (2)
      assert(coordinates.BoundaryType(i,j,k) == 1); //boundary-face
      bool found = false;
      for(auto it = ghost_nodes_outer.begin(); it != ghost_nodes_outer.end(); it++) {
        if(it->ijk[0] == i && it->ijk[1] == j && it->ijk[2] == k) {
        
          // verify that it is next to a symmetry boundary (of the laser)
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
/*
  for(int i=0; i<ghost_tmp.size(); i++) {
    Vec3D& x(coords[ghost_tmp[i][2]][ghost_tmp[i][1]][ghost_tmp[i][0]]);
    int tmp;
    double dtmp = DistanceToSource(image[i], tmp);
    fprintf(stderr,"[%d] (%d,%d,%d)(%e,%e,%e)(%e) --> status = %d, (%e,%e,%e)(%e).\n", mpi_rank, ghost_tmp[i][0], ghost_tmp[i][1],
            ghost_tmp[i][2], x[0], x[1], x[2], ghost_phi[i], status[i], image[i][0], image[i][1], image[i][2], dtmp);
  }
*/
  //----------------------------------------------------------------
  // Step 3: Setup trilinear interpolation coeffs for each mirror image
  //----------------------------------------------------------------
  double tol = eps/5.0; //needed at least for 2D sims, where z-direction has only 1 node
  vector<Int3> image_ijk(image.size(), Int3(INT_MAX)); //lower-left corner of the interpolation box
  vector<Vec3D> image_xi(image.size(), Vec3D(0.0));
  vector<int> prob_nodes; //nodes whose image is outside of the computational domain
  vector<int> far_nodes; //nodes whose image is inside a different subdomain
  level = Level.GetDataPointer(); //level only tags nodes inside subdomain
  for(int n=0; n<image.size(); n++) {
    if(status[n] == -1) {
      image_ijk[n] = Int3((int)image[n][0], (int)image[n][1], (int)image[n][2]);
      assert(level[image_ijk[n][2]][image_ijk[n][1]][image_ijk[n][0]] >= 0);
    }
    else {
      // Here, we try to find the location of the image w.r.t. the primal mesh (defined by nodes and edges). Note that
      // the subdomains have gaps in between, which are covered by the (internal) ghost layer.
      // The image can be located (1) inside the current subdomain or its internal ghost layer, (2) inside another subdomain, 
      // (3) outside the domain (defined by the primal mesh)
      //
      if(image[n][0]<xmin_glob[0]-tol || image[n][0]>=xmax_glob[0]+tol ||
         image[n][1]<xmin_glob[1]-tol || image[n][1]>=xmax_glob[1]+tol ||
         image[n][2]<xmin_glob[2]-tol || image[n][2]>=xmax_glob[2]+tol) { // case (3)
        status[n] = -99; //change its status to -99
        prob_nodes.push_back(n);
        continue; //deal with all the problematic nodes later
      }

      int ihere(INT_MAX), jhere(INT_MAX), khere(INT_MAX);
      for(int i=ii0; i<iimax-1; i++) 
        if(image[n][0]>=coords[kk0][jj0][i][0]-tol && image[n][0]<coords[kk0][jj0][i+1][0]+tol) {ihere = i; break;} 
      for(int j=jj0; j<jjmax-1; j++) 
        if(image[n][1]>=coords[kk0][j][ii0][1]-tol && image[n][1]<coords[kk0][j+1][ii0][1]+tol) {jhere = j; break;} 
      for(int k=kk0; k<kkmax-1; k++) 
        if(image[n][2]>=coords[k][jj0][ii0][2]-tol && image[n][2]<coords[k+1][jj0][ii0][2]+tol) {khere = k; break;} 

      if(ihere==INT_MAX || jhere==INT_MAX || khere==INT_MAX) {//not found inside the subdomain or its ghost layer, case (2)
        status[n] = -98; //change its status to -98
        far_nodes.push_back(n);
//        fprintf(stderr,"[%d] Found a far-node. (%d,%d,%d). image(%e,%e,%e). (%d,%d,%d).\n", mpi_rank, 
//                ghost_tmp[n][0], ghost_tmp[n][1], ghost_tmp[n][2], image[n][0], image[n][1], image[n][2], ihere, jhere, khere);
        continue;  //deal with all the problematic and "far" nodes later
      }

      // case (1)
      image_ijk[n] = Int3(ihere,jhere,khere);
      for(int i=0; i<3; i++)
        image_xi[n][i] = (image[n][i] - coords[khere][jhere][ihere][i])
                       / (coords[khere+1][jhere+1][ihere+1][i] - coords[khere][jhere][ihere][i]);
    }
  }
  Level.RestoreDataPointerToLocalVector();

  // Now, work on the "far" nodes (image is inside another subdomain)
  int total_far_nodes = far_nodes.size();
  MPI_Allreduce(MPI_IN_PLACE, &total_far_nodes, 1, MPI_INT, MPI_SUM, comm);
  if(total_far_nodes>0) {
    for(int proc=0; proc<mpi_size; proc++) { 

      // 1. proc sends far nodes (ijk) and their images (coords) to everyone else
      int local_size = 0;
      if(mpi_rank == proc)
        local_size = far_nodes.size();
      MPI_Bcast(&local_size, 1, MPI_INT, proc, comm);

      if(local_size == 0) //proc does not have "far" nodes
        continue;

      int*    nod = new int[3*local_size];
      double* img = new double[3*local_size];
      if(mpi_rank == proc)
        for(int i = 0; i < far_nodes.size(); i++) {
          int n = far_nodes[i];
          for(int j = 0; j < 3; j++) {
            nod[3*i+j] = ghost_tmp[n][j];
            img[3*i+j] = image[n][j];
          }
        }
      MPI_Bcast(nod, 3*local_size, MPI_INT, proc, comm);
      MPI_Bcast(img, 3*local_size, MPI_DOUBLE, proc, comm);
       
      // 2. everyone else figures out if it owns some of these nodes, and inform proc
      vector<int> far_nodes_owner(local_size, INT_MAX);
      vector<Int3> far_nodes_image_ijk(local_size, INT_MAX);
      if(mpi_rank != proc) {
        for(int n=0; n<local_size; n++) {

/*
          if(proc==17 && nod[3*n]==51 && nod[3*n+1]==88 && nod[3*n+2]==0) {
            fprintf(stderr,"[%d] Working on (51,88,0) from proc 17... Image(%e,%e,%e)\n", mpi_rank,
                    img[3*n],img[3*n+1],img[3*n+2]);
          }
*/
          int ihere(INT_MAX), jhere(INT_MAX), khere(INT_MAX);
          for(int i=ii0; i<iimax-1; i++) 
            if(img[3*n]>=coords[kk0][jj0][i][0]-tol && img[3*n]<coords[kk0][jj0][i+1][0]+tol) {ihere = i; break;} 
          for(int j=jj0; j<jjmax-1; j++) 
            if(img[3*n+1]>=coords[kk0][j][ii0][1]-tol && img[3*n+1]<coords[kk0][j+1][ii0][1]+tol) {jhere = j; break;} 
          for(int k=kk0; k<kkmax-1; k++) 
            if(img[3*n+2]>=coords[k][jj0][ii0][2]-tol && img[3*n+2]<coords[k+1][jj0][ii0][2]+tol) {khere = k; break;} 

/*
          if(proc==17 && nod[3*n]==51 && nod[3*n+1]==88 && nod[3*n+2]==0) {
            fprintf(stderr,"[%d] Working on (51,88,0) from proc 17... (%d,%d,%d)\n", mpi_rank, ihere, jhere, khere);
          }
*/

          if(ihere<INT_MAX && jhere<INT_MAX && khere<INT_MAX) {//inside this subdomain
            far_nodes_owner[n] = mpi_rank;
            far_nodes_image_ijk[n] = Int3(ihere,jhere,khere);
          } 
        }
      }

      MPI_Allreduce(MPI_IN_PLACE, far_nodes_owner.data(), far_nodes_owner.size(), MPI_INT, MPI_MIN, comm);

      // 3. proc digests the received data
      if(mpi_rank == proc) { 
        std::map<int,int> rank_to_id;
        for(int n=0; n<far_nodes_owner.size(); n++) {
          int owner_rank = far_nodes_owner[n];
/*
          if(owner_rank<0 || owner_rank>=mpi_size) {
            fprintf(stderr,"[%d] Can't find an owner of far-node (%d,%d,%d)\n", mpi_rank, 
                    ghost_tmp[far_nodes[n]][0], ghost_tmp[far_nodes[n]][1], ghost_tmp[far_nodes[n]][2]);
          }
*/
          assert(owner_rank>=0 && owner_rank<mpi_size);
          if(rank_to_id.find(owner_rank) == rank_to_id.end()) {
            ebm.ghostNodes2.push_back(vector<Int3>());
            ebm.ghostNodes2_sender.push_back(owner_rank);
            rank_to_id[owner_rank] = ebm.ghostNodes2.size() - 1;
          }
          int id = rank_to_id[owner_rank];
          ebm.ghostNodes2[id].push_back(ghost_tmp[far_nodes[n]]);
        }
      }

      // 4. others build friendsGhostNodes
      if(mpi_rank != proc) {
        int id = -1;
        for(int n=0; n<far_nodes_owner.size(); n++) {
          if(far_nodes_owner[n] == mpi_rank) {//mine!
            if(id==-1) {
              ebm.friendsGhostNodes.push_back(vector<std::pair<Int3, EmbeddedBoundaryFormula> >() );  
              ebm.friendsGhostNodes_receiver.push_back(proc);
              id = ebm.friendsGhostNodes.size() - 1;
            }
            Int3 nod_int3(nod[3*n],nod[3*n+1],nod[3*n+2]);
            ebm.friendsGhostNodes[id].push_back(std::make_pair(nod_int3, EmbeddedBoundaryFormula(eps)));
            Vec3D xi;
            int ihere = far_nodes_image_ijk[n][0];
            int jhere = far_nodes_image_ijk[n][1];
            int khere = far_nodes_image_ijk[n][2];
            for(int i=0; i<3; i++)
              xi[i] = (img[3*n+i] - coords[khere][jhere][ihere][i])
                    / (coords[khere+1][jhere+1][ihere+1][i] - coords[khere][jhere][ihere][i]);
            ebm.friendsGhostNodes[id].back().second.BuildMirroringFormula(nod_int3, far_nodes_image_ijk[n], xi);
          } 
        }
      }

      delete [] nod;
      delete [] img;
    }

    //verification
    int counter = 0;
    for(int i=0; i<ebm.ghostNodes2.size(); i++)
      counter += ebm.ghostNodes2[i].size();
    assert(counter == far_nodes.size());

  }


/*
  for(int i=0; i<ebm.ghostNodes2.size(); i++)
    for(int j=0; j<ebm.ghostNodes2[i].size(); j++)
      fprintf(stderr,"[%d] GN2 (%d, %d, %d) owner = %d.\n", mpi_rank, ebm.ghostNodes2[i][j][0], ebm.ghostNodes2[i][j][1],
              ebm.ghostNodes2[i][j][2], ebm.ghostNodes2_sender[i]);

  for(int i=0; i<ebm.friendsGhostNodes.size(); i++)
    for(int j=0; j<ebm.friendsGhostNodes[i].size(); j++)
      fprintf(stderr,"[%d] FG (%d, %d, %d), friend = %d.\n", mpi_rank, ebm.friendsGhostNodes[i][j].first[0],
              ebm.friendsGhostNodes[i][j].first[1], ebm.friendsGhostNodes[i][j].first[2], ebm.friendsGhostNodes_receiver[i]);
*/


  // Now, work on the problematic nodes (KW believes this should rarely happen, and if it does, it must
  // be a 2D (or 2D-cylindrical) simulation. 
  // We only try to find the closest valid point within the SAME subdomain, instead of doing another round of MPI
  // communications like what we did for "far_nodes". This means we may get slightly different results when running
  // the same case using different mesh partitions --- but the difference should be minimal.
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
        fprintf(stderr,"*** Error: [%d] Cannot find a valid image for ghost node (%d,%d,%d)(%e,%e,%e). "
                       "Original image: (%e,%e,%e).\n", mpi_rank, ghost_tmp[n][0], ghost_tmp[n][1], ghost_tmp[n][2],
                       coords[ghost_tmp[n][2]][ghost_tmp[n][1]][ghost_tmp[n][0]][0],
                       coords[ghost_tmp[n][2]][ghost_tmp[n][1]][ghost_tmp[n][0]][1],
                       coords[ghost_tmp[n][2]][ghost_tmp[n][1]][ghost_tmp[n][0]][2],
                       image[n][0], image[n][1], image[n][2]);
        exit(-1);
      }
      image_ijk[n] = closest[p].first;
      image_xi[n]  = Vec3D(0,0,0);

/*
      fprintf(stderr,"[%d] Prob node (%d %d %d), image(%e,%e,%e), closest (%d %d %d).\n", mpi_rank, ghost_tmp[n][0], 
              ghost_tmp[n][1], ghost_tmp[n][2], image[n][0], image[n][1], image[n][2], image_ijk[n][0], image_ijk[n][1],
              image_ijk[n][2]);
*/
    }
    Phi.RestoreDataPointerToLocalVector();

  }

  coordinates.RestoreDataPointerToLocalVector();
  

  //----------------------------------------------------------------
  // Step 4: Sort the list of ghost nodes based on distance to source
  //         (ghost_phi).
  //----------------------------------------------------------------
  vector<int> ordered(ghost_phi.size());
  std::iota(std::begin(ordered), std::end(ordered), 0); //ordered = {0, 1, 2, 3, ...}
  std::sort(std::begin(ordered), std::end(ordered), 
            [&ghost_phi](int lhs, int rhs) {return ghost_phi[lhs]<ghost_phi[rhs];} );  

  //----------------------------------------------------------------
  // Step 5: Figure out formula for ghost node update
  //----------------------------------------------------------------
  ebm.ghostNodes1.resize(ordered.size() - far_nodes.size(), std::make_pair(Int3(0), EmbeddedBoundaryFormula(eps)));
  int counter = 0;
  for(int o = 0; o < ordered.size(); o++) {
    int n = ordered[o]; 

    if(status[n] == -98) //added to ghostNodes2 already
      continue;

    ebm.ghostNodes1[counter].first = ghost_tmp[n];
    if(status[n] == -1 || status[n] == -99) {
      vector<Int3> node(1,image_ijk[n]);
      vector<double> coeff(1,1.0);
      ebm.ghostNodes1[counter].second.SpecifyFormula(EmbeddedBoundaryFormula::MIRRORING, 
                                                     EmbeddedBoundaryFormula::NODE, node, coeff);
    } else
      ebm.ghostNodes1[counter].second.BuildMirroringFormula(ghost_tmp[n], image_ijk[n], image_xi[n]);

    counter++;
  }
  assert(counter == ordered.size() - far_nodes.size());


  // Allocate space for storing laser radiance at ghost nodes
  ebm.l1.resize(ebm.ghostNodes1.size(), 0.0);
  ebm.l2.resize(ebm.ghostNodes2.size());
  for(int n1=0; n1<ebm.ghostNodes2.size(); n1++)
    ebm.l2[n1].resize(ebm.ghostNodes2[n1].size(), 0.0);
 
}

//--------------------------------------------------------------------------

void
LaserAbsorptionSolver::SetSourceRadiance(double*** l, Vec3D*** coords, double t)
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

    if(l0<lmin) {
      l0 = lmin;
      if(verbose >= OutputData::HIGH)
        print("Warning: source radiance (%e) is lower than the cut-off radiance (%e).\n", l0, lmin);
    }

    for(int n=0; n<queueCounter[0]; n++) //level 0
      l[sortedNodes[n].k][sortedNodes[n].j][sortedNodes[n].i] = l0;
    levelcomm[0]->ExchangeAndInsert(l);

  }
  else if(iod.laser.source_distribution == LaserData::GAUSSIAN){ //The Laser source is Gaussia

    //give the direction of X
    double beamwaist = iod.laser.source_beam_waist;
    double PI = 2.0*acos(0.0);

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

      if(l[k][j][i]<lmin) {
        l[k][j][i] = lmin;
        if(verbose >= OutputData::HIGH)
          fprintf(stderr,"Warning: [%d] source radiance at (%d,%d,%d)(%e) is below cut-off(%e).\n", 
                  mpi_rank, i,j,k, l[k][j][i], lmin);
      }

    }
    levelcomm[0]->ExchangeAndInsert(l);

  }

}

//--------------------------------------------------------------------------

void
LaserAbsorptionSolver::InitializeLaserDomainAndGhostNodes(double*** l, double*** level)
{
  //Assuming level 0 (i.e. nodes wihin the depth of the source) has been updated.

  //---------------------------------------------------------------------
  // Step 1. Specify a value of L for nodes inside the laser domain, layer by layer.
  //---------------------------------------------------------------------
  auto it = sortedNodes.begin() + queueCounter[0];
  for(int lvl = 1; lvl<queueCounter.size(); lvl++) {
    for(int n = 0; n < queueCounter[lvl]; n++) {

      int i(it->i), j(it->j), k(it->k);
      double sum(0);
      int count(0);
      if(level[k][j][i-1]>=0 && level[k][j][i-1]<lvl) {sum += l[k][j][i-1]; count++;}
      if(level[k][j][i+1]>=0 && level[k][j][i+1]<lvl) {sum += l[k][j][i+1]; count++;}
      if(level[k][j-1][i]>=0 && level[k][j-1][i]<lvl) {sum += l[k][j-1][i]; count++;}
      if(level[k][j+1][i]>=0 && level[k][j+1][i]<lvl) {sum += l[k][j+1][i]; count++;}
      if(level[k-1][j][i]>=0 && level[k-1][j][i]<lvl) {sum += l[k-1][j][i]; count++;}
      if(level[k+1][j][i]>=0 && level[k+1][j][i]<lvl) {sum += l[k+1][j][i]; count++;}
      assert(count>0);
      l[k][j][i] = sum/count;

      if(l[k][j][i]<lmin) {
        l[k][j][i] = lmin;
        if(verbose >= OutputData::MEDIUM)
          fprintf(stderr,"Warning: [%d] Initial radiance at (%d,%d,%d)(%e) is below cut-off(%e).\n", 
                  mpi_rank, i,j,k, l[k][j][i], lmin);
      }

      it++;
    }
    levelcomm[lvl]->ExchangeAndInsert(l);
  }


  //---------------------------------------------------------------------
  // Step 2. Update ghost nodes
  //---------------------------------------------------------------------
  UpdateGhostNodes(l);

}

//--------------------------------------------------------------------------

void
LaserAbsorptionSolver::UpdateGhostNodes(double ***l, int custom_max_iter)
{

  bool custom_max_iter_specified = custom_max_iter>0;
  int max_iter = custom_max_iter_specified ? custom_max_iter : iod.laser.max_iter;

  // Store "old" values. (These vectors will be small.)
  for(int n=0; n<ebm.ghostNodes1.size(); n++) {
    int i(ebm.ghostNodes1[n].first[0]), j(ebm.ghostNodes1[n].first[1]), k(ebm.ghostNodes1[n].first[2]);
    ebm.l1[n] = l[k][j][i];
  }
  for(int n1=0; n1<ebm.ghostNodes2.size(); n1++) {
    for(int n=0; n<ebm.ghostNodes2[n1].size(); n++) {
      int i(ebm.ghostNodes2[n1][n][0]), j(ebm.ghostNodes2[n1][n][1]), k(ebm.ghostNodes2[n1][n][2]);
      ebm.l2[n1][n] = l[k][j][i];
    }
  }

  // Update ghost nodes iteratively
  int iter;
  double max_err;
  for(iter = 0; iter < max_iter; iter++) {

    //---------------------------------
    UpdateGhostNodesOneIteration(l);    
    //---------------------------------

    // Find local max error
    max_err = 0.0;
    for(int n=0; n<ebm.ghostNodes1.size(); n++) {
      int i(ebm.ghostNodes1[n].first[0]), j(ebm.ghostNodes1[n].first[1]), k(ebm.ghostNodes1[n].first[2]);
      double lnew = l[k][j][i];
      max_err = std::max(max_err, fabs((ebm.l1[n] - lnew)/lnew));
      ebm.l1[n] = lnew;
    }
    for(int n1=0; n1<ebm.ghostNodes2.size(); n1++)
      for(int n=0; n<ebm.ghostNodes2[n1].size(); n++) {
        int i(ebm.ghostNodes2[n1][n][0]), j(ebm.ghostNodes2[n1][n][1]), k(ebm.ghostNodes2[n1][n][2]);
        double lnew = l[k][j][i];
        max_err = std::max(max_err, fabs((ebm.l2[n1][n] - lnew)/lnew));
        ebm.l2[n1][n] = lnew;
      }

    // Find global max error and compare with tolerance
    MPI_Allreduce(MPI_IN_PLACE, &max_err, 1, MPI_DOUBLE, MPI_MAX, comm);
    //print("It %d: max_err = %e!\n", iter, max_err);
    if(max_err<iod.laser.convergence_tol)
      break;

  }

  if(!custom_max_iter_specified) { //no need to print anything if max_iter is "customized"
    if(iter == max_iter && verbose >= OutputData::HIGH)
      print("Warning: Ghost nodes in the laser solver didn't converge in %d iterations "
            "(err = %e, tol = %e).\n", iter, max_err, iod.laser.convergence_tol);
    else if(verbose >= OutputData::HIGH)
      print("  o Ghost nodes in the laser solver converged in %d iteration(s) "
            "(err = %e, tol = %e).\n", iter+1, max_err, iod.laser.convergence_tol);
  }

}

//--------------------------------------------------------------------------

void
LaserAbsorptionSolver::UpdateGhostNodesOneIteration(double ***l)
{

  // Step 1: Work on ghostNodes1
  for(auto it = ebm.ghostNodes1.begin(); it != ebm.ghostNodes1.end(); it++) {
    int i(it->first[0]), j(it->first[1]), k(it->first[2]);
    l[k][j][i] = it->second.Evaluate(l);

    if(l[k][j][i]<lmin) {
      l[k][j][i] = lmin;
      if(verbose >= OutputData::MEDIUM)
        fprintf(stderr,"Warning: [%d] Radiance at ghost (%d,%d,%d)(%e) is below cut-off(%e).\n", 
                mpi_rank, i,j,k, l[k][j][i], lmin);
    }
//    fprintf(stderr,"[%d] GN1 l[%d][%d][%d] = %e.\n", mpi_rank, k,j,i, l[k][j][i]);
  }

  // Step 2: Work on friendsGhostNodes
  int nSend = ebm.friendsGhostNodes.size();
  vector<MPI_Request> send_requests(nSend, MPI_Request());
  for(int p=0; p<nSend; p++) {
    int receiver = ebm.friendsGhostNodes_receiver[p]; 
    vector<double> buffer(ebm.friendsGhostNodes[p].size());
    for(int i=0; i<ebm.friendsGhostNodes[p].size(); i++) {
      buffer[i] = ebm.friendsGhostNodes[p][i].second.Evaluate(l);

      if(buffer[i]<lmin) {
        buffer[i] = lmin;
        if(verbose >= OutputData::MEDIUM)
          fprintf(stderr,"Warning: [%d] Radiance at friend's ghost (%d,%d,%d)(%e) is below cut-off(%e).\n", 
                  mpi_rank, ebm.friendsGhostNodes[p][i].first[0], ebm.friendsGhostNodes[p][i].first[1],
                  ebm.friendsGhostNodes[p][i].first[2], buffer[i], lmin);
      }
    }
    MPI_Isend(buffer.data(), buffer.size(), MPI_DOUBLE, receiver, mpi_rank, comm, &(send_requests[p]));
  }

  // Step 3: Receiving data for ghostNodes2
  int nRecv = ebm.ghostNodes2.size(); 
  vector<vector<double> > buffers(nRecv, vector<double>());
  vector<MPI_Request> recv_requests(nRecv, MPI_Request());
  for(int p=0; p<nRecv; p++) {
    int sender = ebm.ghostNodes2_sender[p];
    buffers[p].resize(ebm.ghostNodes2[p].size());
    MPI_Irecv(buffers[p].data(), buffers[p].size(), MPI_DOUBLE, sender, sender, comm, &(recv_requests[p]));
  }

  // Step 4: Wait.
  MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE); //not necessary
  MPI_Waitall(recv_requests.size(), recv_requests.data(), MPI_STATUSES_IGNORE);

  // Step 5: Apply data to ghostNodes2
  for(int p=0; p<nRecv; p++) {
    for(int q=0; q<ebm.ghostNodes2[p].size(); q++) {
      int i(ebm.ghostNodes2[p][q][0]), j(ebm.ghostNodes2[p][q][1]), k(ebm.ghostNodes2[p][q][2]);
      l[k][j][i] = buffers[p][q];
      if(l[k][j][i]<lmin) {
        l[k][j][i] = lmin;
        if(verbose >= OutputData::MEDIUM)
          fprintf(stderr,"Warning: [%d] Received radiance at ghost (%d,%d,%d)(%e) is below cut-off(%e).\n", 
                  mpi_rank, i,j,k, l[k][j][i], lmin);
      }
//      fprintf(stderr,"[%d] GN2 l[%d][%d][%d] = %e.\n", mpi_rank, k,j,i, l[k][j][i]);
    }
  }

}

//--------------------------------------------------------------------------

void
LaserAbsorptionSolver::ComputeLaserRadiance(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &L,
                                            const double t, bool initialized)
{
  bool success;

  success = ComputeLaserRadianceMeanFluxMethod(V,ID,L,t,mfm_alpha,sor_relax,initialized);

  if(!success) {// triger the "failsafe" procedure
    // Method 1: gradually decreasing the relaxation_factor in SOR
    for(int iter = 0; iter < 5; iter++) {
      sor_relax /= 2.0;
      success = ComputeLaserRadianceMeanFluxMethod(V,ID,L,t,mfm_alpha,sor_relax,initialized); 
      if(success)
        break;
    }
  }

  if(!success) {
    //Method 2: switch to the "Step method"
    sor_relax = iod.laser.relax_coeff; //restore user input
    mfm_alpha = 1.0;
    for(int iter = 0; iter < 5; iter++) {
      success = ComputeLaserRadianceMeanFluxMethod(V,ID,L,t,mfm_alpha,sor_relax,initialized); 
      if(success)
        break;
      sor_relax /= 2.0;
    }
  }

  if(!success)
    print_error("*** Error: Laser radiance solver failed to converge.\n");
    
}

//--------------------------------------------------------------------------

bool
LaserAbsorptionSolver::ComputeLaserRadianceMeanFluxMethod(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &L,
                                                          const double t, double alpha, double relax_coeff,
                                                          bool initialized)
{
  // This function will use customized communicators. The sub-functions deal with double*** directly
  //
  double*** l     = L.GetDataPointer();
  double*** T     = Temperature.GetDataPointer();
  double*** l0    = L0.GetDataPointer();
  double*** lbk   = Lbk.GetDataPointer();
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
  Vec3D*** dxyz   = (Vec3D***)delta_xyz.GetDataPointer();
  double*** vol   = volume.GetDataPointer();
  double*** id    = ID.GetDataPointer();
  double*** level = Level.GetDataPointer();
  double*** phi   = Phi.GetDataPointer();
  Vec5D*** v      = (Vec5D***)V.GetDataPointer();


  bool success = true;
  CopyValues(l, lbk); //if the solver fails, return the original 

  // Note that when L enters this function, its laser ghost nodes (i.e. ghostNodes1 and 2) do not have values.

  ComputeTemperatureInLaserDomain(v, id, T);

  //If not initialized or source_power changes, set source radiance; initialize laser if needed
  if(!initialized || !source_power_timehistory.empty())
    SetSourceRadiance(l, coords, t);

  //Make sure both real and ghost laser nodes have some meaningful values
  if(!initialized)
    InitializeLaserDomainAndGhostNodes(l, level);
  else
    ApplyStoredGhostNodesRadiance(l);

  //Gauss-Seidel iterations 
  int GSiter = 0;
  double max_error(DBL_MAX), avg_error(DBL_MAX);
  for(GSiter=0; GSiter<iod.laser.max_iter; GSiter++) {

    CopyValues(l, l0); //l0 = l

    //-----------------------------------------------------------------
    RunMeanFluxMethodOneIteration(l, T, coords, dxyz, vol, id, level, phi, 
                                  alpha, relax_coeff);
    //-----------------------------------------------------------------

    ComputeErrorsInLaserDomain(l0, l, max_error, avg_error);
    //print("It %d. max_error = %e, avg_error = %e!\n", GSiter, max_error, avg_error);

    if(max_error<iod.laser.convergence_tol)
      break;

    UpdateGhostNodes(l);

  }

  if(GSiter == iod.laser.max_iter) {
    print("Warning: Laser radiation solver failed to converge in %d iterations. (Error = %e, Tol = %e) Re-trying.\n",
          GSiter, max_error, iod.laser.convergence_tol);
    CopyValues(lbk, l); //l = lbk
    success = false;

  } else if(verbose >= OutputData::MEDIUM)
    print("- Laser radiation solver converged in %d iteration(s)."
          "(Error = %e, Tol = %e).\n", GSiter+1, max_error, iod.laser.convergence_tol);


  // clean up ghost nodes (for outputting) and store ghost nodes values internally
  CleanUpGhostNodes(l, true);


  L.RestoreDataPointerAndInsert(); //data exchanged using levelcomm. But need to "insert" for outputting
  Temperature.RestoreDataPointerToLocalVector();
  L0.RestoreDataPointerToLocalVector();
  Lbk.RestoreDataPointerToLocalVector();
  coordinates.RestoreDataPointerToLocalVector();
  delta_xyz.RestoreDataPointerToLocalVector();
  volume.RestoreDataPointerToLocalVector(); 
  ID.RestoreDataPointerToLocalVector();
  Level.RestoreDataPointerToLocalVector();
  Phi.RestoreDataPointerToLocalVector();
  V.RestoreDataPointerToLocalVector();

  return success;
}


//--------------------------------------------------------------------------

void
LaserAbsorptionSolver::RunMeanFluxMethodOneIteration(double*** l, double*** T, Vec3D*** coords, 
                                                     Vec3D*** dxyz, double*** vol, double*** id, 
                                                     double*** level, double*** phi, 
                                                     double alpha, double relax)
{
  auto it = sortedNodes.begin() + queueCounter[0];
  for(int lvl = 1; lvl<queueCounter.size(); lvl++) {
    for(int n = 0; n < queueCounter[lvl]; n++) { //sortedNodes on lvl
      int i(it->i), j(it->j), k(it->k);
      
      it->fin  = 0.0;
      it->fout = 0.0;

      // go over the 6 direct neighbors to compute flux across cell boundaries
      Vec3D dir; //laser direction at cell interface
      double proj; //projection of dir to the direction from [k][j][i] to neighbor
      Vec3D xinter;

      //1. left
      xinter = coords[k][j][i] - Vec3D(0.5*dxyz[k][j][i][0],0,0);
      source.GetDirection(xinter, dir);
      proj = -dir[0]; //dir*(-1,0,0)
      proj *= dxyz[k][j][i][1]*dxyz[k][j][i][2];
      if(proj<0) {
        it->fin  += alpha*proj*l[k][j][i-1];
        it->fout += (1.0-alpha)*proj;
      } else {
        it->fin  += (1.0-alpha)*proj*l[k][j][i-1];
        it->fout += alpha*proj;
      }
      //2. right
      xinter = coords[k][j][i] + Vec3D(0.5*dxyz[k][j][i][0],0,0);
      source.GetDirection(xinter, dir);
      proj = dir[0]; //dir*(1,0,0)
      proj *= dxyz[k][j][i][1]*dxyz[k][j][i][2];
      if(proj<0) {
        it->fin  += alpha*proj*l[k][j][i+1];
        it->fout += (1.0-alpha)*proj;
      } else {
        it->fin  += (1.0-alpha)*proj*l[k][j][i+1];
        it->fout += alpha*proj;
      }
      //3. bottom 
      xinter = coords[k][j][i] - Vec3D(0,0.5*dxyz[k][j][i][1],0);
      source.GetDirection(xinter, dir);
      proj = -dir[1]; //dir*(0,-1,0)
      proj *= dxyz[k][j][i][0]*dxyz[k][j][i][2];
      if(proj<0) {
        it->fin  += alpha*proj*l[k][j-1][i];
        it->fout += (1.0-alpha)*proj;
      } else {
        it->fin  += (1.0-alpha)*proj*l[k][j-1][i];
        it->fout += alpha*proj;
      }
      //4. top 
      xinter = coords[k][j][i] + Vec3D(0,0.5*dxyz[k][j][i][1],0);
      source.GetDirection(xinter, dir);
      proj = dir[1]; //dir*(0,1,0)
      proj *= dxyz[k][j][i][0]*dxyz[k][j][i][2];
      if(proj<0) {
        it->fin  += alpha*proj*l[k][j+1][i];
        it->fout += (1.0-alpha)*proj;
      } else {
        it->fin  += (1.0-alpha)*proj*l[k][j+1][i];
        it->fout += alpha*proj;
      }
      //5. back 
      xinter = coords[k][j][i] - Vec3D(0,0,0.5*dxyz[k][j][i][2]);
      source.GetDirection(xinter, dir);
      proj = -dir[2]; //dir*(0,0,-1)
      proj *= dxyz[k][j][i][0]*dxyz[k][j][i][1];
      if(proj<0) {
        it->fin  += alpha*proj*l[k-1][j][i];
        it->fout += (1.0-alpha)*proj;
      } else {
        it->fin  += (1.0-alpha)*proj*l[k-1][j][i];
        it->fout += alpha*proj;
      }
      //6. front 
      xinter = coords[k][j][i] + Vec3D(0,0,0.5*dxyz[k][j][i][2]);
      source.GetDirection(xinter, dir);
      proj = dir[2]; //dir*(0,0,1)
      proj *= dxyz[k][j][i][0]*dxyz[k][j][i][1];
      if(proj<0) {
        it->fin  += alpha*proj*l[k+1][j][i];
        it->fout += (1.0-alpha)*proj;
      } else {
        it->fin  += (1.0-alpha)*proj*l[k+1][j][i];
        it->fout += alpha*proj;
      }

      //update l[k][j][i]
      double eta = GetAbsorptionCoefficient(T[k][j][i], id[k][j][i]); //absorption coeff.
      l[k][j][i] = (1.0-relax)*l[k][j][i] 
                 + relax*(-(it->fin)/(vol[k][j][i]*eta + it->fout));

      if(l[k][j][i]<lmin) {
        l[k][j][i] = lmin;
        if(verbose >= OutputData::HIGH)
          fprintf(stderr, "Warning: [%d] Applied cut-off radiance (%e) to (%d,%d,%d) (orig:%e).\n", 
                  mpi_rank, lmin, i,j,k, l[k][j][i]);
      }

      it++;
    }

    levelcomm[lvl]->ExchangeAndInsert(l);
    UpdateGhostNodes(l, 2); //a "soft" update. ok if not converged.
  }

}

//--------------------------------------------------------------------------

void
LaserAbsorptionSolver::ComputeTemperatureInLaserDomain(Vec5D*** v, double*** id, double*** T)
{

  auto it = sortedNodes.begin();
  int myid;
  double e;
  for(int lvl = 0; lvl<queueCounter.size(); lvl++) {
    for(int n = 0; n < queueCounter[lvl]; n++) {
      int i(it->i),j(it->j),k(it->k);
      myid = id[k][j][i];
      e = varFcn[myid]->GetInternalEnergyPerUnitMass(v[k][j][i][0], v[k][j][i][4]);
      T[k][j][i] = varFcn[myid]->GetTemperature(v[k][j][i][0], e);
      it++;
    }
    levelcomm[lvl]->ExchangeAndInsert(T);
  }

}

//--------------------------------------------------------------------------

void
LaserAbsorptionSolver::CopyValues(double*** from, double*** to)
{
  //Copy values in the laser domain and at the ghost nodes
  for(auto it = sortedNodes.begin(); it != sortedNodes.end(); it++)
    to[it->k][it->j][it->i] = from[it->k][it->j][it->i];

  for(auto it = ebm.ghostNodes1.begin(); it != ebm.ghostNodes1.end(); it++)
    to[it->first[2]][it->first[1]][it->first[0]] = from[it->first[2]][it->first[1]][it->first[0]];

  for(int n=0; n<ebm.ghostNodes2.size(); n++)
    for(auto it = ebm.ghostNodes2[n].begin(); it != ebm.ghostNodes2[n].end(); it++)
      to[(*it)[2]][(*it)[1]][(*it)[0]] = from[(*it)[2]][(*it)[1]][(*it)[0]];
}

//--------------------------------------------------------------------------

void
LaserAbsorptionSolver::ComputeErrorsInLaserDomain(double*** lold, double*** lnew, double &max_error,
                                                  double &avg_error)
{

  max_error = 0.0;
  avg_error = 0.0;

  // loop through levels 1, 2, 3, ...
  double loc_error;
  for(int n = queueCounter[0]; n < sortedNodes.size(); n++) {
    int i(sortedNodes[n].i), j(sortedNodes[n].j), k(sortedNodes[n].k);
    assert(lnew[k][j][i]>0); //at least, lmin should have been applied.
    loc_error = fabs((lold[k][j][i] - lnew[k][j][i])/lnew[k][j][i]);
    max_error = std::max(max_error, loc_error);
    avg_error += loc_error;
  } 
  int N = sortedNodes.size() - queueCounter[0];

  // MPI communications
  MPI_Allreduce(MPI_IN_PLACE, &max_error, 1, MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(MPI_IN_PLACE, &avg_error, 1, MPI_DOUBLE, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, &N, 1, MPI_INT, MPI_SUM, comm);
  avg_error /= N;

}

//--------------------------------------------------------------------------

void
LaserAbsorptionSolver::CleanUpGhostNodes(double*** l, bool backup)
{

  // Step 1: Clean up ghostNodes1
  for(auto it = ebm.ghostNodes1.begin(); it != ebm.ghostNodes1.end(); it++) {
    int i(it->first[0]), j(it->first[1]), k(it->first[2]);
    if(backup)
      ebm.l1[it - ebm.ghostNodes1.begin()] = l[k][j][i];
    l[k][j][i] = 0.0;
  }

  // Step 2: Clean up ghostNodes2
  for(int p=0; p<ebm.ghostNodes2.size(); p++)
    for(int q=0; q<ebm.ghostNodes2[p].size(); q++) {
      int i(ebm.ghostNodes2[p][q][0]), j(ebm.ghostNodes2[p][q][1]), k(ebm.ghostNodes2[p][q][2]);
      if(backup)
        ebm.l2[p][q] = l[k][j][i];
      l[k][j][i] = 0.0;
    } 

}

//--------------------------------------------------------------------------

void
LaserAbsorptionSolver::ApplyStoredGhostNodesRadiance(double*** l)
{

  for(auto it = ebm.ghostNodes1.begin(); it != ebm.ghostNodes1.end(); it++) {
    int i(it->first[0]), j(it->first[1]), k(it->first[2]);
    l[k][j][i] = ebm.l1[it - ebm.ghostNodes1.begin()];
  }

  for(int p=0; p<ebm.ghostNodes2.size(); p++)
    for(int q=0; q<ebm.ghostNodes2[p].size(); q++) {
      int i(ebm.ghostNodes2[p][q][0]), j(ebm.ghostNodes2[p][q][1]), k(ebm.ghostNodes2[p][q][2]);
      l[k][j][i] = ebm.l2[p][q];
    } 

}

//--------------------------------------------------------------------------

void
LaserAbsorptionSolver::AddHeatToNavierStokesResidual(SpaceVariable3D &R, SpaceVariable3D &L, SpaceVariable3D &ID,
                                                     SpaceVariable3D *V)
{
  Vec5D***  r  = (Vec5D***) R.GetDataPointer();
  double*** l  = L.GetDataPointer();
  double*** id = ID.GetDataPointer();
  Vec5D***  v  = V ? (Vec5D***)V->GetDataPointer() : NULL;
  double*** T  = Temperature.GetDataPointer();

  if(V)
    ComputeTemperatureInLaserDomain(v, id, T);
  //otherwise, use stored temeprature

  auto it = sortedNodes.begin() + queueCounter[0];
  for(int lvl = 1; lvl<queueCounter.size(); lvl++) {
    for(int n = 0; n < queueCounter[lvl]; n++) { //sortedNodes on lvl
      int i(it->i), j(it->j), k(it->k);
      double eta = GetAbsorptionCoefficient(T[k][j][i], id[k][j][i]); //absorption coeff.
      r[k][j][i][4] += eta*l[k][j][i];
      it++;
    }
    // No need of data exchange on r
  }

  Temperature.RestoreDataPointerToLocalVector();
  if(V) V->RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();
  L.RestoreDataPointerToLocalVector();
  R.RestoreDataPointerToLocalVector(); //although data has been udpated, no need to communicate.
                                       //one subdomain does not need the residual info of another
                                       //subdomain
}

//--------------------------------------------------------------------------







