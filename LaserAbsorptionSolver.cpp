#include<LaserAbsorptionSolver.h>
#include<GeoTools.h>
#include<algorithm> //std::sort
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

  // Subdomain info
  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);

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

  // Calculate distance from each node to source
  CalculateDistanceToSource(Phi); 

  // Sort nodes in scope (Not all the nodes) by queue-level (primary) and dist-to-source (secondary)
  BuildSortedNodeList();

  // Create a custom communicator for each level
  BuildCustomizedCommunicators();

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
LaserAbsorptionSolver::CalculateDistanceToSource(SpaceVariable3D &Phi)
{

  LaserData &laser(iod.laser);

  Vec3D x0(laser.source_center_x, laser.source_center_y, laser.source_center_z);
  Vec3D dir(laser.source_dir_x, laser.source_dir_y, laser.source_dir_z);
  dir /= dir.norm(); //normalize direction

  //Define the boundary of laser scope
  //if laser light is parallel (or diverging), the range of laser is given by user
  //if laser light is focused, the distance from focusing point to source center is the range(determined by focusing angle)
  double laser_range;
  double PI = 2.0*acos(0.0);
  double angle = laser.focusing_angle_degrees/180.0*PI; //in radians
  if(angle>0.0){
    double d2f = laser.source_radius / tan(angle);
    laser_range = std::min(laser.range, d2f);
  }
  else
    laser_range = laser.range;

  print("- Laser absorption solver activated. Focusing angle: %f degrees; Source radius: %f; Range: %f\n",
        laser.focusing_angle_degrees, laser.source_radius, laser_range);


  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
  double*** phi   = Phi.GetDataPointer();
  double source_radius = laser.source_radius;

  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        phi[k][j][i] = -1.0; //default (out of scope)

        double d2p  = (coords[k][j][i] - x0)*dir; //distance to source plane
        if(d2p<0 || d2p>source_radius)
          continue;//out of scope
        
        Vec3D xp = coords[k][j][i] - x0 - d2p*dir; //projection of x on the plane of the source
        double r = xp.norm();
        double xradius = source_radius - d2p*tan(angle);

        if(r>xradius) 
          continue;//out of scope

        phi[k][j][i] = d2p;
        numNodesInScope++;
      }

  coordinates.RestoreDataPointerToLocalVector();
  Phi.RestoreDataPointerAndInsert(); 

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

/*
  //---------------------------------------------------
  // Verification
  double*** l0 = L0.GetDataPointer();
  level = Level.GetDataPointer();
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++)
        l0[k][j][i] = level[k][j][i];
  L0.RestoreDataPointerToLocalVector();

  for(int lvl = 0; lvl<queueCounter.size(); lvl++) {
    levelcomm[lvl]->ExchangeAndInsert(L0);
    l0 = L0.GetDataPointer();
    for(auto it = ghost_nodes_inner.begin(); it != ghost_nodes_inner.end(); it++) {
      int i(it->ijk[0]), j(it->ijk[1]), k(it->ijk[2]);
      if(level[k][j][i]==lvl) 
        assert(l0[k][j][i] == lvl);
    }
    L0.RestoreDataPointerToLocalVector();
  }
  //---------------------------------------------------
*/

}

//--------------------------------------------------------------------------







