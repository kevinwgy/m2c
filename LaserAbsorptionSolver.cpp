#include<LaserAbsorptionSolver.h>
#include<GeoTools.h>
using std::vector;
using std::pair;

extern int verbose;

//--------------------------------------------------------------------------

LaserAbsorptionSolver::LaserAbsorptionSolver(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_, 
                         std::vector<VarFcnBase*> &varFcn_, SpaceVariable3D &coordinates_, 
                         SpaceVariable3D &delta_xyz_, SpaceVariable3D &volume_,
                         std::vector<GhostPoint> &ghost_nodes_inner_, 
                         std::vector<GhostPoint> &ghost_nodes_outer_)
                     : iod(iod_), varFcn(varFcn_), coordinates(coordinates_), delta_xyz(delta_xyz_),
                       volume(volume_), ghost_nodes_inner(ghost_nodes_inner_),
                       ghost_nodes_outer(ghost_nodes_outer_),
                       L0(comm_, &(dm_all_.ghosted1_1dof)),
                       FluxIn(comm_, &(dm_all_.ghosted1_1dof)),
                       FluxOut(comm_, &(dm_all_.ghosted1_1dof)),
                       Phi(comm_, &(dm_all_.ghosted1_1dof)),
                       Tag(comm_, &(dm_all_.ghosted1_1dof))
{
  // Check input parameters
  CheckForInputErrors(iod.laser);

  // Subdomain info
  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);

  // Get absorption coefficient for each material
  int numMaterials = iod.eqs.materials.dataMap.size();
  abs.resize(numMaterials, std::make_pair(0.0, 0.0)); //by default, set coeff = 0
  for (auto it = iod.laser.abs.dataMap.begin(); it != iod.laser.abs.dataMap.end(); it++) {
    if(it->first < 0 || it->first >= numMaterials) {
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


}

//--------------------------------------------------------------------------

LaserAbsorptionSolver::~LaserAbsorptionSolver()
{

}

//--------------------------------------------------------------------------

void
LaserAbsorptionSolver::Destroy()
{
  L0.Destroy();
  FluxIn.Destroy();
  FluxOut.Destroy();
  Phi.Destroy();
  Tag.Destroy();
}

//--------------------------------------------------------------------------

void
LaserAbsorptionSolver::CheckForInputErrors(LaserData &laser)
{
  int error = 0;
  if(laser.source_distribution == LaserAbsorption::GAUSSIAN &&
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
  input.open(filename, fstream::in);
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
LaserAbsorptionSolver::CalculateDistanceToSource(SpaceVariable3D &Phi); 
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
        laser.focusing_angle, laser.source_radius, laser_range);


  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
  double*** phi   = Phi.GetDataPointer();
  double source_radius = laser.source_radius;

  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++) {

        phi[k][j][i] = -1.0; //default (out of scope)

        double d2p = (coords[k][j][i] - x0)*dir;
        if(d2p<0 || d2p>laser_range)
          continue;

        Vec3D xp = coords[k][j][i] - x0 - d2p*dir; //projection of x on the plane of the source
        double xradius = source_radius - d2p*tan(angle);

        if(xp.norm()>xradius) 
          continue;//out of scope

        phi[k][j][i] = d2p;
      }

  coordinates.RestoreDataPointerToLocalVector();
  Phi.RestoreDataPointerAndInsert(); 

}

//--------------------------------------------------------------------------

void
LaserAbsorptionSolver::BuildSortedNodeList()
{







}

//--------------------------------------------------------------------------

void
LaserAbsorptionSolver::ResetTag()
{
  //-1: out of scope (no need to compute intensity); 
  // 0: intensity not computed yet, but need to be computed;  
  // 1: intensity has been computed

  double*** tag = Tag.GetDataPointer();

  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++) {


  for (iSub=0; iSub<dom.getNumLocSub(); iSub++) {
    numNodesInScope[iSub] = 0;
    for (int i = 0; i < nodeTag(iSub).size(); ++i)
      if (phi(iSub)[i] > -epsilon) {
        nodeTag(iSub)[i] = 0;
        numNodesInScope[iSub]++;
      }
  


}

//--------------------------------------------------------------------------

