#include <ProbeOutput.h>
#include <GeoTools.h>

//-------------------------------------------------------------------------

// This constructor is for explicitly specified probe nodes (i.e. not a line)
ProbeOutput::ProbeOutput(MPI_Comm &comm_, OutputData &iod_output_) : 
             comm(comm_), iod_output(iod_output_)
{
  line_number = -1; //not used in this case

  frequency = iod_output.probes.frequency;

  numNodes = iod_output.probes.myNodes.dataMap.size();
  locations.resize(numNodes);
  for (auto it = iod_output.probes.myNodes.dataMap.begin();
       it!= iod_output.probes.myNodes.dataMap.end(); 
       it++) {
    int pid = it->first;
    if(pid<0 || pid>=numNodes) {
      print_error("*** Error: Probe node index (%d) out of range. Should be between 0 and %d.\n",
                  pid, numNodes-1);
      exit_mpi();
    }
    locations[pid] = Vec3D(it->second->locationX, it->second->locationY, it->second->locationZ);

    print("- [Probe] Node %d: Coords = (%e, %e, %e).\n", 
          it->first, locations[pid][0], locations[pid][1], locations[pid][2]);
  }

  for(int i=0; i<Probes::SIZE; i++) {
    file[i] = NULL;
  }

  int spn = strlen(iod_output.prefix) + 1;

  if (iod_output.probes.density[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.density)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_output.probes.density);
    file[Probes::DENSITY] = fopen(filename, "w");
    delete [] filename;
  }

  if (iod_output.probes.velocity_x[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.velocity_x)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_output.probes.velocity_x);
    file[Probes::VELOCITY_X] = fopen(filename, "w");
    delete [] filename;
  }

  if (iod_output.probes.velocity_y[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.velocity_y)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_output.probes.velocity_y);
    file[Probes::VELOCITY_Y] = fopen(filename, "w");
    delete [] filename;
  }

  if (iod_output.probes.velocity_z[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.velocity_z)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_output.probes.velocity_z);
    file[Probes::VELOCITY_Z] = fopen(filename, "w");
    delete [] filename;
  }

  if (iod_output.probes.pressure[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.pressure)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_output.probes.pressure);
    file[Probes::PRESSURE] = fopen(filename, "w");
    delete [] filename;
  }

  if (iod_output.probes.temperature[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.temperature)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_output.probes.temperature);
    file[Probes::TEMPERATURE] = fopen(filename, "w");
    delete [] filename;
  }

  if (iod_output.probes.materialid[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.materialid)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_output.probes.materialid);
    file[Probes::MATERIALID] = fopen(filename, "w");
    delete [] filename;
  }

  if (iod_output.probes.levelset0[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.levelset0)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_output.probes.levelset0);
    file[Probes::LEVELSET0] = fopen(filename, "w");
    delete [] filename;
  }

  if (iod_output.probes.levelset1[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.levelset1)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_output.probes.levelset1);
    file[Probes::LEVELSET1] = fopen(filename, "w");
    delete [] filename;
  }

  if (iod_output.probes.levelset2[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.levelset2)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_output.probes.levelset2);
    file[Probes::LEVELSET2] = fopen(filename, "w");
    delete [] filename;
  }

  if (iod_output.probes.levelset3[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.levelset3)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_output.probes.levelset3);
    file[Probes::LEVELSET3] = fopen(filename, "w");
    delete [] filename;
  }

  if (iod_output.probes.levelset4[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.levelset4)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_output.probes.levelset4);
    file[Probes::LEVELSET4] = fopen(filename, "w");
    delete [] filename;
  }

  for(int i=0; i<Probes::SIZE; i++)
    if(file[i]) { //write header
      for(int iNode = 0; iNode<numNodes; iNode++)
        print(file[i], "## Probe %d: %e, %e, %e\n", iNode, locations[i][0], locations[i][1], locations[i][2]);
      print(file[i], "## Time step  |  Time  |  Solutions at probe nodes (0, 1, 2, etc.)\n");
      fflush(file[i]);
    }

}

//-------------------------------------------------------------------------
// This constructor is for "line plots"
ProbeOutput::ProbeOutput(MPI_Comm &comm_, OutputData &iod_output_, int line_number_) : 
             comm(comm_), iod_output(iod_output_)
{

  line_number = line_number_;

  LinePlot* line(iod_output.linePlots.dataMap[line_number]);

  numNodes = line->numPoints;
  frequency = line->frequency;

  locations.clear();
  if(numNodes > 0 && frequency > 0) {
    double dx = (line->x1 - line->x0)/(line->numPoints - 1);
    double dy = (line->y1 - line->y0)/(line->numPoints - 1);
    double dz = (line->z1 - line->z0)/(line->numPoints - 1);
    for(int i=0; i<line->numPoints; i++)
      locations.push_back(Vec3D(line->x0 + i*dx, line->y0 + i*dy, line->z0 + i*dz));
  }

  for(int i=0; i<Probes::SIZE; i++) {
    file[i] = NULL;
  }
}

//-------------------------------------------------------------------------

ProbeOutput::~ProbeOutput()
{
  for(int i=0; i<Probes::SIZE; i++)
    if(file[i]) fclose(file[i]); 
}

//-------------------------------------------------------------------------

void
ProbeOutput::SetupInterpolation(SpaceVariable3D &coordinates)
{
  if(numNodes <= 0) //nothing to be done
    return;

  ijk.resize(numNodes);
  trilinear_coords.resize(numNodes);

  int found[numNodes];

  // get mesh coordinates
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();

  // find subdomain corners (have to avoid overlapping & 
  // to account for the ghost layer outside physical domain)
  int i0, j0, k0, imax, jmax, kmax;
  int NX, NY, NZ;
  coordinates.GetGhostedCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates.GetGlobalSize(&NX, &NY, &NZ);
  if(imax != NX+1) imax --;
  if(jmax != NY+1) jmax --;
  if(kmax != NZ+1) kmax --;

  Vec3D xyz0, xyzmax;
  xyz0 = coords[k0][j0][i0]; //lower corner of the ghost layer
  xyzmax = coords[kmax-1][jmax-1][imax-1]; 

  for(int iNode=0; iNode<numNodes; iNode++) {

    Vec3D p = locations[iNode];

    if(p[0] < xyz0[0] || p[0] >= xyzmax[0] || p[1] < xyz0[1] || p[1] >= xyzmax[1] ||
       p[2] < xyz0[2] || p[2] >= xyzmax[2]) {//not in this subdomain
      found[iNode] = 0;
      ijk[iNode] = INT_MIN; //default: not in the current subdomain
    }
    else {//in this subdomain
      found[iNode] = 1; 
      for(int i=i0-1; i<imax-1; i++)
        if(p[0] < coords[k0][j0][i+1][0]) {
          ijk[iNode][0] = i;
          trilinear_coords[iNode][0] = (p[0] - coords[k0][j0][i][0]) / 
                                       (coords[k0][j0][i+1][0] - coords[k0][j0][i][0]);
          break;
        }
      for(int j=j0-1; j<jmax-1; j++)
        if(p[1] < coords[k0][j+1][i0][1]) {
          ijk[iNode][1] = j;
          trilinear_coords[iNode][1] = (p[1] - coords[k0][j][i0][1]) / 
                                       (coords[k0][j+1][i0][1] - coords[k0][j][i0][1]);
          break;
        }
      for(int k=k0-1; k<kmax-1; k++)
        if(p[2] < coords[k+1][j0][i0][2]) {
          ijk[iNode][2] = k;
          trilinear_coords[iNode][2] = (p[2] - coords[k][j0][i0][2]) / 
                                       (coords[k+1][j0][i0][2] - coords[k][j0][i0][2]);
          break;
        }
    }
        
  } 

  MPI_Allreduce(MPI_IN_PLACE, found, numNodes, MPI_INT, MPI_SUM, comm);
  for(int iNode = 0; iNode<numNodes; iNode++) {
    if(found[iNode] != 1) {
      print_error("*** Error: Cannot locate probe node %d in the domain (found = %d).\n", iNode, found[iNode]);
      exit_mpi();
    }
  } 
  coordinates.RestoreDataPointerToLocalVector();
}

//-------------------------------------------------------------------------

void 
ProbeOutput::WriteAllSolutionsAlongLine(double time, int time_step, SpaceVariable3D &V, SpaceVariable3D &ID,
                                        std::vector<SpaceVariable3D*> &Phi, bool must_write)
{
  if(numNodes <= 0) //nothing to be done
    return;

  if(frequency<=0)
    return;

  if(time_step % frequency != 0 && !must_write) //should not output at this time step
    return;

  LinePlot* line(iod_output.linePlots.dataMap[line_number]);

  if(line->filename_base[0] == 0) //file not proided
    return;

  double dx = (line->x1 - line->x0)/(line->numPoints - 1);
  double dy = (line->y1 - line->y0)/(line->numPoints - 1);
  double dz = (line->z1 - line->z0)/(line->numPoints - 1);
  double h = sqrt(dx*dx + dy*dy + dz*dz);

  //figure out file name
  int iFrame;
  if(time_step % frequency == 0)
    iFrame = time_step/frequency;
  else //"must_write" must be true
    iFrame = time_step/frequency + 1;

  char full_fname[256];
  if(iFrame<10) 
    sprintf(full_fname, "%s%s_000%d.vtr", iod_output.prefix, line->filename_base, iFrame);
  else if(iFrame<100)
    sprintf(full_fname, "%s%s_00%d.vtr", iod_output.prefix, line->filename_base, iFrame);
  else if(iFrame<1000)
    sprintf(full_fname, "%s%s_0%d.vtr", iod_output.prefix, line->filename_base, iFrame);
  else
    sprintf(full_fname, "%s%s_%d.vtr", iod_output.prefix, line->filename_base, iFrame);

  //open file & write the header
  FILE *file = fopen(full_fname, "w");
  print(file, "## Line: (%e, %e, %e) -> (%e, %e, %e)\n", line->x0, line->y0, line->z0,
        line->x1, line->y1, line->z1);
  print(file, "## Number of points: %d (h = %e)\n", line->numPoints, h);
  print(file, "## Time: %e, Time step: %d.\n", time, time_step);
  print(file, "## Coordinate  |  Density  |  Velocity (Vx,Vy,Vz)  |  Pressure  |  Material ID  |  LevelSet(s)\n");

  //get data
  double***  v  = (double***) V.GetDataPointer();
  double*** id  = (double***)ID.GetDataPointer();
  std::vector<double***> phi;
  for(int i=0; i<Phi.size(); i++)
    phi.push_back((double***)Phi[i]->GetDataPointer());


  //write data to file
  for(int iNode=0; iNode<numNodes; iNode++) {
    double rho = InterpolateSolutionAtProbe(ijk[iNode], trilinear_coords[iNode], v, 5, 0);
    double vx  = InterpolateSolutionAtProbe(ijk[iNode], trilinear_coords[iNode], v, 5, 1);
    double vy  = InterpolateSolutionAtProbe(ijk[iNode], trilinear_coords[iNode], v, 5, 2);
    double vz  = InterpolateSolutionAtProbe(ijk[iNode], trilinear_coords[iNode], v, 5, 3);
    double p   = InterpolateSolutionAtProbe(ijk[iNode], trilinear_coords[iNode], v, 5, 4);
    double myid= InterpolateSolutionAtProbe(ijk[iNode], trilinear_coords[iNode], id, 1, 0);
    print(file, "%12.8e  %12.8e  %12.8e  %12.8e  %12.8e  %12.8e  %12.8e  ", 
                iNode*h, rho, vx, vy, vz, p, myid);
    for(int i=0; i<Phi.size(); i++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], trilinear_coords[iNode], phi[i], 1, 0);
      print(file, "%12.8e  ", sol);
    }
    print(file, "\n");
  }

  fclose(file);

  V.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();
  for(int i=0; i<Phi.size(); i++)
    Phi[i]->RestoreDataPointerToLocalVector();

}

//-------------------------------------------------------------------------

void 
ProbeOutput::WriteSolutionAtProbes(double time, int time_step, SpaceVariable3D &V, SpaceVariable3D &ID,
                                   std::vector<SpaceVariable3D*> &Phi, bool must_write)
{
  if(numNodes <= 0) //nothing to be done
    return;

  if(frequency<=0)
    return;

  if(time_step % frequency != 0 && !must_write) //should not output at this time step
    return;

  double***  v  = (double***) V.GetDataPointer();

  if(file[Probes::DENSITY]) {
    print(file[Probes::DENSITY], "%10d    %12.8e    ", time_step, time);
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], trilinear_coords[iNode], v, 5, 0);
      print(file[Probes::DENSITY], "%12.8e    ", sol);
    }
    print(file[Probes::DENSITY],"\n");
    fflush(file[Probes::DENSITY]);
  }

  if(file[Probes::VELOCITY_X]) {
    print(file[Probes::VELOCITY_X], "%8d    %12.8e    ", time_step, time);
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], trilinear_coords[iNode], v, 5, 1);
      print(file[Probes::VELOCITY_X], "%12.8e    ", sol);
    }
    print(file[Probes::VELOCITY_X],"\n");
    fflush(file[Probes::VELOCITY_X]);
  }

  if(file[Probes::VELOCITY_Y]) {
    print(file[Probes::VELOCITY_Y], "%8d    %12.8e    ", time_step, time);
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], trilinear_coords[iNode], v, 5, 2);
      print(file[Probes::VELOCITY_Y], "%12.8e    ", sol);
    }
    print(file[Probes::VELOCITY_Y],"\n");
    fflush(file[Probes::VELOCITY_Y]);
  }

  if(file[Probes::VELOCITY_Z]) {
    print(file[Probes::VELOCITY_Z], "%8d    %12.8e    ", time_step, time);
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], trilinear_coords[iNode], v, 5, 3);
      print(file[Probes::VELOCITY_Z], "%12.8e    ", sol);
    }
    print(file[Probes::VELOCITY_Z],"\n");
    fflush(file[Probes::VELOCITY_Z]);
  }

  if(file[Probes::PRESSURE]) {
    print(file[Probes::PRESSURE], "%8d    %12.8e    ", time_step, time);
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], trilinear_coords[iNode], v, 5, 4);
      print(file[Probes::PRESSURE], "%12.8e    ", sol);
    }
    print(file[Probes::PRESSURE],"\n");
    fflush(file[Probes::PRESSURE]);
  }

  if(file[Probes::TEMPERATURE]) {
    print("Warning: Unable to write temperature result at probe locations.\n");
  }

  if(file[Probes::MATERIALID]) {
    print(file[Probes::MATERIALID], "%8d    %12.8e    ", time_step, time);
    double*** id  = (double***)ID.GetDataPointer();
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], trilinear_coords[iNode], id, 1, 0);
      print(file[Probes::MATERIALID], "%12.8e    ", sol);
    }
    print(file[Probes::MATERIALID],"\n");
    fflush(file[Probes::MATERIALID]);
    ID.RestoreDataPointerToLocalVector();
  }

  if(file[Probes::LEVELSET0] && Phi.size()>=1) {
    print(file[Probes::LEVELSET0], "%8d    %12.8e    ", time_step, time);
    double*** phi = (double***)Phi[0]->GetDataPointer();
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], trilinear_coords[iNode], phi, 1, 0);
      print(file[Probes::LEVELSET0], "%12.8e    ", sol);
    }
    print(file[Probes::LEVELSET0],"\n");
    fflush(file[Probes::LEVELSET0]);
    Phi[0]->RestoreDataPointerToLocalVector(); //no changes made
  }

  if(file[Probes::LEVELSET1] && Phi.size()>=2) {
    print(file[Probes::LEVELSET1], "%8d    %12.8e    ", time_step, time);
    double*** phi = (double***)Phi[1]->GetDataPointer();
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], trilinear_coords[iNode], phi, 1, 0);
      print(file[Probes::LEVELSET1], "%12.8e    ", sol);
    }
    print(file[Probes::LEVELSET1],"\n");
    fflush(file[Probes::LEVELSET1]);
    Phi[1]->RestoreDataPointerToLocalVector(); //no changes made
  }

  if(file[Probes::LEVELSET2] && Phi.size()>=3) {
    print(file[Probes::LEVELSET2], "%8d    %12.8e    ", time_step, time);
    double*** phi = (double***)Phi[2]->GetDataPointer();
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], trilinear_coords[iNode], phi, 1, 0);
      print(file[Probes::LEVELSET2], "%12.8e    ", sol);
    }
    print(file[Probes::LEVELSET2],"\n");
    fflush(file[Probes::LEVELSET2]);
    Phi[2]->RestoreDataPointerToLocalVector(); //no changes made
  }

  if(file[Probes::LEVELSET3] && Phi.size()>=4) {
    print(file[Probes::LEVELSET3], "%8d    %12.8e    ", time_step, time);
    double*** phi = (double***)Phi[3]->GetDataPointer();
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], trilinear_coords[iNode], phi, 1, 0);
      print(file[Probes::LEVELSET3], "%12.8e    ", sol);
    }
    print(file[Probes::LEVELSET3],"\n");
    fflush(file[Probes::LEVELSET3]);
    Phi[3]->RestoreDataPointerToLocalVector(); //no changes made
  }

  if(file[Probes::LEVELSET4] && Phi.size()>=5) {
    print(file[Probes::LEVELSET4], "%8d    %12.8e    ", time_step, time);
    double*** phi = (double***)Phi[4]->GetDataPointer();
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], trilinear_coords[iNode], phi, 1, 0);
      print(file[Probes::LEVELSET4], "%12.8e    ", sol);
    }
    print(file[Probes::LEVELSET4],"\n");
    fflush(file[Probes::LEVELSET4]);
    Phi[4]->RestoreDataPointerToLocalVector(); //no changes made
  }

  V.RestoreDataPointerToLocalVector();

}

//-------------------------------------------------------------------------

double
ProbeOutput::InterpolateSolutionAtProbe(Int3& ijk, Vec3D &trilinear_coords, double ***v, int dim, int p)
{
  double sol = 0.0;

  int i = ijk[0], j = ijk[1], k = ijk[2];

  if(i!=INT_MIN && j!=INT_MIN && k!=INT_MIN) {//this probe node is in the current subdomain
    double c000 = v[k][j][i*dim+p];
    double c100 = v[k][j][(i+1)*dim+p];
    double c010 = v[k][j+1][i*dim+p];
    double c110 = v[k][j+1][(i+1)*dim+p];
    double c001 = v[k+1][j][i*dim+p];
    double c101 = v[k+1][j][(i+1)*dim+p];
    double c011 = v[k+1][j+1][i*dim+p];
    double c111 = v[k+1][j+1][(i+1)*dim+p];
    sol = GeoTools::TrilinearInterpolation(trilinear_coords, c000, c100, c010, c110, c001, c101, c011, c111);
  }

  MPI_Allreduce(MPI_IN_PLACE, &sol, 1, MPI_DOUBLE, MPI_SUM, comm);
  return sol;
}

//-------------------------------------------------------------------------

