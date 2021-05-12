#include <ProbeOutput.h>
#include <GeoTools.h>

//-------------------------------------------------------------------------

ProbeOutput::ProbeOutput(MPI_Comm &comm_, OutputData &iod_output) : comm(comm_)
{

  frequency = iod_output.probes.frequency;

  int i;
  for (i = 0; i < Probes::MAXNODES; ++i) {
    locations.push_back(Vec3D(iod_output.probes.myNodes[i].locationX,
                         iod_output.probes.myNodes[i].locationY,
                         iod_output.probes.myNodes[i].locationZ));

    if(locations[i][0] < -1e15)
      break;
    else
      print("- [Probe] Node %d: Coords = (%e, %e, %e).\n", 
            i, locations[i][0], locations[i][1], locations[i][2]);
  }
  numNodes = i;

  for(i=0; i<Probes::SIZE; i++) {
    file[i] = NULL;
  }

  int spn = strlen(iod_output.probes.prefix) + 1;

  if (iod_output.probes.density[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.density)];
    sprintf(filename, "%s%s", iod_output.probes.prefix, iod_output.probes.density);
    file[Probes::DENSITY] = fopen(filename, "w");
    delete [] filename;
  }

  if (iod_output.probes.velocity_x[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.velocity_x)];
    sprintf(filename, "%s%s", iod_output.probes.prefix, iod_output.probes.velocity_x);
    file[Probes::VELOCITY_X] = fopen(filename, "w");
    delete [] filename;
  }

  if (iod_output.probes.velocity_y[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.velocity_y)];
    sprintf(filename, "%s%s", iod_output.probes.prefix, iod_output.probes.velocity_y);
    file[Probes::VELOCITY_Y] = fopen(filename, "w");
    delete [] filename;
  }

  if (iod_output.probes.velocity_z[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.velocity_z)];
    sprintf(filename, "%s%s", iod_output.probes.prefix, iod_output.probes.velocity_z);
    file[Probes::VELOCITY_Z] = fopen(filename, "w");
    delete [] filename;
  }

  if (iod_output.probes.pressure[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.pressure)];
    sprintf(filename, "%s%s", iod_output.probes.prefix, iod_output.probes.pressure);
    file[Probes::PRESSURE] = fopen(filename, "w");
    delete [] filename;
  }

  if (iod_output.probes.temperature[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.temperature)];
    sprintf(filename, "%s%s", iod_output.probes.prefix, iod_output.probes.temperature);
    file[Probes::TEMPERATURE] = fopen(filename, "w");
    delete [] filename;
  }

  if (iod_output.probes.materialid[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.materialid)];
    sprintf(filename, "%s%s", iod_output.probes.prefix, iod_output.probes.materialid);
    file[Probes::MATERIALID] = fopen(filename, "w");
    delete [] filename;
  }

  if (iod_output.probes.levelset0[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.levelset0)];
    sprintf(filename, "%s%s", iod_output.probes.prefix, iod_output.probes.levelset0);
    file[Probes::LEVELSET0] = fopen(filename, "w");
    delete [] filename;
  }

  if (iod_output.probes.levelset1[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.levelset1)];
    sprintf(filename, "%s%s", iod_output.probes.prefix, iod_output.probes.levelset1);
    file[Probes::LEVELSET1] = fopen(filename, "w");
    delete [] filename;
  }

  if (iod_output.probes.levelset2[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.levelset2)];
    sprintf(filename, "%s%s", iod_output.probes.prefix, iod_output.probes.levelset2);
    file[Probes::LEVELSET2] = fopen(filename, "w");
    delete [] filename;
  }

  if (iod_output.probes.levelset3[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.levelset3)];
    sprintf(filename, "%s%s", iod_output.probes.prefix, iod_output.probes.levelset3);
    file[Probes::LEVELSET3] = fopen(filename, "w");
    delete [] filename;
  }

  if (iod_output.probes.levelset4[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.levelset4)];
    sprintf(filename, "%s%s", iod_output.probes.prefix, iod_output.probes.levelset4);
    file[Probes::LEVELSET4] = fopen(filename, "w");
    delete [] filename;
  }

  for(i=0; i<Probes::SIZE; i++)
    if(file[i]) { //write header
      for(int iNode = 0; iNode<numNodes; iNode++)
        print(file[i], "## Probe %d: %e, %e, %e\n", locations[i][0], locations[i][1], locations[i][2]);
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
  if(numNodes == 0) //nothing to be done
    return;

  ijk.resize(numNodes);
  trilinear_coords.resize(numNodes);

  int found[numNodes];

  // get mesh coordinates
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();

  // find subdomain corners (not including the ghost layer)
  int i0, j0, k0, imax, jmax, kmax;
  Vec3D xyz0, xyzmax;
  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  xyz0 = coords[k0-1][j0-1][i0-1]; //lower corner of the ghost layer
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
  for(int iNode = 0; iNode<numNodes; iNode++)
    if(found[iNode] != 1) {
      print_error("*** Error: Cannot locate probe node %d in the mesh (found = %d).\n", iNode, found[iNode]);
      exit_mpi();
    }
  
  coordinates.RestoreDataPointerToLocalVector();
}

//-------------------------------------------------------------------------

void 
ProbeOutput::WriteSolutionAtProbes(double time, int time_step, SpaceVariable3D &V, SpaceVariable3D &ID,
                                   std::vector<SpaceVariable3D*> &Phi)
{
  if(numNodes == 0) //nothing to be done
    return;

  if(time_step % frequency != 0) //should not output at this time step
    return;

  double***  v  = (double***) V.GetDataPointer();

  if(file[Probes::DENSITY]) {
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], trilinear_coords[iNode], v, 5, 0);
      print(file[Probes::DENSITY], "%12.8e    ", sol);
    }
    print(file[Probes::DENSITY],"\n");
  }

  if(file[Probes::VELOCITY_X]) {
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], trilinear_coords[iNode], v, 5, 1);
      print(file[Probes::VELOCITY_X], "%12.8e    ", sol);
    }
    print(file[Probes::VELOCITY_X],"\n");
  }

  if(file[Probes::VELOCITY_Y]) {
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], trilinear_coords[iNode], v, 5, 2);
      print(file[Probes::VELOCITY_Y], "%12.8e    ", sol);
    }
    print(file[Probes::VELOCITY_Y],"\n");
  }

  if(file[Probes::VELOCITY_Z]) {
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], trilinear_coords[iNode], v, 5, 3);
      print(file[Probes::VELOCITY_Z], "%12.8e    ", sol);
    }
    print(file[Probes::VELOCITY_Z],"\n");
  }

  if(file[Probes::PRESSURE]) {
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], trilinear_coords[iNode], v, 5, 4);
      print(file[Probes::PRESSURE], "%12.8e    ", sol);
    }
    print(file[Probes::PRESSURE],"\n");
  }

  if(file[Probes::TEMPERATURE]) {
    print("Warning: Unable to write temperature result at probe locations.\n");
  }

  if(file[Probes::MATERIALID]) {
    double*** id  = (double***)ID.GetDataPointer();
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], trilinear_coords[iNode], id, 1, 0);
      print(file[Probes::MATERIALID], "%12.8e    ", sol);
    }
    print(file[Probes::MATERIALID],"\n");
    ID.RestoreDataPointerToLocalVector();
  }

  if(file[Probes::LEVELSET0]) {
    double*** phi = (double***)Phi[0]->GetDataPointer();
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], trilinear_coords[iNode], phi, 1, 0);
      print(file[Probes::LEVELSET0], "%12.8e    ", sol);
    }
    print(file[Probes::LEVELSET0],"\n");
    Phi[0]->RestoreDataPointerToLocalVector(); //no changes made
  }

  if(file[Probes::LEVELSET1]) {
    double*** phi = (double***)Phi[1]->GetDataPointer();
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], trilinear_coords[iNode], phi, 1, 0);
      print(file[Probes::LEVELSET1], "%12.8e    ", sol);
    }
    print(file[Probes::LEVELSET1],"\n");
    Phi[1]->RestoreDataPointerToLocalVector(); //no changes made
  }

  if(file[Probes::LEVELSET2]) {
    double*** phi = (double***)Phi[2]->GetDataPointer();
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], trilinear_coords[iNode], phi, 1, 0);
      print(file[Probes::LEVELSET2], "%12.8e    ", sol);
    }
    print(file[Probes::LEVELSET2],"\n");
    Phi[2]->RestoreDataPointerToLocalVector(); //no changes made
  }

  if(file[Probes::LEVELSET3]) {
    double*** phi = (double***)Phi[3]->GetDataPointer();
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], trilinear_coords[iNode], phi, 1, 0);
      print(file[Probes::LEVELSET3], "%12.8e    ", sol);
    }
    print(file[Probes::LEVELSET3],"\n");
    Phi[3]->RestoreDataPointerToLocalVector(); //no changes made
  }

  if(file[Probes::LEVELSET4]) {
    double*** phi = (double***)Phi[4]->GetDataPointer();
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], trilinear_coords[iNode], phi, 1, 0);
      print(file[Probes::LEVELSET4], "%12.8e    ", sol);
    }
    print(file[Probes::LEVELSET4],"\n");
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
