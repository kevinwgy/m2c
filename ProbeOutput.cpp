/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <ProbeOutput.h>
#include <trilinear_interpolation.h>
using std::pair;
using std::array;

//-------------------------------------------------------------------------

// This constructor is for explicitly specified probe nodes (i.e. not a line)
ProbeOutput::ProbeOutput(MPI_Comm &comm_, OutputData &iod_output_, std::vector<VarFcnBase*> &vf_,
                         IonizationOperator* ion_, HyperelasticityOperator *heo_) : 
             comm(comm_), iod_output(iod_output_), vf(vf_), 
             ion(ion_), heo(heo_)
{
  iFrame = 0;

  line_number = -1; //not used in this case

  frequency = iod_output.probes.frequency;
  frequency_dt = iod_output.probes.frequency_dt;

  last_snapshot_time = -1.0;

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

    print("\n- [Probe] Node %d: Coords = (%e, %e, %e).\n", 
          it->first, locations[pid][0], locations[pid][1], locations[pid][2]);
    mpi_barrier();
  }

  for(int i=0; i<Probes::SIZE; i++) {
    file[i] = NULL;
  }

  int spn = strlen(iod_output.prefix) + 1;

  if (iod_output.probes.density[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.density)];
    snprintf(filename, spn + strlen(iod_output.probes.density), //avoid writing out-of-range
             "%s%s", iod_output.prefix, iod_output.probes.density);
    file[Probes::DENSITY] = fopen(filename, "w");

    if(!file[Probes::DENSITY]) {
      print_error("*** Error: Cannot open file '%s' for output.\n", filename);
      exit_mpi();
    }

    delete [] filename;
  }

  if (iod_output.probes.velocity_x[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.velocity_x)];
    snprintf(filename, spn + strlen(iod_output.probes.velocity_x),
             "%s%s", iod_output.prefix, iod_output.probes.velocity_x);
    file[Probes::VELOCITY_X] = fopen(filename, "w");

    if(!file[Probes::VELOCITY_X]) {
      print_error("*** Error: Cannot open file '%s' for output.\n", filename);
      exit_mpi();
    }

    delete [] filename;
  }

  if (iod_output.probes.velocity_y[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.velocity_y)];
    snprintf(filename, spn + strlen(iod_output.probes.velocity_y),
             "%s%s", iod_output.prefix, iod_output.probes.velocity_y);
    file[Probes::VELOCITY_Y] = fopen(filename, "w");

    if(!file[Probes::VELOCITY_Y]) {
      print_error("*** Error: Cannot open file '%s' for output.\n", filename);
      exit_mpi();
    }

    delete [] filename;
  }

  if (iod_output.probes.velocity_z[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.velocity_z)];
    snprintf(filename, spn + strlen(iod_output.probes.velocity_z),
             "%s%s", iod_output.prefix, iod_output.probes.velocity_z);
    file[Probes::VELOCITY_Z] = fopen(filename, "w");

    if(!file[Probes::VELOCITY_Z]) {
      print_error("*** Error: Cannot open file '%s' for output.\n", filename);
      exit_mpi();
    }

    delete [] filename;
  }

  if (iod_output.probes.pressure[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.pressure)];
    snprintf(filename, spn + strlen(iod_output.probes.pressure),
             "%s%s", iod_output.prefix, iod_output.probes.pressure);
    file[Probes::PRESSURE] = fopen(filename, "w");

    if(!file[Probes::PRESSURE]) {
      print_error("*** Error: Cannot open file '%s' for output.\n", filename);
      exit_mpi();
    }

    delete [] filename;
  }

  if (iod_output.probes.temperature[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.temperature)];
    snprintf(filename, spn + strlen(iod_output.probes.temperature),
             "%s%s", iod_output.prefix, iod_output.probes.temperature);
    file[Probes::TEMPERATURE] = fopen(filename, "w");

    if(!file[Probes::TEMPERATURE]) {
      print_error("*** Error: Cannot open file '%s' for output.\n", filename);
      exit_mpi();
    }

    delete [] filename;
  }

  if (iod_output.probes.delta_temperature[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.delta_temperature)];
    snprintf(filename, spn + strlen(iod_output.probes.delta_temperature),
             "%s%s", iod_output.prefix, iod_output.probes.delta_temperature);
    file[Probes::DELTA_TEMPERATURE] = fopen(filename, "w");

    if(!file[Probes::DELTA_TEMPERATURE]) {
      print_error("*** Error: Cannot open file '%s' for output.\n", filename);
      exit_mpi();
    }

    delete [] filename;
  }

  if (iod_output.probes.materialid[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.materialid)];
    snprintf(filename, spn + strlen(iod_output.probes.materialid),
             "%s%s", iod_output.prefix, iod_output.probes.materialid);
    file[Probes::MATERIALID] = fopen(filename, "w");

    if(!file[Probes::MATERIALID]) {
      print_error("*** Error: Cannot open file '%s' for output.\n", filename);
      exit_mpi();
    }

    delete [] filename;
  }

  if (iod_output.probes.laser_radiance[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.laser_radiance)];
    snprintf(filename, spn + strlen(iod_output.probes.laser_radiance),
             "%s%s", iod_output.prefix, iod_output.probes.laser_radiance);
    file[Probes::LASERRADIANCE] = fopen(filename, "w");

    if(!file[Probes::LASERRADIANCE]) {
      print_error("*** Error: Cannot open file '%s' for output.\n", filename);
      exit_mpi();
    }

    delete [] filename;
  }

  if (iod_output.probes.levelset0[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.levelset0)];
    snprintf(filename, spn + strlen(iod_output.probes.levelset0),
             "%s%s", iod_output.prefix, iod_output.probes.levelset0);
    file[Probes::LEVELSET0] = fopen(filename, "w");

    if(!file[Probes::LEVELSET0]) {
      print_error("*** Error: Cannot open file '%s' for output.\n", filename);
      exit_mpi();
    }

    delete [] filename;
  }

  if (iod_output.probes.levelset1[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.levelset1)];
    snprintf(filename, spn + strlen(iod_output.probes.levelset1),
             "%s%s", iod_output.prefix, iod_output.probes.levelset1);
    file[Probes::LEVELSET1] = fopen(filename, "w");

    if(!file[Probes::LEVELSET1]) {
      print_error("*** Error: Cannot open file '%s' for output.\n", filename);
      exit_mpi();
    }

    delete [] filename;
  }

  if (iod_output.probes.levelset2[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.levelset2)];
    snprintf(filename, spn + strlen(iod_output.probes.levelset2),
             "%s%s", iod_output.prefix, iod_output.probes.levelset2);
    file[Probes::LEVELSET2] = fopen(filename, "w");

    if(!file[Probes::LEVELSET2]) {
      print_error("*** Error: Cannot open file '%s' for output.\n", filename);
      exit_mpi();
    }

    delete [] filename;
  }

  if (iod_output.probes.levelset3[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.levelset3)];
    snprintf(filename, spn + strlen(iod_output.probes.levelset3),
             "%s%s", iod_output.prefix, iod_output.probes.levelset3);
    file[Probes::LEVELSET3] = fopen(filename, "w");

    if(!file[Probes::LEVELSET3]) {
      print_error("*** Error: Cannot open file '%s' for output.\n", filename);
      exit_mpi();
    }

    delete [] filename;
  }

  if (iod_output.probes.levelset4[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.levelset4)];
    snprintf(filename, spn + strlen(iod_output.probes.levelset4),
             "%s%s", iod_output.prefix, iod_output.probes.levelset4);
    file[Probes::LEVELSET4] = fopen(filename, "w");

    if(!file[Probes::LEVELSET4]) {
      print_error("*** Error: Cannot open file '%s' for output.\n", filename);
      exit_mpi();
    }

    delete [] filename;
  }

  if (iod_output.probes.ionization_result[0] != 0) {
    if(!ion) {
      print_error("*** Error: Requested ionization result at probe(s), without specifying an ionization model.\n");
      exit_mpi();
    }
    char *filename = new char[spn + strlen(iod_output.probes.ionization_result)];
    snprintf(filename, spn + strlen(iod_output.probes.ionization_result),
             "%s%s", iod_output.prefix, iod_output.probes.ionization_result);
    file[Probes::IONIZATION] = fopen(filename, "w");

    if(!file[Probes::IONIZATION]) {
      print_error("*** Error: Cannot open file '%s' for output.\n", filename);
      exit_mpi();
    }

    delete [] filename;
  }

  if (iod_output.probes.reference_map[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.probes.reference_map)];
    snprintf(filename, spn + strlen(iod_output.probes.reference_map),
             "%s%s", iod_output.prefix, iod_output.probes.reference_map);
    file[Probes::REFERENCE_MAP] = fopen(filename, "w");

    if(!file[Probes::REFERENCE_MAP]) {
      print_error("*** Error: Cannot open file '%s' for output.\n", filename);
      exit_mpi();
    }

    delete [] filename;
  }

  if (iod_output.probes.principal_elastic_stresses[0] != 0) {
    if(!heo) {
      print_error("*** Error: Requested elasticity result at probe(s), without specifying a hyperelasticity model.\n");
      exit_mpi();
    }
    char *filename = new char[spn + strlen(iod_output.probes.principal_elastic_stresses)];
    snprintf(filename, spn + strlen(iod_output.probes.principal_elastic_stresses),
             "%s%s", iod_output.prefix, iod_output.probes.principal_elastic_stresses);
    file[Probes::PRINCIPAL_ELASTIC_STRESSES] = fopen(filename, "w");

    if(!file[Probes::PRINCIPAL_ELASTIC_STRESSES]) {
      print_error("*** Error: Cannot open file '%s' for output.\n", filename);
      exit_mpi();
    }

    delete [] filename;
  }


  for(int i=0; i<Probes::SIZE; i++)
    if(file[i]) { //write header
      for(int iNode = 0; iNode<numNodes; iNode++)
        print(file[i], "## Probe %d: %e, %e, %e\n", iNode, locations[iNode][0], locations[iNode][1], 
                       locations[iNode][2]);
      print(file[i], "## Time step  |  Time  |  Solutions at probe nodes (0, 1, 2, etc.)\n");
      mpi_barrier();
      fflush(file[i]);
    }

}

//-------------------------------------------------------------------------
// This constructor is for "line plots"
ProbeOutput::ProbeOutput(MPI_Comm &comm_, OutputData &iod_output_, std::vector<VarFcnBase*> &vf_,
                         IonizationOperator* ion_, int line_number_) : 
             comm(comm_), iod_output(iod_output_), vf(vf_), ion(ion_)
{

  iFrame = 0;

  line_number = line_number_;

  LinePlot* line(iod_output.linePlots.dataMap[line_number]);

  numNodes = line->numPoints;
  frequency = line->frequency;
  frequency_dt = line->frequency_dt;
  last_snapshot_time = -1.0;

  print("\n- [Line Plot] (%e, %e, %e) --> (%e, %e, %e), %d points.\n", 
        line->x0, line->y0, line->z0, line->x1, line->y1, line->z1, line->numPoints);
  mpi_barrier();

  locations.clear();
  if(numNodes > 1) {
    double dx = (line->x1 - line->x0)/(line->numPoints - 1);
    double dy = (line->y1 - line->y0)/(line->numPoints - 1);
    double dz = (line->z1 - line->z0)/(line->numPoints - 1);
    for(int i=0; i<line->numPoints; i++)
      locations.push_back(Vec3D(line->x0 + i*dx, line->y0 + i*dy, line->z0 + i*dz));
  } else if(numNodes == 1) {
    print_error("*** Error: Must have more than 1 point for a line plot.\n");
    exit(-1);
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
  ijk_valid.resize(numNodes);

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
      ijk_valid[iNode].first = 0;
      ijk_valid[iNode].second.fill(false);
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


      int i(ijk[iNode][0]),j(ijk[iNode][1]),k(ijk[iNode][2]);

      ijk_valid[iNode].second.fill(true);
      if(coordinates.OutsidePhysicalDomainAndUnpopulated(i,j,k)) //c000
        ijk_valid[iNode].second[0] = false;
      if(coordinates.OutsidePhysicalDomainAndUnpopulated(i+1,j,k)) //c100
        ijk_valid[iNode].second[1] = false;
      if(coordinates.OutsidePhysicalDomainAndUnpopulated(i,j+1,k)) //c010
        ijk_valid[iNode].second[2] = false;
      if(coordinates.OutsidePhysicalDomainAndUnpopulated(i+1,j+1,k)) //c110
        ijk_valid[iNode].second[3] = false;
      if(coordinates.OutsidePhysicalDomainAndUnpopulated(i,j,k+1)) //c001
        ijk_valid[iNode].second[4] = false;
      if(coordinates.OutsidePhysicalDomainAndUnpopulated(i+1,j,k+1)) //c101
        ijk_valid[iNode].second[5] = false;
      if(coordinates.OutsidePhysicalDomainAndUnpopulated(i,j+1,k+1)) //c011
        ijk_valid[iNode].second[6] = false;
      if(coordinates.OutsidePhysicalDomainAndUnpopulated(i+1,j+1,k+1)) //c111
        ijk_valid[iNode].second[7] = false;

      ijk_valid[iNode].first = 0;
      for(auto&& val : ijk_valid[iNode].second)
        if(val)
          ijk_valid[iNode].first++;

      if(ijk_valid[iNode].first==0) {
        fprintf(stdout,"\033[0;31m*** Error: Location of probe node %d is too close to the edges "
                       "of the domain boundary.\033[0m\n", iNode);
        exit(-1);
      }

    }
        
  } 

/*
  for(int iNode = 0; iNode<numNodes; iNode++) {
    if(found[iNode]) {
      fprintf(stdout, "Probe %d: (%d, %d, %d): %e %e %e.\n", iNode, ijk[iNode][0], ijk[iNode][1], ijk[iNode][2],
             trilinear_coords[iNode][0], trilinear_coords[iNode][1], trilinear_coords[iNode][2]);
    }
  }
*/

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
ProbeOutput::WriteAllSolutionsAlongLine(double time, double dt, int time_step, SpaceVariable3D &V, SpaceVariable3D &ID,
                                        std::vector<SpaceVariable3D*> &Phi, SpaceVariable3D *L, SpaceVariable3D *Nu_T, bool force_write)
{
  if(numNodes <= 0)
    return;

  if(!isTimeToWrite(time, dt, time_step, frequency_dt, frequency, last_snapshot_time, force_write))
    return;

  LinePlot* line(iod_output.linePlots.dataMap[line_number]);
  if(line->filename_base[0] == 0)
    return;

  double dx = (line->x1 - line->x0)/(line->numPoints - 1);
  double dy = (line->y1 - line->y0)/(line->numPoints - 1);
  double dz = (line->z1 - line->z0)/(line->numPoints - 1);
  double h = sqrt(dx*dx + dy*dy + dz*dz);

  char full_fname[256];
  if(iFrame<10) 
    snprintf(full_fname, 256, "%s%s_000%d.txt", iod_output.prefix, line->filename_base, iFrame);
  else if(iFrame<100)
    snprintf(full_fname, 256, "%s%s_00%d.txt", iod_output.prefix, line->filename_base, iFrame);
  else if(iFrame<1000)
    snprintf(full_fname, 256, "%s%s_0%d.txt", iod_output.prefix, line->filename_base, iFrame);
  else
    snprintf(full_fname, 256, "%s%s_%d.txt", iod_output.prefix, line->filename_base, iFrame);

  //open file & write the header
  FILE *file = fopen(full_fname, "w");
  if(!file) {
    print_error("*** Error: Cannot open file '%s' for output.\n", full_fname);
    exit_mpi();
  }

  print(file, "## Line: (%e, %e, %e) -> (%e, %e, %e)\n", line->x0, line->y0, line->z0,
        line->x1, line->y1, line->z1);
  print(file, "## Number of points: %d (h = %e)\n", line->numPoints, h);
  print(file, "## Time: %e, Time step: %d.\n", time, time_step);
  if(L) 
    print(file, "## Coordinate  |  Density  |  Velocity (Vx,Vy,Vz)  |  Pressure  |  Temperature  |  Material ID  "
                "|  Laser Radiance  |  LevelSet(s)");
  if(Nu_T) 
    print(file, "## Coordinate  |  Density  |  Velocity (Vx,Vy,Vz)  |  Pressure  |  Temperature  |  Material ID  "
                "|  Eddy Viscosity  |  LevelSet(s)");
  else
    print(file, "## Coordinate  |  Density  |  Velocity (Vx,Vy,Vz)  |  Pressure  |  Temperature  |  Material ID  "
                "|  LevelSet(s)");

  if(ion)
    print(file, "  |  Mean Charge  |  Heavy Particles Density");

  print(file, "\n");
  MPI_Barrier(comm);

  //get data
  double***  v  = (double***) V.GetDataPointer();
  double*** id  = (double***)ID.GetDataPointer();
  double***  l  = L? L->GetDataPointer() : NULL;
  double***  nu_t  = Nu_T? Nu_T->GetDataPointer() : NULL;
  std::vector<double***> phi;
  for(int i=0; i<(int)Phi.size(); i++)
    phi.push_back((double***)Phi[i]->GetDataPointer());


  //write data to file
  for(int iNode=0; iNode<numNodes; iNode++) {
    double rho = InterpolateSolutionAtProbe(ijk[iNode], ijk_valid[iNode], trilinear_coords[iNode], v, 5, 0);
    double vx  = InterpolateSolutionAtProbe(ijk[iNode], ijk_valid[iNode], trilinear_coords[iNode], v, 5, 1);
    double vy  = InterpolateSolutionAtProbe(ijk[iNode], ijk_valid[iNode], trilinear_coords[iNode], v, 5, 2);
    double vz  = InterpolateSolutionAtProbe(ijk[iNode], ijk_valid[iNode], trilinear_coords[iNode], v, 5, 3);
    double p   = InterpolateSolutionAtProbe(ijk[iNode], ijk_valid[iNode], trilinear_coords[iNode], v, 5, 4);
    double T   = CalculateTemperatureAtProbe(ijk[iNode], ijk_valid[iNode], trilinear_coords[iNode], v, id);
    double myid= InterpolateSolutionAtProbe(ijk[iNode], ijk_valid[iNode], trilinear_coords[iNode], id, 1, 0);
    print(file, "%16.8e  %16.8e  %16.8e  %16.8e  %16.8e  %16.8e  %16.8e  %16.8e", 
                iNode*h, rho, vx, vy, vz, p, T, myid);
    if(l) {
      double laser_rad = InterpolateSolutionAtProbe(ijk[iNode], ijk_valid[iNode], trilinear_coords[iNode], l, 1, 0);
      print(file, "%16.8e  ", laser_rad);
    }

    if(nu_t) {
      double eddy_visc  = InterpolateSolutionAtProbe(ijk[iNode], ijk_valid[iNode], trilinear_coords[iNode], nu_t, 1, 0);
      print(file, "%16.8e  ", eddy_visc);
    }
    for(int i=0; i<(int)Phi.size(); i++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], ijk_valid[iNode], trilinear_coords[iNode], phi[i], 1, 0);
      print(file, "%16.8e  ", sol);
    }

    if(ion) {
      Vec3D ion_res = CalculateIonizationAtProbe(ijk[iNode], ijk_valid[iNode], trilinear_coords[iNode], v, id);
      print(file, "%16.8e  %16.8e  ", ion_res[0], ion_res[1]);
    }

    print(file, "\n");
  }
  MPI_Barrier(comm);

  fclose(file);


  V.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();
  if(L) L->RestoreDataPointerToLocalVector();
  if(Nu_T) Nu_T->RestoreDataPointerToLocalVector();
  for(int i=0; i<(int)Phi.size(); i++)
    Phi[i]->RestoreDataPointerToLocalVector();

  iFrame++;
  last_snapshot_time = time;
}

//-------------------------------------------------------------------------

void 
ProbeOutput::WriteSolutionAtProbes(double time, double dt, int time_step, SpaceVariable3D &V, SpaceVariable3D &ID,
                                   std::vector<SpaceVariable3D*> &Phi, SpaceVariable3D *L, 
                                   SpaceVariable3D *Xi, bool force_write)
{

  if(numNodes <= 0) //nothing to be done
    return;

  if(!isTimeToWrite(time,dt,time_step,frequency_dt,frequency,last_snapshot_time,force_write))
    return;

  double***  v  = (double***) V.GetDataPointer();

  if(file[Probes::DENSITY]) {
    print(file[Probes::DENSITY], "%10d    %16.8e    ", time_step, time);
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], ijk_valid[iNode], trilinear_coords[iNode], v, 5, 0);
      print(file[Probes::DENSITY], "%16.8e    ", sol);
    }
    print(file[Probes::DENSITY],"\n");
    mpi_barrier();
    fflush(file[Probes::DENSITY]);
  }

  if(file[Probes::VELOCITY_X]) {
    print(file[Probes::VELOCITY_X], "%8d    %16.8e    ", time_step, time);
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], ijk_valid[iNode], trilinear_coords[iNode], v, 5, 1);
      print(file[Probes::VELOCITY_X], "%16.8e    ", sol);
    }
    print(file[Probes::VELOCITY_X],"\n");
    mpi_barrier();
    fflush(file[Probes::VELOCITY_X]);
  }

  if(file[Probes::VELOCITY_Y]) {
    print(file[Probes::VELOCITY_Y], "%8d    %16.8e    ", time_step, time);
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], ijk_valid[iNode], trilinear_coords[iNode], v, 5, 2);
      print(file[Probes::VELOCITY_Y], "%16.8e    ", sol);
    }
    print(file[Probes::VELOCITY_Y],"\n");
    mpi_barrier();
    fflush(file[Probes::VELOCITY_Y]);
  }

  if(file[Probes::VELOCITY_Z]) {
    print(file[Probes::VELOCITY_Z], "%8d    %16.8e    ", time_step, time);
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], ijk_valid[iNode], trilinear_coords[iNode], v, 5, 3);
      print(file[Probes::VELOCITY_Z], "%16.8e    ", sol);
    }
    print(file[Probes::VELOCITY_Z],"\n");
    mpi_barrier();
    fflush(file[Probes::VELOCITY_Z]);
  }

  if(file[Probes::PRESSURE]) {
    print(file[Probes::PRESSURE], "%8d    %16.8e    ", time_step, time);
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], ijk_valid[iNode], trilinear_coords[iNode], v, 5, 4);
      print(file[Probes::PRESSURE], "%16.8e    ", sol);
    }
    print(file[Probes::PRESSURE],"\n");
    mpi_barrier();
    fflush(file[Probes::PRESSURE]);
  }

  if(file[Probes::TEMPERATURE]) {
    print(file[Probes::TEMPERATURE], "%8d    %16.8e    ", time_step, time);
    double*** id  = (double***)ID.GetDataPointer();
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = CalculateTemperatureAtProbe(ijk[iNode], ijk_valid[iNode], trilinear_coords[iNode], v, id);
      print(file[Probes::TEMPERATURE], "%16.8e    ", sol);
    }
    print(file[Probes::TEMPERATURE],"\n");
    mpi_barrier();
    fflush(file[Probes::TEMPERATURE]);
    ID.RestoreDataPointerToLocalVector();
  }

  if(file[Probes::DELTA_TEMPERATURE]) {
    print(file[Probes::DELTA_TEMPERATURE], "%8d    %16.8e    ", time_step, time);
    double*** id  = (double***)ID.GetDataPointer();
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = CalculateDeltaTemperatureAtProbe(ijk[iNode], ijk_valid[iNode], trilinear_coords[iNode], v, id);
      print(file[Probes::DELTA_TEMPERATURE], "%16.8e    ", sol);
    }
    print(file[Probes::DELTA_TEMPERATURE],"\n");
    mpi_barrier();
    fflush(file[Probes::DELTA_TEMPERATURE]);
    ID.RestoreDataPointerToLocalVector();
  }

  if(file[Probes::MATERIALID]) {
    print(file[Probes::MATERIALID], "%8d    %16.8e    ", time_step, time);
    double*** id  = (double***)ID.GetDataPointer();
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], ijk_valid[iNode], trilinear_coords[iNode], id, 1, 0);
      print(file[Probes::MATERIALID], "%16.8e    ", sol);
    }
    print(file[Probes::MATERIALID],"\n");
    mpi_barrier();
    fflush(file[Probes::MATERIALID]);
    ID.RestoreDataPointerToLocalVector();
  }

  if(file[Probes::LASERRADIANCE]) {
    print(file[Probes::LASERRADIANCE], "%8d    %16.8e    ", time_step, time);
    if(L == NULL) {
      print_error("*** Error: Requested laser radiance probe, but laser source is not specified.\n");
      exit_mpi();
    }
    double*** l  = (double***)L->GetDataPointer();
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], ijk_valid[iNode], trilinear_coords[iNode], l, 1, 0);
      print(file[Probes::LASERRADIANCE], "%16.8e    ", sol);
    }
    print(file[Probes::LASERRADIANCE],"\n");
    mpi_barrier();
    fflush(file[Probes::LASERRADIANCE]);
    L->RestoreDataPointerToLocalVector();
  }

  if(file[Probes::LEVELSET0] && Phi.size()>=1) {
    print(file[Probes::LEVELSET0], "%8d    %16.8e    ", time_step, time);
    double*** phi = (double***)Phi[0]->GetDataPointer();
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], ijk_valid[iNode], trilinear_coords[iNode], phi, 1, 0);
      print(file[Probes::LEVELSET0], "%16.8e    ", sol);
    }
    print(file[Probes::LEVELSET0],"\n");
    mpi_barrier();
    fflush(file[Probes::LEVELSET0]);
    Phi[0]->RestoreDataPointerToLocalVector(); //no changes made
  }

  if(file[Probes::LEVELSET1] && Phi.size()>=2) {
    print(file[Probes::LEVELSET1], "%8d    %16.8e    ", time_step, time);
    double*** phi = (double***)Phi[1]->GetDataPointer();
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], ijk_valid[iNode], trilinear_coords[iNode], phi, 1, 0);
      print(file[Probes::LEVELSET1], "%16.8e    ", sol);
    }
    print(file[Probes::LEVELSET1],"\n");
    mpi_barrier();
    fflush(file[Probes::LEVELSET1]);
    Phi[1]->RestoreDataPointerToLocalVector(); //no changes made
  }

  if(file[Probes::LEVELSET2] && Phi.size()>=3) {
    print(file[Probes::LEVELSET2], "%8d    %16.8e    ", time_step, time);
    double*** phi = (double***)Phi[2]->GetDataPointer();
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], ijk_valid[iNode], trilinear_coords[iNode], phi, 1, 0);
      print(file[Probes::LEVELSET2], "%16.8e    ", sol);
    }
    print(file[Probes::LEVELSET2],"\n");
    mpi_barrier();
    fflush(file[Probes::LEVELSET2]);
    Phi[2]->RestoreDataPointerToLocalVector(); //no changes made
  }

  if(file[Probes::LEVELSET3] && Phi.size()>=4) {
    print(file[Probes::LEVELSET3], "%8d    %16.8e    ", time_step, time);
    double*** phi = (double***)Phi[3]->GetDataPointer();
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], ijk_valid[iNode], trilinear_coords[iNode], phi, 1, 0);
      print(file[Probes::LEVELSET3], "%16.8e    ", sol);
    }
    print(file[Probes::LEVELSET3],"\n");
    mpi_barrier();
    fflush(file[Probes::LEVELSET3]);
    Phi[3]->RestoreDataPointerToLocalVector(); //no changes made
  }

  if(file[Probes::LEVELSET4] && Phi.size()>=5) {
    print(file[Probes::LEVELSET4], "%8d    %16.8e    ", time_step, time);
    double*** phi = (double***)Phi[4]->GetDataPointer();
    for(int iNode=0; iNode<numNodes; iNode++) {
      double sol = InterpolateSolutionAtProbe(ijk[iNode], ijk_valid[iNode], trilinear_coords[iNode], phi, 1, 0);
      print(file[Probes::LEVELSET4], "%16.8e    ", sol);
    }
    print(file[Probes::LEVELSET4],"\n");
    mpi_barrier();
    fflush(file[Probes::LEVELSET4]);
    Phi[4]->RestoreDataPointerToLocalVector(); //no changes made
  }

  if(file[Probes::IONIZATION]) {
    print(file[Probes::IONIZATION], "%8d    %16.8e    ", time_step, time);
    double*** id  = (double***)ID.GetDataPointer();
    for(int iNode=0; iNode<numNodes; iNode++) {
      Vec3D sol = CalculateIonizationAtProbe(ijk[iNode], ijk_valid[iNode], trilinear_coords[iNode], v, id);
      print(file[Probes::IONIZATION], "%16.8e    %16.8e    ", sol[0], sol[1]);
    }
    print(file[Probes::IONIZATION],"\n");
    mpi_barrier();
    fflush(file[Probes::IONIZATION]);
    ID.RestoreDataPointerToLocalVector();
  }

  if(file[Probes::REFERENCE_MAP]) {
    print(file[Probes::REFERENCE_MAP], "%8d    %16.8e    ", time_step, time);
    if(Xi == NULL) {
      print_error("*** Error: Requested reference map probe, but ReferenceMapOperator is not constructed.\n");
      exit_mpi();
    }
    double*** xi  = (double***)Xi->GetDataPointer();
    for(int iNode=0; iNode<numNodes; iNode++) {
      Vec3D sol;
      sol[0] = InterpolateSolutionAtProbe(ijk[iNode], ijk_valid[iNode], trilinear_coords[iNode], xi, 3, 0);
      sol[1] = InterpolateSolutionAtProbe(ijk[iNode], ijk_valid[iNode], trilinear_coords[iNode], xi, 3, 1);
      sol[2] = InterpolateSolutionAtProbe(ijk[iNode], ijk_valid[iNode], trilinear_coords[iNode], xi, 3, 2);
      print(file[Probes::REFERENCE_MAP], "%16.8e    %16.8e    %16.8e    ", sol[0], sol[1], sol[2]);
    }
    print(file[Probes::REFERENCE_MAP],"\n");
    mpi_barrier();
    fflush(file[Probes::REFERENCE_MAP]);
    Xi->RestoreDataPointerToLocalVector();
  }


  if(file[Probes::PRINCIPAL_ELASTIC_STRESSES]) {
    print(file[Probes::PRINCIPAL_ELASTIC_STRESSES], "%8d    %16.8e    ", time_step, time);
    assert(heo && Xi); //this is a just a redundant check
    vector<Vec3D> sol; //to be filled
    heo->ComputePrincipalStressesAtProbes(*Xi, ID, ijk, ijk_valid, trilinear_coords, (Vec5D***)v, sol);
    assert(numNodes == (int)sol.size());
    for(int iNode=0; iNode<numNodes; iNode++) {
      print(file[Probes::PRINCIPAL_ELASTIC_STRESSES], "%16.8e    %16.8e    %16.8e   ",
            sol[iNode][0], sol[iNode][1], sol[iNode][2]);
    }
    print(file[Probes::PRINCIPAL_ELASTIC_STRESSES],"\n");
    mpi_barrier();
    fflush(file[Probes::PRINCIPAL_ELASTIC_STRESSES]);
  }


  MPI_Barrier(comm);

  V.RestoreDataPointerToLocalVector();

  last_snapshot_time = time;

}

//-------------------------------------------------------------------------

double
ProbeOutput::InterpolateSolutionAtProbe(Int3& ijk, pair<int, array<bool,8> >& ijk_valid,
                                        Vec3D &trilinear_coords, double ***v, int dim, int p)
{
  double sol = 0.0;

  int i = ijk[0], j = ijk[1], k = ijk[2];

  if(i!=INT_MIN && j!=INT_MIN && k!=INT_MIN) {//this probe node is in the current subdomain
    double c000 = ijk_valid.second[0] ? v[k][j][i*dim+p]         : 0.0;
    double c100 = ijk_valid.second[1] ? v[k][j][(i+1)*dim+p]     : 0.0;
    double c010 = ijk_valid.second[2] ? v[k][j+1][i*dim+p]       : 0.0;
    double c110 = ijk_valid.second[3] ? v[k][j+1][(i+1)*dim+p]   : 0.0;
    double c001 = ijk_valid.second[4] ? v[k+1][j][i*dim+p]       : 0.0;
    double c101 = ijk_valid.second[5] ? v[k+1][j][(i+1)*dim+p]   : 0.0;
    double c011 = ijk_valid.second[6] ? v[k+1][j+1][i*dim+p]     : 0.0;
    double c111 = ijk_valid.second[7] ? v[k+1][j+1][(i+1)*dim+p] : 0.0;

    if(ijk_valid.first<8) {//fill invalid slots with average value
      double c_avg = (c000+c100+c010+c110+c001+c101+c011+c111)/ijk_valid.first;
      if(!ijk_valid.second[0])  c000 = c_avg;
      if(!ijk_valid.second[1])  c100 = c_avg;
      if(!ijk_valid.second[2])  c010 = c_avg;
      if(!ijk_valid.second[3])  c110 = c_avg;
      if(!ijk_valid.second[4])  c001 = c_avg;
      if(!ijk_valid.second[5])  c101 = c_avg;
      if(!ijk_valid.second[6])  c011 = c_avg;
      if(!ijk_valid.second[7])  c111 = c_avg;
    }

    sol = MathTools::trilinear_interpolation(trilinear_coords, c000, c100, c010, c110, c001, c101, c011, c111);
  }

  MPI_Allreduce(MPI_IN_PLACE, &sol, 1, MPI_DOUBLE, MPI_SUM, comm);
  return sol;
}

//-------------------------------------------------------------------------

double
ProbeOutput::CalculateTemperatureAtProbe(Int3& ijk, pair<int, array<bool,8> >& ijk_valid,
                                         Vec3D &trilinear_coords, double ***v, double ***id)
{
  double sol = 0.0;

  int i = ijk[0], j = ijk[1], k = ijk[2];
  int dim = 5;
  double rho,p,e;
  int myid;

  if(i!=INT_MIN && j!=INT_MIN && k!=INT_MIN) {//this probe node is in the current subdomain

    // c000
    double c000 = 0.0;
    if(ijk_valid.second[0]) {
      myid = id[k][j][i]; 
      if(vf[myid]->type != VarFcnBase::HOMOGENEOUS_INCOMPRESSIBLE) {
        rho  =  v[k][j][i*dim];
        p    =  v[k][j][i*dim+4];
        e    = vf[myid]->GetInternalEnergyPerUnitMass(rho,p);
        c000 = vf[myid]->GetTemperature(rho,e);
      }
    }

    // c100
    double c100 = 0.0;
    if(ijk_valid.second[1]) {
      myid = id[k][j][i+1];
      if(vf[myid]->type != VarFcnBase::HOMOGENEOUS_INCOMPRESSIBLE) {
        rho  =  v[k][j][(i+1)*dim];
        p    =  v[k][j][(i+1)*dim+4];
        e    = vf[myid]->GetInternalEnergyPerUnitMass(rho,p);
        c100 = vf[myid]->GetTemperature(rho,e);
      }
    }

    // c010
    double c010 = 0.0;
    if(ijk_valid.second[2]) {
      myid = id[k][j+1][i];
      if(vf[myid]->type != VarFcnBase::HOMOGENEOUS_INCOMPRESSIBLE) {
        rho  =  v[k][j+1][i*dim];
        p    =  v[k][j+1][i*dim+4];
        e    = vf[myid]->GetInternalEnergyPerUnitMass(rho,p);
        c010 = vf[myid]->GetTemperature(rho,e);
      }
    }

    // c110
    double c110 = 0.0;
    if(ijk_valid.second[3]) {
      myid = id[k][j+1][i+1];
      if(vf[myid]->type != VarFcnBase::HOMOGENEOUS_INCOMPRESSIBLE) {
        rho  =  v[k][j+1][(i+1)*dim];
        p    =  v[k][j+1][(i+1)*dim+4];
        e    = vf[myid]->GetInternalEnergyPerUnitMass(rho,p);
        c110 = vf[myid]->GetTemperature(rho,e);
      }
    }

    // c001
    double c001 = 0.0;
    if(ijk_valid.second[4]) {
      myid = id[k+1][j][i];
      if(vf[myid]->type != VarFcnBase::HOMOGENEOUS_INCOMPRESSIBLE) {
        rho  =  v[k+1][j][i*dim];
        p    =  v[k+1][j][i*dim+4];
        e    = vf[myid]->GetInternalEnergyPerUnitMass(rho,p);
        c001 = vf[myid]->GetTemperature(rho,e);
      }
    }

    // c101
    double c101 = 0.0;
    if(ijk_valid.second[5]) {
      myid = id[k+1][j][i+1];
      if(vf[myid]->type != VarFcnBase::HOMOGENEOUS_INCOMPRESSIBLE) {
        rho  =  v[k+1][j][(i+1)*dim];
        p    =  v[k+1][j][(i+1)*dim+4];
        e    = vf[myid]->GetInternalEnergyPerUnitMass(rho,p);
        c101 = vf[myid]->GetTemperature(rho,e);
      }
    }

    // c011
    double c011 = 0.0;
    if(ijk_valid.second[6]) {
      myid = id[k+1][j+1][i];
      if(vf[myid]->type != VarFcnBase::HOMOGENEOUS_INCOMPRESSIBLE) {
        rho  =  v[k+1][j+1][i*dim];
        p    =  v[k+1][j+1][i*dim+4];
        e    = vf[myid]->GetInternalEnergyPerUnitMass(rho,p);
        c011 = vf[myid]->GetTemperature(rho,e);
      }
    }

    // c111
    double c111 = 0.0;
    if(ijk_valid.second[7]) {
      myid = id[k+1][j+1][i+1];
      if(vf[myid]->type != VarFcnBase::HOMOGENEOUS_INCOMPRESSIBLE) {
        rho  =  v[k+1][j+1][(i+1)*dim];
        p    =  v[k+1][j+1][(i+1)*dim+4];
        e    = vf[myid]->GetInternalEnergyPerUnitMass(rho,p);
        c111 = vf[myid]->GetTemperature(rho,e);
      }
    }

    if(ijk_valid.first<8) {//fill invalid slots with average value
      double c_avg = (c000+c100+c010+c110+c001+c101+c011+c111)/ijk_valid.first;
      if(!ijk_valid.second[0])  c000 = c_avg;
      if(!ijk_valid.second[1])  c100 = c_avg;
      if(!ijk_valid.second[2])  c010 = c_avg;
      if(!ijk_valid.second[3])  c110 = c_avg;
      if(!ijk_valid.second[4])  c001 = c_avg;
      if(!ijk_valid.second[5])  c101 = c_avg;
      if(!ijk_valid.second[6])  c011 = c_avg;
      if(!ijk_valid.second[7])  c111 = c_avg;
    }

    sol = MathTools::trilinear_interpolation(trilinear_coords, c000, c100, c010, c110, c001, c101, c011, c111);
  }

  MPI_Allreduce(MPI_IN_PLACE, &sol, 1, MPI_DOUBLE, MPI_SUM, comm);
  return sol;
}

//-------------------------------------------------------------------------

double
ProbeOutput::CalculateDeltaTemperatureAtProbe(Int3& ijk, pair<int, array<bool,8> >& ijk_valid, 
                                              Vec3D &trilinear_coords, double ***v, double ***id)
{
  double sol = 0.0;

  int i = ijk[0], j = ijk[1], k = ijk[2];
  int dim = 5;
  double rho,p,e;
  int myid;

  if(i!=INT_MIN && j!=INT_MIN && k!=INT_MIN) {//this probe node is in the current subdomain

    // c000
    double c000 = 0.0;
    if(ijk_valid.second[0]) {
      myid = id[k][j][i]; 
      rho  =  v[k][j][i*dim];
      p    =  v[k][j][i*dim+4];
      e    = vf[myid]->GetInternalEnergyPerUnitMass(rho,p);
      c000 = vf[myid]->GetTemperature(rho,e) - vf[myid]->GetReferenceTemperature();
    }

    // c100
    double c100 = 0.0;
    if(ijk_valid.second[1]) {
      myid = id[k][j][i+1];
      rho  =  v[k][j][(i+1)*dim];
      p    =  v[k][j][(i+1)*dim+4];
      e    = vf[myid]->GetInternalEnergyPerUnitMass(rho,p);
      c100 = vf[myid]->GetTemperature(rho,e) - vf[myid]->GetReferenceTemperature();
    }

    // c010
    double c010 = 0.0;
    if(ijk_valid.second[2]) { 
      myid = id[k][j+1][i];
      rho  =  v[k][j+1][i*dim];
      p    =  v[k][j+1][i*dim+4];
      e    = vf[myid]->GetInternalEnergyPerUnitMass(rho,p);
      c010 = vf[myid]->GetTemperature(rho,e) - vf[myid]->GetReferenceTemperature();
    }

    // c110
    double c110 = 0.0;
    if(ijk_valid.second[3]) {
      myid = id[k][j+1][i+1];
      rho  =  v[k][j+1][(i+1)*dim];
      p    =  v[k][j+1][(i+1)*dim+4];
      e    = vf[myid]->GetInternalEnergyPerUnitMass(rho,p);
      c110 = vf[myid]->GetTemperature(rho,e) - vf[myid]->GetReferenceTemperature();
    }

    // c001
    double c001 = 0.0;
    if(ijk_valid.second[4]) {
      myid = id[k+1][j][i];
      rho  =  v[k+1][j][i*dim];
      p    =  v[k+1][j][i*dim+4];
      e    = vf[myid]->GetInternalEnergyPerUnitMass(rho,p);
      c001 = vf[myid]->GetTemperature(rho,e) - vf[myid]->GetReferenceTemperature();
    }

    // c101
    double c101 = 0.0;
    if(ijk_valid.second[5]) {
      myid = id[k+1][j][i+1];
      rho  =  v[k+1][j][(i+1)*dim];
      p    =  v[k+1][j][(i+1)*dim+4];
      e    = vf[myid]->GetInternalEnergyPerUnitMass(rho,p);
      c101 = vf[myid]->GetTemperature(rho,e) - vf[myid]->GetReferenceTemperature();
    }

    // c011
    double c011 = 0.0;
    if(ijk_valid.second[6]) {
      myid = id[k+1][j+1][i];
      rho  =  v[k+1][j+1][i*dim];
      p    =  v[k+1][j+1][i*dim+4];
      e    = vf[myid]->GetInternalEnergyPerUnitMass(rho,p);
      c011 = vf[myid]->GetTemperature(rho,e) - vf[myid]->GetReferenceTemperature();
    }

    // c111
    double c111 = 0.0;
    if(ijk_valid.second[7]) {
      myid = id[k+1][j+1][i+1];
      rho  =  v[k+1][j+1][(i+1)*dim];
      p    =  v[k+1][j+1][(i+1)*dim+4];
      e    = vf[myid]->GetInternalEnergyPerUnitMass(rho,p);
      c111 = vf[myid]->GetTemperature(rho,e) - vf[myid]->GetReferenceTemperature();
    }


    if(ijk_valid.first<8) {//fill invalid slots with average value
      double c_avg = (c000+c100+c010+c110+c001+c101+c011+c111)/ijk_valid.first;
      if(!ijk_valid.second[0])  c000 = c_avg;
      if(!ijk_valid.second[1])  c100 = c_avg;
      if(!ijk_valid.second[2])  c010 = c_avg;
      if(!ijk_valid.second[3])  c110 = c_avg;
      if(!ijk_valid.second[4])  c001 = c_avg;
      if(!ijk_valid.second[5])  c101 = c_avg;
      if(!ijk_valid.second[6])  c011 = c_avg;
      if(!ijk_valid.second[7])  c111 = c_avg;
    }

    sol = MathTools::trilinear_interpolation(trilinear_coords, c000, c100, c010, c110, c001, c101, c011, c111);
  }

  MPI_Allreduce(MPI_IN_PLACE, &sol, 1, MPI_DOUBLE, MPI_SUM, comm);
  return sol;
}

//-------------------------------------------------------------------------

Vec3D
ProbeOutput::CalculateIonizationAtProbe(Int3& ijk, pair<int, array<bool,8> >& ijk_valid,
                                        Vec3D &trilinear_coords, double ***v, double ***id)
{
  Vec3D sol(0.0);

  int i = ijk[0], j = ijk[1], k = ijk[2];
  int dim = 5;
  double rho,p;
  int myid;

  if(i!=INT_MIN && j!=INT_MIN && k!=INT_MIN) {//this probe node is in the current subdomain

    // c000
    Vec3D c000 = 0.0;
    if(ijk_valid.second[0]) {
      myid = id[k][j][i]; 
      rho  =  v[k][j][i*dim];
      p    =  v[k][j][i*dim+4];
      c000 = ion->ComputeIonizationAtOnePoint(myid, rho, p);
    }

    // c100
    Vec3D c100 = 0.0;
    if(ijk_valid.second[1]) {
      myid = id[k][j][i+1];
      rho  =  v[k][j][(i+1)*dim];
      p    =  v[k][j][(i+1)*dim+4];
      c100 = ion->ComputeIonizationAtOnePoint(myid, rho, p);
    }

    // c010
    Vec3D c010 = 0.0;
    if(ijk_valid.second[2]) {
      myid = id[k][j+1][i];
      rho  =  v[k][j+1][i*dim];
      p    =  v[k][j+1][i*dim+4];
      c010 = ion->ComputeIonizationAtOnePoint(myid, rho, p);
    }

    // c110
    Vec3D c110 = 0.0;
    if(ijk_valid.second[3]) {
      myid = id[k][j+1][i+1];
      rho  =  v[k][j+1][(i+1)*dim];
      p    =  v[k][j+1][(i+1)*dim+4];
      c110 = ion->ComputeIonizationAtOnePoint(myid, rho, p);
    }

    // c001
    Vec3D c001 = 0.0;
    if(ijk_valid.second[4]) {
      myid = id[k+1][j][i];
      rho  =  v[k+1][j][i*dim];
      p    =  v[k+1][j][i*dim+4];
      c001 = ion->ComputeIonizationAtOnePoint(myid, rho, p);
    }

    // c101
    Vec3D c101 = 0.0;
    if(ijk_valid.second[5]) {
      myid = id[k+1][j][i+1];
      rho  =  v[k+1][j][(i+1)*dim];
      p    =  v[k+1][j][(i+1)*dim+4];
      c101 = ion->ComputeIonizationAtOnePoint(myid, rho, p);
    }

    // c011
    Vec3D c011 = 0.0;
    if(ijk_valid.second[6]) {
      myid = id[k+1][j+1][i];
      rho  =  v[k+1][j+1][i*dim];
      p    =  v[k+1][j+1][i*dim+4];
      c011 = ion->ComputeIonizationAtOnePoint(myid, rho, p);
    }

    // c111
    Vec3D c111 = 0.0;
    if(ijk_valid.second[7]) {
      myid = id[k+1][j+1][i+1];
      rho  =  v[k+1][j+1][(i+1)*dim];
      p    =  v[k+1][j+1][(i+1)*dim+4];
      c111 = ion->ComputeIonizationAtOnePoint(myid, rho, p);
    }


    if(ijk_valid.first<8) {//fill invalid slots with average value
      Vec3D c_avg = (c000+c100+c010+c110+c001+c101+c011+c111)/ijk_valid.first;
      if(!ijk_valid.second[0])  c000 = c_avg;
      if(!ijk_valid.second[1])  c100 = c_avg;
      if(!ijk_valid.second[2])  c010 = c_avg;
      if(!ijk_valid.second[3])  c110 = c_avg;
      if(!ijk_valid.second[4])  c001 = c_avg;
      if(!ijk_valid.second[5])  c101 = c_avg;
      if(!ijk_valid.second[6])  c011 = c_avg;
      if(!ijk_valid.second[7])  c111 = c_avg;
    }

    sol = MathTools::trilinear_interpolation(trilinear_coords, c000, c100, c010, c110, c001, c101, c011, c111);
  }

  MPI_Allreduce(MPI_IN_PLACE, (double*)sol, 3, MPI_DOUBLE, MPI_SUM, comm);
  return sol;
}

//-------------------------------------------------------------------------

//-------------------------------------------------------------------------

//-------------------------------------------------------------------------

