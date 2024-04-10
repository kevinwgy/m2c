#include <EnergyIntegrationOutput.h>

using namespace std;

//---------------------------------------------------------------------------------------------------------------------------

// This constructor is for integration in a specified region
EnergyIntegrationOutput::EnergyIntegrationOutput(MPI_Comm &comm_, IoData &iod, OutputData &iod_output_, 
                                                 MeshData &iod_mesh_, EquationsData &iod_eqs_,
                                                 LaserAbsorptionSolver* laser_,
                                                 std::vector<VarFcnBase*> &vf_, SpaceVariable3D& coordinates_,
                                                 SpaceVariable3D& delta_xyz_, SpaceVariable3D& cell_volume_) : 
                                                 comm(comm_), iod_output(iod_output_), 
                                                 iod_mesh(iod_mesh_), iod_eqs(iod_eqs_),laser(laser_),vf(vf_),
                                                 coordinates(coordinates_),
                                                 delta_xyz(delta_xyz_),
                                                 cell_volume(cell_volume_) 
{

   frequency = iod_output.energy_integration.frequency;
   frequency_dt = iod_output.energy_integration.frequency_dt;

   last_snapshot_time = -1.0;

   numMaterials = iod.eqs.materials.dataMap.size() + 1; //an extra one for "ghost/inactive"

   for(int i=0; i<EnergyIntegrationData::SIZE; i++) {
    file[i] = NULL;
   }

   int spn = strlen(iod_output.prefix) + 1;

   if (iod_output.energy_integration.volume[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.energy_integration.volume)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_output.energy_integration.volume);
    file[EnergyIntegrationData::VOLUME] = fopen(filename, "w");
    delete [] filename;
   }

   if (iod_output.energy_integration.mass[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.energy_integration.mass)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_output.energy_integration.mass);
    file[EnergyIntegrationData::MASS] = fopen(filename, "w");
    delete [] filename;
   }

   if (iod_output.energy_integration.total_energy[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.energy_integration.total_energy)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_output.energy_integration.total_energy);
    file[EnergyIntegrationData::TOTAL_ENERGY] = fopen(filename, "w");
    delete [] filename;
   }

   if (iod_output.energy_integration.total_enthalpy[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.energy_integration.total_enthalpy)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_output.energy_integration.total_enthalpy);
    file[EnergyIntegrationData::TOTAL_ENTHALPY] = fopen(filename, "w");
    delete [] filename;
   }

   if (iod_output.energy_integration.kinetic_energy[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.energy_integration.kinetic_energy)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_output.energy_integration.kinetic_energy);
    file[EnergyIntegrationData::KINETIC_ENERGY] = fopen(filename, "w");
    delete [] filename;
   }

   if (iod_output.energy_integration.internal_energy[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.energy_integration.internal_energy)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_output.energy_integration.internal_energy);
    file[EnergyIntegrationData::INTERNAL_ENERGY] = fopen(filename, "w");
    delete [] filename;
   }
  
   if (iod_output.energy_integration.potential_energy[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.energy_integration.potential_energy)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_output.energy_integration.potential_energy);
    file[EnergyIntegrationData::POTENTIAL_ENERGY] = fopen(filename, "w");
    delete [] filename;
   }

   if (iod_output.energy_integration.laser_radiation[0] != 0) {
    char *filename = new char[spn + strlen(iod_output.energy_integration.laser_radiation)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_output.energy_integration.laser_radiation);
    file[EnergyIntegrationData::LASER_RADIATION] = fopen(filename, "w");
    delete [] filename;
   }

   for(int i=0; i<EnergyIntegrationData::SIZE; i++)
    if(file[i]) { //write header
      print(file[i], "## Number of physical materials: %d (0 - %d); ID for ghost/inactive cells: %d.\n",
                numMaterials-1, numMaterials-2, numMaterials-1);
      print(file[i], "## Time step  |  Time  ");
      for(int j=0; j<numMaterials; j++) {
         print(file[i], "|  Solutions(Material %d)  ", j);
      }     
      print(file[i],"   |  Sum(including ghost/inactive)\n");
      fflush(file[i]);
    }

}

//----------------------------------------------------------------------------------------------------

EnergyIntegrationOutput::~EnergyIntegrationOutput()
{
   for(int i=0; i<EnergyIntegrationData::SIZE; i++)
    if(file[i]) fclose(file[i]);
}

//------------------------------------------------------------------------------------------------------
void EnergyIntegrationOutput::WriteSolutionOfIntegrationEnergy(double time, double dt, int time_step, SpaceVariable3D &V, SpaceVariable3D &ID,
                                        SpaceVariable3D* L, bool force_write)
{
 
  if(!isTimeToWrite(time,dt,time_step,frequency_dt,frequency,last_snapshot_time,force_write))
    return;


  if(file[EnergyIntegrationData::VOLUME]) {
    print(file[EnergyIntegrationData::VOLUME], "%10d    %16.8e    ", time_step, time);
    double vol[numMaterials];
    double sum = 0.0;
    IntegrateVolume(ID,vol);
    for(int i=0; i<numMaterials; i++) {
      print(file[EnergyIntegrationData::VOLUME], "%16.8e  ", vol[i]);
      sum += vol[i];
    } 
    print(file[EnergyIntegrationData::VOLUME], "%16.8e\n", sum);
    fflush(file[EnergyIntegrationData::VOLUME]);
  }

  if(file[EnergyIntegrationData::MASS]) {
    print(file[EnergyIntegrationData::MASS], "%10d    %16.8e    ", time_step, time);
    double mass[numMaterials];
    double sum = 0.0;
    IntegrateMass(V,ID,mass);
    for(int i=0; i<numMaterials; i++) {
      print(file[EnergyIntegrationData::MASS], "%16.8e  ", mass[i]);
      sum += mass[i];
    }
    print(file[EnergyIntegrationData::MASS], "%16.8e\n", sum);
    fflush(file[EnergyIntegrationData::MASS]);
  } 

  if(file[EnergyIntegrationData::TOTAL_ENERGY]) {
    print(file[EnergyIntegrationData::TOTAL_ENERGY], "%10d    %16.8e    ", time_step, time);
    double E[numMaterials];
    double sum = 0.0;
    IntegrateTotalEnergy(V,ID,E);
    for(int i=0; i<numMaterials; i++) {
      print(file[EnergyIntegrationData::TOTAL_ENERGY], "%16.8e  ", E[i]);
      sum += E[i];
    }
    print(file[EnergyIntegrationData::TOTAL_ENERGY], "%16.8e\n", sum);
    fflush(file[EnergyIntegrationData::TOTAL_ENERGY]);
  }

  if(file[EnergyIntegrationData::TOTAL_ENTHALPY]) {
    print(file[EnergyIntegrationData::TOTAL_ENTHALPY], "%10d    %16.8e    ", time_step, time);
    double H[numMaterials];
    double sum = 0.0;
    IntegrateTotalEnthalpy(V,ID,H);
    for(int i=0; i<numMaterials; i++) {
      print(file[EnergyIntegrationData::TOTAL_ENTHALPY], "%16.8e  ", H[i]);
      sum += H[i];
    }
    print(file[EnergyIntegrationData::TOTAL_ENTHALPY], "%16.8e\n", sum);
    fflush(file[EnergyIntegrationData::TOTAL_ENTHALPY]);
  }

  if(file[EnergyIntegrationData::KINETIC_ENERGY]) {
    print(file[EnergyIntegrationData::KINETIC_ENERGY], "%10d    %16.8e    ", time_step, time);
    double kinetic[numMaterials];
    double sum = 0.0;
    IntegrateKineticEnergy(V,ID,kinetic);
    for(int i=0; i<numMaterials; i++) {
      print(file[EnergyIntegrationData::KINETIC_ENERGY], "%16.8e  ", kinetic[i]);
      sum += kinetic[i];
    }
    print(file[EnergyIntegrationData::KINETIC_ENERGY], "%16.8e\n", sum);
    fflush(file[EnergyIntegrationData::KINETIC_ENERGY]);
  }

  if(file[EnergyIntegrationData::INTERNAL_ENERGY]) {
    print(file[EnergyIntegrationData::INTERNAL_ENERGY], "%10d    %16.8e    ", time_step, time);
    double internal[numMaterials];
    double sum = 0.0;
    IntegrateInternalEnergy(V,ID,internal);
    for(int i=0; i<numMaterials; i++) {
      print(file[EnergyIntegrationData::INTERNAL_ENERGY], "%16.8e  ", internal[i]);
      sum += internal[i];
    }
    print(file[EnergyIntegrationData::INTERNAL_ENERGY], "%16.8e\n", sum);
    fflush(file[EnergyIntegrationData::INTERNAL_ENERGY]);
  }

  if(file[EnergyIntegrationData::POTENTIAL_ENERGY]) {
    print(file[EnergyIntegrationData::POTENTIAL_ENERGY], "%10d    %16.8e    ", time_step, time);
    double potential[numMaterials];
    double sum = 0.0;
    IntegratePotentialEnergy(V,ID,potential);
    for(int i=0; i<numMaterials; i++) {
      print(file[EnergyIntegrationData::POTENTIAL_ENERGY], "%16.8e  ", potential[i]);
      sum += potential[i];
    }
    print(file[EnergyIntegrationData::POTENTIAL_ENERGY], "%16.8e\n", sum);
    fflush(file[EnergyIntegrationData::POTENTIAL_ENERGY]);
  }

  if(file[EnergyIntegrationData::LASER_RADIATION]) {
    print(file[EnergyIntegrationData::LASER_RADIATION], "%10d    %16.8e    ", time_step, time);
    if(L == NULL) {
      print_error("*** Error: Requested laser radiation integration, but laser source is not specified.\n");
      exit_mpi();
    }
    double radiation[numMaterials];
    double sum = 0.0;
    IntegrateLaserRadiation(V,ID,L,radiation);
    for(int i=0; i<numMaterials; i++) {
      print(file[EnergyIntegrationData::LASER_RADIATION], "%16.8e  ", radiation[i]);
      sum += radiation[i];
    }
    print(file[EnergyIntegrationData::LASER_RADIATION], "%16.8e\n", sum);
    fflush(file[EnergyIntegrationData::LASER_RADIATION]);
  }

  last_snapshot_time = time;
    
}   
    
//--------------------------------------------------------------------------------------------------------------

void EnergyIntegrationOutput::IntegrateVolume(SpaceVariable3D &ID, double* vol)     
{
  double*** id   = (double***)ID.GetDataPointer();
  double*** cell = (double***)cell_volume.GetDataPointer();
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
  Vec3D*** dxyz = (Vec3D***)delta_xyz.GetDataPointer();

  int i0, j0, k0, imax, jmax, kmax;
  ID.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);

  //initialize volumes to 0
  for(int i=0; i<numMaterials; i++)
      vol[i] = 0.0;
 
  // loop through the true domain (excluding the ghost layer)
  double PI = acos(0.0)*2.0;
  int myid;
  for(int k=k0; k<kmax; k++) {
    for(int j=j0; j<jmax; j++) {
      for(int i=i0; i<imax; i++) {
        myid = id[k][j][i];
        if(myid<0 || myid>=numMaterials) {
          fprintf(stderr,"*** Error: Detected an unrecognized material id (%d)\n",
                  myid);
          exit(-1);
        }
        if(coords[k][j][i][0] >= iod_output.energy_integration.x_min &&
           coords[k][j][i][0] <= iod_output.energy_integration.x_max)
          if(coords[k][j][i][1] >= iod_output.energy_integration.y_min &&
             coords[k][j][i][1] <= iod_output.energy_integration.y_max)
            if(coords[k][j][i][2] >= iod_output.energy_integration.z_min &&
               coords[k][j][i][2] <= iod_output.energy_integration.z_max){
               if(iod_mesh.type == MeshData::SPHERICAL){
                 double scale = PI*4.0*coords[k][j][i][0]*coords[k][j][i][0]/dxyz[k][j][i][2]/dxyz[k][j][i][1];
                 vol[myid] += cell[k][j][i]*scale;
               }
               else if(iod_mesh.type == MeshData::CYLINDRICAL){
                 double scalor = PI*2.0*coords[k][j][i][1]/dxyz[k][j][i][2];
                 vol[myid] += cell[k][j][i]*scalor;
               }
               else
                 vol[myid] += cell[k][j][i];
            }
      }
    }
  }

  //MPI communication
  MPI_Allreduce(MPI_IN_PLACE, vol, numMaterials, MPI_DOUBLE, MPI_SUM, comm);

  ID.RestoreDataPointerToLocalVector();
  cell_volume.RestoreDataPointerToLocalVector();
  coordinates.RestoreDataPointerToLocalVector();
  delta_xyz.RestoreDataPointerToLocalVector();
 
}


//--------------------------------------------------------------------------------------------------------------

void EnergyIntegrationOutput::IntegrateMass(SpaceVariable3D &V, SpaceVariable3D &ID, double* mass)
{
  Vec5D***  v  = (Vec5D***) V.GetDataPointer();
  double*** id   = (double***)ID.GetDataPointer();
  double*** cell = (double***)cell_volume.GetDataPointer();
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
  Vec3D*** dxyz = (Vec3D***)delta_xyz.GetDataPointer();  

  int i0, j0, k0, imax, jmax, kmax;
  ID.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);

  for(int i=0; i<numMaterials; i++)
      mass[i] = 0.0;

  double PI = acos(0.0)*2.0;
  int myid;
  for(int k=k0; k<kmax; k++) {
    for(int j=j0; j<jmax; j++) {
      for(int i=i0; i<imax; i++) {
        myid = id[k][j][i];
        if(myid<0 || myid>=numMaterials) {
          fprintf(stderr,"*** Error: Detected an unrecognized material id (%d)\n",
                  myid);
          exit(-1);
        }
        if(coords[k][j][i][0] >= iod_output.energy_integration.x_min &&
           coords[k][j][i][0] <= iod_output.energy_integration.x_max)
          if(coords[k][j][i][1] >= iod_output.energy_integration.y_min &&
             coords[k][j][i][1] <= iod_output.energy_integration.y_max)
            if(coords[k][j][i][2] >= iod_output.energy_integration.z_min &&
               coords[k][j][i][2] <= iod_output.energy_integration.z_max){
               if(iod_mesh.type == MeshData::SPHERICAL){
                 double scale = PI*4.0*coords[k][j][i][0]*coords[k][j][i][0]/dxyz[k][j][i][2]/dxyz[k][j][i][1];
                 mass[myid] += v[k][j][i][0]*cell[k][j][i]*scale;
               }
               else if(iod_mesh.type == MeshData::CYLINDRICAL){
                 double scalor = PI*2.0*coords[k][j][i][1]/dxyz[k][j][i][2];
                 mass[myid] += v[k][j][i][0]*cell[k][j][i]*scalor;
               }
               else
                 mass[myid] += v[k][j][i][0]*cell[k][j][i];
            }
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, mass, numMaterials, MPI_DOUBLE, MPI_SUM, comm);

  V.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();
  cell_volume.RestoreDataPointerToLocalVector();
  coordinates.RestoreDataPointerToLocalVector();
  delta_xyz.RestoreDataPointerToLocalVector();

}

//----------------------------------------------------------------------------------------------------------------

void EnergyIntegrationOutput::IntegrateTotalEnergy(SpaceVariable3D &V, SpaceVariable3D &ID, double* E)
{
  Vec5D***  v  = (Vec5D***) V.GetDataPointer();
  double*** id   = (double***)ID.GetDataPointer();
  double*** cell = (double***)cell_volume.GetDataPointer();
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
  Vec3D*** dxyz = (Vec3D***)delta_xyz.GetDataPointer();

  int i0, j0, k0, imax, jmax, kmax;
  ID.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);

  for(int i=0; i<numMaterials; i++)
      E[i] = 0.0;

  double PI = acos(0.0)*2.0;
  int myid;
  for(int k=k0; k<kmax; k++) {
    for(int j=j0; j<jmax; j++) {
      for(int i=i0; i<imax; i++) {
        myid = id[k][j][i];
        if(myid<0 || myid>=numMaterials) {
          fprintf(stderr,"*** Error: Detected an unrecognized material id (%d)\n",
                  myid);
          exit(-1);
        }
        if(coords[k][j][i][0] >= iod_output.energy_integration.x_min &&
           coords[k][j][i][0] <= iod_output.energy_integration.x_max)
          if(coords[k][j][i][1] >= iod_output.energy_integration.y_min &&
             coords[k][j][i][1] <= iod_output.energy_integration.y_max)
            if(coords[k][j][i][2] >= iod_output.energy_integration.z_min &&
               coords[k][j][i][2] <= iod_output.energy_integration.z_max){
               double vsquare = v[k][j][i][1]*v[k][j][i][1] + v[k][j][i][2]*v[k][j][i][2]
                              + v[k][j][i][3]*v[k][j][i][3];
               double e = vf[myid]->GetInternalEnergyPerUnitMass(v[k][j][i][0], v[k][j][i][4]);
               if(iod_mesh.type == MeshData::SPHERICAL){
                 double scale = PI*4.0*coords[k][j][i][0]*coords[k][j][i][0]/dxyz[k][j][i][2]/dxyz[k][j][i][1];
                 E[myid] += v[k][j][i][0]*(0.5*vsquare+e)*cell[k][j][i]*scale;
               }
               else if(iod_mesh.type == MeshData::CYLINDRICAL){
                 double scalor = PI*2.0*coords[k][j][i][1]/dxyz[k][j][i][2];
                 E[myid] += v[k][j][i][0]*(0.5*vsquare+e)*cell[k][j][i]*scalor;
               }
               else
                 E[myid] += v[k][j][i][0]*(0.5*vsquare+e)*cell[k][j][i];
            }
      }
    }
  }   


  MPI_Allreduce(MPI_IN_PLACE, E, numMaterials, MPI_DOUBLE, MPI_SUM, comm);

  V.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();
  cell_volume.RestoreDataPointerToLocalVector();
  coordinates.RestoreDataPointerToLocalVector();
  delta_xyz.RestoreDataPointerToLocalVector();

}

//-----------------------------------------------------------------------------------------------------------------

void EnergyIntegrationOutput::IntegrateTotalEnthalpy(SpaceVariable3D &V, SpaceVariable3D &ID, double* H)
{
  Vec5D***  v  = (Vec5D***) V.GetDataPointer();
  double*** id   = (double***)ID.GetDataPointer();
  double*** cell = (double***)cell_volume.GetDataPointer();
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
  Vec3D*** dxyz = (Vec3D***)delta_xyz.GetDataPointer();

  int i0, j0, k0, imax, jmax, kmax;
  ID.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);

  for(int i=0; i<numMaterials; i++)
      H[i] = 0.0;

  double PI = acos(0.0)*2.0;
  int myid;
  for(int k=k0; k<kmax; k++) {
    for(int j=j0; j<jmax; j++) {
      for(int i=i0; i<imax; i++) {
        myid = id[k][j][i];
        if(myid<0 || myid>=numMaterials) {
          fprintf(stderr,"*** Error: Detected an unrecognized material id (%d)\n",
                  myid);
          exit(-1);
        }
        if(coords[k][j][i][0] >= iod_output.energy_integration.x_min &&
           coords[k][j][i][0] <= iod_output.energy_integration.x_max)
          if(coords[k][j][i][1] >= iod_output.energy_integration.y_min &&
             coords[k][j][i][1] <= iod_output.energy_integration.y_max)
            if(coords[k][j][i][2] >= iod_output.energy_integration.z_min &&
               coords[k][j][i][2] <= iod_output.energy_integration.z_max){
               double vsquare = v[k][j][i][1]*v[k][j][i][1] + v[k][j][i][2]*v[k][j][i][2]
                              + v[k][j][i][3]*v[k][j][i][3];
               double e = vf[myid]->GetInternalEnergyPerUnitMass(v[k][j][i][0], v[k][j][i][4]);
               if(iod_mesh.type == MeshData::SPHERICAL){
                 double scale = PI*4.0*coords[k][j][i][0]*coords[k][j][i][0]/dxyz[k][j][i][2]/dxyz[k][j][i][1];
                 H[myid] += (v[k][j][i][0]*(0.5*vsquare+e)+abs(v[k][j][i][4]))*cell[k][j][i]*scale;
               }
               else if(iod_mesh.type == MeshData::CYLINDRICAL){
                 double scalor = PI*2.0*coords[k][j][i][1]/dxyz[k][j][i][2];
                 H[myid] += (v[k][j][i][0]*(0.5*vsquare+e)+abs(v[k][j][i][4]))*cell[k][j][i]*scalor;
               }
               else
                 H[myid] += (v[k][j][i][0]*(0.5*vsquare+e)+abs(v[k][j][i][4]))*cell[k][j][i];
            }
      }
    }
  }

 
  MPI_Allreduce(MPI_IN_PLACE, H, numMaterials, MPI_DOUBLE, MPI_SUM, comm);
 
  V.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();
  cell_volume.RestoreDataPointerToLocalVector();
  coordinates.RestoreDataPointerToLocalVector();
  delta_xyz.RestoreDataPointerToLocalVector();

}

//-----------------------------------------------------------------------------------------------------------------

void EnergyIntegrationOutput::IntegrateKineticEnergy(SpaceVariable3D &V, SpaceVariable3D &ID, double* kinetic)
{
  Vec5D***  v  = (Vec5D***) V.GetDataPointer();
  double*** id   = (double***)ID.GetDataPointer();
  double*** cell = (double***)cell_volume.GetDataPointer();
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
  Vec3D*** dxyz = (Vec3D***)delta_xyz.GetDataPointer();

  int i0, j0, k0, imax, jmax, kmax;
  ID.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);

  for(int i=0; i<numMaterials; i++)
      kinetic[i] = 0.0;

  double PI = acos(0.0)*2.0;
  int myid;
  for(int k=k0; k<kmax; k++) {
    for(int j=j0; j<jmax; j++) {
      for(int i=i0; i<imax; i++) {
        myid = id[k][j][i];
        if(myid<0 || myid>=numMaterials) {
          fprintf(stderr,"*** Error: Detected an unrecognized material id (%d)\n",
                  myid);
          exit(-1);
        }
        if(coords[k][j][i][0] >= iod_output.energy_integration.x_min &&
           coords[k][j][i][0] <= iod_output.energy_integration.x_max)
          if(coords[k][j][i][1] >= iod_output.energy_integration.y_min &&
             coords[k][j][i][1] <= iod_output.energy_integration.y_max)
            if(coords[k][j][i][2] >= iod_output.energy_integration.z_min &&
               coords[k][j][i][2] <= iod_output.energy_integration.z_max){
               double vsquare = v[k][j][i][1]*v[k][j][i][1] + v[k][j][i][2]*v[k][j][i][2]
                              + v[k][j][i][3]*v[k][j][i][3];
               if(iod_mesh.type == MeshData::SPHERICAL){
                 double scale = PI*4.0*coords[k][j][i][0]*coords[k][j][i][0]/dxyz[k][j][i][2]/dxyz[k][j][i][1];
                 kinetic[myid] += 0.5*v[k][j][i][0]*vsquare*cell[k][j][i]*scale;
               }
               else if(iod_mesh.type == MeshData::CYLINDRICAL){
                 double scalor = PI*2.0*coords[k][j][i][1]/dxyz[k][j][i][2];
                 kinetic[myid] += 0.5*v[k][j][i][0]*vsquare*cell[k][j][i]*scalor;
               }
               else
                 kinetic[myid] += 0.5*v[k][j][i][0]*vsquare*cell[k][j][i];
            }
      }
    }
  }


  MPI_Allreduce(MPI_IN_PLACE, kinetic, numMaterials, MPI_DOUBLE, MPI_SUM, comm);

  V.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();
  cell_volume.RestoreDataPointerToLocalVector();
  coordinates.RestoreDataPointerToLocalVector();
  delta_xyz.RestoreDataPointerToLocalVector();

}

//-----------------------------------------------------------------------------------------------------------------

void EnergyIntegrationOutput::IntegrateInternalEnergy(SpaceVariable3D &V, SpaceVariable3D &ID, double* internal)
{
  Vec5D***  v  = (Vec5D***) V.GetDataPointer();
  double*** id   = (double***)ID.GetDataPointer();
  double*** cell = (double***)cell_volume.GetDataPointer();
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
  Vec3D*** dxyz = (Vec3D***)delta_xyz.GetDataPointer();

  int i0, j0, k0, imax, jmax, kmax;
  ID.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);

  for(int i=0; i<numMaterials; i++)
      internal[i] = 0.0;

  double PI = acos(0.0)*2.0;
  int myid;
  for(int k=k0; k<kmax; k++) {
    for(int j=j0; j<jmax; j++) {
      for(int i=i0; i<imax; i++) {
        myid = id[k][j][i];
        if(myid<0 || myid>=numMaterials) {
          fprintf(stderr,"*** Error: Detected an unrecognized material id (%d)\n",
                  myid);
          exit(-1);
        }
        if(coords[k][j][i][0] >= iod_output.energy_integration.x_min &&
           coords[k][j][i][0] <= iod_output.energy_integration.x_max)
          if(coords[k][j][i][1] >= iod_output.energy_integration.y_min &&
             coords[k][j][i][1] <= iod_output.energy_integration.y_max)
            if(coords[k][j][i][2] >= iod_output.energy_integration.z_min &&
               coords[k][j][i][2] <= iod_output.energy_integration.z_max){
               double e = vf[myid]->GetInternalEnergyPerUnitMass(v[k][j][i][0], v[k][j][i][4]);
               if(iod_mesh.type == MeshData::SPHERICAL){
                 double scale = PI*4.0*coords[k][j][i][0]*coords[k][j][i][0]/dxyz[k][j][i][2]/dxyz[k][j][i][1];
                 internal[myid] += v[k][j][i][0]*e*cell[k][j][i]*scale;
               }
               else if(iod_mesh.type == MeshData::CYLINDRICAL){
                 double scalor = PI*2.0*coords[k][j][i][1]/dxyz[k][j][i][2];
                 internal[myid] += v[k][j][i][0]*e*cell[k][j][i]*scalor;
               }
               else
                 internal[myid] += v[k][j][i][0]*e*cell[k][j][i];
            }
      }
    }
  }  

 
  MPI_Allreduce(MPI_IN_PLACE, internal, numMaterials, MPI_DOUBLE, MPI_SUM, comm);
 
  V.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();
  cell_volume.RestoreDataPointerToLocalVector();
  coordinates.RestoreDataPointerToLocalVector();
  delta_xyz.RestoreDataPointerToLocalVector();

}
    
 //-----------------------------------------------------------------------------------------------------------------------------------------
 
void EnergyIntegrationOutput::IntegratePotentialEnergy(SpaceVariable3D &V, SpaceVariable3D &ID, double* potential)
{
  Vec5D***  v  = (Vec5D***) V.GetDataPointer();
  double*** id   = (double***)ID.GetDataPointer();
  double*** cell = (double***)cell_volume.GetDataPointer();
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
  Vec3D*** dxyz = (Vec3D***)delta_xyz.GetDataPointer();

  int i0, j0, k0, imax, jmax, kmax;
  ID.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);

  for(int i=0; i<numMaterials; i++)
      potential[i] = 0.0;

  double PI = acos(0.0)*2.0;
  int myid;
  for(int k=k0; k<kmax; k++) {
    for(int j=j0; j<jmax; j++) {
      for(int i=i0; i<imax; i++) {
        myid = id[k][j][i];
        if(myid<0 || myid>=numMaterials) {
          fprintf(stderr,"*** Error: Detected an unrecognized material id (%d)\n",
                  myid);
          exit(-1);
        }
        if(coords[k][j][i][0] >= iod_output.energy_integration.x_min &&
           coords[k][j][i][0] <= iod_output.energy_integration.x_max)
          if(coords[k][j][i][1] >= iod_output.energy_integration.y_min &&
             coords[k][j][i][1] <= iod_output.energy_integration.y_max)
            if(coords[k][j][i][2] >= iod_output.energy_integration.z_min &&
               coords[k][j][i][2] <= iod_output.energy_integration.z_max){
               if(iod_mesh.type == MeshData::SPHERICAL){
                 double scale = PI*4.0*coords[k][j][i][0]*coords[k][j][i][0]/dxyz[k][j][i][2]/dxyz[k][j][i][1];
                 potential[myid] += abs(v[k][j][i][4])*cell[k][j][i]*scale;
               }
               else if(iod_mesh.type == MeshData::CYLINDRICAL){
                 double scalor = PI*2.0*coords[k][j][i][1]/dxyz[k][j][i][2];
                 potential[myid] += abs(v[k][j][i][4])*cell[k][j][i]*scalor;
               }
               else
                 potential[myid] += abs(v[k][j][i][4])*cell[k][j][i];
            }
      }
    }
  } 


  MPI_Allreduce(MPI_IN_PLACE, potential, numMaterials, MPI_DOUBLE, MPI_SUM, comm);

  V.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();
  cell_volume.RestoreDataPointerToLocalVector();
  coordinates.RestoreDataPointerToLocalVector();
  delta_xyz.RestoreDataPointerToLocalVector();

}

//-----------------------------------------------------------------------------------------------------------------------------------------------
void EnergyIntegrationOutput::IntegrateLaserRadiation(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D *L, double* radiation)
{
  Vec5D***  v  = (Vec5D***) V.GetDataPointer();
  double*** id   = (double***)ID.GetDataPointer();
  double*** cell = (double***)cell_volume.GetDataPointer();
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
  Vec3D*** dxyz = (Vec3D***)delta_xyz.GetDataPointer();
  double*** l  = (double***)L->GetDataPointer();

  int i0, j0, k0, imax, jmax, kmax;
  ID.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);

  for(int i=0; i<numMaterials; i++)
      radiation[i] = 0.0;

  double PI = acos(0.0)*2.0;
  int myid;
  double e, myT, eta;
  for(int k=k0; k<kmax; k++) {
    for(int j=j0; j<jmax; j++) {
      for(int i=i0; i<imax; i++) {
        myid = id[k][j][i];
        if(myid<0 || myid>=numMaterials) {
          fprintf(stderr,"*** Error: Detected an unrecognized material id (%d)\n",
                  myid);
          exit(-1);
        }
        if(coords[k][j][i][0] >= iod_output.energy_integration.x_min &&
           coords[k][j][i][0] <= iod_output.energy_integration.x_max)
          if(coords[k][j][i][1] >= iod_output.energy_integration.y_min &&
             coords[k][j][i][1] <= iod_output.energy_integration.y_max)
            if(coords[k][j][i][2] >= iod_output.energy_integration.z_min &&
               coords[k][j][i][2] <= iod_output.energy_integration.z_max){
               if(iod_mesh.type == MeshData::SPHERICAL){
                 double scale = PI*4.0*coords[k][j][i][0]*coords[k][j][i][0]/dxyz[k][j][i][2]/dxyz[k][j][i][1];
                 e = vf[myid]->GetInternalEnergyPerUnitMass(v[k][j][i][0], v[k][j][i][4]);
                 myT = vf[myid]->GetTemperature(v[k][j][i][0], e);
                 eta = laser->GetAbsorptionCoefficient(myT, myid);
                 radiation[myid] += eta*l[k][j][i]*cell[k][j][i]*scale;
               }
               else if(iod_mesh.type == MeshData::CYLINDRICAL){
                 double scalor = PI*2.0*coords[k][j][i][1]/dxyz[k][j][i][2];
                 e = vf[myid]->GetInternalEnergyPerUnitMass(v[k][j][i][0], v[k][j][i][4]);
                 myT = vf[myid]->GetTemperature(v[k][j][i][0], e);
                 eta = laser->GetAbsorptionCoefficient(myT, myid);
                 radiation[myid] += eta*l[k][j][i]*cell[k][j][i]*scalor;
               }
               else{
                 e = vf[myid]->GetInternalEnergyPerUnitMass(v[k][j][i][0], v[k][j][i][4]);
                 myT = vf[myid]->GetTemperature(v[k][j][i][0], e);
                 eta = laser->GetAbsorptionCoefficient(myT, myid);
                 radiation[myid] += eta*l[k][j][i]*cell[k][j][i];
               }
            }
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, radiation, numMaterials, MPI_DOUBLE, MPI_SUM, comm);

  V.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();
  cell_volume.RestoreDataPointerToLocalVector();
  coordinates.RestoreDataPointerToLocalVector();
  delta_xyz.RestoreDataPointerToLocalVector();
  L->RestoreDataPointerToLocalVector();


}
