/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <MaterialVolumeOutput.h>
#include <IoData.h>

//-------------------------------------------------------------------------

MaterialVolumeOutput::MaterialVolumeOutput(MPI_Comm &comm_, IoData &iod, 
                                           SpaceVariable3D& cell_volume_) :
                      comm(comm_), cell_volume(cell_volume_)
{
  frequency = iod.output.materialVolumes.frequency;
  frequency_dt = iod.output.materialVolumes.frequency_dt;

  last_snapshot_time = -1.0;

  numMaterials = iod.eqs.materials.dataMap.size() + 1; //an extra one for "ghost/inactive" 

  // open file (all processors will open the same file, but only proc 0 will write to it)
  int spn = strlen(iod.output.prefix) + 1;
  if(iod.output.materialVolumes.filename[0] != 0) {
    char *full_fname = new char[spn + strlen(iod.output.materialVolumes.filename)];
    sprintf(full_fname, "%s%s", iod.output.prefix, iod.output.materialVolumes.filename);
    file = fopen(full_fname, "w");

    if(!file) {
      print_error("*** Error: Cannot open file '%s' for output.\n", full_fname);
      exit_mpi();
    }

    print(file, "## Number of physical materials: %d (0 - %d); ID for ghost/inactive cells: %d.\n",
                numMaterials-1, numMaterials-2, numMaterials-1);
    print(file, "## Total volume of computational domain: %e.\n",
                (iod.mesh.xmax-iod.mesh.x0)*(iod.mesh.ymax-iod.mesh.y0)*
                (iod.mesh.zmax-iod.mesh.z0));
    print(file, "## The last column is the sum over all material subdomains "
                "(including ghost/inactive).\n");
    mpi_barrier();
    fflush(file);
    delete [] full_fname;
  } 
  else
    file = NULL;

}

//-------------------------------------------------------------------------

MaterialVolumeOutput::~MaterialVolumeOutput()
{
  if(file)
    fclose(file);
}

//-------------------------------------------------------------------------

void
MaterialVolumeOutput::WriteSolution(double time, double dt, int time_step, SpaceVariable3D& ID, 
                                    bool force_write)
{

  if(file == NULL)
    return;

  if(!isTimeToWrite(time,dt,time_step,frequency_dt,frequency,last_snapshot_time,force_write))
    return;

  double vol[numMaterials];
  ComputeMaterialVolumes(ID, vol);

  print(file, "%8d  %16.8e  ", time_step, time);
  double sum = 0.0;
  for(int i=0; i<numMaterials; i++) {
    print(file, "%16.8e  ", vol[i]);
    sum += vol[i];
  }
  print(file, "%16.8e\n", sum);
  mpi_barrier();
  fflush(file);

  last_snapshot_time = time;

}

//-------------------------------------------------------------------------

void
MaterialVolumeOutput::ComputeMaterialVolumes(SpaceVariable3D& ID, double* vol)
{

  double*** id   = (double***)ID.GetDataPointer();
  double*** cell = (double***)cell_volume.GetDataPointer();

  int i0, j0, k0, imax, jmax, kmax;
  ID.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);

  //initialize volumes to 0
  for(int i=0; i<numMaterials; i++)
    vol[i] = 0.0;

  // loop through the true domain (excluding the ghost layer)
  int myid;
  for(int k=k0; k<kmax; k++) {
    for(int j=j0; j<jmax; j++) {
      for(int i=i0; i<imax; i++) {
        myid = id[k][j][i]; 
        if(myid<0 || myid>=numMaterials) {
          fprintf(stdout,"*** Error: Detected an unrecognized material id (%d)\n", 
                  myid);
          exit(-1);
        }
        vol[myid] += cell[k][j][i];
      }
    }
  }

  //MPI communication
  MPI_Allreduce(MPI_IN_PLACE, vol, numMaterials, MPI_DOUBLE, MPI_SUM, comm);

  ID.RestoreDataPointerToLocalVector();
  cell_volume.RestoreDataPointerToLocalVector();

}

//-------------------------------------------------------------------------






