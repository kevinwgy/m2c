#include <MaterialVolumeOutput.h>
#include <IoData.h>

//-------------------------------------------------------------------------

MaterialVolumeOutput::MaterialVolumeOutput(MPI_Comm &comm_, IoData &iod, 
                                           SpaceVariable3D& cell_volume_) :
                      comm(comm_), cell_volume(cell_volume_)
{
  frequency = iod.output.materialVolumes.frequency;

  numMaterials = iod.eqs.materials.dataMap.size() + 1; //an extra one for "ghost/inactive" 

  // open file (all processors will open the same file, but only proc 0 will write to it)
  int spn = strlen(iod.output.prefix) + 1;
  if(iod.output.materialVolumes.filename[0] != 0) {
    char *full_fname = new char[spn + strlen(iod.output.materialVolumes.filename)];
    sprintf(full_fname, "%s%s", iod.output.prefix, iod.output.materialVolumes.filename);
    file = fopen(full_fname, "w");
    print(file, "## Number of physical materials: %d (0 - %d); ID for ghost/inactive cells: %d.\n",
                numMaterials-1, numMaterials-2, numMaterials-1);
    print(file, "## Total volume of computational domain: %e.\n",
                (iod.mesh.xmax-iod.mesh.x0)*(iod.mesh.ymax-iod.mesh.y0)*
                (iod.mesh.zmax-iod.mesh.z0));
    print(file, "## The last column is the sum over all material subdomains "
                "(including ghost/inactive).\n");
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
MaterialVolumeOutput::WriteSolution(double time, int time_step, SpaceVariable3D& ID, 
                                    bool must_write)
{

  if(frequency <=0 || file == NULL)
    return;

  if(time_step % frequency != 0 && !must_write) //should not output at this time step
    return;

  double vol[numMaterials];
  ComputeMaterialVolumes(ID, vol);

  print(file, "%8d  %12.8e  ", time_step, time);
  double sum = 0.0;
  for(int i=0; i<numMaterials; i++) {
    print(file, "%12.8e  ", vol[i]);
    sum += vol[i];
  }
  print(file, "%12.8e\n", sum);
  fflush(file);
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
          fprintf(stderr,"*** Error: Detected an unrecognized material id (%d)\n", 
                  myid);
          exit_mpi();
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






