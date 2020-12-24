#include <Utils.h>
#include <Vector5D.h>
#include <Output.h>

//--------------------------------------------------------------------------

Output::Output(MPI_Comm &comm_, DataManagers2D &dms, IoData &iod_, VarFcnBase &vf_) : comm(comm_), 
    iod(iod_), vf(vf_),
    scalar(comm_, &(dms.ghosted1_1dof)),
    vector3(comm_, &(dms.ghosted1_3dof))
{
  iFrame = 0;

  last_snapshot_time = -1.0;

  char f1[256];
  sprintf(f1, "%s%s.pvd", iod.output.prefix, iod.output.solution_filename_base);
  pvdfile  = fopen(f1,"w");
  print(pvdfile, "<?xml version=\"1.0\"?>\n");
  print(pvdfile, "<VTKFile type=\"Collection\" version=\"0.1\"\n");
  print(pvdfile, "byte_order=\"LittleEndian\">\n");
  print(pvdfile, "  <Collection>\n");
}

//--------------------------------------------------------------------------

Output::~Output()
{
  if(pvdfile) fclose(pvdfile);
  scalar.Destroy();
  vector3.Destroy();
}

//--------------------------------------------------------------------------

void Output::InitializeOutput(SpaceVariable2D &coordinates)
{
  scalar.StoreMeshCoordinates(coordinates);
  vector3.StoreMeshCoordinates(coordinates);
}

//--------------------------------------------------------------------------

bool Output::ToWriteSolutionSnapshot(double time, double dt, int time_step)
{
  //! First check frequency_dt. If it is not specified, use frequency
  if(iod.output.frequency_dt > 0) {
    if(time - last_snapshot_time >= iod.output.frequency_dt - 0.01*dt/*a small tolerance*/)
      return true;
    else
      return false;
  } 
  else { //!< use frequency
    if((iod.output.frequency > 0) && (time_step % iod.output.frequency == 0))
      return true;
    else
      return false;
  } 
}

//--------------------------------------------------------------------------

void Output::WriteSolutionSnapshot(double time, int time_step, SpaceVariable2D &V)
{
  //! Define vtr file name
  char full_fname[256];
  char fname[256];
  if(iFrame<10) {
    sprintf(fname, "%s_000%d.vtr", 
            iod.output.solution_filename_base, iFrame);
    sprintf(full_fname, "%s%s_000%d.vtr", iod.output.prefix,
            iod.output.solution_filename_base, iFrame); 
  }
  else if(iFrame<100){
    sprintf(fname, "%s_00%d.vtr", 
            iod.output.solution_filename_base, iFrame);
    sprintf(full_fname, "%s%s_00%d.vtr", iod.output.prefix,
            iod.output.solution_filename_base, iFrame); 
  }
  else if(iFrame<1000){
    sprintf(fname, "%s_0%d.vtr", 
            iod.output.solution_filename_base, iFrame);
    sprintf(full_fname, "%s%s_0%d.vtr", iod.output.prefix,
            iod.output.solution_filename_base, iFrame); 
  }
  else {
    sprintf(fname, "%s_%d.vtr", 
            iod.output.solution_filename_base, iFrame);
    sprintf(full_fname, "%s%s_%d.vtr", iod.output.prefix,
            iod.output.solution_filename_base, iFrame); 
  }


  // Write solution snapshot
  PetscViewer viewer;
  PetscViewerVTKOpen(PETSC_COMM_WORLD, full_fname, FILE_MODE_WRITE, &viewer); 

  Vec5D**  v  = (Vec5D**) V.GetDataPointer();

  int i0, j0, imax, jmax;
  V.GetCornerIndices(&i0, &j0, &imax, &jmax);

  if(iod.output.density==OutputData::ON) {
    double** s  = (double**) scalar.GetDataPointer();
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++)
        s[j][i] = v[j][i][0];
    scalar.RestoreDataPointerAndInsert();
    PetscObjectSetName((PetscObject)(scalar.GetRefToGlobalVec()), "density");
    VecView(scalar.GetRefToGlobalVec(), viewer);
  }

  if(iod.output.velocity==OutputData::ON) {
    Vec3D** v3  = (Vec3D**) vector3.GetDataPointer();
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++)
        for(int k=0; k<3; k++)
          v3[j][i][k] = v[j][i][1+k];
    vector3.RestoreDataPointerAndInsert();
    PetscObjectSetName((PetscObject)(vector3.GetRefToGlobalVec()), "velocity");
    VecView(vector3.GetRefToGlobalVec(), viewer);
  }

  if(iod.output.pressure==OutputData::ON) {
    double** s  = (double**) scalar.GetDataPointer();
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++)
        s[j][i] = v[j][i][4];
    scalar.RestoreDataPointerAndInsert();
    PetscObjectSetName((PetscObject)(scalar.GetRefToGlobalVec()), "pressure");
    VecView(scalar.GetRefToGlobalVec(), viewer);
  }

  // Add a line to the pvd file to record the new solutio snapshot
  print(pvdfile, "  <DataSet timestep=\"%e\" file=\"%s\"/>\n", time, fname);

  // clean up
  PetscViewerDestroy(&viewer);
  V.RestoreDataPointerToLocalVector(); //no changes made to V.

  // bookkeeping
  iFrame++;
  last_snapshot_time = time;
   
}

//--------------------------------------------------------------------------

void Output::FinalizeOutput()
{
  print(pvdfile, "  </Collection>\n");
  print(pvdfile, "</VTKFile>\n");

  fclose(pvdfile); pvdfile = NULL;
}

//--------------------------------------------------------------------------

















