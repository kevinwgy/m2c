#include <Utils.h>
#include <Vector5D.h>
#include <Output.h>

//--------------------------------------------------------------------------

Output::Output(MPI_Comm &comm_, DataManagers3D &dms, IoData &iod_, vector<VarFcnBase*> &vf_) : comm(comm_), 
    iod(iod_), vf(vf_),
    scalar(comm_, &(dms.ghosted1_1dof)),
    vector3(comm_, &(dms.ghosted1_3dof)),
    probe_output(comm_, iod_.output)
{
  iFrame = 0;

  last_snapshot_time = -1.0;

  char f1[256];
  sprintf(f1, "%s%s.pvd", iod.output.prefix, iod.output.solution_filename_base);

  pvdfile  = fopen(f1,"w");
  if(!pvdfile) {
    print_error("*** Error: Cannot open file '%s%s.pvd' for output.\n", iod.output.prefix, iod.output.solution_filename_base);
    exit_mpi();
  }

  print(pvdfile, "<?xml version=\"1.0\"?>\n");
  print(pvdfile, "<VTKFile type=\"Collection\" version=\"0.1\"\n");
  print(pvdfile, "byte_order=\"LittleEndian\">\n");
  print(pvdfile, "  <Collection>\n");

  print(pvdfile, "  </Collection>\n");
  print(pvdfile, "</VTKFile>\n");

  fclose(pvdfile); pvdfile = NULL;

}

//--------------------------------------------------------------------------

Output::~Output()
{
  if(pvdfile) fclose(pvdfile);
}

//--------------------------------------------------------------------------

void Output::InitializeOutput(SpaceVariable3D &coordinates)
{
  scalar.StoreMeshCoordinates(coordinates);
  vector3.StoreMeshCoordinates(coordinates);
  probe_output.SetupInterpolation(coordinates);
}

//--------------------------------------------------------------------------

bool Output::ToWriteSolutionSnapshot(double time, double dt, int time_step)
{
  //! First check frequency_dt. If it is not specified, use frequency
  if(iod.output.frequency_dt > 0) {
    if(floor(time/iod.output.frequency_dt) != floor((time-dt)/iod.output.frequency_dt))
//    if(time - last_snapshot_time >= iod.output.frequency_dt - 0.01*dt/*a small tolerance*/)
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

void Output::WriteSolutionSnapshot(double time, int time_step, SpaceVariable3D &V, 
                                   SpaceVariable3D &ID, std::vector<SpaceVariable3D*> &Phi)
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

  //Open file
  PetscViewer viewer;
  int code = PetscViewerVTKOpen(PETSC_COMM_WORLD, full_fname, FILE_MODE_WRITE, &viewer); 
  if(code) {
    print_error("*** Error: Cannot open file '%s' for output. (code: %d)\n", full_fname, code);
    exit_mpi();
  }

  // Write solution snapshot
  Vec5D***  v  = (Vec5D***) V.GetDataPointer();

  int i0, j0, k0, imax, jmax, kmax;
  V.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);

  if(iod.output.density==OutputData::ON) {
    double*** s  = (double***) scalar.GetDataPointer();
    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++)
          s[k][j][i] = v[k][j][i][0];
    scalar.RestoreDataPointerAndInsert();
    PetscObjectSetName((PetscObject)(scalar.GetRefToGlobalVec()), "density");
    VecView(scalar.GetRefToGlobalVec(), viewer);
  }

  if(iod.output.velocity==OutputData::ON) {
    Vec3D*** v3  = (Vec3D***) vector3.GetDataPointer();
    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++)
          for(int p=0; p<3; p++)
            v3[k][j][i][p] = v[k][j][i][1+p];
    vector3.RestoreDataPointerAndInsert();
    PetscObjectSetName((PetscObject)(vector3.GetRefToGlobalVec()), "velocity");
    VecView(vector3.GetRefToGlobalVec(), viewer);
  }

  if(iod.output.pressure==OutputData::ON) {
    double*** s  = (double***) scalar.GetDataPointer();
    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++)
          s[k][j][i] = v[k][j][i][4];
    scalar.RestoreDataPointerAndInsert();
    PetscObjectSetName((PetscObject)(scalar.GetRefToGlobalVec()), "pressure");
    VecView(scalar.GetRefToGlobalVec(), viewer);
  }

  if(iod.output.internal_energy==OutputData::ON) {
    double*** s  = (double***) scalar.GetDataPointer();
    double*** id  = (double***)ID.GetDataPointer();
    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++)
          s[k][j][i] = vf[(int)id[k][j][i]]->GetInternalEnergyPerUnitMass(v[k][j][i][0], v[k][j][i][4]);
    scalar.RestoreDataPointerAndInsert();
    ID.RestoreDataPointerToLocalVector();
    PetscObjectSetName((PetscObject)(scalar.GetRefToGlobalVec()), "internal_energy");
    VecView(scalar.GetRefToGlobalVec(), viewer);
  }

  if(iod.output.materialid==OutputData::ON) {
    
// TODO(KW): not sure why this does not work...
//    PetscObjectSetName((PetscObject)(ID.GetRefToGlobalVec()), "materialid"); //adding the name directly to ID.
//    VecView(ID.GetRefToGlobalVec(), viewer);

    double*** s  = (double***) scalar.GetDataPointer();
    double*** id  = (double***)ID.GetDataPointer();
    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++)
          s[k][j][i] = id[k][j][i];
    scalar.RestoreDataPointerAndInsert();
    ID.RestoreDataPointerToLocalVector();
    PetscObjectSetName((PetscObject)(scalar.GetRefToGlobalVec()), "materialid");
    VecView(scalar.GetRefToGlobalVec(), viewer);

  }

  for(auto it = iod.schemes.ls.dataMap.begin(); it != iod.schemes.ls.dataMap.end(); it++) {
    if(it->first >= OutputData::MAXLS) {
      print_error("*** Error: Not able to output level set %d (id must be less than %d).\n", it->first, OutputData::MAXLS);
      exit_mpi();
    }
    if(iod.output.levelset[it->first]==OutputData::ON) {
      char word[12];
      sprintf(word, "levelset%d", it->first);
      PetscObjectSetName((PetscObject)(Phi[it->first]->GetRefToGlobalVec()), word); //adding the name directly to Phi[i].
      VecView(Phi[it->first]->GetRefToGlobalVec(), viewer);
    }
  }

  // Add a line to the pvd file to record the new solutio snapshot
  char f1[256];
  sprintf(f1, "%s%s.pvd", iod.output.prefix, iod.output.solution_filename_base);
  pvdfile  = fopen(f1,"r+");
  if(!pvdfile) {
    print_error("*** Error: Cannot open file '%s%s.pvd' for output.\n", iod.output.prefix, iod.output.solution_filename_base);
    exit_mpi();
  }
  fseek(pvdfile, -27, SEEK_END); //!< overwrite the previous end of file script
  print(pvdfile, "  <DataSet timestep=\"%e\" file=\"%s\"/>\n", time, fname);
  print(pvdfile, "  </Collection>\n");
  print(pvdfile, "</VTKFile>\n");
  fclose(pvdfile); pvdfile = NULL;

  // clean up
  PetscViewerDestroy(&viewer);
  V.RestoreDataPointerToLocalVector(); //no changes made to V.

  // bookkeeping
  iFrame++;
  last_snapshot_time = time;
   
  //print("\033[0;36m- Wrote solution at %e to %s.\033[0m\n", time, fname);
  print("- Wrote solution at %e to %s.\n", time, fname);
}

//--------------------------------------------------------------------------

void Output::WriteSolutionAtProbes(double time, int time_step, SpaceVariable3D &V, 
                                   SpaceVariable3D &ID, vector<SpaceVariable3D*> &Phi)
{
  probe_output.WriteSolutionAtProbes(time, time_step, V, ID, Phi);
}

//--------------------------------------------------------------------------

void Output::FinalizeOutput()
{
  scalar.Destroy();
  vector3.Destroy(); 
}

//--------------------------------------------------------------------------

















