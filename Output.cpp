/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <Utils.h>
#include <Vector5D.h>
#include <Output.h>
#include <float.h> //DBL_MAX

//--------------------------------------------------------------------------

Output::Output(MPI_Comm &comm_, DataManagers3D &dms, IoData &iod_, GlobalMeshInfo &global_mesh_,
               vector<VarFcnBase*> &vf_, 
               SpaceVariable3D &cell_volume, IonizationOperator* ion_) : 
    comm(comm_), 
    iod(iod_), global_mesh(global_mesh_), vf(vf_),
    scalar(comm_, &(dms.ghosted1_1dof)),
    vector3(comm_, &(dms.ghosted1_3dof)),
    probe_output(comm_, iod_.output, vf_, ion_),
    matvol_output(comm_, iod_, cell_volume),
    ion(ion_),
    terminal(comm_, iod_.terminal_visualization, global_mesh_, vf_, ion_)
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

  // setup line plots
  int numLines = iod.output.linePlots.dataMap.size();
  line_outputs.resize(numLines, NULL);
  for(auto it = iod.output.linePlots.dataMap.begin(); it != iod.output.linePlots.dataMap.end(); it++) {
    int line_number = it->first;
    if(line_number<0 || line_number>=numLines) {
      print_error("*** Error: Detected error in line output. Line number = %d (should be between 0 and %d)\n",
                  line_number, numLines-1); 
      exit_mpi();
    }
    line_outputs[line_number] = new ProbeOutput(comm, iod.output, vf, ion, line_number);
  }

  // setup plane plots
  int numPlanes = iod.output.planePlots.dataMap.size();
  plane_outputs.resize(numPlanes, NULL);
  for(auto&& plane : iod.output.planePlots.dataMap) {
    int plane_number = plane.first;
    if(plane_number<0 || plane_number>=numPlanes) {
      print_error("*** Error: Detected error in plane output. Plane number = %d (should be between 0 and %d)\n",
                  plane_number, numPlanes-1); 
      exit_mpi();
    }
    plane_outputs[plane_number] = new PlaneOutput(comm, iod.output, *plane.second, vf, global_mesh, ion);
  }

  // check ionization requests
  if(iod.output.ionization_output_requested() && !ion) {
    print_error("*** Error: User requested ionization output(s) without specifying an ionization model.\n");
    exit_mpi();
  }

}

//--------------------------------------------------------------------------

Output::~Output()
{
  if(pvdfile) fclose(pvdfile);
  for(int i=0; i<(int)line_outputs.size(); i++)
    if(line_outputs[i]) delete line_outputs[i];
  for(int i=0; i<(int)plane_outputs.size(); i++)
    if(plane_outputs[i]) delete plane_outputs[i];
}

//--------------------------------------------------------------------------

void Output::InitializeOutput(SpaceVariable3D &coordinates)
{
  scalar.StoreMeshCoordinates(coordinates);
  vector3.StoreMeshCoordinates(coordinates);

  probe_output.SetupInterpolation(coordinates);

  for(int i=0; i<(int)line_outputs.size(); i++)
    line_outputs[i]->SetupInterpolation(coordinates);

  for(auto&& plane : plane_outputs)
    plane->InitializeOutput(coordinates);

  if(iod.output.mesh_filename[0] != 0)
    OutputMeshInformation(coordinates);

  if(iod.output.mesh_partition[0] != 0)
    OutputMeshPartition();
}

//--------------------------------------------------------------------------

void Output::OutputSolutions(double time, double dt, int time_step, SpaceVariable3D &V, 
                             SpaceVariable3D &ID, std::vector<SpaceVariable3D*> &Phi, 
                             SpaceVariable3D *L, SpaceVariable3D *Xi, bool force_write)
{

  //write solution snapshot
  if(isTimeToWrite(time, dt, time_step, iod.output.frequency_dt, iod.output.frequency, 
     last_snapshot_time, force_write))
    WriteSolutionSnapshot(time, time_step, V, ID, Phi, L, Xi);

  //write solutions at probes
  probe_output.WriteSolutionAtProbes(time, dt, time_step, V, ID, Phi, L, force_write);

  //write solutions along lines
  for(int i=0; i<(int)line_outputs.size(); i++)
    line_outputs[i]->WriteAllSolutionsAlongLine(time, dt, time_step, V, ID, Phi, L, force_write);

  //write solutions on planes
  for(auto&& plane : plane_outputs)
    plane->WriteSolutionOnPlane(time, dt, time_step, V, ID, Phi, L, force_write);

  //write material volumes
  matvol_output.WriteSolution(time, dt, time_step, ID, force_write);

  //write terminal visualization 
  terminal.PrintSolutionSnapshot(time, dt, time_step, V, ID, Phi, L, force_write);

}
//--------------------------------------------------------------------------

void Output::WriteSolutionSnapshot(double time, [[maybe_unused]] int time_step, SpaceVariable3D &V, 
                                   SpaceVariable3D &ID, std::vector<SpaceVariable3D*> &Phi,
                                   SpaceVariable3D *L, SpaceVariable3D *Xi)
{
  //! Post-processing
  if(iod.output.ionization_output_requested())
    ion->ComputeIonization(V,ID);

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
  double*** id = (double***)ID.GetDataPointer();

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


  //KW: Looks like "VecView" may have a problem depending on the order. I had to move "molar_fractions" up here.
  bool molar_fractions_requested = false;
  for(int j=0; j<OutputData::MAXSPECIES; j++)
    if(iod.output.molar_fractions[j]==OutputData::ON) {
      molar_fractions_requested = true;
      break;
    }

  if(molar_fractions_requested) {
    std::map<int, SpaceVariable3D*>& AlphaRJ(ion->GetReferenceToAlphaRJ());
    for(int j=0; j<OutputData::MAXSPECIES; j++) {
      if(iod.output.molar_fractions[j] != OutputData::ON)
        continue;

      auto it = AlphaRJ.find(j);
      if(it == AlphaRJ.end()) {
        print_error("*** Error: User requested output of molar fractions for species %d, but it does not exist.\n", j);
        exit_mpi();
      }

      SpaceVariable3D* AlphaR = it->second;
      char word[40];
      sprintf(word, "molar_fractions_%d", j); 
      PetscObjectSetName((PetscObject)(AlphaR->GetRefToGlobalVec()), word);
      VecView(AlphaR->GetRefToGlobalVec(), viewer);
    } 
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
    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++)
          s[k][j][i] = vf[(int)id[k][j][i]]->GetInternalEnergyPerUnitMass(v[k][j][i][0], v[k][j][i][4]);
    scalar.RestoreDataPointerAndInsert();
    PetscObjectSetName((PetscObject)(scalar.GetRefToGlobalVec()), "internal_energy");
    VecView(scalar.GetRefToGlobalVec(), viewer);
  }

  if(iod.output.delta_internal_energy==OutputData::ON) {
    double*** s  = (double***) scalar.GetDataPointer();
    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++)
          s[k][j][i] = vf[(int)id[k][j][i]]->GetInternalEnergyPerUnitMass(v[k][j][i][0], v[k][j][i][4])
                     - vf[(int)id[k][j][i]]->GetReferenceInternalEnergyPerUnitMass();
    scalar.RestoreDataPointerAndInsert();
    PetscObjectSetName((PetscObject)(scalar.GetRefToGlobalVec()), "delta_internal_energy");
    VecView(scalar.GetRefToGlobalVec(), viewer);
  }

  if(iod.output.materialid==OutputData::ON) {
    
// TODO(KW): not sure why this does not work...
//    PetscObjectSetName((PetscObject)(ID.GetRefToGlobalVec()), "materialid"); //adding the name directly to ID.
//    VecView(ID.GetRefToGlobalVec(), viewer);

    double*** s  = (double***) scalar.GetDataPointer();
    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++)
          s[k][j][i] = id[k][j][i];
    scalar.RestoreDataPointerAndInsert();
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


  if(iod.output.temperature==OutputData::ON) {
    double*** s  = (double***) scalar.GetDataPointer();
    double e;
    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++) {
          e = vf[(int)id[k][j][i]]->GetInternalEnergyPerUnitMass(v[k][j][i][0], v[k][j][i][4]);
          s[k][j][i] = vf[(int)id[k][j][i]]->GetTemperature(v[k][j][i][0], e);
        }
    scalar.RestoreDataPointerAndInsert();
    PetscObjectSetName((PetscObject)(scalar.GetRefToGlobalVec()), "temperature");
    VecView(scalar.GetRefToGlobalVec(), viewer);
  }


  if(iod.output.delta_temperature==OutputData::ON) {
    double*** s  = (double***) scalar.GetDataPointer();
    double e;
    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++) {
          e = vf[(int)id[k][j][i]]->GetInternalEnergyPerUnitMass(v[k][j][i][0], v[k][j][i][4]);
          s[k][j][i] = vf[(int)id[k][j][i]]->GetTemperature(v[k][j][i][0], e)
                     - vf[(int)id[k][j][i]]->GetReferenceTemperature();
        }
    scalar.RestoreDataPointerAndInsert();
    PetscObjectSetName((PetscObject)(scalar.GetRefToGlobalVec()), "delta_temperature");
    VecView(scalar.GetRefToGlobalVec(), viewer);
  }


  if(iod.output.laser_radiance==OutputData::ON) {
    if(L == NULL) {
      print_error("*** Error: Cannot output laser radiance. Solver is not activated.\n");
      exit_mpi();
    }
    PetscObjectSetName((PetscObject)(L->GetRefToGlobalVec()), "laser_radiance");
    VecView(L->GetRefToGlobalVec(), viewer);
  }


  if(iod.output.mean_charge==OutputData::ON) {
    SpaceVariable3D& Zav(ion->GetReferenceToZav());
    PetscObjectSetName((PetscObject)(Zav.GetRefToGlobalVec()), "mean_charge_number");
    VecView(Zav.GetRefToGlobalVec(), viewer);
  }

  if(iod.output.heavy_particles_density==OutputData::ON) {
    SpaceVariable3D& Nh(ion->GetReferenceToNh());
    PetscObjectSetName((PetscObject)(Nh.GetRefToGlobalVec()), "heavy_particles_density");
    VecView(Nh.GetRefToGlobalVec(), viewer);
  }

  if(iod.output.electron_density==OutputData::ON) {
    SpaceVariable3D& Ne(ion->GetReferenceToNe());
    PetscObjectSetName((PetscObject)(Ne.GetRefToGlobalVec()), "electron_density");
    VecView(Ne.GetRefToGlobalVec(), viewer);
  }

  if(iod.output.electron_density==OutputData::ON) {
    SpaceVariable3D& Ne(ion->GetReferenceToNe());
    PetscObjectSetName((PetscObject)(Ne.GetRefToGlobalVec()), "electron_density");
    VecView(Ne.GetRefToGlobalVec(), viewer);
  }

  if(iod.output.reference_map==OutputData::ON) {
    if(Xi == NULL) {
      print_error("*** Error: Cannot output reference map. Solver is not activated.\n");
      exit_mpi();
    }
    PetscObjectSetName((PetscObject)(Xi->GetRefToGlobalVec()), "reference_map");
    VecView(Xi->GetRefToGlobalVec(), viewer);
  }


  MPI_Barrier(comm); //this might be needed to avoid file corruption (incomplete output)

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
  ID.RestoreDataPointerToLocalVector();

  // bookkeeping
  iFrame++;
  last_snapshot_time = time;
   
  //print("\033[0;36m- Wrote solution at %e to %s.\033[0m\n", time, fname);
  print("- Wrote solution at %e to %s.\n", time, fname);
}

//--------------------------------------------------------------------------

void Output::OutputMeshInformation(SpaceVariable3D& coordinates)
{
  if(iod.output.mesh_filename[0] == 0)
    return; //nothing to do

  char fname[256];
  sprintf(fname, "%s%s", iod.output.prefix, iod.output.mesh_filename);
  FILE* file = fopen(fname, "w");

  int ii0, jj0, kk0, iimax, jjmax, kkmax;
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);

  int NX, NY, NZ; 
  coordinates.GetGlobalSize(&NX, &NY, &NZ);

  int nGhost = coordinates.NumGhostLayers();

  print(file, "## Number of Cells/Nodes (Excluding Ghost Layer(s)): NX = %d, NY = %d, NZ = %d.\n",
        NX, NY, NZ); 
  print(file, "## Number of Ghost Layers: %d\n", nGhost);
  print(file, "## Index  |  x  |  y  |  z\n");

  vector<double> x,y,z;
  x.resize(NX+2*nGhost, -DBL_MAX);
  y.resize(NY+2*nGhost, -DBL_MAX);
  z.resize(NZ+2*nGhost, -DBL_MAX);

  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();

  for(int i=ii0; i<iimax; i++) {
    if((i==ii0 && i>=0) || (i==iimax-1 && i<NX)) //internal ghost layer
      continue;
    x[i+nGhost] = coords[kk0][jj0][i][0];
  }
  for(int j=jj0; j<jjmax; j++) {
    if((j==jj0 && j>=0) || (j==jjmax-1 && j<NY)) //internal ghost layer
      continue;
    y[j+nGhost] = coords[kk0][j][ii0][1];
  }
  for(int k=kk0; k<kkmax; k++) {
    if((k==kk0 && k>=0) || (k==kkmax-1 && k<NZ)) //internal ghost layer
      continue;
    z[k+nGhost] = coords[k][jj0][ii0][2];
  }

  // Collect data
  MPI_Allreduce(MPI_IN_PLACE, x.data(), x.size(), MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(MPI_IN_PLACE, y.data(), y.size(), MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(MPI_IN_PLACE, z.data(), z.size(), MPI_DOUBLE, MPI_MAX, comm);

  for(int i=0; i<std::max(std::max((int)x.size(),(int)y.size()),(int)z.size()); i++) {
    print(file,"%8d\t", i-nGhost);
    if(i<(int)x.size())
      print(file,"%16.8e\t", x[i]);
    else
      print(file,"                \t");
    if(i<(int)y.size())
      print(file,"%16.8e\t", y[i]);
    else
      print(file,"                \t");
    if(i<(int)z.size())
      print(file,"%16.8e", z[i]);
    else
      print(file,"                ");
    print(file,"\n");
  }

  coordinates.RestoreDataPointerToLocalVector();
  fclose(file);
}

//--------------------------------------------------------------------------

void Output::OutputMeshPartition()
{
  if(iod.output.mesh_partition[0] == 0)
    return; //nothing to do

  char fname[256];
  sprintf(fname, "%s%s.vtr", iod.output.prefix, iod.output.mesh_partition);

  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);

  int i0, j0, k0, imax, jmax, kmax;
  scalar.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);

  double*** s  = (double***) scalar.GetDataPointer();
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++)
        s[k][j][i] = mpi_rank;
  scalar.RestoreDataPointerAndInsert();

  scalar.WriteToVTRFile(fname, "partition");
  
  MPI_Barrier(comm); //this might be needed to avoid file corruption (incomplete output)

  scalar.SetConstantValue(0.0); //clean up the internal variable (not necessary)
}

//--------------------------------------------------------------------------

void Output::FinalizeOutput()
{
  scalar.Destroy();
  vector3.Destroy(); 
}

//--------------------------------------------------------------------------

















