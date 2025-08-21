/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <SpaceVariable.h>
#include <GlobalMeshInfo.h>
#include <petscviewer.h>
//#include <petscviewerhdf5.h>
#include <Utils.h>
#include <algorithm> //min_element, max_element
#include <cfloat>
using std::vector;
extern int INACTIVE_MATERIAL_ID;

//---------------------------------------------------------
// DataManagers3D
//---------------------------------------------------------
// static member variables
/*
DM DataManagers3D::ghosted1_1dof; //ghosted"1" --> stencil width is 1
DM DataManagers3D::ghosted1_2dof;
DM DataManagers3D::ghosted1_3dof;
DM DataManagers3D::ghosted1_4dof;
DM DataManagers3D::ghosted1_5dof;
DM DataManagers3D::ghosted1_6dof;
DM DataManagers3D::ghosted1_9dof;

DM DataManagers3D::ghosted2_1dof;
DM DataManagers3D::ghosted2_3dof;
*/
//---------------------------------------------------------

DataManagers3D::DataManagers3D()
{

}

//---------------------------------------------------------

DataManagers3D::DataManagers3D(MPI_Comm comm, int NX, int NY, int NZ)
{
  CreateAllDataManagers(comm, NX, NY, NZ);
}

//---------------------------------------------------------

DataManagers3D::~DataManagers3D()
{
 //DMDestroy needs to be called before PetscFinalize()!
}

//---------------------------------------------------------

int DataManagers3D::CreateAllDataManagers(MPI_Comm comm, int NX, int NY, int NZ)
{
  int nProcX, nProcY, nProcZ; //All DM's should use the same domain partition

  auto ierr = DMDACreate3d(comm, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                           DMDA_STENCIL_BOX,
                           NX, NY, NZ,
                           PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                           1/*dof*/, 1/*stencil width*/, 
                           NULL, NULL, NULL,
                           &ghosted1_1dof);
  CHKERRQ(ierr);
  DMSetFromOptions(ghosted1_1dof);
  DMSetUp(ghosted1_1dof);

  DMDAGetInfo(ghosted1_1dof, NULL, NULL, NULL, NULL, &nProcX, &nProcY, &nProcZ, NULL, NULL, NULL,
              NULL, NULL, NULL); 

  ierr = DMDACreate3d(comm, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                      DMDA_STENCIL_BOX,
                      NX, NY, NZ,
                      PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                      2/*dof*/, 1/*stencil width*/, 
                      NULL, NULL, NULL,
                      &ghosted1_2dof);
  CHKERRQ(ierr);
  DMSetFromOptions(ghosted1_2dof);
  DMSetUp(ghosted1_2dof);

  ierr = DMDACreate3d(comm, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                      DMDA_STENCIL_BOX,
                      NX, NY, NZ,
                      PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                      3/*dof*/, 1/*stencil width*/, 
                      NULL, NULL, NULL,
                      &ghosted1_3dof);
  CHKERRQ(ierr);
  DMSetFromOptions(ghosted1_3dof);
  DMSetUp(ghosted1_3dof);

  ierr = DMDACreate3d(comm, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                      DMDA_STENCIL_BOX,
                      NX, NY, NZ,
                      PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                      4/*dof*/, 1/*stencil width*/, 
                      NULL, NULL, NULL,
                      &ghosted1_4dof);
  CHKERRQ(ierr);
  DMSetFromOptions(ghosted1_4dof);
  DMSetUp(ghosted1_4dof);

  ierr = DMDACreate3d(comm, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                      DMDA_STENCIL_BOX,
                      NX, NY, NZ,
                      PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                      5/*dof*/, 1/*stencil width*/, 
                      NULL, NULL, NULL,
                      &ghosted1_5dof);
  CHKERRQ(ierr);
  DMSetFromOptions(ghosted1_5dof);
  DMSetUp(ghosted1_5dof);

  ierr = DMDACreate3d(comm, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                      DMDA_STENCIL_BOX,
                      NX, NY, NZ,
                      PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                      6/*dof*/, 1/*stencil width*/, 
                      NULL, NULL, NULL,
                      &ghosted1_6dof);
  CHKERRQ(ierr);
  DMSetFromOptions(ghosted1_6dof);
  DMSetUp(ghosted1_6dof);

  ierr = DMDACreate3d(comm, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                      DMDA_STENCIL_BOX,
                      NX, NY, NZ,
                      PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                      9/*dof*/, 1/*stencil width*/, 
                      NULL, NULL, NULL,
                      &ghosted1_9dof);
  CHKERRQ(ierr);
  DMSetFromOptions(ghosted1_9dof);
  DMSetUp(ghosted1_9dof);

  ierr = DMDACreate3d(comm, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                      DMDA_STENCIL_BOX,
                      NX, NY, NZ,
                      PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                      1/*dof*/, 2/*stencil width*/, 
                      NULL, NULL, NULL,
                      &ghosted2_1dof);
  CHKERRQ(ierr);
  DMSetFromOptions(ghosted2_1dof);
  DMSetUp(ghosted2_1dof);

  ierr = DMDACreate3d(comm, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                      DMDA_STENCIL_BOX,
                      NX, NY, NZ,
                      PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                      3/*dof*/, 2/*stencil width*/, 
                      NULL, NULL, NULL,
                      &ghosted2_3dof);
  CHKERRQ(ierr);
  DMSetFromOptions(ghosted2_3dof);
  DMSetUp(ghosted2_3dof);

  return 0;
}

//---------------------------------------------------------

void DataManagers3D::DestroyAllDataManagers()
{
  DMDestroy(&ghosted1_1dof);
  DMDestroy(&ghosted1_2dof);
  DMDestroy(&ghosted1_3dof);
  DMDestroy(&ghosted1_4dof);
  DMDestroy(&ghosted1_5dof);
  DMDestroy(&ghosted1_6dof);
  DMDestroy(&ghosted1_9dof);

  DMDestroy(&ghosted2_1dof);
  DMDestroy(&ghosted2_3dof);
}

//---------------------------------------------------------
// SpaceVariable3D
//---------------------------------------------------------

SpaceVariable3D::SpaceVariable3D() : comm(NULL), dm(NULL), globalVec(), localVec()
{
  array = NULL;
}

//---------------------------------------------------------

SpaceVariable3D::SpaceVariable3D(MPI_Comm &comm_, DM *dm_) : globalVec(), localVec()
{
  Setup(comm_, dm_);
}
  
//---------------------------------------------------------

void SpaceVariable3D::Setup(MPI_Comm &comm_, DM *dm_)
{
  comm = &comm_;
  dm = dm_;

  DMCreateGlobalVector(*dm, &globalVec); 
  VecSet(globalVec, 0.0);

  DMCreateLocalVector(*dm, &localVec);
  VecSet(localVec, 0.0);

  array = NULL;

  DMBoundaryType bx, by, bz;

  DMDAGetInfo(*dm, NULL, &NX, &NY, &NZ, &nProcX, &nProcY, &nProcZ, &dof, &ghost_width, 
              &bx, &by, &bz, NULL);
  
  if(bx==DM_BOUNDARY_NONE && by==DM_BOUNDARY_NONE && bz==DM_BOUNDARY_NONE) {
    ghosted = false;
    ghost_width = 0;
  } else if(bx==DM_BOUNDARY_GHOSTED && by==DM_BOUNDARY_GHOSTED && bz==DM_BOUNDARY_GHOSTED) {
    ghosted = true;
  } else {
    PetscPrintf(*comm, "*** Error: Unsupported ghost type.\n");
    MPI_Abort(*comm, 1); 
  }

  DMDAGetCorners(*dm, &i0, &j0, &k0, &nx, &ny, &nz);
  imax = i0 + nx;
  jmax = j0 + ny;
  kmax = k0 + nz;

  DMDAGetGhostCorners(*dm, &ghost_i0, &ghost_j0, &ghost_k0, &ghost_nx, &ghost_ny, &ghost_nz);
  ghost_imax = ghost_i0 + ghost_nx;
  ghost_jmax = ghost_j0 + ghost_ny;
  ghost_kmax = ghost_k0 + ghost_nz;

  internal_ghost_i0   = ghost_i0<0    ?   i0 : ghost_i0;
  internal_ghost_imax = ghost_imax>NX ? imax : ghost_imax;
  internal_ghost_j0   = ghost_j0<0    ?   j0 : ghost_j0;
  internal_ghost_jmax = ghost_jmax>NY ? jmax : ghost_jmax;
  internal_ghost_k0   = ghost_k0<0    ?   k0 : ghost_k0;
  internal_ghost_kmax = ghost_kmax>NZ ? kmax : ghost_kmax;

  numNodes0 = (long long)nx*ny*nz;
  numNodes1 = (long long)(internal_ghost_imax-internal_ghost_i0)
            * (internal_ghost_jmax-internal_ghost_j0)
            * (internal_ghost_kmax-internal_ghost_k0);
  numNodes2 = (long long)ghost_nx*ghost_ny*ghost_nz;
}

//---------------------------------------------------------

SpaceVariable3D::~SpaceVariable3D() 
{
  if(dm) //"Setup(...)" has been called.
    Destroy(); //Actually, the vector should have already been destroyed.
}

//---------------------------------------------------------

double*** SpaceVariable3D::GetDataPointer()
{
  if(!dm) return NULL;

  DMDAVecGetArray(*dm, localVec, &array);
  return array;
}

//---------------------------------------------------------

void SpaceVariable3D::RestoreDataPointerAndInsert()
{
  if(!dm)
    return;

  RestoreDataPointerToLocalVector();
  DMLocalToGlobal(*dm, localVec, INSERT_VALUES, globalVec);

  // sync local to global
  DMGlobalToLocalBegin(*dm, globalVec, INSERT_VALUES, localVec);
  DMGlobalToLocalEnd(*dm, globalVec, INSERT_VALUES, localVec);
}

//---------------------------------------------------------

void SpaceVariable3D::RestoreDataPointerAndAdd()
{
  if(!dm)
    return;

  RestoreDataPointerToLocalVector();
  DMLocalToGlobal(*dm, localVec, ADD_VALUES, globalVec);

  // sync local to global
  DMGlobalToLocalBegin(*dm, globalVec, INSERT_VALUES, localVec);
  DMGlobalToLocalEnd(*dm, globalVec, INSERT_VALUES, localVec);
}

//---------------------------------------------------------

void SpaceVariable3D::RestoreDataPointerToLocalVector()
{
  if(!dm)
    return;

  DMDAVecRestoreArray(*dm, localVec, &array);  
}

//---------------------------------------------------------

void SpaceVariable3D::SyncLocalToGlobal()
{
  if(!dm)
    return;

  // sync local to global
  DMGlobalToLocalBegin(*dm, globalVec, INSERT_VALUES, localVec);
  DMGlobalToLocalEnd(*dm, globalVec, INSERT_VALUES, localVec);
}

//---------------------------------------------------------

void SpaceVariable3D::Destroy()
{
  if(!dm)
    return;

  VecDestroy(&globalVec);
  VecDestroy(&localVec);
}

//---------------------------------------------------------

void SpaceVariable3D::StoreMeshCoordinates(SpaceVariable3D &coordinates)
{
  if(!dm)
    return;

  DMSetCoordinateDim(*dm, 3/*3D*/);
  DMSetCoordinates(*dm, coordinates.globalVec);
}

//---------------------------------------------------------

void SpaceVariable3D::WriteToVTRFile(const char *filename, const char *varname)
{
  if(!dm)
    return;

  if(varname)
    SetOutputVariableName(varname);

  PetscViewer viewer(nullptr);
  PetscViewerVTKOpen(PetscObjectComm((PetscObject)*dm), filename, FILE_MODE_WRITE, &viewer);
  VecView(globalVec, viewer);
  PetscViewerDestroy(&viewer);
  MPI_Barrier(*comm); 
}

//---------------------------------------------------------

void SpaceVariable3D::WriteToHDF5File([[maybe_unused]] const char *filename, const char *varname)
{
  if(!dm)
    return;

  if(varname)
    SetOutputVariableName(varname);

  PetscViewer viewer(nullptr);
  //PetscViewerHDF5Open(PetscObjectComm((PetscObject)*dm), filename, FILE_MODE_WRITE, &viewer); //NOT AVAILABLE?
  VecView(globalVec, viewer);
  PetscViewerDestroy(&viewer);
  MPI_Barrier(*comm); 
}

//---------------------------------------------------------

void SpaceVariable3D::WriteToCGNSFile(const char *filename, const char *varname)
{
  if(!dm)
    return;

  if(varname)
    SetOutputVariableName(varname);

  PetscViewer viewer(nullptr);
  PetscViewerCreate(*comm, &viewer);
  //PetscViewerSetType(viewer, PETSCVIEWERCGNS); //NOT AVAILABLE?
  PetscViewerFileSetMode(viewer, FILE_MODE_WRITE);
  PetscViewerFileSetName(viewer, filename);

  VecView(globalVec, viewer);

  PetscViewerDestroy(&viewer);
  MPI_Barrier(*comm);
}

//---------------------------------------------------------

void SpaceVariable3D::WriteToMatlabFile(const char *filename, const char *varname)
{
  if(!dm)
    return;

  if(varname)
    SetOutputVariableName(varname);

  PetscViewer viewer;
  int code = PetscViewerASCIIOpen(*comm, filename, &viewer);
  if(code) {
    print_error("*** Error: Cannot open file '%s' for output. (code: %d)\n", filename, code);
    exit_mpi();
  }

  PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
  VecView(globalVec, viewer);

  PetscViewerDestroy(&viewer);
  MPI_Barrier(*comm); 
}

//---------------------------------------------------------

void SpaceVariable3D::SetOutputVariableName(const char *name)
{
  if(!dm)
    return;

  PetscObjectSetName((PetscObject)(globalVec), name);
}

//---------------------------------------------------------

void SpaceVariable3D::AXPlusB(double a, double b, bool workOnGhost)
{
  if(!dm)
    return;

  double*** v = GetDataPointer();
  int myi0, myj0, myk0, myimax, myjmax, mykmax;

  if(workOnGhost)
    GetGhostedCornerIndices(&myi0, &myj0, &myk0, &myimax, &myjmax, &mykmax);
  else
    GetCornerIndices(&myi0, &myj0, &myk0, &myimax, &myjmax, &mykmax);

  for(int k=myk0; k<mykmax; k++)
    for(int j=myj0; j<myjmax; j++)
      for(int i=myi0; i<myimax; i++)
        for(int p=0; p<dof; p++)
          v[k][j][i*dof+p] = a*v[k][j][i*dof+p] + b;

  RestoreDataPointerAndInsert();
}

//---------------------------------------------------------

void SpaceVariable3D::AXPlusBY(double a, double b, SpaceVariable3D &y, bool workOnGhost)
{
  if(dof != y.NumDOF()) {
    print_error("*** Error: Vector operation (AXPlusBY) failed due to inconsistent sizes (%d vs. %d)\n", dof, y.NumDOF());
    exit_mpi();
  }

  double*** v  = GetDataPointer();
  double*** v2 = y.GetDataPointer();

  int myi0, myj0, myk0, myimax, myjmax, mykmax;

  if(workOnGhost)
    GetGhostedCornerIndices(&myi0, &myj0, &myk0, &myimax, &myjmax, &mykmax);
  else
    GetCornerIndices(&myi0, &myj0, &myk0, &myimax, &myjmax, &mykmax);

  for(int k=myk0; k<mykmax; k++)
    for(int j=myj0; j<myjmax; j++)
      for(int i=myi0; i<myimax; i++)
        for(int p=0; p<dof; p++)
          v[k][j][i*dof+p] = a*v[k][j][i*dof+p] + b*v2[k][j][i*dof+p];

  RestoreDataPointerAndInsert();
  y.RestoreDataPointerToLocalVector(); //no changes
}

//---------------------------------------------------------
// takes in two vectors that identify the indices of X and Y that are involved in this operation
void SpaceVariable3D::AXPlusBY(double a, double b, SpaceVariable3D &y, std::vector<int>& Xindices,
                               std::vector<int>& Yindices, bool workOnGhost)
{
  //check for common errors
  if(Xindices.size() != Yindices.size()) {
    print_error("*** Error: Vector operation (AXPlusBY) failed due to inconsistent number of indices (%d vs. %d)\n", 
                Xindices.size(), Yindices.size());
    exit_mpi();
  }
  if(*std::min_element(Xindices.begin(), Xindices.end())<0 ||  
     *std::max_element(Xindices.begin(), Xindices.end())>=dof) {
    print_error("*** Error: Vector operation (AXPlusBY) failed due to incorrect X index (min:%d, max:%d).\n", 
                *std::min_element(Xindices.begin(), Xindices.end()),
                *std::max_element(Xindices.begin(), Xindices.end()));
    exit_mpi();
  }
  if(*std::min_element(Yindices.begin(), Yindices.end())<0 ||  
     *std::max_element(Yindices.begin(), Yindices.end())>=y.NumDOF()) {
    print_error("*** Error: Vector operation (AXPlusBY) failed due to incorrect Y index (min:%d, max:%d).\n", 
                *std::min_element(Yindices.begin(), Yindices.end()),
                *std::max_element(Yindices.begin(), Yindices.end()));
    exit_mpi();
  }

  double*** v  = GetDataPointer();
  double*** v2 = y.GetDataPointer();

  int myi0, myj0, myk0, myimax, myjmax, mykmax;

  if(workOnGhost)
    GetGhostedCornerIndices(&myi0, &myj0, &myk0, &myimax, &myjmax, &mykmax);
  else
    GetCornerIndices(&myi0, &myj0, &myk0, &myimax, &myjmax, &mykmax);

  int px, py;
  for(int k=myk0; k<mykmax; k++)
    for(int j=myj0; j<myjmax; j++)
      for(int i=myi0; i<myimax; i++)
        for(int id=0; id<(int)Xindices.size(); id++) {
          px = Xindices[id];
          py = Yindices[id];
          v[k][j][i*dof+px] = a*v[k][j][i*dof+px] + b*v2[k][j][i*dof+py];
        }

  RestoreDataPointerAndInsert();
  y.RestoreDataPointerToLocalVector(); //no changes
}

//---------------------------------------------------------

void SpaceVariable3D::SetConstantValue(double a, bool workOnGhost)
{
  double*** v  = GetDataPointer();

  int myi0, myj0, myk0, myimax, myjmax, mykmax;

  if(workOnGhost)
    GetGhostedCornerIndices(&myi0, &myj0, &myk0, &myimax, &myjmax, &mykmax);
  else
    GetCornerIndices(&myi0, &myj0, &myk0, &myimax, &myjmax, &mykmax);

  for(int k=myk0; k<mykmax; k++)
    for(int j=myj0; j<myjmax; j++)
      for(int i=myi0; i<myimax; i++)
        for(int p=0; p<dof; p++)
          v[k][j][i*dof+p] = a;

  RestoreDataPointerAndInsert();
}

//---------------------------------------------------------

double
SpaceVariable3D::CalculateGlobalMin(int mydof, bool workOnGhost)
{
  assert(mydof<dof && mydof>=0);

  double global_min = DBL_MAX;

  double*** v  = GetDataPointer();

  int myi0, myj0, myk0, myimax, myjmax, mykmax;

  if(workOnGhost)
    GetGhostedCornerIndices(&myi0, &myj0, &myk0, &myimax, &myjmax, &mykmax);
  else
    GetCornerIndices(&myi0, &myj0, &myk0, &myimax, &myjmax, &mykmax);

  for(int k=myk0; k<mykmax; k++)
    for(int j=myj0; j<myjmax; j++)
      for(int i=myi0; i<myimax; i++)
        global_min = std::min(global_min, v[k][j][i*dof+mydof]); 

  MPI_Allreduce(MPI_IN_PLACE, &global_min, 1, MPI_DOUBLE, MPI_MIN, *comm);

  RestoreDataPointerToLocalVector(); 

  return global_min;  
}

//---------------------------------------------------------

double
SpaceVariable3D::CalculateGlobalMax(int mydof, bool workOnGhost)
{
  assert(mydof<dof && mydof>=0);

  double global_max = -DBL_MAX;

  double*** v  = GetDataPointer();

  int myi0, myj0, myk0, myimax, myjmax, mykmax;

  if(workOnGhost)
    GetGhostedCornerIndices(&myi0, &myj0, &myk0, &myimax, &myjmax, &mykmax);
  else
    GetCornerIndices(&myi0, &myj0, &myk0, &myimax, &myjmax, &mykmax);

  for(int k=myk0; k<mykmax; k++)
    for(int j=myj0; j<myjmax; j++)
      for(int i=myi0; i<myimax; i++)
        global_max = std::max(global_max, v[k][j][i*dof+mydof]); 

  MPI_Allreduce(MPI_IN_PLACE, &global_max, 1, MPI_DOUBLE, MPI_MAX, *comm);

  RestoreDataPointerToLocalVector(); 

  return global_max;  
}

//---------------------------------------------------------

double
SpaceVariable3D::CalculateVectorOneNorm()
{
  double norm(0.0);
  VecNorm(globalVec, NORM_1, &norm);
  return norm;
}

//---------------------------------------------------------

double
SpaceVariable3D::CalculateVectorTwoNorm()
{
  double norm(0.0);
  VecNorm(globalVec, NORM_2, &norm);
  return norm;
}

//---------------------------------------------------------

double
SpaceVariable3D::CalculateVectorInfNorm()
{
  double norm(0.0);
  VecNorm(globalVec, NORM_INFINITY, &norm);
  return norm;
}

//---------------------------------------------------------

void
SpaceVariable3D::CalculateFunctionNormsConRec(SpaceVariable3D &volume, vector<double> &norm1_dofs,
                                              vector<double> &norm2_dofs, vector<double> &norminf_dofs)
{
  norm1_dofs.assign(dof,0.0);
  norm2_dofs.assign(dof,0.0);
  norminf_dofs.assign(dof,0.0);

  double*** v   = GetDataPointer();
  double*** vol = volume.GetDataPointer();

  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++)
        for(int p=0; p<dof; p++) {
          norm1_dofs[p] += fabs(v[k][j][i*dof+p])*vol[k][j][i];
          norm2_dofs[p] += v[k][j][i*dof+p]*v[k][j][i*dof+p]*vol[k][j][i];
          norminf_dofs[p] = std::max(norminf_dofs[p], fabs(v[k][j][i*dof+p]));
        }

  MPI_Allreduce(MPI_IN_PLACE, norm1_dofs.data(), dof, MPI_DOUBLE, MPI_SUM, *comm);
  MPI_Allreduce(MPI_IN_PLACE, norm2_dofs.data(), dof, MPI_DOUBLE, MPI_SUM, *comm);
  MPI_Allreduce(MPI_IN_PLACE, norminf_dofs.data(), dof, MPI_DOUBLE, MPI_MAX, *comm);

  for(auto&& norm2 : norm2_dofs)
    norm2 = sqrt(norm2);
 
  RestoreDataPointerToLocalVector();
  volume.RestoreDataPointerToLocalVector();
}

//---------------------------------------------------------

void
SpaceVariable3D::CalculateFunctionNormsConRec(SpaceVariable3D &ID, SpaceVariable3D &volume, vector<double> &norm1_dofs,
                                              vector<double> &norm2_dofs, vector<double> &norminf_dofs)
{
  norm1_dofs.assign(dof,0.0);
  norm2_dofs.assign(dof,0.0);
  norminf_dofs.assign(dof,0.0);

  double*** v   = GetDataPointer();
  double*** vol = volume.GetDataPointer();
  double*** id  = ID.GetDataPointer();

  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        if(id[k][j][i] == INACTIVE_MATERIAL_ID)
          continue;
        for(int p=0; p<dof; p++) {
          norm1_dofs[p] += fabs(v[k][j][i*dof+p])*vol[k][j][i];
          norm2_dofs[p] += v[k][j][i*dof+p]*v[k][j][i*dof+p]*vol[k][j][i];
          norminf_dofs[p] = std::max(norminf_dofs[p], fabs(v[k][j][i*dof+p]));
        }
      }

  MPI_Allreduce(MPI_IN_PLACE, norm1_dofs.data(), dof, MPI_DOUBLE, MPI_SUM, *comm);
  MPI_Allreduce(MPI_IN_PLACE, norm2_dofs.data(), dof, MPI_DOUBLE, MPI_SUM, *comm);
  MPI_Allreduce(MPI_IN_PLACE, norminf_dofs.data(), dof, MPI_DOUBLE, MPI_MAX, *comm);

  for(auto&& norm2 : norm2_dofs)
    norm2 = sqrt(norm2);
 
  RestoreDataPointerToLocalVector();
  volume.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();
}

//---------------------------------------------------------

void
SpaceVariable3D::CalculateFunctionNormsConRec(GlobalMeshInfo &global_mesh, vector<double> &norm1_dofs,
                                              vector<double> &norm2_dofs, vector<double> &norminf_dofs)
{
  norm1_dofs.assign(dof,0.0);
  norm2_dofs.assign(dof,0.0);
  norminf_dofs.assign(dof,0.0);

  double*** v   = GetDataPointer();

  double dz, dydz, dxdydz;
  for(int k=k0; k<kmax; k++) {
    dz = global_mesh.GetDz(k);
    for(int j=j0; j<jmax; j++) {
      dydz = dz*global_mesh.GetDy(j);
      for(int i=i0; i<imax; i++) {
        dxdydz = dydz*global_mesh.GetDx(i);
        for(int p=0; p<dof; p++) {
          norm1_dofs[p] += fabs(v[k][j][i*dof+p])*dxdydz;
          norm2_dofs[p] += v[k][j][i*dof+p]*v[k][j][i*dof+p]*dxdydz;
          norminf_dofs[p] = std::max(norminf_dofs[p], fabs(v[k][j][i*dof+p]));
        }
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, norm1_dofs.data(), dof, MPI_DOUBLE, MPI_SUM, *comm);
  MPI_Allreduce(MPI_IN_PLACE, norm2_dofs.data(), dof, MPI_DOUBLE, MPI_SUM, *comm);
  MPI_Allreduce(MPI_IN_PLACE, norminf_dofs.data(), dof, MPI_DOUBLE, MPI_MAX, *comm);

  for(auto&& norm2 : norm2_dofs)
    norm2 = sqrt(norm2);
 
  RestoreDataPointerToLocalVector();
}

//---------------------------------------------------------

void
SpaceVariable3D::CalculateFunctionNormsConRec(SpaceVariable3D &ID,
                                              GlobalMeshInfo &global_mesh, vector<double> &norm1_dofs,
                                              vector<double> &norm2_dofs, vector<double> &norminf_dofs)
{
  norm1_dofs.assign(dof,0.0);
  norm2_dofs.assign(dof,0.0);
  norminf_dofs.assign(dof,0.0);

  double*** v  = GetDataPointer();
  double*** id = ID.GetDataPointer();

  double dz, dydz, dxdydz;
  for(int k=k0; k<kmax; k++) {
    dz = global_mesh.GetDz(k);
    for(int j=j0; j<jmax; j++) {
      dydz = dz*global_mesh.GetDy(j);
      for(int i=i0; i<imax; i++) {
        if(id[k][j][i] == INACTIVE_MATERIAL_ID)
          continue;
        dxdydz = dydz*global_mesh.GetDx(i);
        for(int p=0; p<dof; p++) {
          norm1_dofs[p] += fabs(v[k][j][i*dof+p])*dxdydz;
          norm2_dofs[p] += v[k][j][i*dof+p]*v[k][j][i*dof+p]*dxdydz;
          norminf_dofs[p] = std::max(norminf_dofs[p], fabs(v[k][j][i*dof+p]));
        }
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, norm1_dofs.data(), dof, MPI_DOUBLE, MPI_SUM, *comm);
  MPI_Allreduce(MPI_IN_PLACE, norm2_dofs.data(), dof, MPI_DOUBLE, MPI_SUM, *comm);
  MPI_Allreduce(MPI_IN_PLACE, norminf_dofs.data(), dof, MPI_DOUBLE, MPI_MAX, *comm);

  for(auto&& norm2 : norm2_dofs)
    norm2 = sqrt(norm2);
 
  RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();
}

//---------------------------------------------------------


//---------------------------------------------------------



//---------------------------------------------------------



