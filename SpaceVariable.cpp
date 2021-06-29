#include <SpaceVariable.h>
#include <Utils.h>
#include <bits/stdc++.h> //min_element, max_element

//---------------------------------------------------------
// DataManagers3D
//---------------------------------------------------------
// static member variables
DM DataManagers3D::ghosted1_1dof; //ghosted"1" --> stencil width is 1
DM DataManagers3D::ghosted1_2dof;
DM DataManagers3D::ghosted1_3dof;
DM DataManagers3D::ghosted1_5dof;

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
                      5/*dof*/, 1/*stencil width*/, 
                      NULL, NULL, NULL,
                      &ghosted1_5dof);
  CHKERRQ(ierr);
  DMSetFromOptions(ghosted1_5dof);
  DMSetUp(ghosted1_5dof);

  return 0;
}

//---------------------------------------------------------

void DataManagers3D::DestroyAllDataManagers()
{
  DMDestroy(&ghosted1_1dof);
  DMDestroy(&ghosted1_2dof);
  DMDestroy(&ghosted1_3dof);
  DMDestroy(&ghosted1_5dof);
}

//---------------------------------------------------------
// SpaceVariable3D
//---------------------------------------------------------

SpaceVariable3D::SpaceVariable3D(MPI_Comm &comm_, DM *dm_) : comm(comm_), globalVec(), localVec()
{
  dm = dm_;

  auto ierr = DMCreateGlobalVector(*dm, &globalVec); 
  //CHKERRQ(ierr);
  VecSet(globalVec, 0.0);

  ierr = DMCreateLocalVector(*dm, &localVec);
  //CHKERRQ(ierr);
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
    PetscPrintf(comm, "*** Error: Unsupported ghost type.\n");
    MPI_Abort(comm, 1); 
  }

  DMDAGetCorners(*dm, &i0, &j0, &k0, &nx, &ny, &nz);
  imax = i0 + nx;
  jmax = j0 + ny;
  kmax = k0 + nz;

  DMDAGetGhostCorners(*dm, &ghost_i0, &ghost_j0, &ghost_k0, &ghost_nx, &ghost_ny, &ghost_nz);
  ghost_imax = ghost_i0 + ghost_nx;
  ghost_jmax = ghost_j0 + ghost_ny;
  ghost_kmax = ghost_k0 + ghost_nz;

}

//---------------------------------------------------------

SpaceVariable3D::~SpaceVariable3D() 
{
  Destroy();
}

//---------------------------------------------------------

double*** SpaceVariable3D::GetDataPointer()
{
  auto ierr = DMDAVecGetArray(*dm, localVec, &array);
  //CHKERRQ(ierr);
  return array;
}

//---------------------------------------------------------

void SpaceVariable3D::RestoreDataPointerAndInsert()
{
  RestoreDataPointerToLocalVector();
  auto ierr = DMLocalToGlobal(*dm, localVec, INSERT_VALUES, globalVec);
  //CHKERRQ(ierr);

  // sync local to global
  DMGlobalToLocalBegin(*dm, globalVec, INSERT_VALUES, localVec);
  DMGlobalToLocalEnd(*dm, globalVec, INSERT_VALUES, localVec);
}

//---------------------------------------------------------

void SpaceVariable3D::RestoreDataPointerAndAdd()
{
  RestoreDataPointerToLocalVector();
  auto ierr = DMLocalToGlobal(*dm, localVec, ADD_VALUES, globalVec);
  //CHKERRQ(ierr);

  // sync local to global
  DMGlobalToLocalBegin(*dm, globalVec, INSERT_VALUES, localVec);
  DMGlobalToLocalEnd(*dm, globalVec, INSERT_VALUES, localVec);
}

//---------------------------------------------------------

void SpaceVariable3D::RestoreDataPointerToLocalVector()
{
  auto ierr = DMDAVecRestoreArray(*dm, localVec, &array);  
  //CHKERRQ(ierr);
}

//---------------------------------------------------------

void SpaceVariable3D::Destroy()
{
  VecDestroy(&globalVec);
  VecDestroy(&localVec);
}

//---------------------------------------------------------

void SpaceVariable3D::StoreMeshCoordinates(SpaceVariable3D &coordinates)
{
  auto ierr = DMSetCoordinateDim(*dm, 3/*3D*/);
  //CHKERRQ(ierr);
  ierr = DMSetCoordinates(*dm, coordinates.globalVec);
  //CHKERRQ(ierr);
}

//---------------------------------------------------------

void SpaceVariable3D::WriteToVTRFile(const char *filename)
{
  PetscViewer viewer;
  PetscViewerVTKOpen(PetscObjectComm((PetscObject)*dm), filename, FILE_MODE_WRITE, &viewer);
  VecView(globalVec, viewer);
  PetscViewerDestroy(&viewer);
}

//---------------------------------------------------------

void SpaceVariable3D::AXPlusB(double a, double b, bool workOnGhost)
{
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
void SpaceVariable3D::AXPlusBY(double a, double b, SpaceVariable3D &y, std::vector<int> Xindices,
                               std::vector<int> Yindices, bool workOnGhost)
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
        for(int id=0; id<Xindices.size(); id++) {
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

  if(!workOnGhost)
    RestoreDataPointerAndInsert(); 
  else
    RestoreDataPointerToLocalVector(); //no need to communicate because the ghost region has been set to a const
}

//---------------------------------------------------------

