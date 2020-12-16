#include <SpaceVariable.h>

//---------------------------------------------------------
// DataManagers2D
//---------------------------------------------------------
// static member variables
DM DataManagers2D::ghosted1_1dof; //ghosted"1" --> stencil width is 1
DM DataManagers2D::ghosted1_2dof;
DM DataManagers2D::ghosted1_3dof;
DM DataManagers2D::ghosted1_4dof;
DM DataManagers2D::ghosted1_5dof;

//---------------------------------------------------------

DataManagers2D::DataManagers2D()
{

}

//---------------------------------------------------------

DataManagers2D::DataManagers2D(MPI_Comm comm, int NX, int NY)
{
  CreateAllDataManagers(comm, NX, NY);
}

//---------------------------------------------------------

DataManagers2D::~DataManagers2D()
{
 //DMDestroy needs to be called before PetscFinalize()!
}

//---------------------------------------------------------

int DataManagers2D::CreateAllDataManagers(MPI_Comm comm, int NX, int NY)
{
  int nProcX, nProcY; //All DM's should use the same domain partition

  auto ierr = DMDACreate2d(comm, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX,
                           NX, NY, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL,
                           &ghosted1_1dof);
  CHKERRQ(ierr);
  DMSetFromOptions(ghosted1_1dof);
  DMSetUp(ghosted1_1dof);

  DMDAGetInfo(ghosted1_1dof, NULL, NULL, NULL, NULL, &nProcX, &nProcY, NULL, NULL, NULL, NULL,
              NULL, NULL, NULL); 

  ierr = DMDACreate2d(comm, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX,
                      NX, NY, nProcX, nProcY, 2, 1, NULL, NULL,
                      &ghosted1_2dof);
  CHKERRQ(ierr);
  DMSetFromOptions(ghosted1_2dof);
  DMSetUp(ghosted1_2dof);

  ierr = DMDACreate2d(comm, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX,
                      NX, NY, nProcX, nProcY, 3, 1, NULL, NULL,
                      &ghosted1_3dof);
  CHKERRQ(ierr);
  DMSetFromOptions(ghosted1_3dof);
  DMSetUp(ghosted1_3dof);

  ierr = DMDACreate2d(comm, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX,
                      NX, NY, nProcX, nProcY, 4, 1, NULL, NULL,
                      &ghosted1_4dof);
  CHKERRQ(ierr);
  DMSetFromOptions(ghosted1_4dof);
  DMSetUp(ghosted1_4dof);

  ierr = DMDACreate2d(comm, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX,
                      NX, NY, nProcX, nProcY, 5, 1, NULL, NULL,
                      &ghosted1_5dof);
  CHKERRQ(ierr);
  DMSetFromOptions(ghosted1_5dof);
  DMSetUp(ghosted1_5dof);

  return 0;
}

//---------------------------------------------------------

void DataManagers2D::DestroyAllDataManagers()
{
  DMDestroy(&ghosted1_1dof);
  DMDestroy(&ghosted1_2dof);
  DMDestroy(&ghosted1_3dof);
  DMDestroy(&ghosted1_4dof);
  DMDestroy(&ghosted1_5dof);
}

//---------------------------------------------------------
// SpaceVariable2D
//---------------------------------------------------------

SpaceVariable2D::SpaceVariable2D(MPI_Comm &comm_, DM *dm_) : comm(comm_), globalVec(), localVec()
{
  dm = dm_;

  auto ierr = DMCreateGlobalVector(*dm, &globalVec); 
  //CHKERRQ(ierr);
  VecSet(globalVec, 0.0);

  ierr = DMCreateLocalVector(*dm, &localVec);
  //CHKERRQ(ierr);
  VecSet(localVec, 0.0);

  array = NULL;

  DMBoundaryType bx, by;

  DMDAGetInfo(*dm, NULL, &NX, &NY, NULL, &nProcX, &nProcY, NULL, &dof, &ghost_width, 
              &bx, &by, NULL, NULL);
  
  if(bx==DM_BOUNDARY_NONE && by==DM_BOUNDARY_NONE) {
    ghosted = false;
    ghost_width = 0;
  } else if(bx==DM_BOUNDARY_GHOSTED && by==DM_BOUNDARY_GHOSTED) {
    ghosted = true;
  } else {
    PetscPrintf(comm, "ERROR: Unsupported ghost type.\n");
    MPI_Abort(comm, 1); 
  }

  DMDAGetCorners(*dm, &i0, &j0, NULL, &nx, &ny, NULL);
  imax = i0 + nx;
  jmax = j0 + ny;

  DMDAGetGhostCorners(*dm, &ghost_i0, &ghost_j0, NULL, &ghost_nx, &ghost_ny, NULL);
  ghost_imax = ghost_i0 + ghost_nx;
  ghost_jmax = ghost_j0 + ghost_ny;

}

//---------------------------------------------------------

SpaceVariable2D::~SpaceVariable2D() 
{
  Destroy();
}

//---------------------------------------------------------

double** SpaceVariable2D::GetDataPointer()
{
  DMGlobalToLocalBegin(*dm, globalVec, INSERT_VALUES, localVec);
  DMGlobalToLocalEnd(*dm, globalVec, INSERT_VALUES, localVec);
  auto ierr = DMDAVecGetArray(*dm, localVec, &array);
  //CHKERRQ(ierr);
  return array;
}

//---------------------------------------------------------

void SpaceVariable2D::RestoreDataPointerAndInsert()
{
  RestoreDataPointerToLocalVector();
  auto ierr = DMLocalToGlobal(*dm, localVec, INSERT_VALUES, globalVec);
  //CHKERRQ(ierr);
}

//---------------------------------------------------------

void SpaceVariable2D::RestoreDataPointerAndAdd()
{
  RestoreDataPointerToLocalVector();
  auto ierr = DMLocalToGlobal(*dm, localVec, ADD_VALUES, globalVec);
  //CHKERRQ(ierr);
}

//---------------------------------------------------------

void SpaceVariable2D::RestoreDataPointerToLocalVector()
{
  auto ierr = DMDAVecRestoreArray(*dm, localVec, &array);  
  //CHKERRQ(ierr);
}

//---------------------------------------------------------

void SpaceVariable2D::Destroy()
{
  VecDestroy(&globalVec);
  VecDestroy(&localVec);
}

//---------------------------------------------------------

void SpaceVariable2D::StoreMeshCoordinates(SpaceVariable2D &coordinates)
{
  auto ierr = DMSetCoordinateDim(*dm, 2/*2D*/);
  //CHKERRQ(ierr);
  ierr = DMSetCoordinates(*dm, coordinates.globalVec);
  //CHKERRQ(ierr);
}

//---------------------------------------------------------

void SpaceVariable2D::WriteToVTRFile(const char *filename)
{
  PetscViewer viewer;
  PetscViewerVTKOpen(PetscObjectComm((PetscObject)*dm), filename, FILE_MODE_WRITE, &viewer);
  VecView(globalVec, viewer);
  PetscViewerDestroy(&viewer);
}


