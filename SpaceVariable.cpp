#include <SpaceVariable.h>

//---------------------------------------------------------
// DataManagers2D
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

void DataManagers2D::CreateAllDataManagers(MPI_Comm comm, int NX, int NY)
{
  int nProcX, nProcY; //All DM's should use the same domain partition

  auto ierr = DMDACreate2D(comm, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_STENCIL_BOX,
                           NX, NY, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL,
                           &ghosted_1dof);
  CHKERRXX(ierr);

  DMDAGetInfo(ghosted_1dof, NULL, NULL, NULL, NULL, &nProcX, &nProcY, NULL, NULL, NULL, NULL,
              NULL, NULL, NULL); 

  ierr = DMDACreate2D(comm, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_STENCIL_BOX,
                      NX, NY, nProcX, nProcY, PETSC_DECIDE, 2, 1, NULL, NULL,
                      &ghosted_2dof);
  CHKERRXX(ierr);

  ierr = DMDACreate2D(comm, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_STENCIL_BOX,
                      NX, NY, nProcX, nProcY, PETSC_DECIDE, 3, 1, NULL, NULL,
                      &ghosted_3dof);
  CHKERRXX(ierr);

  ierr = DMDACreate2D(comm, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_STENCIL_BOX,
                      NX, NY, nProcX, nProcY, PETSC_DECIDE, 4, 1, NULL, NULL,
                      &ghosted_4dof);
  CHKERRXX(ierr);

  ierr = DMDACreate2D(comm, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_STENCIL_BOX,
                      NX, NY, nProcX, nProcY, PETSC_DECIDE, 5, 1, NULL, NULL,
                      &ghosted_5dof);
  CHKERRXX(ierr);
}

//---------------------------------------------------------

void DestroyAllDataManagers()
{
  DMDestroy(&ghosted_1dof);
  DMDestroy(&ghosted_2dof);
  DMDestroy(&ghosted_3dof);
  DMDestroy(&ghosted_4dof);
  DMDestroy(&ghosted_5dof);
}

//---------------------------------------------------------
// SpaceVariable2D
//---------------------------------------------------------

SpaceVariable2D::SpaceVariable2D(DM *dm_) : comm(comm_), globalVec(), localVec()
{
  dm = dm_;

  auto ierr = DMCreateGlobalVector(*dm, &globalVec); 
  CHKERRXX(ierr);
  VecSet(globalVec, 0.0);

  ierr = DMCreateLocalVector(*dm, &localVec);
  CHKERRXX(ierr);
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
    PETSCABORT(comm, 1); 
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

}

//---------------------------------------------------------

void* SpaceVariable2D::GetDataPointer()
{
  DMGlobalToLocalBegin(*dm, globalVec, INSERT_VALUES, localVec);
  DMGlobalToLocalEnd(*dm, globalVec, INSERT_VALUES, localVec);
  auto ierr = DMDAVecGetArray(*dm, localVec, array);
  CHKERRXX(ierr);
  return array;
}

//---------------------------------------------------------

void SpaceVariable2D::RestoreDataPointer()
{
  auto ierr = DMDAVecRestoreArray(*dm, localVec, array);  
  CHKERRXX(ierr);
}

//---------------------------------------------------------

void SpaceVariable2D::AssembleInsert()
{
  auto ierr = DMDALocalToGlobal(*dm, localVec, INSERT_VALUES, globalVec);
  CHKERRXX(ierr);
}

//---------------------------------------------------------

void SpaceVariable2D::AssembleAdd()
{
  auto ierr = DMDALocalToGlobal(*dm, localVec, ADD_VALUES, globalVec);
  CHKERRXX(ierr);
}

//---------------------------------------------------------

void SpaceVariable2D::Destroy()
{
  VecDestroy(&globalVec);
  VecDestroy(&localVec);
}

//---------------------------------------------------------
















