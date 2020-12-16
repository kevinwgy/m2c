#include <SpaceOperator.h>
#include <Vector2D.h>
#include <Vector5D.h>
using std::cout;
using std::endl;

//-----------------------------------------------------

SpaceOperator::SpaceOperator(MPI_Comm &comm_, DataManagers2D &dm_all_, IoData &iod_,
                             VarFcnBase &vf_) 
  : comm(comm_), dm_all(dm_all_), iod(iod_), vf(vf_), 
    coordinates(comm_, &(dm_all_.ghosted1_2dof)),
    delta_xy(comm_, &(dm_all_.ghosted1_2dof)),
    face_rt(comm_, &(dm_all_.ghosted1_2dof)),
    volume(comm_, &(dm_all_.ghosted1_1dof))
{
  SetupMesh();
}

//-----------------------------------------------------

SpaceOperator::~SpaceOperator()
{

}

//-----------------------------------------------------

void SpaceOperator::Destroy()
{
  coordinates.Destroy();
  delta_xy.Destroy();
  face_rt.Destroy();
  volume.Destroy();
}

//-----------------------------------------------------

void SpaceOperator::SetupMesh()
{
  // Setup nodal coordinates
  if(true)
    SetupNodalCoordinatesUniformRectangularDomain();
  // TODO: add more choices later


  // Compute mesh information
  Vec2D** coords = (Vec2D**)coordinates.GetDataPointer(); 
  Vec2D**  dxy   = (Vec2D**)delta_xy.GetDataPointer();
  Vec2D**  frt   = (Vec2D**)face_rt.GetDataPointer();
  double** vol   = (double**)volume.GetDataPointer();

  int i0, j0, imax, jmax, ii0, jj0;
  coordinates.GetCornerIndices(&i0, &j0, &imax, &jmax);
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, NULL, NULL);

  // Calculate dx and dy in the forward difference fashion. 
  // Include ghost nodes outside the bottom and the left domain boundaries
  for(int j=jj0; j<jmax; j++)
    for(int i=ii0; i<imax; i++) {
      dxy[j][i][0] = coords[j][i+1][0] - coords[j][i][0]; 
      dxy[j][i][1] = coords[j+1][i][1] - coords[j][i][1]; 
    }

  // Calculate the area/length of right and top faces
  // Calculate the volume/area of node-centered control volumes ("cells")
  // Include ghost nodes outside the bottom and the left domain boundaries,
  // but not the ghost node at the bottom-left corner.
  for(int j=jj0; j<jmax; j++)
    for(int i=ii0; i<imax; i++) {
      if(j!=-1)
        frt[j][i][0]/*right face*/ = 0.5*(dxy[j-1][i][1] + dxy[j][i][1]);
      if(i!=-1)
        frt[j][i][1]/*top face*/   = 0.5*(dxy[j][i-1][0] + dxy[j][i][0]); 
      if(i!=-1 && j!=-1)
        vol[j][i]   /*area of cv*/ = frt[j][i][0]*frt[j][i][1];
    }

  coordinates.RestoreDataPointerToLocalVector(); //no changes have been made
  delta_xy.RestoreDataPointerAndInsert();
  face_rt.RestoreDataPointerAndInsert();
  volume.RestoreDataPointerAndInsert();
}

//-----------------------------------------------------

void SpaceOperator::SetupNodalCoordinatesUniformRectangularDomain()
{
  int i0, j0, imax, jmax, NX, NY;
  coordinates.GetCornerIndices(&i0, &j0, &imax, &jmax);
  coordinates.GetGlobalSize(&NX, &NY);

  double dx = (iod.mesh.xmax - iod.mesh.x0)/NX;
  double dy = (iod.mesh.ymax - iod.mesh.y0)/NY;

  // get array to edit
  Vec2D** coords = (Vec2D**)coordinates.GetDataPointer();

  // Fill the actual subdomain, w/o ghost layer
  for(int j=j0; j<jmax; j++)
    for(int i=i0; i<imax; i++) {
      coords[j][i][0] = iod.mesh.x0 + 0.5*dx + i*dx; 
      coords[j][i][1] = iod.mesh.y0 + 0.5*dy + j*dy; 
    } 

  // restore array
  coordinates.RestoreDataPointerAndInsert(); //update localVec and globalVec;

  // Populate the ghost boundary layer
  PopulateGhostBoundaryNodalCoordinates();
}

//-----------------------------------------------------

void SpaceOperator::PopulateGhostBoundaryNodalCoordinates()
{
  Vec2D** v = (Vec2D**) coordinates.GetDataPointer();

  int ii0, jj0, nnx, nny, NX, NY;
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, NULL, NULL);
  coordinates.GetGhostedSize(&nnx, &nny);
  coordinates.GetGlobalSize(&NX, &NY);

  if(ii0 == -1) {//Left
    for(int j=jj0; j<jj0+nny; j++) {
      v[j][ii0][0] = 2.0*iod.mesh.x0 - v[j][ii0+1][0];
      v[j][ii0][1] = v[j][ii0+1][1];
    }
  }
  if(ii0+nnx == NX+1) {//Right
    for(int j=jj0; j<jj0+nny; j++) {
      v[j][ii0+nnx-1][0] = 2.0*iod.mesh.xmax - v[j][ii0+nnx-2][0];
      v[j][ii0+nnx-1][1] = v[j][ii0+nnx-2][1];
    }
  }
  if(jj0 == -1) {//Bottom
    for(int i=ii0; i<ii0+nnx; i++) {
      v[jj0][i][0] = v[jj0+1][i][0];
      v[jj0][i][1] = 2.0*iod.mesh.y0 - v[jj0+1][i][1];
    }
  }
  if(jj0+nny == NY+1) {//Top
    for(int i=ii0; i<ii0+nnx; i++) {
      v[jj0+nny-1][i][0] = v[jj0+nny-2][i][0];
      v[jj0+nny-1][i][1] = 2.0*iod.mesh.ymax - v[jj0+nny-2][i][1];
    }
  }
  if(ii0 == -1 && jj0 == -1) {//Bottom-Left Corner
    v[jj0][ii0][0] = v[jj0+1][ii0][0];
    v[jj0][ii0][1] = v[jj0][ii0+1][1];
  }
  if(ii0 == -1 && jj0+nny == NY+1) {//Top-Left Corner
    v[jj0+nny-1][ii0][0] = v[jj0+nny-2][ii0][0];
    v[jj0+nny-1][ii0][1] = v[jj0+nny-1][ii0+1][1];
  }
  if(ii0+nnx == NX+1 && jj0 == -1) {//Bottom-Right Corner
    v[jj0][ii0+nnx-1][0] = v[jj0+1][ii0+nnx-1][0];
    v[jj0][ii0+nnx-1][1] = v[jj0][ii0+nnx-2][1];
  }
  if(ii0+nnx == NX+1 && jj0+nny == NY+1) {//Top-Right Corner
    v[jj0+nny-1][ii0+nnx-1][0] = v[jj0+nny-2][ii0+nnx-1][0];
    v[jj0+nny-1][ii0+nnx-1][1] = v[jj0+nny-1][ii0+nnx-2][1];
  }

  coordinates.RestoreDataPointerAndInsert();
}

//-----------------------------------------------------
// U and V should have 5 DOFs per node, like in 3D. This allows
// us to re-use the same functions for 3D. Just set z-velocity/momentum
// to 0.
void SpaceOperator::ConservativeToPrimitive(SpaceVariable2D &U, SpaceVariable2D &V, bool workOnGhost)
{
  Vec5D** u = (Vec5D**) U.GetDataPointer();
  Vec5D** v = (Vec5D**) V.GetDataPointer();

  int i0, j0, imax, jmax;
  if(workOnGhost)
    U.GetGhostedCornerIndices(&i0, &j0, &imax, &jmax);
  else
    U.GetCornerIndices(&i0, &j0, &imax, &jmax);

  for(int j=j0; j<jmax; j++)
    for(int i=i0; i<imax; i++)
      vf.conservativeToPrimitive((double*)u[j][i], (double*)v[j][i]); 

  U.RestoreDataPointerToLocalVector(); //no changes made
  V.RestoreDataPointerAndInsert();
}

//-----------------------------------------------------
// U and V should have 5 DOFs per node, like in 3D. This allows
// us to re-use the same functions for 3D. Just set z-velocity/momentum
// to 0.
void SpaceOperator::PrimitiveToConservative(SpaceVariable2D &V, SpaceVariable2D &U, bool workOnGhost)
{
  Vec5D** v = (Vec5D**) V.GetDataPointer();
  Vec5D** u = (Vec5D**) U.GetDataPointer();

  int i0, j0, imax, jmax;
  if(workOnGhost)
    V.GetGhostedCornerIndices(&i0, &j0, &imax, &jmax);
  else
    V.GetCornerIndices(&i0, &j0, &imax, &jmax);

  for(int j=j0; j<jmax; j++)
    for(int i=i0; i<imax; i++)
      vf.primitiveToConservative((double*)v[j][i], (double*)u[j][i]); 

  V.RestoreDataPointerToLocalVector(); //no changes made
  U.RestoreDataPointerAndInsert();
}

//-----------------------------------------------------

void SpaceOperator::SetInitialCondition(SpaceVariable2D &V) //apply IC within the real domain
{
  Vec5D** v = (Vec5D**) V.GetDataPointer();

  int i0, j0, imax, jmax;
  V.GetCornerIndices(&i0, &j0, &imax, &jmax);

  // First, apply the farfield state
  for(int j=j0; j<jmax; j++)
    for(int i=i0; i<imax; i++) {
      v[j][i][0] = iod.bc.farfield.density;
      v[j][i][1] = iod.bc.farfield.velocity_x;
      v[j][i][2] = iod.bc.farfield.velocity_y;
      v[j][i][3] = iod.bc.farfield.velocity_z;
      v[j][i][4] = iod.bc.farfield.pressure;
    }

  // Second, apply user-specified function


  // Apply boundary condition to populate ghost nodes
}

//-----------------------------------------------------

void SpaceOperator::SetBoundaryConditions(SpaceVariable2D &U)
{
  //TODO
}

//-----------------------------------------------------

void SpaceOperator::ComputeAdvectionFluxes(SpaceVariable2D &U, SpaceVariable2D &F)
{
  //TODO
}

//-----------------------------------------------------

