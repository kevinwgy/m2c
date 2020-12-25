#include <Utils.h>
#include <SpaceOperator.h>
#include <Vector2D.h>
#include <Vector5D.h>
#include <algorithm> //std::upper_bound
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
  
  coordinates.GetCornerIndices(&i0, &j0, &imax, &jmax);
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, &iimax, &jjmax);

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
  //! Setup nodal coordinates
  if(true)
    SetupNodalCoordinatesUniformRectangularDomain();
  //! TODO: add more choices later


  //! Compute mesh information
  Vec2D** coords = (Vec2D**)coordinates.GetDataPointer(); 
  Vec2D**  dxy   = (Vec2D**)delta_xy.GetDataPointer();
  Vec2D**  frt   = (Vec2D**)face_rt.GetDataPointer();
  double** vol   = (double**)volume.GetDataPointer();

  //! Calculate dx and dy in the forward difference fashion. 
  //! Include ghost nodes outside the bottom and the left domain boundaries
  for(int j=jj0; j<jmax; j++)
    for(int i=ii0; i<imax; i++) {
      dxy[j][i][0] = coords[j][i+1][0] - coords[j][i][0]; 
      dxy[j][i][1] = coords[j+1][i][1] - coords[j][i][1]; 
    }

  /** Calculate the area/length of right and top faces
   *  Calculate the volume/area of node-centered control volumes ("cells")
   *  Include ghost nodes outside the bottom and the left domain boundaries,
   *  but not the ghost node at the bottom-left corner.
   */
  for(int j=jj0; j<jmax; j++)
    for(int i=ii0; i<imax; i++) {
      if(j!=jj0)
        frt[j][i][0]/*right face*/ = 0.5*(dxy[j-1][i][1] + dxy[j][i][1]);
      if(i!=ii0)
        frt[j][i][1]/*top face*/   = 0.5*(dxy[j][i-1][0] + dxy[j][i][0]); 
      if(i!=ii0 && j!=jj0)
        vol[j][i]   /*area of cv*/ = frt[j][i][0]*frt[j][i][1];
    }

  coordinates.RestoreDataPointerToLocalVector(); //!< no changes have been made
  delta_xy.RestoreDataPointerAndInsert();
  face_rt.RestoreDataPointerAndInsert();
  volume.RestoreDataPointerAndInsert();
}

//-----------------------------------------------------

void SpaceOperator::SetupNodalCoordinatesUniformRectangularDomain()
{
  int NX, NY;
  coordinates.GetGlobalSize(&NX, &NY);

  double dx = (iod.mesh.xmax - iod.mesh.x0)/NX;
  double dy = (iod.mesh.ymax - iod.mesh.y0)/NY;

  //! get array to edit
  Vec2D** coords = (Vec2D**)coordinates.GetDataPointer();

  //! Fill the actual subdomain, w/o ghost layer
  for(int j=j0; j<jmax; j++)
    for(int i=i0; i<imax; i++) {
      coords[j][i][0] = iod.mesh.x0 + 0.5*dx + i*dx; 
      coords[j][i][1] = iod.mesh.y0 + 0.5*dy + j*dy; 
    } 

  //! restore array
  coordinates.RestoreDataPointerAndInsert(); //update localVec and globalVec;

  //! Populate the ghost boundary layer
  PopulateGhostBoundaryNodalCoordinates();
}

//-----------------------------------------------------

void SpaceOperator::PopulateGhostBoundaryNodalCoordinates()
{
  Vec2D** v = (Vec2D**) coordinates.GetDataPointer();

  int nnx, nny, NX, NY;
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
/** 
 * U and V should have 5 DOFs per node, like in 3D. This allows
 * us to re-use the same functions for 3D. Just set z-velocity/momentum
 * to 0.
 */
void SpaceOperator::ConservativeToPrimitive(SpaceVariable2D &U, SpaceVariable2D &V, bool workOnGhost)
{
  Vec5D** u = (Vec5D**) U.GetDataPointer();
  Vec5D** v = (Vec5D**) V.GetDataPointer();

  int myi0, myj0, myimax, myjmax;
  if(workOnGhost)
    U.GetGhostedCornerIndices(&myi0, &myj0, &myimax, &myjmax);
  else
    U.GetCornerIndices(&myi0, &myj0, &myimax, &myjmax);

  for(int j=myj0; j<myjmax; j++)
    for(int i=myi0; i<myimax; i++)
      vf.conservativeToPrimitive((double*)u[j][i], (double*)v[j][i]); 

  U.RestoreDataPointerToLocalVector(); //no changes made
  V.RestoreDataPointerAndInsert();
}

//-----------------------------------------------------
/** U and V should have 5 DOFs per node, like in 3D. This allows
 *  us to re-use the same functions for 3D. Just set z-velocity/momentum
 *  to 0.
 */
void SpaceOperator::PrimitiveToConservative(SpaceVariable2D &V, SpaceVariable2D &U, bool workOnGhost)
{
  Vec5D** v = (Vec5D**) V.GetDataPointer();
  Vec5D** u = (Vec5D**) U.GetDataPointer();

  int myi0, myj0, myimax, myjmax;
  if(workOnGhost)
    U.GetGhostedCornerIndices(&myi0, &myj0, &myimax, &myjmax);
  else
    U.GetCornerIndices(&myi0, &myj0, &myimax, &myjmax);

  for(int j=myj0; j<myjmax; j++)
    for(int i=myi0; i<myimax; i++)
      vf.primitiveToConservative((double*)v[j][i], (double*)u[j][i]); 

  V.RestoreDataPointerToLocalVector(); //no changes made
  U.RestoreDataPointerAndInsert();
}

//-----------------------------------------------------

void SpaceOperator::SetInitialCondition(SpaceVariable2D &V) //apply IC within the real domain
{
  Vec5D** v = (Vec5D**) V.GetDataPointer();

  //! First, apply the inlet (i.e. farfield) state
  for(int j=jj0; j<jjmax; j++)
    for(int i=ii0; i<iimax; i++) {
      v[j][i][0] = iod.bc.inlet.density;
      v[j][i][1] = iod.bc.inlet.velocity_x;
      v[j][i][2] = iod.bc.inlet.velocity_y;
      v[j][i][3] = iod.bc.inlet.velocity_z;
      v[j][i][4] = iod.bc.inlet.pressure;
    }

  //! Second, apply user-specified function
  if(iod.ic.type != IcData::NONE) {

    //! Get coordinates
    Vec2D** coords = (Vec2D**)coordinates.GetDataPointer();
    Vec2D   x0(iod.ic.x0[0], iod.ic.x0[1]); //!< 3D -> 2D
    Vec2D   dir(iod.ic.dir[0], iod.ic.dir[1]); //!< 3D -> 2D
    dir /= dir.norm();

    if(iod.ic.type == IcData::PLANAR) {
      print("- Applying user-specified initial condition (with planar symmetry).\n");
 
      double x;
      int n = iod.ic.user_data[IcData::COORDINATE].size(); //!< number of data points provided by user

      int k0, k1;    
      double a0, a1;
//      cout << "i0 = " << i0 << ", imax = " << imax << ", j0 = " << j0 << ", jmax = " << jmax << endl;
//      cout << "ii0 = " << ii0 << ", iimax = " << iimax << ", jj0 = " << jj0 << ", jjmax = " << jjmax << endl;
//      cout << "n = " << n << endl;
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++) {

          x = (coords[j][i] - x0)*dir; //!< projection onto the 1D axis
//          cout << "coords: " << coords[j][i][0] << ", " << coords[j][i][1] << "; x = " << x << endl;
          if(x<0 || x>iod.ic.user_data[IcData::COORDINATE][n-1])
            continue;
 
          //! Find the first 1D coordinate greater than x
          auto upper_it = std::upper_bound(iod.ic.user_data[IcData::COORDINATE].begin(),
                                           iod.ic.user_data[IcData::COORDINATE].end(),
                                           x); 
          k1 = (int)(upper_it - iod.ic.user_data[IcData::COORDINATE].begin());

          if(k1==0) // exactly the first node in 1D
            k1 = 1;

          k0 = k1 - 1;

          //! calculate interpolation weights
          a0 = (iod.ic.user_data[IcData::COORDINATE][k1] - x) /
               (iod.ic.user_data[IcData::COORDINATE][k1] - iod.ic.user_data[IcData::COORDINATE][k0]);
          a1 = 1.0 - a0;

//          cout << "k0 = " << k0 << ", k1 = " << k1 << ", a0 = " << a0 << ", a1 = " << a1 << endl;
//          cout << "coord_k1:" << iod.ic.user_data[IcData::COORDINATE][k1] << ", coord_k0:" << iod.ic.user_data[IcData::COORDINATE][k0] << endl;

          //! specify i.c. on node (cell center)
          v[j][i][0] = a0*iod.ic.user_data[IcData::DENSITY][k0] + a1*iod.ic.user_data[IcData::DENSITY][k1];
          v[j][i][1] = (a0*iod.ic.user_data[IcData::VELOCITY][k0] + a1*iod.ic.user_data[IcData::VELOCITY][k1])*dir[0];
          v[j][i][2] = (a0*iod.ic.user_data[IcData::VELOCITY][k0] + a1*iod.ic.user_data[IcData::VELOCITY][k1])*dir[1];
          v[j][i][3] = 0.0; //To be updated for 3D
          v[j][i][4] = a0*iod.ic.user_data[IcData::PRESSURE][k0] + a1*iod.ic.user_data[IcData::PRESSURE][k1];
        }
    } 
    else if (iod.ic.type == IcData::CYLINDRICAL) {
      print_error("ERROR: Cannot handle cylindrical i.c. at the moment.\n");
      exit_mpi();
    } 
    else if (iod.ic.type == IcData::SPHERICAL) {
      print_error("ERROR: Cannot handle spherical i.c. at the moment.\n");
      exit_mpi();
    }

    coordinates.RestoreDataPointerToLocalVector(); //!< data was not changed.
  }

  V.RestoreDataPointerAndInsert();

  //! Apply boundary condition to populate ghost nodes
  ApplyBoundaryConditions(V);   

}

//-----------------------------------------------------
//! Apply boundary conditions by populating the ghost cells.
void SpaceOperator::ApplyBoundaryConditions(SpaceVariable2D &V)
{
  Vec5D** v = (Vec5D**) V.GetDataPointer();

  int NX, NY;
  V.GetGlobalSize(&NX, &NY);
//  cout << "NX = " << NX << ", NY = " << NY << endl;

  //! Left boundary
  if(ii0==-1) { 
    switch (iod.mesh.bc_x0) {
      case MeshData::INLET :
        for(int j=j0; j<jmax; j++) {
          v[j][ii0][0] = iod.bc.inlet.density;
          v[j][ii0][1] = iod.bc.inlet.velocity_x;
          v[j][ii0][2] = iod.bc.inlet.velocity_y;
          v[j][ii0][3] = iod.bc.inlet.velocity_z;
          v[j][ii0][4] = iod.bc.inlet.pressure;
        }
        break;
      case MeshData::OUTLET :
        for(int j=j0; j<jmax; j++) {
          v[j][ii0][0] = iod.bc.outlet.density;
          v[j][ii0][1] = iod.bc.outlet.velocity_x;
          v[j][ii0][2] = iod.bc.outlet.velocity_y;
          v[j][ii0][3] = iod.bc.outlet.velocity_z;
          v[j][ii0][4] = iod.bc.outlet.pressure;
        }
        break; 
      case MeshData::WALL :
      case MeshData::SYMMETRY :
        for(int j=j0; j<jmax; j++) {
          v[j][ii0][0] =      v[j][ii0+1][0];
          v[j][ii0][1] = -1.0*v[j][ii0+1][1];
          v[j][ii0][2] =      v[j][ii0+1][2];
          v[j][ii0][3] =      v[j][ii0+1][3]; 
          v[j][ii0][4] =      v[j][ii0+1][4];
        }
        break;
      default :
        print_error("ERROR: Boundary condition at x=x0 cannot be specified!\n");
        exit_mpi();
    }
  }

  //! Right boundary
  if(iimax==NX+1) { 
    switch (iod.mesh.bc_xmax) {
      case MeshData::INLET :
        for(int j=j0; j<jmax; j++) {
          v[j][iimax-1][0] = iod.bc.inlet.density;
          v[j][iimax-1][1] = iod.bc.inlet.velocity_x;
          v[j][iimax-1][2] = iod.bc.inlet.velocity_y;
          v[j][iimax-1][3] = iod.bc.inlet.velocity_z;
          v[j][iimax-1][4] = iod.bc.inlet.pressure;
        }
        break;
      case MeshData::OUTLET :
        for(int j=j0; j<jmax; j++) {
          v[j][iimax-1][0] = iod.bc.outlet.density;
          v[j][iimax-1][1] = iod.bc.outlet.velocity_x;
          v[j][iimax-1][2] = iod.bc.outlet.velocity_y;
          v[j][iimax-1][3] = iod.bc.outlet.velocity_z;
          v[j][iimax-1][4] = iod.bc.outlet.pressure;
        }
        break; 
      case MeshData::WALL :
      case MeshData::SYMMETRY :
        for(int j=j0; j<jmax; j++) {
          v[j][iimax-1][0] =      v[j][iimax-2][0];
          v[j][iimax-1][1] = -1.0*v[j][iimax-2][1];
          v[j][iimax-1][2] =      v[j][iimax-2][2];
          v[j][iimax-1][3] =      v[j][iimax-2][3]; 
          v[j][iimax-1][4] =      v[j][iimax-2][4];
        }
        break;
      default :
        print_error("ERROR: Boundary condition at x=xmax cannot be specified!\n");
        exit_mpi();
    }
  }

  //! Bottom boundary
  if(jj0==-1) { 
    switch (iod.mesh.bc_y0) {
      case MeshData::INLET :
        for(int i=i0; i<imax; i++) {
          v[jj0][i][0] = iod.bc.inlet.density;
          v[jj0][i][1] = iod.bc.inlet.velocity_x;
          v[jj0][i][2] = iod.bc.inlet.velocity_y;
          v[jj0][i][3] = iod.bc.inlet.velocity_z;
          v[jj0][i][4] = iod.bc.inlet.pressure;
        }
        break;
      case MeshData::OUTLET :
        for(int i=i0; i<imax; i++) {
          v[jj0][i][0] = iod.bc.outlet.density;
          v[jj0][i][1] = iod.bc.outlet.velocity_x;
          v[jj0][i][2] = iod.bc.outlet.velocity_y;
          v[jj0][i][3] = iod.bc.outlet.velocity_z;
          v[jj0][i][4] = iod.bc.outlet.pressure;
        }
        break; 
      case MeshData::WALL :
      case MeshData::SYMMETRY :
        for(int i=i0; i<imax; i++) {
          v[jj0][i][0] =      v[jj0+1][i][0];
          v[jj0][i][1] =      v[jj0+1][i][1];
          v[jj0][i][2] = -1.0*v[jj0+1][i][2];
          v[jj0][i][3] =      v[jj0+1][i][3]; 
          v[jj0][i][4] =      v[jj0+1][i][4];
        }
        break;
      default :
        print_error("ERROR: Boundary condition at y=y0 cannot be specified!\n");
        exit_mpi();
    }
  }

  //! Bottom boundary
  if(jjmax==NY+1) { 
    switch (iod.mesh.bc_ymax) {
      case MeshData::INLET :
        for(int i=i0; i<imax; i++) {
          v[jjmax-1][i][0] = iod.bc.inlet.density;
          v[jjmax-1][i][1] = iod.bc.inlet.velocity_x;
          v[jjmax-1][i][2] = iod.bc.inlet.velocity_y;
          v[jjmax-1][i][3] = iod.bc.inlet.velocity_z;
          v[jjmax-1][i][4] = iod.bc.inlet.pressure;
        }
        break;
      case MeshData::OUTLET :
        for(int i=i0; i<imax; i++) {
          v[jjmax-1][i][0] = iod.bc.outlet.density;
          v[jjmax-1][i][1] = iod.bc.outlet.velocity_x;
          v[jjmax-1][i][2] = iod.bc.outlet.velocity_y;
          v[jjmax-1][i][3] = iod.bc.outlet.velocity_z;
          v[jjmax-1][i][4] = iod.bc.outlet.pressure;
        }
        break; 
      case MeshData::WALL :
      case MeshData::SYMMETRY :
        for(int i=i0; i<imax; i++) {
          v[jjmax-1][i][0] =      v[jjmax-2][i][0];
          v[jjmax-1][i][1] =      v[jjmax-2][i][1];
          v[jjmax-1][i][2] = -1.0*v[jjmax-2][i][2];
          v[jjmax-1][i][3] =      v[jjmax-2][i][3]; 
          v[jjmax-1][i][4] =      v[jjmax-2][i][4];
        }
        break;
      default :
        print_error("ERROR: Boundary condition at y=ymax cannot be specified!\n");
        exit_mpi();
    }
  }

  V.RestoreDataPointerAndInsert();
}

//-----------------------------------------------------

void SpaceOperator::ComputeTimeStepSize(SpaceVariable2D &V, double &dt, double &cfl)
{
  if(iod.ts.timestep > 0)
    dt = iod.ts.timestep;
  else {
    print_error("ERROR: Unable to calculate time step size.\n");
    exit_mpi();
  }
}

//-----------------------------------------------------

void SpaceOperator::ComputeAdvectionFluxes(SpaceVariable2D &U, SpaceVariable2D &F)
{
  //TODO
}

//-----------------------------------------------------

