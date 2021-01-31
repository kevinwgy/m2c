#include <Utils.h>
#include <SpaceOperator.h>
#include <Vector3D.h>
#include <Vector5D.h>
#include <algorithm> //std::upper_bound
#include <cfloat> //DBL_MAX
using std::cout;
using std::endl;
using std::max;
using std::min;

//-----------------------------------------------------

SpaceOperator::SpaceOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_,
                             VarFcnBase &varFcn_, FluxFcnBase &fluxFcn_) 
  : comm(comm_), dm_all(dm_all_), iod(iod_), varFcn(varFcn_), fluxFcn(fluxFcn_),
    coordinates(comm_, &(dm_all_.ghosted1_3dof)),
    delta_xyz(comm_, &(dm_all_.ghosted1_3dof)),
    volume(comm_, &(dm_all_.ghosted1_1dof)),
    rec(comm_, dm_all_, iod_, coordinates, delta_xyz),
    Vl(comm_, &(dm_all_.ghosted1_5dof)),
    Vr(comm_, &(dm_all_.ghosted1_5dof)),
    Vb(comm_, &(dm_all_.ghosted1_5dof)),
    Vt(comm_, &(dm_all_.ghosted1_5dof)),
    Vk(comm_, &(dm_all_.ghosted1_5dof)),
    Vf(comm_, &(dm_all_.ghosted1_5dof))
{
  
  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);

  SetupMesh();

  rec.Setup(); //this function requires mesh info (dxyz)
  
}

//-----------------------------------------------------

SpaceOperator::~SpaceOperator()
{

}

//-----------------------------------------------------

void SpaceOperator::Destroy()
{
  rec.Destroy();
  coordinates.Destroy();
  delta_xyz.Destroy();
  volume.Destroy();
  Vl.Destroy();
  Vr.Destroy();
  Vb.Destroy();
  Vt.Destroy();
  Vk.Destroy();
  Vf.Destroy();
}

//-----------------------------------------------------

void SpaceOperator::SetupMesh()
{
  //! Setup coordinates of cell centers and dx, dy, dz
  if(true)
    SetupMeshUniformRectangularDomain();
  //! TODO: add more choices later


  //! Compute mesh information
  Vec3D***  dxyz = (Vec3D***)delta_xyz.GetDataPointer();
  double*** vol  = (double***)volume.GetDataPointer();


  /** Calculate the volume/area of node-centered control volumes ("cells")
   *  Include ghost cells. 
   */
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++) {
        vol[k][j][i] /*volume of cv*/ = dxyz[k][j][i][0]*dxyz[k][j][i][1]*dxyz[k][j][i][2];
      }

  delta_xyz.RestoreDataPointerAndInsert();
  volume.RestoreDataPointerAndInsert();
}

//-----------------------------------------------------

void SpaceOperator::SetupMeshUniformRectangularDomain()
{
  int NX, NY, NZ;
  coordinates.GetGlobalSize(&NX, &NY, &NZ);

  double dx = (iod.mesh.xmax - iod.mesh.x0)/NX;
  double dy = (iod.mesh.ymax - iod.mesh.y0)/NY;
  double dz = (iod.mesh.zmax - iod.mesh.z0)/NZ;

  //! get array to edit
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
  Vec3D*** dxyz   = (Vec3D***)delta_xyz.GetDataPointer();

  //! Fill the actual subdomain, w/o ghost cells 
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        coords[k][j][i][0] = iod.mesh.x0 + 0.5*dx + i*dx; 
        coords[k][j][i][1] = iod.mesh.y0 + 0.5*dy + j*dy; 
        coords[k][j][i][2] = iod.mesh.z0 + 0.5*dz + k*dz; 
        dxyz[k][j][i][0] = dx;
        dxyz[k][j][i][1] = dy;
        dxyz[k][j][i][2] = dz;
      } 

  //! restore array
  coordinates.RestoreDataPointerAndInsert(); //update localVec and globalVec;
  delta_xyz.RestoreDataPointerAndInsert(); //update localVec and globalVec;

  //! Populate the ghost cells (coordinates, dx, dy, dz)
  PopulateGhostBoundaryCoordinates();
}

//-----------------------------------------------------
/** Populate the coordinates, dx, dy, and dz of ghost cells */
void SpaceOperator::PopulateGhostBoundaryCoordinates()
{
  Vec3D*** v    = (Vec3D***) coordinates.GetDataPointer();
  Vec3D*** dxyz = (Vec3D***) delta_xyz.GetDataPointer();

  int nnx, nny, nnz, NX, NY, NZ;
  coordinates.GetGhostedSize(&nnx, &nny, &nnz);
  coordinates.GetGlobalSize(&NX, &NY, &NZ);

  // capture the mesh info of the corners 
  double v0[3], v1[3];
  double dxyz0[3], dxyz1[3];
  for(int p=0; p<3; p++) {
    v0[p]    = v[k0][j0][i0][p] - dxyz[k0][j0][i0][p];
    v1[p]    = v[kmax][jmax][imax][p] + dxyz[kmax][jmax][imax][p];
    dxyz0[p] = dxyz[k0][j0][i0][p];
    dxyz1[p] = dxyz[kmax][jmax][imax][p];
  }

  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++) {

        if(k!=-1 && k!=NZ+1 && j!=-1 && j!=NY+1 && i!=-1 && i!=NX+1)
          continue; //not in the ghost layer of the physical domain

        Vec3D& X  = v[k][j][i];
        Vec3D& dX = dxyz[k][j][i];

        bool xdone = false, ydone = false, zdone = false;

        if(i==-1) {
          X[0]  = v0[0];
          dX[0] = dxyz0[0];
          xdone = true;
        }

        if(i==NX+1) {
          X[0]  = v1[0];
          dX[0] = dxyz1[0];
          xdone = true;
        }

        if(j==-1) {
          X[1]  = v0[1];
          dX[1] = dxyz0[1];
          ydone = true;
        }

        if(j==NY+1) {
          X[1]  = v1[1];
          dX[1] = dxyz1[1];
          ydone = true;
        }

        if(k==-1) {
          X[2]  = v0[2];
          dX[2] = dxyz0[2];
          zdone = true;
        }

        if(k==NZ+1) {
          X[2]  = v1[2];
          dX[2] = dxyz1[2];
          zdone = true;
        }

        if(!xdone) {
          X[0]  = v[k0][j0][i][0];    //x[i]
          dX[0] = dxyz[k0][j0][i][0]; //dx[i]
        }

        if(!ydone) {
          X[1]  = v[k0][j][i0][1];    //y[j]
          dX[1] = dxyz[k0][j][i0][1]; //dy[j]
        }

        if(!zdone) {
          X[2]  = v[k][j0][i0][2];    //z[k]
          dX[2] = dxyz[k][j0][i0][2]; //dz[k]
        }

      }
  
  coordinates.RestoreDataPointerAndInsert();
  delta_xyz.RestoreDataPointerAndInsert();
}

//-----------------------------------------------------

void SpaceOperator::ConservativeToPrimitive(SpaceVariable3D &U, SpaceVariable3D &V, bool workOnGhost)
{
  Vec5D*** u = (Vec5D***) U.GetDataPointer();
  Vec5D*** v = (Vec5D***) V.GetDataPointer();

  int myi0, myj0, myk0, myimax, myjmax, mykmax;
  if(workOnGhost)
    U.GetGhostedCornerIndices(&myi0, &myj0, &myk0, &myimax, &myjmax, &mykmax);
  else
    U.GetCornerIndices(&myi0, &myj0, &myk0, &myimax, &myjmax, &mykmax);

  for(int k=myk0; k<mykmax; k++)
    for(int j=myj0; j<myjmax; j++)
      for(int i=myi0; i<myimax; i++)
        varFcn.ConservativeToPrimitive((double*)u[k][j][i], (double*)v[k][j][i]); 

  U.RestoreDataPointerToLocalVector(); //no changes made
  V.RestoreDataPointerAndInsert();
}

//-----------------------------------------------------

void SpaceOperator::PrimitiveToConservative(SpaceVariable3D &V, SpaceVariable3D &U, bool workOnGhost)
{
  Vec5D*** v = (Vec5D***) V.GetDataPointer();
  Vec5D*** u = (Vec5D***) U.GetDataPointer();

  int myi0, myj0, myk0, myimax, myjmax, mykmax;
  if(workOnGhost)
    U.GetGhostedCornerIndices(&myi0, &myj0, &myk0, &myimax, &myjmax, &mykmax);
  else
    U.GetCornerIndices(&myi0, &myj0, &myk0, &myimax, &myjmax, &mykmax);

  for(int k=myk0; k<mykmax; k++)
    for(int j=myj0; j<myjmax; j++)
      for(int i=myi0; i<myimax; i++)
        varFcn.PrimitiveToConservative((double*)v[k][j][i], (double*)u[k][j][i]); 

  V.RestoreDataPointerToLocalVector(); //no changes made
  U.RestoreDataPointerAndInsert();
}

//-----------------------------------------------------

void SpaceOperator::SetInitialCondition(SpaceVariable3D &V) //apply IC within the real domain
{
  Vec5D*** v = (Vec5D***) V.GetDataPointer();

  //! First, apply the inlet (i.e. farfield) state
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++) {
        v[k][j][i][0] = iod.bc.inlet.density;
        v[k][j][i][1] = iod.bc.inlet.velocity_x;
        v[k][j][i][2] = iod.bc.inlet.velocity_y;
        v[k][j][i][3] = iod.bc.inlet.velocity_z;
        v[k][j][i][4] = iod.bc.inlet.pressure;
      }

  //! Second, apply user-specified function
  if(iod.ic.type != IcData::NONE) {

    //! Get coordinates
    Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
    Vec3D    x0(iod.ic.x0[0], iod.ic.x0[1], iod.ic.x0[2]); 
    Vec3D    dir(iod.ic.dir[0], iod.ic.dir[1], iod.ic.dir[2]); 
    dir /= dir.norm();

    if(iod.ic.type == IcData::PLANAR) {
      print("- Applying user-specified initial condition (with planar symmetry).\n");
 
      double x;
      int n = iod.ic.user_data[IcData::COORDINATE].size(); //!< number of data points provided by user

      int t0, t1;    
      double a0, a1;
//      cout << "i0 = " << i0 << ", imax = " << imax << ", j0 = " << j0 << ", jmax = " << jmax << endl;
//      cout << "ii0 = " << ii0 << ", iimax = " << iimax << ", jj0 = " << jj0 << ", jjmax = " << jjmax << endl;
//      cout << "n = " << n << endl;
      for(int k=k0; k<kmax; k++)
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {

            x = (coords[k][j][i] - x0)*dir; //!< projection onto the 1D axis
//            cout << "coords: " << coords[j][i][0] << ", " << coords[j][i][1] << "; x = " << x << endl;
            if(x<0 || x>iod.ic.user_data[IcData::COORDINATE][n-1])
              continue;
 
            //! Find the first 1D coordinate greater than x
            auto upper_it = std::upper_bound(iod.ic.user_data[IcData::COORDINATE].begin(),
                                             iod.ic.user_data[IcData::COORDINATE].end(),
                                             x); 
            t1 = (int)(upper_it - iod.ic.user_data[IcData::COORDINATE].begin());

            if(t1==0) // exactly the first node in 1D
              t1 = 1;

            t0 = t1 - 1;

            //! calculate interpolation weights
            a0 = (iod.ic.user_data[IcData::COORDINATE][t1] - x) /
                 (iod.ic.user_data[IcData::COORDINATE][t1] - iod.ic.user_data[IcData::COORDINATE][t0]);
            a1 = 1.0 - a0;

//            cout << "t0 = " << t0 << ", t1 = " << t1 << ", a0 = " << a0 << ", a1 = " << a1 << endl;
//            cout << "coord_t1:" << iod.ic.user_data[IcData::COORDINATE][t1] << ", coord_t0:" << iod.ic.user_data[IcData::COORDINATE][t0] << endl;

            //! specify i.c. on node (cell center)
            v[k][j][i][0] =  a0*iod.ic.user_data[IcData::DENSITY][t0]  + a1*iod.ic.user_data[IcData::DENSITY][t1];
            v[k][j][i][1] = (a0*iod.ic.user_data[IcData::VELOCITY][t0] + a1*iod.ic.user_data[IcData::VELOCITY][t1])*dir[0];
            v[k][j][i][2] = (a0*iod.ic.user_data[IcData::VELOCITY][t0] + a1*iod.ic.user_data[IcData::VELOCITY][t1])*dir[1];
            v[k][j][i][3] = (a0*iod.ic.user_data[IcData::VELOCITY][t0] + a1*iod.ic.user_data[IcData::VELOCITY][t1])*dir[2];
            v[k][j][i][4] =  a0*iod.ic.user_data[IcData::PRESSURE][t0] + a1*iod.ic.user_data[IcData::PRESSURE][t1];
          }
    } 
    else if (iod.ic.type == IcData::CYLINDRICAL) {
      print_error("Error: Cannot handle cylindrical i.c. at the moment.\n");
      exit_mpi();
    } 
    else if (iod.ic.type == IcData::SPHERICAL) {
      print_error("Error: Cannot handle spherical i.c. at the moment.\n");
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
void SpaceOperator::ApplyBoundaryConditions(SpaceVariable3D &V)
{
  Vec5D*** v = (Vec5D***) V.GetDataPointer();

  int NX, NY, NZ;
  V.GetGlobalSize(&NX, &NY, &NZ);
//  cout << "NX = " << NX << ", NY = " << NY << endl;

  //! Left boundary
  if(ii0==-1) { 
    switch (iod.mesh.bc_x0) {
      case MeshData::INLET :
        for(int k=k0; k<kmax; k++)
          for(int j=j0; j<jmax; j++) {
            v[k][j][ii0][0] = iod.bc.inlet.density;
            v[k][j][ii0][1] = iod.bc.inlet.velocity_x;
            v[k][j][ii0][2] = iod.bc.inlet.velocity_y;
            v[k][j][ii0][3] = iod.bc.inlet.velocity_z;
            v[k][j][ii0][4] = iod.bc.inlet.pressure;
          }
        break;
      case MeshData::OUTLET :
        for(int k=k0; k<kmax; k++)
          for(int j=j0; j<jmax; j++) {
            v[k][j][ii0][0] = iod.bc.outlet.density;
            v[k][j][ii0][1] = iod.bc.outlet.velocity_x;
            v[k][j][ii0][2] = iod.bc.outlet.velocity_y;
            v[k][j][ii0][3] = iod.bc.outlet.velocity_z;
            v[k][j][ii0][4] = iod.bc.outlet.pressure;
          }
        break; 
      case MeshData::WALL :
      case MeshData::SYMMETRY :
        for(int k=k0; k<kmax; k++)
          for(int j=j0; j<jmax; j++) {
            v[k][j][ii0][0] =      v[k][j][ii0+1][0];
            v[k][j][ii0][1] = -1.0*v[k][j][ii0+1][1];
            v[k][j][ii0][2] =      v[k][j][ii0+1][2];
            v[k][j][ii0][3] =      v[k][j][ii0+1][3]; 
            v[k][j][ii0][4] =      v[k][j][ii0+1][4];
          }
        break;
      default :
        print_error("Error: Boundary condition at x=x0 cannot be specified!\n");
        exit_mpi();
    }
  }

  //! Right boundary
  if(iimax==NX+1) { 
    switch (iod.mesh.bc_xmax) {
      case MeshData::INLET :
        for(int k=k0; k<kmax; k++)
          for(int j=j0; j<jmax; j++) {
            v[k][j][iimax-1][0] = iod.bc.inlet.density;
            v[k][j][iimax-1][1] = iod.bc.inlet.velocity_x;
            v[k][j][iimax-1][2] = iod.bc.inlet.velocity_y;
            v[k][j][iimax-1][3] = iod.bc.inlet.velocity_z;
            v[k][j][iimax-1][4] = iod.bc.inlet.pressure;
          }
        break;
      case MeshData::OUTLET :
        for(int k=k0; k<kmax; k++)
          for(int j=j0; j<jmax; j++) {
            v[k][j][iimax-1][0] = iod.bc.outlet.density;
            v[k][j][iimax-1][1] = iod.bc.outlet.velocity_x;
            v[k][j][iimax-1][2] = iod.bc.outlet.velocity_y;
            v[k][j][iimax-1][3] = iod.bc.outlet.velocity_z;
            v[k][j][iimax-1][4] = iod.bc.outlet.pressure;
          }
        break; 
      case MeshData::WALL :
      case MeshData::SYMMETRY :
        for(int k=k0; k<kmax; k++)
          for(int j=j0; j<jmax; j++) {
            v[k][j][iimax-1][0] =      v[k][j][iimax-2][0];
            v[k][j][iimax-1][1] = -1.0*v[k][j][iimax-2][1];
            v[k][j][iimax-1][2] =      v[k][j][iimax-2][2];
            v[k][j][iimax-1][3] =      v[k][j][iimax-2][3]; 
            v[k][j][iimax-1][4] =      v[k][j][iimax-2][4];
          }
        break;
      default :
        print_error("Error: Boundary condition at x=xmax cannot be specified!\n");
        exit_mpi();
    }
  }

  //! Bottom boundary
  if(jj0==-1) { 
    switch (iod.mesh.bc_y0) {
      case MeshData::INLET :
        for(int k=k0; k<kmax; k++)
          for(int i=i0; i<imax; i++) {
            v[k][jj0][i][0] = iod.bc.inlet.density;
            v[k][jj0][i][1] = iod.bc.inlet.velocity_x;
            v[k][jj0][i][2] = iod.bc.inlet.velocity_y;
            v[k][jj0][i][3] = iod.bc.inlet.velocity_z;
            v[k][jj0][i][4] = iod.bc.inlet.pressure;
          }
        break;
      case MeshData::OUTLET :
        for(int k=k0; k<kmax; k++)
          for(int i=i0; i<imax; i++) {
            v[k][jj0][i][0] = iod.bc.outlet.density;
            v[k][jj0][i][1] = iod.bc.outlet.velocity_x;
            v[k][jj0][i][2] = iod.bc.outlet.velocity_y;
            v[k][jj0][i][3] = iod.bc.outlet.velocity_z;
            v[k][jj0][i][4] = iod.bc.outlet.pressure;
          }
        break; 
      case MeshData::WALL :
      case MeshData::SYMMETRY :
        for(int k=k0; k<kmax; k++)
          for(int i=i0; i<imax; i++) {
            v[k][jj0][i][0] =      v[k][jj0+1][i][0];
            v[k][jj0][i][1] =      v[k][jj0+1][i][1];
            v[k][jj0][i][2] = -1.0*v[k][jj0+1][i][2];
            v[k][jj0][i][3] =      v[k][jj0+1][i][3]; 
            v[k][jj0][i][4] =      v[k][jj0+1][i][4];
          }
        break;
      default :
        print_error("Error: Boundary condition at y=y0 cannot be specified!\n");
        exit_mpi();
    }
  }

  //! Bottom boundary
  if(jjmax==NY+1) { 
    switch (iod.mesh.bc_ymax) {
      case MeshData::INLET :
        for(int k=k0; k<kmax; k++)
          for(int i=i0; i<imax; i++) {
            v[k][jjmax-1][i][0] = iod.bc.inlet.density;
            v[k][jjmax-1][i][1] = iod.bc.inlet.velocity_x;
            v[k][jjmax-1][i][2] = iod.bc.inlet.velocity_y;
            v[k][jjmax-1][i][3] = iod.bc.inlet.velocity_z;
            v[k][jjmax-1][i][4] = iod.bc.inlet.pressure;
          }
        break;
      case MeshData::OUTLET :
        for(int k=k0; k<kmax; k++)
          for(int i=i0; i<imax; i++) {
            v[k][jjmax-1][i][0] = iod.bc.outlet.density;
            v[k][jjmax-1][i][1] = iod.bc.outlet.velocity_x;
            v[k][jjmax-1][i][2] = iod.bc.outlet.velocity_y;
            v[k][jjmax-1][i][3] = iod.bc.outlet.velocity_z;
            v[k][jjmax-1][i][4] = iod.bc.outlet.pressure;
          }
        break; 
      case MeshData::WALL :
      case MeshData::SYMMETRY :
        for(int k=k0; k<kmax; k++)
          for(int i=i0; i<imax; i++) {
            v[k][jjmax-1][i][0] =      v[k][jjmax-2][i][0];
            v[k][jjmax-1][i][1] =      v[k][jjmax-2][i][1];
            v[k][jjmax-1][i][2] = -1.0*v[k][jjmax-2][i][2];
            v[k][jjmax-1][i][3] =      v[k][jjmax-2][i][3]; 
            v[k][jjmax-1][i][4] =      v[k][jjmax-2][i][4];
          }
        break;
      default :
        print_error("Error: Boundary condition at y=ymax cannot be specified!\n");
        exit_mpi();
    }
  }

  //! Back boundary (z min)
  if(kk0==-1) { 
    switch (iod.mesh.bc_z0) {
      case MeshData::INLET :
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {
            v[kk0][j][i][0] = iod.bc.inlet.density;
            v[kk0][j][i][1] = iod.bc.inlet.velocity_x;
            v[kk0][j][i][2] = iod.bc.inlet.velocity_y;
            v[kk0][j][i][3] = iod.bc.inlet.velocity_z;
            v[kk0][j][i][4] = iod.bc.inlet.pressure;
          }
        break;
      case MeshData::OUTLET :
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {
            v[kk0][j][i][0] = iod.bc.outlet.density;
            v[kk0][j][i][1] = iod.bc.outlet.velocity_x;
            v[kk0][j][i][2] = iod.bc.outlet.velocity_y;
            v[kk0][j][i][3] = iod.bc.outlet.velocity_z;
            v[kk0][j][i][4] = iod.bc.outlet.pressure;
          }
        break; 
      case MeshData::WALL :
      case MeshData::SYMMETRY :
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {
            v[kk0][j][i][0] =      v[kk0+1][j][i][0];
            v[kk0][j][i][1] =      v[kk0+1][j][i][1];
            v[kk0][j][i][2] =      v[kk0+1][j][i][2];
            v[kk0][j][i][3] = -1.0*v[kk0+1][j][i][3]; 
            v[kk0][j][i][4] =      v[kk0+1][j][i][4];
          }
        break;
      default :
        print_error("Error: Boundary condition at z=z0 cannot be specified!\n");
        exit_mpi();
    }
  }

  //! Front boundary (z max)
  if(kkmax==NZ+1) { 
    switch (iod.mesh.bc_zmax) {
      case MeshData::INLET :
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {
            v[kkmax-1][j][i][0] = iod.bc.inlet.density;
            v[kkmax-1][j][i][1] = iod.bc.inlet.velocity_x;
            v[kkmax-1][j][i][2] = iod.bc.inlet.velocity_y;
            v[kkmax-1][j][i][3] = iod.bc.inlet.velocity_z;
            v[kkmax-1][j][i][4] = iod.bc.inlet.pressure;
          }
        break;
      case MeshData::OUTLET :
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {
            v[kkmax-1][j][i][0] = iod.bc.outlet.density;
            v[kkmax-1][j][i][1] = iod.bc.outlet.velocity_x;
            v[kkmax-1][j][i][2] = iod.bc.outlet.velocity_y;
            v[kkmax-1][j][i][3] = iod.bc.outlet.velocity_z;
            v[kkmax-1][j][i][4] = iod.bc.outlet.pressure;
          }
        break; 
      case MeshData::WALL :
      case MeshData::SYMMETRY :
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {
            v[kkmax-1][j][i][0] =      v[kkmax-2][j][i][0];
            v[kkmax-1][j][i][1] =      v[kkmax-2][j][i][1];
            v[kkmax-1][j][i][2] =      v[kkmax-2][j][i][2];
            v[kkmax-1][j][i][3] = -1.0*v[kkmax-2][j][i][3]; 
            v[kkmax-1][j][i][4] =      v[kkmax-2][j][i][4];
          }
        break;
      default :
        print_error("Error: Boundary condition at z=zmax cannot be specified!\n");
        exit_mpi();
    }
  }

  V.RestoreDataPointerAndInsert();
}

//-----------------------------------------------------

void SpaceOperator::FindExtremeValuesOfFlowVariables(SpaceVariable3D &V,
                        double *Vmin, double *Vmax, double &cmin, double &cmax,
                        double &Machmax, double &char_speed_max,
                        double &dx_over_char_speed_min)
{
  Vec5D*** v    = (Vec5D***)V.GetDataPointer();
  Vec3D*** dxyz = (Vec3D***)delta_xyz.GetDataPointer();

  for(int i=0; i<5; i++) {
    Vmin[i] = DBL_MAX; //max. double precision number
    Vmax[i] = -DBL_MAX;
  }
  cmin = DBL_MAX;
  cmax = Machmax = char_speed_max = -DBL_MAX;
  dx_over_char_speed_min = DBL_MAX;

  // Loop through the real domain (excluding the ghost layer)
  double c, mach, lam_f, lam_g, lam_h;
  for(int k=k0; k<kmax; k++) {
    for(int j=j0; j<jmax; j++) {
      for(int i=i0; i<imax; i++) {

        for(int p=0; p<5; p++) {
          Vmin[p] = min(Vmin[p], v[k][j][i][p]);
          Vmax[p] = max(Vmax[p], v[k][j][i][p]);
        } 

        c = varFcn.ComputeSoundSpeed(v[k][j][i][0]/*rho*/, varFcn.GetInternalEnergyPerUnitMass(v[k][j][i][0],v[k][j][i][4])/*e*/);
        cmin = min(cmin, c);
        cmax = max(cmax, c);
        mach = varFcn.ComputeMachNumber(v[k][j][i]); 
        Machmax = max(Machmax, mach); 

        fluxFcn.EvaluateMaxEigenvalues(v[k][j][i], lam_f, lam_g, lam_h);
        char_speed_max = max(max(max(char_speed_max, lam_f), lam_g), lam_h);

        dx_over_char_speed_min = min(dx_over_char_speed_min, 
                                     min(dxyz[k][j][i][0]/lam_f, 
                                         min(dxyz[k][j][i][1]/lam_g, dxyz[k][j][i][2]/lam_h) ) );
      }
    }
  }

  MPI_Allreduce(Vmin, Vmin, 5, MPI_DOUBLE, MPI_MIN, comm);
  MPI_Allreduce(Vmax, Vmax, 5, MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(&cmin, &cmin, 1, MPI_DOUBLE, MPI_MIN, comm);
  MPI_Allreduce(&cmax, &cmax, 1, MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(&Machmax, &Machmax, 1, MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(&char_speed_max, &char_speed_max, 1, MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(&dx_over_char_speed_min, &dx_over_char_speed_min, 1, MPI_DOUBLE, MPI_MIN, comm);

  V.RestoreDataPointerToLocalVector(); 
  delta_xyz.RestoreDataPointerToLocalVector(); 
}

//-----------------------------------------------------

void SpaceOperator::ComputeTimeStepSize(SpaceVariable3D &V, double &dt, double &cfl)
{
  double Vmin[5], Vmax[5], cmin, cmax, Machmax, char_speed_max, dx_over_char_speed_min; 
  FindExtremeValuesOfFlowVariables(V, Vmin, Vmax, cmin, cmax, Machmax, char_speed_max, dx_over_char_speed_min);

  if(iod.output.verbose == OutputData::ON)
    print("  - Maximum values: rho = %e, p = %e, c = %e, Mach = %e, char. speed = %e.\n", 
          Vmax[0], Vmax[4], cmax, Machmax, char_speed_max);

  if(iod.ts.timestep > 0) {
    dt = iod.ts.timestep;
    cfl = dt/dx_over_char_speed_min;
  } else {//apply the CFL number
    cfl = iod.ts.cfl;      
    dt = cfl*dx_over_char_speed_min;
  }

}

//-----------------------------------------------------

void SpaceOperator::ComputeAdvectionFluxes(SpaceVariable3D &V, SpaceVariable3D &F)
{
  //------------------------------------
  // Reconstruction w/ slope limiters.
  //------------------------------------
  rec.Reconstruct(V, Vl, Vr, Vb, Vt, Vk, Vf);

  Vec5D*** v  = (Vec5D***) V.GetDataPointer();
  Vec5D*** vl = (Vec5D***) Vl.GetDataPointer();
  Vec5D*** vr = (Vec5D***) Vr.GetDataPointer();
  Vec5D*** vb = (Vec5D***) Vb.GetDataPointer();
  Vec5D*** vt = (Vec5D***) Vt.GetDataPointer();
  Vec5D*** vk = (Vec5D***) Vk.GetDataPointer();
  Vec5D*** vf = (Vec5D***) Vf.GetDataPointer();
  Vec5D*** f  = (Vec5D***) F.GetDataPointer();

  //------------------------------------
  // Clip pressure and density for the reconstructed state
  // Verify hyperbolicity (i.e. c^2 > 0).
  //------------------------------------
  int nClipped = 0;
  bool error = false;
  for(int k=kk0; k<kkmax; k++) {
    for(int j=jj0; j<jjmax; j++) {
      for(int i=ii0; i<iimax; i++) {
        nClipped += (int)varFcn.ClipDensityAndPressure(vl[k][j][i]);
        nClipped += (int)varFcn.ClipDensityAndPressure(vr[k][j][i]);
        nClipped += (int)varFcn.ClipDensityAndPressure(vb[k][j][i]);
        nClipped += (int)varFcn.ClipDensityAndPressure(vt[k][j][i]);
        nClipped += (int)varFcn.ClipDensityAndPressure(vk[k][j][i]);
        nClipped += (int)varFcn.ClipDensityAndPressure(vf[k][j][i]);
         
        error = varFcn.CheckState(vl[k][j][i]) || varFcn.CheckState(vr[k][j][i]) || 
                varFcn.CheckState(vb[k][j][i]) || varFcn.CheckState(vt[k][j][i]) ||
                varFcn.CheckState(vk[k][j][i]) || varFcn.CheckState(vf[k][j][i]);

        if(error) {
          print_error("Error: Reconstructed state at (%d,%d,%d) violates hyperbolicity.\n", i,j,k);
          exit_mpi();
        } 
      }
    }
  }
  MPI_Allreduce(&nClipped, &nClipped, 1, MPI_INT, MPI_SUM, comm);
  if(nClipped)
    print("Warning: Clipped pressure and/or density in %d reconstructed states.\n", nClipped);
  
  //------------------------------------
  // Compute fluxes
  //------------------------------------
  int nDOF = 5;
  Vec5D localflux;
  Vec3D*** dxyz = (Vec3D***)delta_xyz.GetDataPointer();

  // Initialize F to 0
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++) 
      for(int i=ii0; i<iimax; i++)
          f[k][j][i] = 0.0; //setting f[k][j][i][0] = ... = f[k][j][i][4] = 0.0;

  // Loop through the domain interior, and the right and top ghost layers. For each cell, calculate the
  // numerical flux across the left and lower cell boundaries/interfaces
  for(int k=k0; k<kkmax; k++) {
    for(int j=j0; j<jjmax; j++) {
      for(int i=i0; i<iimax; i++) {

        //calculate flux function F_{i-1/2,j,k}
        fluxFcn.ComputeNumericalFluxAtCellInterface(0/*F*/, vr[k][j][i-1]/*Vm*/, vl[k][j][i]/*Vp*/, localflux);
        localflux *= dxyz[k][j][i][1]*dxyz[k][j][i][2];
        f[k][j][i-1] += localflux;
        f[k][j][i]   -= localflux;  // the scheme is conservative 

        //calculate flux function G_{i,j-1/2,k}
        fluxFcn.ComputeNumericalFluxAtCellInterface(1/*G*/, vt[k][j-1][i]/*Vm*/, vb[k][j][i]/*Vp*/, localflux);
        localflux *= dxyz[k][j][i][0]*dxyz[k][j][i][2];
        f[k][j-1][i] += localflux;
        f[k][j][i]   -= localflux;  // the scheme is conservative 

        //calculate flux function H_{i,j,k-1/2}
        fluxFcn.ComputeNumericalFluxAtCellInterface(2/*H*/, vf[k-1][j][i]/*Vm*/, vk[k][j][i]/*Vp*/, localflux);
        localflux *= dxyz[k][j][i][0]*dxyz[k][j][i][1];
        f[k-1][j][i] += localflux;
        f[k][j][i]   -= localflux;  // the scheme is conservative 
      
      }
    }
  }
        
  //------------------------------------
  // Restore Spatial Variables
  //------------------------------------
  delta_xyz.RestoreDataPointerToLocalVector(); //no changes
  V.RestoreDataPointerToLocalVector(); 
  Vl.RestoreDataPointerToLocalVector(); 
  Vr.RestoreDataPointerToLocalVector(); 
  Vb.RestoreDataPointerToLocalVector(); 
  Vt.RestoreDataPointerToLocalVector(); 
  Vk.RestoreDataPointerToLocalVector(); 
  Vf.RestoreDataPointerToLocalVector(); 

  F.RestoreDataPointerToLocalVector(); //NOTE: although F has been updated, there is no need of 
                                       //      cross-subdomain communications. So, no need to 
                                       //      update the global vec.
}

//-----------------------------------------------------

void SpaceOperator::ComputeResidual(SpaceVariable3D &V, SpaceVariable3D &R)
{
  ComputeAdvectionFluxes(V,R);

  // -------------------------------------------------
  // multiply flux by -1, and divide by cell volume (for cells within the actual domain)
  // -------------------------------------------------
  Vec5D***    r = (Vec5D***) R.GetDataPointer();
  double*** vol = (double***)volume.GetDataPointer();

  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++) 
      for(int i=i0; i<imax; i++)
        r[k][j][i] /= -vol[k][j][i];

  // restore spatial variables
  R.RestoreDataPointerToLocalVector(); //NOTE: although R has been updated, there is no need of 
                                       //      cross-subdomain communications. So, no need to 
                                       //      update the global vec.
  volume.RestoreDataPointerToLocalVector();
}


//-----------------------------------------------------

int SpaceOperator::ClipDensityAndPressure(SpaceVariable3D &V, bool workOnGhost, bool checkState)
{

  Vec5D*** v = (Vec5D***) V.GetDataPointer();

  int myi0, myj0, myk0, myimax, myjmax, mykmax;
  if(workOnGhost)
    V.GetGhostedCornerIndices(&myi0, &myj0, &myk0, &myimax, &myjmax, &mykmax);
  else
    V.GetCornerIndices(&myi0, &myj0, &myk0, &myimax, &myjmax, &mykmax);

  int nClipped = 0;
  for(int k=myk0; k<mykmax; k++) {
    for(int j=myj0; j<myjmax; j++) {
      for(int i=myi0; i<myimax; i++) {
        nClipped += (int)varFcn.ClipDensityAndPressure(v[k][j][i]);

        if(checkState) {
          if(varFcn.CheckState(v[k][j][i])) {
            print_error("Error: State variables at (%d,%d,%d) violate hyperbolicity.\n", i,j,k);
            exit_mpi();
          }
        }
      }
    }
  }

  MPI_Allreduce(&nClipped, &nClipped, 1, MPI_INT, MPI_SUM, comm);
  if(nClipped)
    print("Warning: Clipped pressure and/or density in %d cells.\n", nClipped);

  V.RestoreDataPointerAndInsert();

  return nClipped;
}  

//-----------------------------------------------------











