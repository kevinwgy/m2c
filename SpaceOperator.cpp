#include <Utils.h>
#include <SpaceOperator.h>
#include <Vector3D.h>
#include <Vector5D.h>
#include <GeoTools.h>
#include <algorithm> //std::upper_bound
#include <cfloat> //DBL_MAX
using std::cout;
using std::endl;
using std::max;
using std::min;
using std::map;
using namespace GeoTools;

//-----------------------------------------------------

SpaceOperator::SpaceOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_,
                             vector<VarFcnBase*> &varFcn_, FluxFcnBase &fluxFcn_,
                             ExactRiemannSolverBase &riemann_,
                             vector<double> &x, vector<double> &y, vector<double> &z,
                             vector<double> &dx, vector<double> &dy, vector<double> &dz) 
  : comm(comm_), iod(iod_), varFcn(varFcn_), fluxFcn(fluxFcn_), riemann(riemann_),
    coordinates(comm_, &(dm_all_.ghosted1_3dof)),
    delta_xyz(comm_, &(dm_all_.ghosted1_3dof)),
    volume(comm_, &(dm_all_.ghosted1_1dof)),
    rec(comm_, dm_all_, iod_.schemes.ns.rec, coordinates, delta_xyz),
    Vl(comm_, &(dm_all_.ghosted1_5dof)),
    Vr(comm_, &(dm_all_.ghosted1_5dof)),
    Vb(comm_, &(dm_all_.ghosted1_5dof)),
    Vt(comm_, &(dm_all_.ghosted1_5dof)),
    Vk(comm_, &(dm_all_.ghosted1_5dof)),
    Vf(comm_, &(dm_all_.ghosted1_5dof))
{
  
  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);

  SetupMesh(x,y,z,dx,dy,dz);

  CreateGhostNodeLists(); //create ghost_nodes_inner and ghost_nodes_outer

  rec.Setup(&ghost_nodes_inner, &ghost_nodes_outer); //this function requires mesh info (dxyz)
  
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

void SpaceOperator::SetupMesh(vector<double> &x, vector<double> &y, vector<double> &z,
                              vector<double> &dx, vector<double> &dy, vector<double> &dz)
{
  //! Setup coordinates of cell centers and dx, dy, dz
  int NX, NY, NZ;
  coordinates.GetGlobalSize(&NX, &NY, &NZ);

  //! get array to edit
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
  Vec3D*** dxyz   = (Vec3D***)delta_xyz.GetDataPointer();

  //! Fill the actual subdomain, w/o ghost cells 
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        coords[k][j][i][0] = x[i];
        coords[k][j][i][1] = y[j];
        coords[k][j][i][2] = z[k];
        dxyz[k][j][i][0] = dx[i];
        dxyz[k][j][i][1] = dy[j];
        dxyz[k][j][i][2] = dz[k];
      } 

  //! restore array
  coordinates.RestoreDataPointerAndInsert(); //update localVec and globalVec;
  delta_xyz.RestoreDataPointerAndInsert(); //update localVec and globalVec;

  //! Populate the ghost cells (coordinates, dx, dy, dz)
  PopulateGhostBoundaryCoordinates();


  //(Obsolete. Can be used for debugging purpose though.) 
  //SetupMeshUniformRectangularDomain();


  //! Compute mesh information
  dxyz = (Vec3D***)delta_xyz.GetDataPointer();
  double*** vol  = (double***)volume.GetDataPointer();


  /** Calculate the volume/area of node-centered control volumes ("cells")
   *  Include ghost cells. 
   */
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++) {
        vol[k][j][i] /*volume of cv*/ = dxyz[k][j][i][0]*dxyz[k][j][i][1]*dxyz[k][j][i][2];
//        fprintf(stderr,"(%d,%d,%d), dx = %e, dy = %e, dz = %e, vol = %e.\n", i,j,k, dxyz[k][j][i][0], dxyz[k][j][i][1], dxyz[k][j][i][2], vol[k][j][i]);
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
    v1[p]    = v[kmax-1][jmax-1][imax-1][p] + dxyz[kmax-1][jmax-1][imax-1][p];
    dxyz0[p] = dxyz[k0][j0][i0][p];
    dxyz1[p] = dxyz[kmax-1][jmax-1][imax-1][p];
  }

  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++) {

        if(k!=-1 && k!=NZ && j!=-1 && j!=NY && i!=-1 && i!=NX)
          continue; //not in the ghost layer of the physical domain

        Vec3D& X  = v[k][j][i];
        Vec3D& dX = dxyz[k][j][i];

        bool xdone = false, ydone = false, zdone = false;

        if(i==-1) {
          X[0]  = v0[0];
          dX[0] = dxyz0[0];
          xdone = true;
        }

        if(i==NX) {
          X[0]  = v1[0];
          dX[0] = dxyz1[0];
          xdone = true;
        }

        if(j==-1) {
          X[1]  = v0[1];
          dX[1] = dxyz0[1];
          ydone = true;
        }

        if(j==NY) {
          X[1]  = v1[1];
          dX[1] = dxyz1[1];
          ydone = true;
        }

        if(k==-1) {
          X[2]  = v0[2];
          dX[2] = dxyz0[2];
          zdone = true;
        }

        if(k==NZ) {
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

void SpaceOperator::CreateGhostNodeLists()
{
  ghost_nodes_inner.clear();
  ghost_nodes_outer.clear();

  int NX, NY, NZ;
  coordinates.GetGlobalSize(&NX, &NY, &NZ);

  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();

  Int3 image;
  Vec3D proj(0.0), out_normal(0.0);
  MeshData::BcType bcType = MeshData::NONE;
  int counter;
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++) {

        if(k>=k0 && k<kmax && j>=j0 && j<jmax && i>=i0 && i<imax)
          continue; //interior of the subdomain

        if(coordinates.OutsidePhysicalDomain(i,j,k)) {//outside physical domain

          //determine the image, the projection point and boundary condition
          image      = 0;
          proj       = 0.0; 
          out_normal = 0.0;
          counter    = 0; 
          bcType     = MeshData::NONE;

          if(i<0)        {image[0] = -i-1;
                          proj[0]  = iod.mesh.x0;      out_normal[0] = -1.0; 
                          bcType   = iod.mesh.bc_x0;   counter++;}
          else if(i>=NX) {image[0] = NX+(i-NX)-1;
                          proj[0]  = iod.mesh.xmax;    out_normal[0] =  1.0; 
                          bcType   = iod.mesh.bc_xmax; counter++;}
          else           {image[0] = i;
                          proj[0]  = coords[k][j][i][0];}
                     

          if(j<0)        {image[1] = -j-1;
                          proj[1]  = iod.mesh.y0;      out_normal[1] = -1.0; 
                          bcType   = iod.mesh.bc_y0;   counter++;}
          else if(j>=NY) {image[1] = NY+(j-NY)-1;
                          proj[1]  = iod.mesh.ymax;    out_normal[1] =  1.0; 
                          bcType   = iod.mesh.bc_ymax; counter++;}
          else           {image[1] = j;
                          proj[1]  = coords[k][j][i][1];}
         

          if(k<0)        {image[2] = -k-1;
                          proj[2]  = iod.mesh.z0;      out_normal[2] = -1.0;
                          bcType   = iod.mesh.bc_z0;   counter++;}
          else if(k>=NZ) {image[2] = NZ+(k-NZ)-1;
                          proj[2]  = iod.mesh.zmax;    out_normal[2] =  1.0; 
                          bcType   = iod.mesh.bc_zmax; counter++;}
          else           {image[2] = k;
                          proj[2]  = coords[k][j][i][2];}
         
          out_normal /= out_normal.norm();

          assert(counter<=3 && counter>0);

          if(counter == 1)
            ghost_nodes_outer.push_back(GhostPoint(Int3(i,j,k), image, GhostPoint::FACE,
                                        proj, out_normal, (int)bcType));
          else if(counter == 2)
            ghost_nodes_outer.push_back(GhostPoint(Int3(i,j,k), image, GhostPoint::EDGE,
                                        proj, out_normal, 0));
          else
            ghost_nodes_outer.push_back(GhostPoint(Int3(i,j,k), image, GhostPoint::VERTEX,
                                        proj, out_normal, 0));

        } 
        else //inside physical domain
          ghost_nodes_inner.push_back(GhostPoint(i,j,k));

      }

  int nInner = ghost_nodes_inner.size();
  int nOuter = ghost_nodes_outer.size();
  MPI_Allreduce(MPI_IN_PLACE, &nInner, 1, MPI_INT, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, &nOuter, 1, MPI_INT, MPI_SUM, comm);
  print("  Number of ghost nodes inside computational domain (overlapping between subdomains): %d\n",
        nInner);
  print("  Number of ghost nodes outside computational domain: %d\n",
        nOuter);
  print("\n");


/*
  for(int i=0; i<ghost_nodes_outer.size(); i++)
    fprintf(stderr,"Ghost %d: (%d, %d, %d) | Image: (%d, %d, %d) | ProjType = %d | BcType = %d | Proj: (%e, %e, %e), (%e, %e, %e)\n", i, ghost_nodes_outer[i].ijk[0], ghost_nodes_outer[i].ijk[1], ghost_nodes_outer[i].ijk[2], ghost_nodes_outer[i].image_ijk[0], ghost_nodes_outer[i].image_ijk[1], ghost_nodes_outer[i].image_ijk[2], ghost_nodes_outer[i].type_projection, ghost_nodes_outer[i].bcType, ghost_nodes_outer[i].boundary_projection[0], ghost_nodes_outer[i].boundary_projection[1], ghost_nodes_outer[i].boundary_projection[2], ghost_nodes_outer[i].outward_normal[0], ghost_nodes_outer[i].outward_normal[1], ghost_nodes_outer[i].outward_normal[2]);
*/

  coordinates.RestoreDataPointerToLocalVector();
}

//-----------------------------------------------------

void SpaceOperator::ConservativeToPrimitive(SpaceVariable3D &U, SpaceVariable3D &ID, SpaceVariable3D &V,
                                            bool workOnGhost)
{
  Vec5D*** u = (Vec5D***) U.GetDataPointer();
  Vec5D*** v = (Vec5D***) V.GetDataPointer();
  double*** id = (double***) ID.GetDataPointer();

  int myi0, myj0, myk0, myimax, myjmax, mykmax;
  if(workOnGhost)
    U.GetGhostedCornerIndices(&myi0, &myj0, &myk0, &myimax, &myjmax, &mykmax);
  else
    U.GetCornerIndices(&myi0, &myj0, &myk0, &myimax, &myjmax, &mykmax);

  for(int k=myk0; k<mykmax; k++)
    for(int j=myj0; j<myjmax; j++)
      for(int i=myi0; i<myimax; i++)
        varFcn[id[k][j][i]]->ConservativeToPrimitive((double*)u[k][j][i], (double*)v[k][j][i]); 

  U.RestoreDataPointerToLocalVector(); //no changes made
  V.RestoreDataPointerAndInsert();
  ID.RestoreDataPointerToLocalVector(); //no changes made
}

//-----------------------------------------------------

void SpaceOperator::PrimitiveToConservative(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &U, 
                                            bool workOnGhost)
{
  Vec5D*** v = (Vec5D***) V.GetDataPointer();
  Vec5D*** u = (Vec5D***) U.GetDataPointer();
  double*** id = (double***) ID.GetDataPointer();

  int myi0, myj0, myk0, myimax, myjmax, mykmax;
  if(workOnGhost)
    U.GetGhostedCornerIndices(&myi0, &myj0, &myk0, &myimax, &myjmax, &mykmax);
  else
    U.GetCornerIndices(&myi0, &myj0, &myk0, &myimax, &myjmax, &mykmax);

  for(int k=myk0; k<mykmax; k++)
    for(int j=myj0; j<myjmax; j++)
      for(int i=myi0; i<myimax; i++)
        varFcn[id[k][j][i]]->PrimitiveToConservative((double*)v[k][j][i], (double*)u[k][j][i]); 

  V.RestoreDataPointerToLocalVector(); //no changes made
  U.RestoreDataPointerAndInsert();
  ID.RestoreDataPointerToLocalVector(); //no changes made
}

//-----------------------------------------------------

int SpaceOperator::ClipDensityAndPressure(SpaceVariable3D &V, SpaceVariable3D &ID, 
                                          bool workOnGhost, bool checkState)
{

  Vec5D*** v = (Vec5D***) V.GetDataPointer();
  double*** id = (double***) ID.GetDataPointer();
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();

  int myi0, myj0, myk0, myimax, myjmax, mykmax;
  if(workOnGhost)
    V.GetGhostedCornerIndices(&myi0, &myj0, &myk0, &myimax, &myjmax, &mykmax);
  else
    V.GetCornerIndices(&myi0, &myj0, &myk0, &myimax, &myjmax, &mykmax);

  int nClipped = 0;
  for(int k=myk0; k<mykmax; k++) {
    for(int j=myj0; j<myjmax; j++) {
      for(int i=myi0; i<myimax; i++) {

        nClipped += (int)varFcn[id[k][j][i]]->ClipDensityAndPressure(v[k][j][i]);

        if(checkState) {
          if(varFcn[id[k][j][i]]->CheckState(v[k][j][i])) {
            fprintf(stderr, "\033[0;31m*** Error: State variables at (%e,%e,%e) violate hyperbolicity." 
                    " matid = %d.\n\033[0m", coords[k][j][i][0],coords[k][j][i][1],coords[k][j][i][2], (int)id[k][j][i]);
            fprintf(stderr, "\033[0;31mv[%d(i),%d(j),%d(k)] = [%e, %e, %e, %e, %e]\n\033[0m", 
                    i,j,k, v[k][j][i][0], v[k][j][i][1], v[k][j][i][2], v[k][j][i][3], v[k][j][i][4]);
            exit(-1);
          }
        }
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &nClipped, 1, MPI_INT, MPI_SUM, comm);
  if(nClipped)
    print("Warning: Clipped pressure and/or density in %d cells.\n", nClipped);


  ID.RestoreDataPointerToLocalVector(); //no changes made
  V.RestoreDataPointerAndInsert();

  coordinates.RestoreDataPointerToLocalVector(); //no changes

  return nClipped;
}  

//-----------------------------------------------------

//apply IC within the real domain
void SpaceOperator::SetInitialCondition(SpaceVariable3D &V, SpaceVariable3D &ID) 
{
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();

  Vec5D*** v = (Vec5D***) V.GetDataPointer();
  double*** id = (double***) ID.GetDataPointer();

  //! 1. apply the inlet (i.e. farfield) state
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++) {
        v[k][j][i][0] = iod.bc.inlet.density;
        v[k][j][i][1] = iod.bc.inlet.velocity_x;
        v[k][j][i][2] = iod.bc.inlet.velocity_y;
        v[k][j][i][3] = iod.bc.inlet.velocity_z;
        v[k][j][i][4] = iod.bc.inlet.pressure;
        id[k][j][i]   = iod.bc.inlet.materialid;
      }



  //! 2. apply user-specified function
  if(iod.ic.type != IcData::NONE) {

    //! Get coordinates
    Vec3D    x0(iod.ic.x0[0], iod.ic.x0[1], iod.ic.x0[2]); 

    if (iod.ic.type == IcData::PLANAR || iod.ic.type == IcData::CYLINDRICAL) {

      if(iod.ic.type == IcData::PLANAR)
        print("- Applying initial condition specified in %s (planar).\n", 
              iod.ic.user_specified_ic);
      else
        print("- Applying initial condition specified in %s (with cylindrical symmetry).\n", 
              iod.ic.user_specified_ic);
 
      Vec3D dir(iod.ic.dir[0], iod.ic.dir[1], iod.ic.dir[2]); 
      dir /= dir.norm();

      double x = 0.0;
      int n = iod.ic.user_data[IcData::COORDINATE].size(); //!< number of data points provided by user (in the axial dir)

      double r = 0.0;
      int nrad = iod.ic.user_data2[IcData::COORDINATE].size(); //!< number of data points in the radial dir.

      int t0, t1;    
      double a0, a1;

      for(int k=k0; k<kmax; k++)
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {

            x = (coords[k][j][i] - x0)*dir; //!< projection onto the 1D axis
//            cout << "coords: " << coords[j][i][0] << ", " << coords[j][i][1] << "; x = " << x << endl;
            if(x<0 || x>iod.ic.user_data[IcData::COORDINATE][n-1])
              continue;

            if(nrad>0) {
              r = (coords[k][j][i] - x0 - x*dir).norm();
              if(r>iod.ic.user_data2[IcData::COORDINATE][nrad-1])
                continue;
            }
 
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

            //! specify i.c. on node (cell center)
            if(iod.ic.specified[IcData::DENSITY])
              v[k][j][i][0] = a0*iod.ic.user_data[IcData::DENSITY][t0]  
                            + a1*iod.ic.user_data[IcData::DENSITY][t1];
            if(iod.ic.specified[IcData::VELOCITY]) {
              v[k][j][i][1] = (a0*iod.ic.user_data[IcData::VELOCITY][t0] 
                            +  a1*iod.ic.user_data[IcData::VELOCITY][t1])*dir[0];
              v[k][j][i][2] = (a0*iod.ic.user_data[IcData::VELOCITY][t0]
                            +  a1*iod.ic.user_data[IcData::VELOCITY][t1])*dir[1];
              v[k][j][i][3] = (a0*iod.ic.user_data[IcData::VELOCITY][t0] 
                            +  a1*iod.ic.user_data[IcData::VELOCITY][t1])*dir[2];
            }
            if(iod.ic.specified[IcData::PRESSURE]) 
              v[k][j][i][4] = a0*iod.ic.user_data[IcData::PRESSURE][t0]
                            + a1*iod.ic.user_data[IcData::PRESSURE][t1];
            if(iod.ic.specified[IcData::MATERIALID])
              id[k][j][i]   = std::round(a0*iod.ic.user_data[IcData::MATERIALID][t0] 
                                       + a1*iod.ic.user_data[IcData::MATERIALID][t1]);

            //! apply radial variation (if provided by user)
            if(nrad>0) { 
 
              //! Find the first radial coordinate greater than r
              auto upper_it = std::upper_bound(iod.ic.user_data2[IcData::COORDINATE].begin(),
                                               iod.ic.user_data2[IcData::COORDINATE].end(),
                                               r); 
              t1 = (int)(upper_it - iod.ic.user_data2[IcData::COORDINATE].begin());

              if(t1==0) // exactly the first node
                t1 = 1;

              t0 = t1 - 1;

              //! calculate interpolation weights
              a0 = (iod.ic.user_data2[IcData::COORDINATE][t1] - r) /
                   (iod.ic.user_data2[IcData::COORDINATE][t1] - iod.ic.user_data2[IcData::COORDINATE][t0]);
              a1 = 1.0 - a0;

              if(iod.ic.specified[IcData::DENSITY])
                v[k][j][i][0] *= a0*iod.ic.user_data2[IcData::DENSITY][t0] 
                               + a1*iod.ic.user_data2[IcData::DENSITY][t1];
              if(iod.ic.specified[IcData::VELOCITY])
                for(int p=1; p<=3; p++)
                  v[k][j][i][p] *= a0*iod.ic.user_data2[IcData::VELOCITY][t0] 
                                 + a1*iod.ic.user_data2[IcData::VELOCITY][t1];
              if(iod.ic.specified[IcData::PRESSURE])
                v[k][j][i][4] *= a0*iod.ic.user_data2[IcData::PRESSURE][t0]
                               + a1*iod.ic.user_data2[IcData::PRESSURE][t1];
            }
          }

    } 

    else if (iod.ic.type == IcData::SPHERICAL) {

      print("- Applying initial condition in %s (with spherical symmetry).\n", 
            iod.ic.user_specified_ic);
 
      double x;
      int n = iod.ic.user_data[IcData::COORDINATE].size(); //!< number of data points provided by user
      Vec3D dir;

      int t0, t1;    
      double a0, a1;

      for(int k=k0; k<kmax; k++)
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {

            dir = coords[k][j][i] - x0;
            x = dir.norm();
            dir /= x;
   
//            cout << "coords: " << coords[j][i][0] << ", " << coords[j][i][1] << "; x = " << x << endl;
            if(x>iod.ic.user_data[IcData::COORDINATE][n-1])
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

            //! specify i.c. on node (cell center)
            if(iod.ic.specified[IcData::DENSITY])
              v[k][j][i][0] = a0*iod.ic.user_data[IcData::DENSITY][t0]
                            + a1*iod.ic.user_data[IcData::DENSITY][t1];
            if(iod.ic.specified[IcData::VELOCITY]) {
              v[k][j][i][1] = (a0*iod.ic.user_data[IcData::VELOCITY][t0]
                            +  a1*iod.ic.user_data[IcData::VELOCITY][t1])*dir[0];
              v[k][j][i][2] = (a0*iod.ic.user_data[IcData::VELOCITY][t0] 
                            +  a1*iod.ic.user_data[IcData::VELOCITY][t1])*dir[1];
              v[k][j][i][3] = (a0*iod.ic.user_data[IcData::VELOCITY][t0]
                            +  a1*iod.ic.user_data[IcData::VELOCITY][t1])*dir[2];
            }
            if(iod.ic.specified[IcData::PRESSURE])
              v[k][j][i][4] = a0*iod.ic.user_data[IcData::PRESSURE][t0]
                            + a1*iod.ic.user_data[IcData::PRESSURE][t1];
            if(iod.ic.specified[IcData::MATERIALID])
              id[k][j][i]   = std::round(a0*iod.ic.user_data[IcData::MATERIALID][t0] 
                                       + a1*iod.ic.user_data[IcData::MATERIALID][t1]);
          }

    }
  }



  //! 3. apply i.c. based on geometric objects (planes, cylinder-cones, spheres)
  MultiInitialConditionsData &ic(iod.ic.multiInitialConditions);

  // planes
  for(auto it=ic.planeMap.dataMap.begin(); it!=ic.planeMap.dataMap.end(); it++) {

    print("- Applying initial condition on one side of a plane (material id: %d).\n", 
          it->second->initialConditions.materialid);
    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
    Vec3D dir(it->second->nx, it->second->ny, it->second->nz);
    dir /= dir.norm();
    double dist;

    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++) {
          dist = (coords[k][j][i]-x0)*dir;
          if(dist>0) {
            v[k][j][i][0] = it->second->initialConditions.density;
            v[k][j][i][1] = it->second->initialConditions.velocity_x;
            v[k][j][i][2] = it->second->initialConditions.velocity_y;
            v[k][j][i][3] = it->second->initialConditions.velocity_z;
            v[k][j][i][4] = it->second->initialConditions.pressure;
            id[k][j][i]   = it->second->initialConditions.materialid;
          }
        }
  }

  // cylinder-cone
  for(auto it=ic.cylinderconeMap.dataMap.begin(); it!=ic.cylinderconeMap.dataMap.end(); it++) {

    print("- Applying initial condition within a cylinder-cone (material id: %d).\n",
          it->second->initialConditions.materialid);
    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
    Vec3D dir(it->second->nx, it->second->ny, it->second->nz);
    dir /= dir.norm();

    double L = it->second->L; //cylinder height
    double R = it->second->r; //cylinder radius
    double tan_alpha = tan(it->second->opening_angle_degrees/180.0*acos(-1.0));//opening angle
    double Hmax = R/tan_alpha;
    double H = min(it->second->cone_height, Hmax); //cone's height

    double x, r;
    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++) {
          x = (coords[k][j][i]-x0)*dir;
          r = (coords[k][j][i] - x0 - x*dir).norm();
          if( (x>0 && x<L && r<R) || (x>=L && x<L+H && r<(L+Hmax-x)*tan_alpha) ) {//inside
            v[k][j][i][0] = it->second->initialConditions.density;
            v[k][j][i][1] = it->second->initialConditions.velocity_x;
            v[k][j][i][2] = it->second->initialConditions.velocity_y;
            v[k][j][i][3] = it->second->initialConditions.velocity_z;
            v[k][j][i][4] = it->second->initialConditions.pressure;
            id[k][j][i]   = it->second->initialConditions.materialid;
          }
        }
  }


  // spheres
  for(auto it=ic.sphereMap.dataMap.begin(); it!=ic.sphereMap.dataMap.end(); it++) {

    print("- Applying initial condition within a sphere (material id: %d).\n",
          it->second->initialConditions.materialid);
    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
    double dist;
    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++) {
          dist = (coords[k][j][i]-x0).norm() - it->second->radius;
          if (dist<0) {
            v[k][j][i][0] = it->second->initialConditions.density;
            v[k][j][i][1] = it->second->initialConditions.velocity_x;
            v[k][j][i][2] = it->second->initialConditions.velocity_y;
            v[k][j][i][3] = it->second->initialConditions.velocity_z;
            v[k][j][i][4] = it->second->initialConditions.pressure;
            id[k][j][i]   = it->second->initialConditions.materialid;
          }
        }
  }



  V.RestoreDataPointerAndInsert();
  ID.RestoreDataPointerAndInsert();
  coordinates.RestoreDataPointerToLocalVector(); //!< data was not changed.

  //! Apply boundary condition to populate ghost nodes (no need to do this for ID)
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
        print_error("*** Error: Boundary condition at x=x0 cannot be specified!\n");
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
        print_error("*** Error: Boundary condition at x=xmax cannot be specified!\n");
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
        print_error("*** Error: Boundary condition at y=y0 cannot be specified!\n");
        exit_mpi();
    }
  }

  //! Top boundary
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
        print_error("*** Error: Boundary condition at y=ymax cannot be specified!\n");
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
        print_error("*** Error: Boundary condition at z=z0 cannot be specified!\n");
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
        print_error("*** Error: Boundary condition at z=zmax cannot be specified!\n");
        exit_mpi();
    }
  }

  ApplyBoundaryConditionsGeometricEntities(v);

  V.RestoreDataPointerAndInsert();
}

//-----------------------------------------------------

void
SpaceOperator::ApplyBoundaryConditionsGeometricEntities(Vec5D*** v)
{

  map<int, DiskData* >& disks(iod.bc.multiBoundaryConditions.diskMap.dataMap);
  map<int, RectangleData* >& rectangles(iod.bc.multiBoundaryConditions.rectangleMap.dataMap);

  if(!disks.size() && !rectangles.size())
    return;

  int NX, NY, NZ;
  coordinates.GetGlobalSize(&NX, &NY, &NZ);  

  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
  

  if(ii0==-1) { 
    if (iod.mesh.bc_x0 == MeshData::INLET || iod.mesh.bc_x0 == MeshData::OUTLET || iod.mesh.bc_x0 == MeshData::WALL) {

      vector<DiskData* > mydisks;
      for(auto it=disks.begin(); it!=disks.end(); it++)
        if(it->second->cen_x == iod.mesh.x0) {
          Vec3D n(it->second->normal_x, it->second->normal_y, it->second->normal_z);
          if(fabs(n[0])/n.norm()>1-1e-8)
            mydisks.push_back(it->second); 
        }
      vector<RectangleData* > myrects;
      for(auto it=rectangles.begin(); it!=rectangles.end(); it++)       
        if(it->second->cen_x == iod.mesh.x0) {
          Vec3D n(it->second->normal_x, it->second->normal_y, it->second->normal_z);
          if(fabs(n[0])/n.norm()>1-1e-8)
            myrects.push_back(it->second); 
         }
      if(mydisks.size() || myrects.size()) {

        for(int k=k0; k<kmax; k++)
          for(int j=j0; j<jmax; j++) {

            for(int p=0; p<mydisks.size(); p++) {
              if(IsPointInDisk(coords[k][j][ii0][1], coords[k][j][ii0][2],
                               mydisks[p]->cen_y, mydisks[p]->cen_z, mydisks[p]->radius)){
                v[k][j][ii0][0] = mydisks[p]->state.density;
                v[k][j][ii0][1] = mydisks[p]->state.velocity_x;
                v[k][j][ii0][2] = mydisks[p]->state.velocity_y;
                v[k][j][ii0][3] = mydisks[p]->state.velocity_z;
                v[k][j][ii0][4] = mydisks[p]->state.pressure;
              }
            }

            for(int p=0; p<myrects.size(); p++) {
              if(IsPointInRectangle(coords[k][j][ii0][1], coords[k][j][ii0][2],
                                    myrects[p]->cen_y, myrects[p]->cen_z, myrects[p]->a, myrects[p]->b)){
                v[k][j][ii0][0] = myrects[p]->state.density;
                v[k][j][ii0][1] = myrects[p]->state.velocity_x;
                v[k][j][ii0][2] = myrects[p]->state.velocity_y;
                v[k][j][ii0][3] = myrects[p]->state.velocity_z;
                v[k][j][ii0][4] = myrects[p]->state.pressure;
              }
            }

          }
      }

    }
  }

  if(iimax==NX+1) { 
    if (iod.mesh.bc_xmax == MeshData::INLET || iod.mesh.bc_xmax == MeshData::OUTLET || iod.mesh.bc_xmax == MeshData::WALL) {

      vector<DiskData* > mydisks;
      for(auto it=disks.begin(); it!=disks.end(); it++)
        if(it->second->cen_x == iod.mesh.xmax) {
          Vec3D n(it->second->normal_x, it->second->normal_y, it->second->normal_z);
          if(fabs(n[0])/n.norm()>1-1e-8)
            mydisks.push_back(it->second); 
        } 
      vector<RectangleData* > myrects;
      for(auto it=rectangles.begin(); it!=rectangles.end(); it++)       
        if(it->second->cen_x == iod.mesh.xmax) {
          Vec3D n(it->second->normal_x, it->second->normal_y, it->second->normal_z);
          if(fabs(n[0])/n.norm()>1-1e-8)
            myrects.push_back(it->second); 
         }

      if(mydisks.size() || myrects.size()) {

        for(int k=k0; k<kmax; k++)
          for(int j=j0; j<jmax; j++) {

            for(int p=0; p<mydisks.size(); p++) {
              if(IsPointInDisk(coords[k][j][iimax-1][1], coords[k][j][iimax-1][2],
                               mydisks[p]->cen_y, mydisks[p]->cen_z, mydisks[p]->radius)){
                v[k][j][iimax-1][0] = mydisks[p]->state.density;
                v[k][j][iimax-1][1] = mydisks[p]->state.velocity_x;
                v[k][j][iimax-1][2] = mydisks[p]->state.velocity_y;
                v[k][j][iimax-1][3] = mydisks[p]->state.velocity_z;
                v[k][j][iimax-1][4] = mydisks[p]->state.pressure;
              }
            }

            for(int p=0; p<myrects.size(); p++) {
              if(IsPointInRectangle(coords[k][j][iimax-1][1], coords[k][j][iimax-1][2],
                                    myrects[p]->cen_y, myrects[p]->cen_z, myrects[p]->a, myrects[p]->b)){
                v[k][j][iimax-1][0] = myrects[p]->state.density;
                v[k][j][iimax-1][1] = myrects[p]->state.velocity_x;
                v[k][j][iimax-1][2] = myrects[p]->state.velocity_y;
                v[k][j][iimax-1][3] = myrects[p]->state.velocity_z;
                v[k][j][iimax-1][4] = myrects[p]->state.pressure;
              }
            }

          }
      }

    }
  }

  
  if(jj0==-1) { 
    if (iod.mesh.bc_y0 == MeshData::INLET || iod.mesh.bc_y0 == MeshData::OUTLET || iod.mesh.bc_y0 == MeshData::WALL) {

      vector<DiskData* > mydisks;
      for(auto it=disks.begin(); it!=disks.end(); it++)
        if(it->second->cen_y == iod.mesh.y0) {
          Vec3D n(it->second->normal_x, it->second->normal_y, it->second->normal_z);
          if(fabs(n[1])/n.norm()>1-1e-8)
            mydisks.push_back(it->second); 
        }
      vector<RectangleData* > myrects;
      for(auto it=rectangles.begin(); it!=rectangles.end(); it++)       
        if(it->second->cen_y == iod.mesh.y0) {
          Vec3D n(it->second->normal_x, it->second->normal_y, it->second->normal_z);
          if(fabs(n[1])/n.norm()>1-1e-8)
            myrects.push_back(it->second); 
        }
      if(mydisks.size() || myrects.size()) {

        for(int k=k0; k<kmax; k++)
          for(int i=i0; i<imax; i++) {

            for(int p=0; p<mydisks.size(); p++) {
              if(IsPointInDisk(coords[k][jj0][i][2], coords[k][jj0][i][0],
                               mydisks[p]->cen_z, mydisks[p]->cen_x, mydisks[p]->radius)){
                v[k][jj0][i][0] = mydisks[p]->state.density;
                v[k][jj0][i][1] = mydisks[p]->state.velocity_x;
                v[k][jj0][i][2] = mydisks[p]->state.velocity_y;
                v[k][jj0][i][3] = mydisks[p]->state.velocity_z;
                v[k][jj0][i][4] = mydisks[p]->state.pressure;
              }
            }

            for(int p=0; p<myrects.size(); p++) {
              if(IsPointInRectangle(coords[k][jj0][i][2], coords[k][jj0][i][0],
                                    myrects[p]->cen_z, myrects[p]->cen_x, myrects[p]->a, myrects[p]->b)){
                v[k][jj0][i][0] = myrects[p]->state.density;
                v[k][jj0][i][1] = myrects[p]->state.velocity_x;
                v[k][jj0][i][2] = myrects[p]->state.velocity_y;
                v[k][jj0][i][3] = myrects[p]->state.velocity_z;
                v[k][jj0][i][4] = myrects[p]->state.pressure;
              }
            }

          }
      }

    }
  }

  if(jjmax==NY+1) { 
    if (iod.mesh.bc_ymax == MeshData::INLET || iod.mesh.bc_ymax == MeshData::OUTLET || iod.mesh.bc_ymax == MeshData::WALL) {

      vector<DiskData* > mydisks;
      for(auto it=disks.begin(); it!=disks.end(); it++)
        if(it->second->cen_y == iod.mesh.ymax) {
          Vec3D n(it->second->normal_x, it->second->normal_y, it->second->normal_z);
          if(fabs(n[1])/n.norm()>1-1e-8)
            mydisks.push_back(it->second); 
        }
      vector<RectangleData* > myrects;
      for(auto it=rectangles.begin(); it!=rectangles.end(); it++)       
        if(it->second->cen_y == iod.mesh.ymax) {
          Vec3D n(it->second->normal_x, it->second->normal_y, it->second->normal_z);
          if(fabs(n[1])/n.norm()>1-1e-8)
            myrects.push_back(it->second); 
        }

      if(mydisks.size() || myrects.size()) {

        for(int k=k0; k<kmax; k++)
          for(int i=i0; i<imax; i++) {

            for(int p=0; p<mydisks.size(); p++) {
              if(IsPointInDisk(coords[k][jjmax-1][i][2], coords[k][jjmax-1][i][0],
                               mydisks[p]->cen_z, mydisks[p]->cen_x, mydisks[p]->radius)){
                v[k][jjmax-1][i][0] = mydisks[p]->state.density;
                v[k][jjmax-1][i][1] = mydisks[p]->state.velocity_x;
                v[k][jjmax-1][i][2] = mydisks[p]->state.velocity_y;
                v[k][jjmax-1][i][3] = mydisks[p]->state.velocity_z;
                v[k][jjmax-1][i][4] = mydisks[p]->state.pressure;
              }
            }

            for(int p=0; p<myrects.size(); p++) {
              if(IsPointInRectangle(coords[k][jjmax-1][i][2], coords[k][jjmax-1][i][0],
                                    myrects[p]->cen_z, myrects[p]->cen_x, myrects[p]->a, myrects[p]->b)){
                v[k][jjmax-1][i][0] = myrects[p]->state.density;
                v[k][jjmax-1][i][1] = myrects[p]->state.velocity_x;
                v[k][jjmax-1][i][2] = myrects[p]->state.velocity_y;
                v[k][jjmax-1][i][3] = myrects[p]->state.velocity_z;
                v[k][jjmax-1][i][4] = myrects[p]->state.pressure;
              }
            }

          }
      }

    }
  }

  
  if(kk0==-1) { 
    if (iod.mesh.bc_z0 == MeshData::INLET || iod.mesh.bc_z0 == MeshData::OUTLET || iod.mesh.bc_z0 == MeshData::WALL) {

      vector<DiskData* > mydisks;
      for(auto it=disks.begin(); it!=disks.end(); it++)
        if(it->second->cen_z == iod.mesh.z0) {
          Vec3D n(it->second->normal_x, it->second->normal_y, it->second->normal_z);
          if(fabs(n[2])/n.norm()>1-1e-8)
            mydisks.push_back(it->second); 
        }
      vector<RectangleData* > myrects;
      for(auto it=rectangles.begin(); it!=rectangles.end(); it++)       
        if(it->second->cen_z == iod.mesh.z0) {
          Vec3D n(it->second->normal_x, it->second->normal_y, it->second->normal_z);
          if(fabs(n[2])/n.norm()>1-1e-8)
            myrects.push_back(it->second); 
        }

      if(mydisks.size() || myrects.size()) {

        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {

            for(int p=0; p<mydisks.size(); p++) {
              if(IsPointInDisk(coords[kk0][j][i][0], coords[kk0][j][i][1],
                               mydisks[p]->cen_x, mydisks[p]->cen_y, mydisks[p]->radius)){
                v[kk0][j][i][0] = mydisks[p]->state.density;
                v[kk0][j][i][1] = mydisks[p]->state.velocity_x;
                v[kk0][j][i][2] = mydisks[p]->state.velocity_y;
                v[kk0][j][i][3] = mydisks[p]->state.velocity_z;
                v[kk0][j][i][4] = mydisks[p]->state.pressure;
              }
            }

            for(int p=0; p<myrects.size(); p++) {
              if(IsPointInRectangle(coords[kk0][j][i][0], coords[kk0][j][i][1],
                                    myrects[p]->cen_x, myrects[p]->cen_y, myrects[p]->a, myrects[p]->b)){
                v[kk0][j][i][0] = myrects[p]->state.density;
                v[kk0][j][i][1] = myrects[p]->state.velocity_x;
                v[kk0][j][i][2] = myrects[p]->state.velocity_y;
                v[kk0][j][i][3] = myrects[p]->state.velocity_z;
                v[kk0][j][i][4] = myrects[p]->state.pressure;
              }
            }

          }
      }

    }
  }

  if(kkmax==NZ+1) { 
    if (iod.mesh.bc_zmax == MeshData::INLET || iod.mesh.bc_zmax == MeshData::OUTLET || iod.mesh.bc_zmax == MeshData::WALL) {

      vector<DiskData* > mydisks;
      for(auto it=disks.begin(); it!=disks.end(); it++)
        if(it->second->cen_z == iod.mesh.zmax) {
          Vec3D n(it->second->normal_x, it->second->normal_y, it->second->normal_z);
          if(fabs(n[2])/n.norm()>1-1e-8)
            mydisks.push_back(it->second); 
        }

      vector<RectangleData* > myrects;
      for(auto it=rectangles.begin(); it!=rectangles.end(); it++)       
        if(it->second->cen_z == iod.mesh.zmax) {
          Vec3D n(it->second->normal_x, it->second->normal_y, it->second->normal_z);
          if(fabs(n[2])/n.norm()>1-1e-8)
            myrects.push_back(it->second); 
        }

      if(mydisks.size() || myrects.size()) {

        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {

            for(int p=0; p<mydisks.size(); p++) {
              if(IsPointInDisk(coords[kkmax-1][j][i][0], coords[kkmax-1][j][i][1],
                               mydisks[p]->cen_x, mydisks[p]->cen_y, mydisks[p]->radius)){
                v[kkmax-1][j][i][0] = mydisks[p]->state.density;
                v[kkmax-1][j][i][1] = mydisks[p]->state.velocity_x;
                v[kkmax-1][j][i][2] = mydisks[p]->state.velocity_y;
                v[kkmax-1][j][i][3] = mydisks[p]->state.velocity_z;
                v[kkmax-1][j][i][4] = mydisks[p]->state.pressure;
              }
            }

            for(int p=0; p<myrects.size(); p++) {
              if(IsPointInRectangle(coords[kkmax-1][j][i][0], coords[kkmax-1][j][i][1],
                                    myrects[p]->cen_x, myrects[p]->cen_y, myrects[p]->a, myrects[p]->b)){
                v[kkmax-1][j][i][0] = myrects[p]->state.density;
                v[kkmax-1][j][i][1] = myrects[p]->state.velocity_x;
                v[kkmax-1][j][i][2] = myrects[p]->state.velocity_y;
                v[kkmax-1][j][i][3] = myrects[p]->state.velocity_z;
                v[kkmax-1][j][i][4] = myrects[p]->state.pressure;
              }
            }

          }
      }

    }
  }


  coordinates.RestoreDataPointerToLocalVector(); //no changes
}

//-----------------------------------------------------

void SpaceOperator::FindExtremeValuesOfFlowVariables(SpaceVariable3D &V, SpaceVariable3D &ID,
                        double *Vmin, double *Vmax, double &cmin, double &cmax,
                        double &Machmax, double &char_speed_max,
                        double &dx_over_char_speed_min)
{
  Vec5D*** v    = (Vec5D***)V.GetDataPointer();
  Vec3D*** dxyz = (Vec3D***)delta_xyz.GetDataPointer();
  double*** id = (double***) ID.GetDataPointer();

  for(int i=0; i<5; i++) {
    Vmin[i] = DBL_MAX; //max. double precision number
    Vmax[i] = -DBL_MAX;
  }
  cmin = DBL_MAX;
  cmax = Machmax = char_speed_max = -DBL_MAX;
  dx_over_char_speed_min = DBL_MAX;

  // Loop through the real domain (excluding the ghost layer)
  double c, mach, lam_f, lam_g, lam_h;
  int myid;
  for(int k=k0; k<kmax; k++) {
    for(int j=j0; j<jmax; j++) {
      for(int i=i0; i<imax; i++) {

        for(int p=0; p<5; p++) {
          Vmin[p] = min(Vmin[p], v[k][j][i][p]);
          Vmax[p] = max(Vmax[p], v[k][j][i][p]);
        } 

        myid = id[k][j][i];

        c = varFcn[myid]->ComputeSoundSpeedSquare(v[k][j][i][0]/*rho*/, 
                            varFcn[myid]->GetInternalEnergyPerUnitMass(v[k][j][i][0],v[k][j][i][4])/*e*/);

        if(c<0) {
          fprintf(stderr,"*** Error: c^2 (square of sound speed) = %e in SpaceOperator. V = %e, %e, %e, %e, %e, ID = %d.\n",
                  c, v[k][j][i][0], v[k][j][i][1], v[k][j][i][2], v[k][j][i][3], v[k][j][i][4], myid);
          exit_mpi();
        } else
          c = sqrt(c);

        cmin = min(cmin, c);
        cmax = max(cmax, c);
        mach = varFcn[myid]->ComputeMachNumber(v[k][j][i]); 
        Machmax = max(Machmax, mach); 

        fluxFcn.EvaluateMaxEigenvalues(v[k][j][i], myid, lam_f, lam_g, lam_h);
        char_speed_max = max(max(max(char_speed_max, lam_f), lam_g), lam_h);

        dx_over_char_speed_min = min(dx_over_char_speed_min, 
                                     min(dxyz[k][j][i][0]/lam_f, 
                                         min(dxyz[k][j][i][1]/lam_g, dxyz[k][j][i][2]/lam_h) ) );
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, Vmin, 5, MPI_DOUBLE, MPI_MIN, comm);
  MPI_Allreduce(MPI_IN_PLACE, Vmax, 5, MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(MPI_IN_PLACE, &cmin, 1, MPI_DOUBLE, MPI_MIN, comm);
  MPI_Allreduce(MPI_IN_PLACE, &cmax, 1, MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(MPI_IN_PLACE, &Machmax, 1, MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(MPI_IN_PLACE, &char_speed_max, 1, MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(MPI_IN_PLACE, &dx_over_char_speed_min, 1, MPI_DOUBLE, MPI_MIN, comm);

  V.RestoreDataPointerToLocalVector(); 
  delta_xyz.RestoreDataPointerToLocalVector(); 
  ID.RestoreDataPointerToLocalVector();
}

//-----------------------------------------------------

void SpaceOperator::ComputeTimeStepSize(SpaceVariable3D &V, SpaceVariable3D &ID, double &dt, double &cfl)
{
  double Vmin[5], Vmax[5], cmin, cmax, Machmax, char_speed_max, dx_over_char_speed_min; 
  FindExtremeValuesOfFlowVariables(V, ID, Vmin, Vmax, cmin, cmax, Machmax, char_speed_max, dx_over_char_speed_min);

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

void SpaceOperator::ComputeAdvectionFluxes(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &F,
                                           RiemannSolutions *riemann_solutions)
{
  //------------------------------------
  // Preparation: Delete previous riemann_solutions
  //------------------------------------
  if(riemann_solutions)
    riemann_solutions->Clear();

  //------------------------------------
  // Reconstruction w/ slope limiters.
  //------------------------------------
  rec.Reconstruct(V, Vl, Vr, Vb, Vt, Vk, Vf);

  //------------------------------------
  // Check reconstructed states (clip & check)
  //------------------------------------
  CheckReconstructedStates(V, Vl, Vr, Vb, Vt, Vk, Vf, ID);

  //------------------------------------
  // Extract data
  //------------------------------------
  Vec5D*** v  = (Vec5D***) V.GetDataPointer();
  Vec5D*** vl = (Vec5D***) Vl.GetDataPointer();
  Vec5D*** vr = (Vec5D***) Vr.GetDataPointer();
  Vec5D*** vb = (Vec5D***) Vb.GetDataPointer();
  Vec5D*** vt = (Vec5D***) Vt.GetDataPointer();
  Vec5D*** vk = (Vec5D***) Vk.GetDataPointer();
  Vec5D*** vf = (Vec5D***) Vf.GetDataPointer();
  Vec5D*** f  = (Vec5D***) F.GetDataPointer();

  double*** id = (double***) ID.GetDataPointer();

  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer(); //for debugging
 
  //------------------------------------
  // Compute fluxes
  //------------------------------------
  Vec5D localflux1, localflux2;
  Vec3D*** dxyz = (Vec3D***)delta_xyz.GetDataPointer();

  // Initialize F to 0
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++) 
      for(int i=ii0; i<iimax; i++)
          f[k][j][i] = 0.0; //setting f[k][j][i][0] = ... = f[k][j][i][4] = 0.0;


  int myid;
  int neighborid = -1;
  int midid = -1;
  double Vmid[5];
  double Vsm[5], Vsp[5];
  double area = 0.0;
  Int3 ind;

  // Loop through the domain interior, and the right, top, and front ghost layers. For each cell, calculate the
  // numerical flux across the left, lower, and back cell boundaries/interfaces
  for(int k=k0; k<kkmax; k++) {
    for(int j=j0; j<jjmax; j++) {
      for(int i=i0; i<iimax; i++) {

        myid = id[k][j][i];

/*

        if(k==0 && j==1 && i==569) {

          fprintf(stderr,"(569,1,0): id = %d, V = %e %e %e %e %e\n", myid, v[k][j][i][0], v[k][j][i][1], v[k][j][i][2], v[k][j][i][3], v[k][j][i][4]);


        }
*/


        //*****************************************
        //calculate flux function F_{i-1/2,j,k}
        //*****************************************
        if(k!=kkmax-1 && j!=jjmax-1) {
 
          neighborid = id[k][j][i-1];

          if(neighborid==myid) {

            fluxFcn.ComputeNumericalFluxAtCellInterface(0/*F*/, vr[k][j][i-1]/*Vm*/, vl[k][j][i]/*Vp*/, myid, localflux1);
            localflux2 = localflux1;

          } else {//material interface


            if(iod.multiphase.recon == MultiPhaseData::CONSTANT)//switch back to constant reconstruction (i.e. v)
              riemann.ComputeRiemannSolution(0/*F*/, v[k][j][i-1], neighborid, v[k][j][i], myid, Vmid, midid, Vsm, Vsp);
            else//linear reconstruction w/ limitor
              riemann.ComputeRiemannSolution(0/*F*/, vr[k][j][i-1], neighborid, vl[k][j][i], myid, Vmid, midid, Vsm, Vsp);


            if(riemann_solutions) {//store Riemann solution for "phase-change update"
              ind[0] = k; ind[1] = j; ind[2] = i;
              riemann_solutions->left[ind] = std::make_pair((Vec5D)Vsm, neighborid); 
              ind[2] = i-1;
              riemann_solutions->right[ind] = std::make_pair((Vec5D)Vsp, myid); 
            }


            if(iod.multiphase.flux == MultiPhaseData::EXACT) { //Godunov-type flux
              fluxFcn.EvaluateFluxFunction_F(Vmid, midid, localflux1);
              localflux2 = localflux1;
            } else {//Numerical flux function

              if(iod.multiphase.recon == MultiPhaseData::CONSTANT) {//switch back to constant reconstruction (i.e. v)
                fluxFcn.ComputeNumericalFluxAtCellInterface(0/*F*/, v[k][j][i-1]/*Vm*/, Vsm/*Vp*/, neighborid, localflux1);
                fluxFcn.ComputeNumericalFluxAtCellInterface(0/*F*/, Vsp/*Vm*/, v[k][j][i]/*Vp*/, myid, localflux2);
              } else {//linear reconstruction w/ limiter
                fluxFcn.ComputeNumericalFluxAtCellInterface(0/*F*/, vr[k][j][i-1]/*Vm*/, Vsm/*Vp*/, neighborid, localflux1);
                fluxFcn.ComputeNumericalFluxAtCellInterface(0/*F*/, Vsp/*Vm*/, vl[k][j][i]/*Vp*/, myid, localflux2);
              }

            }

          }

          area = dxyz[k][j][i][1]*dxyz[k][j][i][2];
          f[k][j][i-1] += localflux1*area;
          f[k][j][i]   -= localflux2*area;

        }


        //*****************************************
        //calculate flux function G_{i,j-1/2,k}
        //*****************************************
        if(k!=kkmax-1 && i!=iimax-1) {

          neighborid = id[k][j-1][i];

          if(neighborid==myid) {

            fluxFcn.ComputeNumericalFluxAtCellInterface(1/*G*/, vt[k][j-1][i]/*Vm*/, vb[k][j][i]/*Vp*/, myid, localflux1);
            localflux2 = localflux1;

          } else {//material interface

            if(iod.multiphase.recon == MultiPhaseData::CONSTANT)//switch back to constant reconstruction (i.e. v)
              riemann.ComputeRiemannSolution(1/*G*/, v[k][j-1][i], neighborid, v[k][j][i], myid, Vmid, midid, Vsm, Vsp);
            else
              riemann.ComputeRiemannSolution(1/*G*/, vt[k][j-1][i], neighborid, vb[k][j][i], myid, Vmid, midid, Vsm, Vsp);


            if(riemann_solutions) {//store Riemann solution for "phase-change update"
              ind[0] = k; ind[1] = j; ind[2] = i;
              riemann_solutions->bottom[ind] = std::make_pair((Vec5D)Vsm, neighborid); 
              ind[1] = j-1;
              riemann_solutions->top[ind] = std::make_pair((Vec5D)Vsp, myid); 
            }


            if(iod.multiphase.flux == MultiPhaseData::EXACT) { //Godunov-type flux
              fluxFcn.EvaluateFluxFunction_G(Vmid, midid, localflux1);
              localflux2 = localflux1;
            } else {//Numerical flux function

              if(iod.multiphase.recon == MultiPhaseData::CONSTANT) {//switch back to constant reconstruction (i.e. v)
                fluxFcn.ComputeNumericalFluxAtCellInterface(1/*G*/, v[k][j-1][i]/*Vm*/, Vsm/*Vp*/, neighborid, localflux1);
                fluxFcn.ComputeNumericalFluxAtCellInterface(1/*G*/, Vsp/*Vm*/, v[k][j][i]/*Vp*/, myid, localflux2);
              } else {
                fluxFcn.ComputeNumericalFluxAtCellInterface(1/*G*/, vt[k][j-1][i]/*Vm*/, Vsm/*Vp*/, neighborid, localflux1);
                fluxFcn.ComputeNumericalFluxAtCellInterface(1/*G*/, Vsp/*Vm*/, vb[k][j][i]/*Vp*/, myid, localflux2);
              }

            }

          }

          area = dxyz[k][j][i][0]*dxyz[k][j][i][2];
          f[k][j-1][i] += localflux1*area;
          f[k][j][i]   -= localflux2*area;
        }


        //*****************************************
        //calculate flux function H_{i,j,k-1/2}
        //*****************************************
        if(j!=jjmax-1 && i!=iimax-1) {

          neighborid = id[k-1][j][i];

          if(neighborid==myid) {

            fluxFcn.ComputeNumericalFluxAtCellInterface(2/*H*/, vf[k-1][j][i]/*Vm*/, vk[k][j][i]/*Vp*/, myid, localflux1);
            localflux2 = localflux1;

          } else {//material interface

            if(iod.multiphase.recon == MultiPhaseData::CONSTANT) //switch back to constant reconstruction (i.e. v)
              riemann.ComputeRiemannSolution(2/*H*/, v[k-1][j][i], neighborid, v[k][j][i], myid, Vmid, midid, Vsm, Vsp);
            else
              riemann.ComputeRiemannSolution(2/*H*/, vf[k-1][j][i], neighborid, vk[k][j][i], myid, Vmid, midid, Vsm, Vsp);


            if(riemann_solutions) {//store Riemann solution for "phase-change update"
              ind[0] = k; ind[1] = j; ind[2] = i;
              riemann_solutions->back[ind] = std::make_pair((Vec5D)Vsm, neighborid); 
              ind[0] = k-1;
              riemann_solutions->front[ind] = std::make_pair((Vec5D)Vsp, myid); 
            }


            if(iod.multiphase.flux == MultiPhaseData::EXACT) { //Godunov-type flux
              fluxFcn.EvaluateFluxFunction_H(Vmid, midid, localflux1);
              localflux2 = localflux1;
            } else {//Numerical flux function

              if(iod.multiphase.recon == MultiPhaseData::CONSTANT) {//switch back to constant reconstruction (i.e. v)
                fluxFcn.ComputeNumericalFluxAtCellInterface(2/*H*/, v[k-1][j][i]/*Vm*/, Vsm/*Vp*/, neighborid, localflux1);
                fluxFcn.ComputeNumericalFluxAtCellInterface(2/*H*/, Vsp/*Vm*/, v[k][j][i]/*Vp*/, myid, localflux2);
              } else {
                fluxFcn.ComputeNumericalFluxAtCellInterface(2/*H*/, vf[k-1][j][i]/*Vm*/, Vsm/*Vp*/, neighborid, localflux1);
                fluxFcn.ComputeNumericalFluxAtCellInterface(2/*H*/, Vsp/*Vm*/, vk[k][j][i]/*Vp*/, myid, localflux2);
              }

            }

          }

          area = dxyz[k][j][i][0]*dxyz[k][j][i][1];
          f[k-1][j][i] += localflux1*area;
          f[k][j][i]   -= localflux2*area;
        }
      }
    }
  }
        
  //------------------------------------
  // Restore Spatial Variables
  //------------------------------------
  delta_xyz.RestoreDataPointerToLocalVector(); //no changes
  ID.RestoreDataPointerToLocalVector(); //no changes
  coordinates.RestoreDataPointerToLocalVector(); //no changes

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

void 
SpaceOperator::CheckReconstructedStates(SpaceVariable3D &V,
                                        SpaceVariable3D &Vl, SpaceVariable3D &Vr, SpaceVariable3D &Vb,
                                        SpaceVariable3D &Vt, SpaceVariable3D &Vk, SpaceVariable3D &Vf,
                                        SpaceVariable3D &ID)
{

  Vec5D*** v  = (Vec5D***) V.GetDataPointer();
  Vec5D*** vl = (Vec5D***) Vl.GetDataPointer();
  Vec5D*** vr = (Vec5D***) Vr.GetDataPointer();
  Vec5D*** vb = (Vec5D***) Vb.GetDataPointer();
  Vec5D*** vt = (Vec5D***) Vt.GetDataPointer();
  Vec5D*** vk = (Vec5D***) Vk.GetDataPointer();
  Vec5D*** vf = (Vec5D***) Vf.GetDataPointer();

  double*** id = (double***) ID.GetDataPointer();

  //------------------------------------
  // Clip pressure and density for the reconstructed state
  // Verify hyperbolicity (i.e. c^2 > 0).
  //------------------------------------
  int nClipped = 0;
  bool error = false;
  int boundary;
  int myid;
  for(int k=kk0; k<kkmax; k++) {
    for(int j=jj0; j<jjmax; j++) {
      for(int i=ii0; i<iimax; i++) {

        boundary = 0;
        if(k==kk0 || k==kkmax-1) boundary++;
        if(j==jj0 || j==jjmax-1) boundary++;
        if(i==ii0 || i==iimax-1) boundary++;
        if(boundary>=2) //not needed
          continue;

        myid = id[k][j][i];

        if(boundary==0) {//interior
          nClipped += (int)varFcn[myid]->ClipDensityAndPressure(vl[k][j][i]);
          nClipped += (int)varFcn[myid]->ClipDensityAndPressure(vr[k][j][i]);
          nClipped += (int)varFcn[myid]->ClipDensityAndPressure(vb[k][j][i]);
          nClipped += (int)varFcn[myid]->ClipDensityAndPressure(vt[k][j][i]);
          nClipped += (int)varFcn[myid]->ClipDensityAndPressure(vk[k][j][i]);
          nClipped += (int)varFcn[myid]->ClipDensityAndPressure(vf[k][j][i]);
         
          error = varFcn[myid]->CheckState(vl[k][j][i]) || varFcn[myid]->CheckState(vr[k][j][i]) || 
                  varFcn[myid]->CheckState(vb[k][j][i]) || varFcn[myid]->CheckState(vt[k][j][i]) ||
                  varFcn[myid]->CheckState(vk[k][j][i]) || varFcn[myid]->CheckState(vf[k][j][i]);
        } else {//boundary face
          if(i==ii0) {
            nClipped += (int)varFcn[myid]->ClipDensityAndPressure(vr[k][j][i]);
            error = error || varFcn[myid]->CheckState(vr[k][j][i]);
          } else if (i==iimax-1) {
            nClipped += (int)varFcn[myid]->ClipDensityAndPressure(vl[k][j][i]);
            error = error || varFcn[myid]->CheckState(vl[k][j][i]);
          } else if (j==jj0) {
            nClipped += (int)varFcn[myid]->ClipDensityAndPressure(vt[k][j][i]);
            error = error || varFcn[myid]->CheckState(vt[k][j][i]);
          } else if (j==jjmax-1) {
            nClipped += (int)varFcn[myid]->ClipDensityAndPressure(vb[k][j][i]);
            error = error || varFcn[myid]->CheckState(vb[k][j][i]);
          } else if (k==kk0) {
            nClipped += (int)varFcn[myid]->ClipDensityAndPressure(vf[k][j][i]);
            error = error || varFcn[myid]->CheckState(vf[k][j][i]);
          } else if (k==kkmax-1) {
            nClipped += (int)varFcn[myid]->ClipDensityAndPressure(vk[k][j][i]);
            error = error || varFcn[myid]->CheckState(vk[k][j][i]);
          } 
        }

        if(error) {
          print_error("*** Error: Reconstructed state at (%d,%d,%d) violates hyperbolicity. matid = %d.\n", i,j,k, myid);
          fprintf(stderr, "v[%d,%d,%d]  = [%e, %e, %e, %e, %e]\n", i,j,k, v[k][j][i][0], v[k][j][i][1], v[k][j][i][2], v[k][j][i][3], v[k][j][i][4]);
          fprintf(stderr, "vl[%d,%d,%d] = [%e, %e, %e, %e, %e]\n", i,j,k, vl[k][j][i][0], vl[k][j][i][1], vl[k][j][i][2], vl[k][j][i][3], vl[k][j][i][4]);
          fprintf(stderr, "vr[%d,%d,%d] = [%e, %e, %e, %e, %e]\n", i,j,k, vr[k][j][i][0], vr[k][j][i][1], vr[k][j][i][2], vr[k][j][i][3], vr[k][j][i][4]);
          fprintf(stderr, "vb[%d,%d,%d] = [%e, %e, %e, %e, %e]\n", i,j,k, vb[k][j][i][0], vb[k][j][i][1], vb[k][j][i][2], vb[k][j][i][3], vb[k][j][i][4]);
          fprintf(stderr, "vt[%d,%d,%d] = [%e, %e, %e, %e, %e]\n", i,j,k, vt[k][j][i][0], vt[k][j][i][1], vt[k][j][i][2], vt[k][j][i][3], vt[k][j][i][4]);
          fprintf(stderr, "vk[%d,%d,%d] = [%e, %e, %e, %e, %e]\n", i,j,k, vk[k][j][i][0], vk[k][j][i][1], vk[k][j][i][2], vk[k][j][i][3], vk[k][j][i][4]);
          fprintf(stderr, "vf[%d,%d,%d] = [%e, %e, %e, %e, %e]\n", i,j,k, vf[k][j][i][0], vf[k][j][i][1], vf[k][j][i][2], vf[k][j][i][3], vf[k][j][i][4]);
          exit_mpi();
        } 
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &nClipped, 1, MPI_INT, MPI_SUM, comm);
  if(nClipped)
    print("Warning: Clipped pressure and/or density in %d reconstructed states.\n", nClipped);
 
  Vl.RestoreDataPointerToLocalVector(); //no need to communicate
  Vr.RestoreDataPointerToLocalVector(); 
  Vb.RestoreDataPointerToLocalVector(); 
  Vt.RestoreDataPointerToLocalVector(); 
  Vk.RestoreDataPointerToLocalVector(); 
  Vf.RestoreDataPointerToLocalVector(); 
  ID.RestoreDataPointerToLocalVector();
}

//-----------------------------------------------------

void SpaceOperator::ComputeResidual(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &R,
                                    RiemannSolutions *riemann_solutions)
{
  ComputeAdvectionFluxes(V, ID, R, riemann_solutions);

  // -------------------------------------------------
  // multiply flux by -1, and divide by cell volume (for cells within the actual domain)
  // -------------------------------------------------
  Vec5D***    r = (Vec5D***) R.GetDataPointer();
  double*** vol = (double***)volume.GetDataPointer();

  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++) 
      for(int i=i0; i<imax; i++) {
        r[k][j][i] /= -vol[k][j][i];


/*

        if(k==0 && j==1 && i==569) {

          fprintf(stderr,"(569,1,0): r = %e %e %e %e %e\n", r[k][j][i][0], r[k][j][i][1], r[k][j][i][2], r[k][j][i][3], r[k][j][i][4]);

        }
*/

      }

  // restore spatial variables
  R.RestoreDataPointerToLocalVector(); //NOTE: although R has been updated, there is no need of 
                                       //      cross-subdomain communications. So, no need to 
                                       //      update the global vec.
  volume.RestoreDataPointerToLocalVector();
}


//-----------------------------------------------------









