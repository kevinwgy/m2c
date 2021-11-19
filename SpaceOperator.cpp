#include <Utils.h>
#include <SpaceOperator.h>
#include <Vector3D.h>
#include <Vector5D.h>
#include <GeoTools.h>
#include <DistancePointToSpheroid.h>
#include <algorithm> //std::upper_bound
#include <cfloat> //DBL_MAX
#include <KDTree.h>
#include <rbf_interp_2d.hpp>
using std::cout;
using std::endl;
using std::max;
using std::min;
using std::map;
using namespace GeoTools;

extern int verbose;
//-----------------------------------------------------

SpaceOperator::SpaceOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_,
                             vector<VarFcnBase*> &varFcn_, FluxFcnBase &fluxFcn_,
                             ExactRiemannSolverBase &riemann_,
                             vector<double> &x, vector<double> &y, vector<double> &z,
                             vector<double> &dx, vector<double> &dy, vector<double> &dz,
                             bool screenout) 
  : comm(comm_), dm_all(dm_all_),
    iod(iod_), varFcn(varFcn_), fluxFcn(fluxFcn_), riemann(riemann_),
    coordinates(comm_, &(dm_all_.ghosted1_3dof)),
    delta_xyz(comm_, &(dm_all_.ghosted1_3dof)),
    volume(comm_, &(dm_all_.ghosted1_1dof)),
    rec(comm_, dm_all_, iod_.schemes.ns.rec, coordinates, delta_xyz, &varFcn, &fluxFcn),
    Vl(comm_, &(dm_all_.ghosted1_5dof)),
    Vr(comm_, &(dm_all_.ghosted1_5dof)),
    Vb(comm_, &(dm_all_.ghosted1_5dof)),
    Vt(comm_, &(dm_all_.ghosted1_5dof)),
    Vk(comm_, &(dm_all_.ghosted1_5dof)),
    Vf(comm_, &(dm_all_.ghosted1_5dof)),
    Utmp(comm_, &(dm_all_.ghosted1_5dof)),
    symm(NULL), visco(NULL), smooth(NULL)
{
  
  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);
  coordinates.GetGlobalSize(&NX, &NY, &NZ);

  SetupMesh(x,y,z,dx,dy,dz);

  CreateGhostNodeLists(screenout); //create ghost_nodes_inner and ghost_nodes_outer

  rec.Setup(&ghost_nodes_inner, &ghost_nodes_outer); //this function requires mesh info (dxyz)
  
  if(iod.mesh.type == MeshData::SPHERICAL || iod.mesh.type == MeshData::CYLINDRICAL)
    symm = new SymmetryOperator(comm, dm_all, iod.mesh, varFcn, coordinates, delta_xyz, volume);

  if(iod.schemes.ns.smooth.type != SmoothingData::NONE)
    smooth = new SmoothingOperator(comm, dm_all, iod.schemes.ns.smooth, coordinates, delta_xyz, volume);
    
}

//-----------------------------------------------------

SpaceOperator::~SpaceOperator()
{
  if(symm) delete symm;
  if(visco) delete visco;
  if(smooth) delete smooth;
}

//-----------------------------------------------------

void SpaceOperator::Destroy()
{
  rec.Destroy();

  if(symm) symm->Destroy();

  if(visco) visco->Destroy();

  if(smooth) smooth->Destroy();

  coordinates.Destroy();
  delta_xyz.Destroy();
  volume.Destroy();
  Vl.Destroy();
  Vr.Destroy();
  Vb.Destroy();
  Vt.Destroy();
  Vk.Destroy();
  Vf.Destroy();
  Utmp.Destroy();
}

//-----------------------------------------------------

void SpaceOperator::SetupMesh(vector<double> &x, vector<double> &y, vector<double> &z,
                              vector<double> &dx, vector<double> &dy, vector<double> &dz)
{
  //! Setup coordinates of cell centers and dx, dy, dz

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

  int nnx, nny, nnz;
  coordinates.GetGhostedSize(&nnx, &nny, &nnz);

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

/** Re-set ghost nodes */
void SpaceOperator::ResetGhostLayer(double* xminus, double* xplus, double* yminus, double* yplus,
                        double* zminus, double* zplus, double* dxminus, double* dxplus, double* dyminus, 
                        double* dyplus, double* dzminus, double* dzplus)
{
  Vec3D*** v    = (Vec3D***) coordinates.GetDataPointer();
  Vec3D*** dxyz = (Vec3D***) delta_xyz.GetDataPointer();

  int nnx, nny, nnz;
  coordinates.GetGhostedSize(&nnx, &nny, &nnz);

  // capture the mesh info of the corners 
  double v0[3], v1[3];
  double dxyz0[3], dxyz1[3];
  for(int p=0; p<3; p++) {
    v0[p]    = v[k0][j0][i0][p] - dxyz[k0][j0][i0][p];
    v1[p]    = v[kmax-1][jmax-1][imax-1][p] + dxyz[kmax-1][jmax-1][imax-1][p];
    dxyz0[p] = dxyz[k0][j0][i0][p];
    dxyz1[p] = dxyz[kmax-1][jmax-1][imax-1][p];
  }

  //-----------------------------------------------------------------
  // RESET
  if(xminus) {v0[0] = *xminus;} if(xplus) {v1[0] = *xplus;}
  if(yminus) {v0[1] = *yminus;} if(yplus) {v1[1] = *yplus;}
  if(zminus) {v0[2] = *zminus;} if(zplus) {v1[2] = *zplus;}
  if(dxminus) {dxyz0[0] = *dxminus;} if(dxplus) {dxyz1[0] = *dxplus;}
  if(dyminus) {dxyz0[1] = *dyminus;} if(dyplus) {dxyz1[1] = *dyplus;}
  if(dzminus) {dxyz0[2] = *dzminus;} if(dzplus) {dxyz1[2] = *dzplus;}
  //-----------------------------------------------------------------


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


  //! Compute mesh information
  dxyz = (Vec3D***)delta_xyz.GetDataPointer();
  double*** vol  = (double***)volume.GetDataPointer();

  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++)
        vol[k][j][i] /*volume of cv*/ = dxyz[k][j][i][0]*dxyz[k][j][i][1]*dxyz[k][j][i][2];

  delta_xyz.RestoreDataPointerAndInsert();
  volume.RestoreDataPointerAndInsert();


  CreateGhostNodeLists(false); //create ghost_nodes_inner and ghost_nodes_outer

}


//-----------------------------------------------------

void SpaceOperator::CreateGhostNodeLists(bool screenout)
{
  ghost_nodes_inner.clear();
  ghost_nodes_outer.clear();

  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();

  Int3 image;
  Vec3D proj(0.0), out_normal(0.0);
  MeshData::BcType bcType = MeshData::NONE;
  GhostPoint::Side side = GhostPoint::UNDEFINED;
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
          side       = GhostPoint::UNDEFINED;

          if(i<0)        {image[0] = -i-1;
                          proj[0]  = iod.mesh.x0;      out_normal[0] = -1.0; 
                          bcType   = iod.mesh.bc_x0;   side = GhostPoint::LEFT;     counter++;}
          else if(i>=NX) {image[0] = NX+(i-NX)-1;
                          proj[0]  = iod.mesh.xmax;    out_normal[0] =  1.0; 
                          bcType   = iod.mesh.bc_xmax; side = GhostPoint::RIGHT;    counter++;}
          else           {image[0] = i;
                          proj[0]  = coords[k][j][i][0];}
                     

          if(j<0)        {image[1] = -j-1;
                          proj[1]  = iod.mesh.y0;      out_normal[1] = -1.0; 
                          bcType   = iod.mesh.bc_y0;   side = GhostPoint::BOTTOM;   counter++;}
          else if(j>=NY) {image[1] = NY+(j-NY)-1;
                          proj[1]  = iod.mesh.ymax;    out_normal[1] =  1.0; 
                          bcType   = iod.mesh.bc_ymax; side = GhostPoint::TOP;      counter++;}
          else           {image[1] = j;
                          proj[1]  = coords[k][j][i][1];}
         

          if(k<0)        {image[2] = -k-1;
                          proj[2]  = iod.mesh.z0;      out_normal[2] = -1.0;
                          bcType   = iod.mesh.bc_z0;   side = GhostPoint::BACK;     counter++;}
          else if(k>=NZ) {image[2] = NZ+(k-NZ)-1;
                          proj[2]  = iod.mesh.zmax;    out_normal[2] =  1.0; 
                          bcType   = iod.mesh.bc_zmax; side = GhostPoint::FRONT;    counter++;}
          else           {image[2] = k;
                          proj[2]  = coords[k][j][i][2];}
         
          out_normal /= out_normal.norm();

          assert(counter<=3 && counter>0);

          if(counter == 1)
            ghost_nodes_outer.push_back(GhostPoint(Int3(i,j,k), image, GhostPoint::FACE,
                                        proj, out_normal, (int)bcType, side));
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
  if(screenout) {
    print(comm,"- Number of ghost nodes inside computational domain (overlapping between subdomains): %d\n",
          nInner);
    print(comm,"- Number of ghost nodes outside computational domain: %d\n",
          nOuter);
    print(comm,"\n");
  }

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
  if(nClipped && verbose>0) {
    print("Warning: Clipped pressure and/or density in %d cells.\n", nClipped);
    V.RestoreDataPointerAndInsert();
  } else
    V.RestoreDataPointerToLocalVector();

  ID.RestoreDataPointerToLocalVector(); //no changes made

  coordinates.RestoreDataPointerToLocalVector(); //no changes

  return nClipped;
}  

//-----------------------------------------------------
//assign interpolator and gradien calculator (pointers) to the viscosity operator
void SpaceOperator::SetupViscosityOperator(InterpolatorBase *interpolator_, GradientCalculatorBase *grad_)
{
  bool hasViscosity = false;
  for(auto it = iod.eqs.materials.dataMap.begin(); it != iod.eqs.materials.dataMap.end(); it++) {
    if(it->second->viscosity.type != ViscosityModelData::NONE) {
      hasViscosity = true;
      break; //initialize visco if any material has viscosity
    }
  }

  if(hasViscosity) {
    assert(interpolator_); //make sure it is not NULL
    assert(grad_);
    visco = new ViscosityOperator(comm, dm_all, iod.eqs, varFcn, coordinates, delta_xyz,
                                  *interpolator_, *grad_);
  }
}

//-----------------------------------------------------

//apply IC within the real domain
void SpaceOperator::SetInitialCondition(SpaceVariable3D &V, SpaceVariable3D &ID) 
{
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();

  Vec5D*** v = (Vec5D***) V.GetDataPointer();
  double*** id = (double***) ID.GetDataPointer();

  //! 1. apply the inlet (i.e. farfield) state
  if(iod.bc.inlet.materialid != 0) {
    print_error("*** Error: Material at the inlet boundary should have MaterialID = 0. (Found %d instead.)\n", 
                iod.bc.inlet.materialid);
    exit_mpi();
  }
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
        print("- Applying the initial condition specified in %s (planar).\n", 
              iod.ic.user_specified_ic);
      else
        print("- Applying the initial condition specified in %s (with cylindrical symmetry).\n", 
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

    else if (iod.ic.type == IcData::GENERALCYLINDRICAL) {

      Vec3D dir(iod.ic.dir[0], iod.ic.dir[1], iod.ic.dir[2]); 
      dir /= dir.norm();

      int N = iod.ic.user_data[IcData::COORDINATE].size(); //!< number of data points provided by user

      print("- Applying the initial condition specified in %s (%d points, cylindrical symmetry).\n", 
            iod.ic.user_specified_ic, N);

      //Store sample points in a KDTree (K = 2)
      Vec2D *xy = new Vec2D[N]; 
      PointIn2D *p = new PointIn2D[N];
      for(int i=0; i<N; i++) {
        xy[i] = Vec2D(iod.ic.user_data[IcData::COORDINATE][i], iod.ic.user_data[IcData::RADIALCOORDINATE][i]);
        p[i]  = PointIn2D(i, xy[i]);
      }
      KDTree<PointIn2D,2/*dim*/> tree(N, p);
      print("    o Constructed a KDTree (K=2) to store the data points.\n");

      int numPoints = 15; //this is the number of points for interpolation
      int maxCand = numPoints*20;
      PointIn2D candidates[maxCand]; 


      // tentative cutoff distance --- will be adjusted
      double cutoff = sqrt((iod.ic.xmax[0]-iod.ic.xmin[0])*(iod.ic.xmax[1]-iod.ic.xmin[1])/N);

      Vec2D pnode;
      for(int k=k0; k<kmax; k++)
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {

            //3D -> 2D
            pnode[0] = (coords[k][j][i] - x0)*dir; //!< projection onto the 1D axis
            if(pnode[0]<iod.ic.xmin[0] || pnode[0]>iod.ic.xmax[0])
              continue;
            pnode[1] = (coords[k][j][i] - x0 - pnode[0]*dir).norm();
            if(pnode[1]<iod.ic.xmin[1] || pnode[1]>iod.ic.xmax[1])
              continue;
 
            //RBF interpolation w/ KDTree
            int nFound = 0, counter = 0;
            while(nFound<numPoints || nFound>maxCand) {
              if(++counter>2000) {
                fprintf(stderr,"*** Error: Cannot find required number of sample points for "
                               "interpolation (by RBF) after 2000 iterations. "
                               "Coord(3D):%e %e %e  | Coord(2D):%e %e. "
                               "Candidates: %d, cutoff = %e.\n", 
                               coords[k][j][i][0], coords[k][j][i][1], coords[k][j][i][2],
                               pnode[0], pnode[1], nFound, cutoff);
                exit(-1);
              }
              nFound = tree.findCandidatesWithin(pnode, candidates, maxCand, cutoff);
              if(nFound==0)              cutoff *= 4.0;
              else if(nFound<numPoints)  cutoff *= 1.5*sqrt((double)numPoints/(double)nFound);
              else if(nFound>maxCand)    cutoff /= 1.5*sqrt((double)nFound/(double)maxCand);
            }

            //figure out the actual points for interpolation (numPoints)
            vector<pair<double,int> > dist2node;
            for(int i=0; i<nFound; i++) 
              dist2node.push_back(std::make_pair((candidates[i].x-pnode).norm(), candidates[i].id));

            if(nFound>numPoints)
              sort(dist2node.begin(), dist2node.end());         
            dist2node.resize(numPoints);

            //prepare to interpolate
            double xd[2*numPoints];
            for(int i=0; i<numPoints; i++) {
              xd[2*i]   = xy[dist2node[i].second][0];
              xd[2*i+1] = xy[dist2node[i].second][1];
            }
            double r0; //smaller than maximum separation, larger than typical separation
            r0 = dist2node.front().first + dist2node.back().first;
            double fd[numPoints];
            double *rbf_weight, *interp;

            if(iod.ic.specified[IcData::DENSITY]) {
              for(int i=0; i<numPoints; i++)
                fd[i] = iod.ic.user_data[IcData::DENSITY][dist2node[i].second];
              rbf_weight = MathTools::rbf_weight(2, numPoints, xd, r0,
                                                 MathTools::phi1, //multiquadric RBF
                                                 fd);
              interp = MathTools::rbf_interp(2, numPoints, xd, r0,
                                             MathTools::phi1, rbf_weight,
                                             1, pnode);
              v[k][j][i][0] = interp[0];

              delete [] rbf_weight;
              delete [] interp;
            }

            if(iod.ic.specified[IcData::VELOCITY] || iod.ic.specified[IcData::RADIALVELOCITY])
              v[k][j][i][1] = v[k][j][i][2] = v[k][j][i][3] = 0.0;

            if(iod.ic.specified[IcData::VELOCITY]) {
              for(int i=0; i<numPoints; i++)
                fd[i] = iod.ic.user_data[IcData::VELOCITY][dist2node[i].second];
              rbf_weight = MathTools::rbf_weight(2, numPoints, xd, r0,
                                                 MathTools::phi1, //multiquadric RBF
                                                 fd);
              interp = MathTools::rbf_interp(2, numPoints, xd, r0,
                                             MathTools::phi1, rbf_weight,
                                             1, pnode);
              v[k][j][i][1] = interp[0]*dir[0]; 
              v[k][j][i][2] = interp[0]*dir[1]; 
              v[k][j][i][3] = interp[0]*dir[2];

              delete [] rbf_weight;
              delete [] interp;
            }

            if(iod.ic.specified[IcData::RADIALVELOCITY]) {
              for(int i=0; i<numPoints; i++)
                fd[i] = iod.ic.user_data[IcData::RADIALVELOCITY][dist2node[i].second];
              rbf_weight = MathTools::rbf_weight(2, numPoints, xd, r0,
                                                 MathTools::phi1, //multiquadric RBF
                                                 fd);
              interp = MathTools::rbf_interp(2, numPoints, xd, r0,
                                             MathTools::phi1, rbf_weight,
                                             1, pnode);
              Vec3D dir2 = coords[k][j][i] - x0 - pnode[0]*dir;
              if(dir2.norm()>0)
                dir2 /= dir2.norm();
              v[k][j][i][1] += interp[0]*dir2[0]; 
              v[k][j][i][2] += interp[0]*dir2[1]; 
              v[k][j][i][3] += interp[0]*dir2[2];

              delete [] rbf_weight;
              delete [] interp;
            }

            if(iod.ic.specified[IcData::PRESSURE]) {
              for(int i=0; i<numPoints; i++)
                fd[i] = iod.ic.user_data[IcData::PRESSURE][dist2node[i].second];
              rbf_weight = MathTools::rbf_weight(2, numPoints, xd, r0,
                                                 MathTools::phi1, //multiquadric RBF
                                                 fd);
              interp = MathTools::rbf_interp(2, numPoints, xd, r0,
                                             MathTools::phi1, rbf_weight,
                                             1, pnode);
              v[k][j][i][4] = interp[0];

              delete [] rbf_weight;
              delete [] interp;
            }

            if(iod.ic.specified[IcData::MATERIALID]) {
              for(int i=0; i<numPoints; i++)
                fd[i] = iod.ic.user_data[IcData::MATERIALID][dist2node[i].second];
              rbf_weight = MathTools::rbf_weight(2, numPoints, xd, r0,
                                                 MathTools::phi1, //multiquadric RBF
                                                 fd);
              interp = MathTools::rbf_interp(2, numPoints, xd, r0,
                                             MathTools::phi1, rbf_weight,
                                             1, pnode);
              id[k][j][i] = std::round(interp[0]);

              delete [] rbf_weight;
              delete [] interp;
            }

          }

      delete [] xy;
      delete [] p;

    }

    else if (iod.ic.type == IcData::SPHERICAL) {

      print("- Applying the initial condition specified in %s (with spherical symmetry).\n", 
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

    print("- Applying initial condition on one side of a plane (material id: %d).\n\n", 
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

    print("- Applying initial condition within a cylinder-cone (material id: %d).\n\n",
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


  // cylinder with spherical cap(s)
  for(auto it=ic.cylindersphereMap.dataMap.begin(); it!=ic.cylindersphereMap.dataMap.end(); it++) {

    print("- Applying initial condition within a cylinder-sphere (material id: %d).\n\n",
          it->second->initialConditions.materialid);
    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
    Vec3D dir(it->second->nx, it->second->ny, it->second->nz);
    dir /= dir.norm();

    double L = it->second->L; //cylinder height
    double R = it->second->r; //cylinder & sphere's radius
    double Lhalf = 0.5*L;
    bool front_cap = (it->second->front_cap == CylinderSphereData::On);
    bool back_cap = (it->second->back_cap == CylinderSphereData::On);

    Vec3D xf = x0 + Lhalf*dir;
    Vec3D xb = x0 - Lhalf*dir;
    double x, r;
    bool inside;
    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++) {
          inside = false;
          x = (coords[k][j][i]-x0)*dir;
          r = (coords[k][j][i] - x0 - x*dir).norm();
          if(x>-Lhalf && x<Lhalf && r<R)
            inside = true;
          else if(front_cap && x>=Lhalf && (coords[k][j][i]-xf).norm() < R)
            inside = true;
          else if(back_cap && x<=-Lhalf && (coords[k][j][i]-xb).norm() < R)
            inside = true;

          if(inside) {
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

    print("- Applying initial condition within a sphere (material id: %d).\n\n",
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


  // spheroids 
  for(auto it=ic.spheroidMap.dataMap.begin(); it!=ic.spheroidMap.dataMap.end(); it++) {

    print("- Applying initial condition within a spheroid (material id: %d).\n\n",
          it->second->initialConditions.materialid);
    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
    Vec3D axis(it->second->axis_x, it->second->axis_y, it->second->axis_z);
    GeoTools::DistanceFromPointToSpheroid distCal(x0, axis, it->second->length, it->second->diameter);
    double dist;
    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++) {
          dist = distCal.Calculate(coords[k][j][i]); //>0 outside the spheroid
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
  ClipDensityAndPressure(V, ID);
  ApplyBoundaryConditions(V);   

}

//-----------------------------------------------------
//! Apply boundary conditions by populating the ghost cells.
void SpaceOperator::ApplyBoundaryConditions(SpaceVariable3D &V)
{
  Vec5D*** v = (Vec5D***) V.GetDataPointer();

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
SpaceOperator::ApplySmoothingFilter(double time, double dt, int time_step, SpaceVariable3D &V, SpaceVariable3D &ID)
{
  if(!smooth)
    return; //nothing to do

  if(!isTimeToWrite(time, dt, time_step, iod.schemes.ns.smooth.frequency_dt, 
                    iod.schemes.ns.smooth.frequency, 0.0, false))
    return; //nothing to do (not the right time)
  
  print("- Applying a smoothing filter to the N-S state variables. (Iterations: %d)\n", 
        iod.schemes.ns.smooth.iteration);

  for(int iter = 0; iter < iod.schemes.ns.smooth.iteration; iter++) {
    PrimitiveToConservative(V, ID, Utmp, false);
    smooth->ApplySmoothingFilter(Utmp,&ID);
    ConservativeToPrimitive(Utmp, ID, V, false);
    ClipDensityAndPressure(V,ID);
    ApplyBoundaryConditions(V);
  } 
}

//-----------------------------------------------------

void
SpaceOperator::ApplyBoundaryConditionsGeometricEntities(Vec5D*** v)
{

  map<int, DiskData* >& disks(iod.bc.multiBoundaryConditions.diskMap.dataMap);
  map<int, RectangleData* >& rectangles(iod.bc.multiBoundaryConditions.rectangleMap.dataMap);

  if(!disks.size() && !rectangles.size())
    return;

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

  if(verbose>1)
    print("- Maximum values: rho = %e, p = %e, c = %e, Mach = %e, char. speed = %e.\n", 
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
                                           RiemannSolutions *riemann_solutions, vector<int> *ls_mat_id, 
                                           vector<SpaceVariable3D*> *Phi)
{
  //------------------------------------
  // Preparation: Delete previous riemann_solutions
  //------------------------------------
  if(riemann_solutions)
    riemann_solutions->Clear();

  //------------------------------------
  // Reconstruction w/ slope limiters.
  //------------------------------------
  rec.Reconstruct(V, Vl, Vr, Vb, Vt, Vk, Vf, &ID);

  //------------------------------------
  // Check reconstructed states (clip & check)
  //------------------------------------
  CheckReconstructedStates(V, Vl, Vr, Vb, Vt, Vk, Vf, ID); //already checked in Reconstructor::Reconstruct
                                                           //but here we also do clipping

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

  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer(); 
  Vec3D*** dxyz   = (Vec3D***)delta_xyz.GetDataPointer(); 
 
  //------------------------------------
  // Extract level set gradient data
  //------------------------------------
  vector<double***> phi;
  if(Phi && ls_mat_id) {
    phi.resize(Phi->size(), NULL);
    for(int i=0; i<phi.size(); i++)
      phi[i] = (*Phi)[i]->GetDataPointer();
  }
 
  //------------------------------------
  // Compute fluxes
  //------------------------------------
  Vec5D localflux1, localflux2;

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

  // Count the riemann solver errors. Note that when this occurs, (1) it may be that the Riemann problem does NOT
  // have a solution, and (2) the solver returns an "approximate" solution (i.e. code may keep running)
  int riemann_errors = 0; 
  int err;

  // Loop through the domain interior, and the right, top, and front ghost layers. For each cell, calculate the
  // numerical flux across the left, lower, and back cell boundaries/interfaces
  for(int k=k0; k<kkmax; k++) {
    for(int j=j0; j<jjmax; j++) {
      for(int i=i0; i<iimax; i++) {

        myid = id[k][j][i];

        //*****************************************
        //calculate flux function F_{i-1/2,j,k}
        //*****************************************
        if(k!=kkmax-1 && j!=jjmax-1) {
 
          neighborid = id[k][j][i-1];

          if(neighborid==myid) {

            fluxFcn.ComputeNumericalFluxAtCellInterface(0/*F*/, vr[k][j][i-1]/*Vm*/, vl[k][j][i]/*Vp*/, myid, localflux1);
            localflux2 = localflux1;

          } else {//material interface

            // determine the axis/direction of the 1D Riemann problem
            Vec3D dir;
            if(iod.multiphase.riemann_normal == MultiPhaseData::MESH) //use mesh-based normal
              dir = Vec3D(1.0,0.0,0.0);
            else { //LEVEL_SET or AVERAGE

              // calculate level-set gradient at cell interface
              dir = CalculateGradPhiAtCellInterface(0/*i-1/2*/,i,j,k,coords,dxyz,myid,neighborid,ls_mat_id,&phi);

              if(iod.multiphase.riemann_normal == MultiPhaseData::AVERAGE) {
                dir[0] += 1.0;
                dir /= dir.norm();
              }
            }

            //Solve 1D Riemann problem
            if(iod.multiphase.recon == MultiPhaseData::CONSTANT)//switch back to constant reconstruction (i.e. v)
              err = riemann.ComputeRiemannSolution(dir, v[k][j][i-1], neighborid, v[k][j][i], myid, Vmid, midid, Vsm, Vsp);
            else//linear reconstruction w/ limitor
              err = riemann.ComputeRiemannSolution(dir, vr[k][j][i-1], neighborid, vl[k][j][i], myid, Vmid, midid, Vsm, Vsp);

            if(err) 
              riemann_errors++;

            //Clip Riemann solution and check it
            varFcn[neighborid]->ClipDensityAndPressure(Vsm);
            varFcn[neighborid]->CheckState(Vsm);
            varFcn[myid]->ClipDensityAndPressure(Vsp);
            varFcn[myid]->CheckState(Vsp);

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

            // determine the axis/direction of the 1D Riemann problem
            Vec3D dir;
            if(iod.multiphase.riemann_normal == MultiPhaseData::MESH) //use mesh-based normal
              dir = Vec3D(0.0,1.0,0.0);
            else { //LEVEL_SET or AVERAGE

              // calculate level-set gradient at cell interface
              dir = CalculateGradPhiAtCellInterface(1/*j-1/2*/,i,j,k,coords,dxyz,myid,neighborid,ls_mat_id,&phi);

              if(iod.multiphase.riemann_normal == MultiPhaseData::AVERAGE) {
                dir[1] += 1.0;
                dir /= dir.norm();
              }
            }

            //Solve 1D Riemann problem
            if(iod.multiphase.recon == MultiPhaseData::CONSTANT)//switch back to constant reconstruction (i.e. v)
              err = riemann.ComputeRiemannSolution(dir, v[k][j-1][i], neighborid, v[k][j][i], myid, Vmid, midid, Vsm, Vsp);
            else
              err = riemann.ComputeRiemannSolution(dir, vt[k][j-1][i], neighborid, vb[k][j][i], myid, Vmid, midid, Vsm, Vsp);

            if(err) 
              riemann_errors++;

            //Clip Riemann solution and check it
            varFcn[neighborid]->ClipDensityAndPressure(Vsm);
            varFcn[neighborid]->CheckState(Vsm);
            varFcn[myid]->ClipDensityAndPressure(Vsp);
            varFcn[myid]->CheckState(Vsp);


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

            // determine the axis/direction of the 1D Riemann problem
            Vec3D dir;
            if(iod.multiphase.riemann_normal == MultiPhaseData::MESH) //use mesh-based normal
              dir = Vec3D(0.0,0.0,1.0);
            else {

              // calculate level-set gradient at cell interface
              dir = CalculateGradPhiAtCellInterface(2/*k-1/2*/,i,j,k,coords,dxyz,myid,neighborid,ls_mat_id,&phi);

              if(iod.multiphase.riemann_normal == MultiPhaseData::AVERAGE) {
                dir[2] += 1.0;
                dir /= dir.norm();
              }
            }

            //Solve 1D Riemann problem
            if(iod.multiphase.recon == MultiPhaseData::CONSTANT) //switch back to constant reconstruction (i.e. v)
              err = riemann.ComputeRiemannSolution(dir, v[k-1][j][i], neighborid, v[k][j][i], myid, Vmid, midid, Vsm, Vsp);
            else
              err = riemann.ComputeRiemannSolution(dir, vf[k-1][j][i], neighborid, vk[k][j][i], myid, Vmid, midid, Vsm, Vsp);

            if(err) 
              riemann_errors++;

            //Clip Riemann solution and check it
            varFcn[neighborid]->ClipDensityAndPressure(Vsm);
            varFcn[neighborid]->CheckState(Vsm);
            varFcn[myid]->ClipDensityAndPressure(Vsp);
            varFcn[myid]->CheckState(Vsp);


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
        
  
  MPI_Allreduce(MPI_IN_PLACE, &riemann_errors, 1, MPI_INT, MPI_SUM, comm);
  if(riemann_errors>0) 
    print("Warning: Riemann solver failed to find a bracketing interval on %d edges.\n", riemann_errors);


  //------------------------------------
  // Restore Spatial Variables
  //------------------------------------

  if(Phi && ls_mat_id) {
    for(int i=0; i<Phi->size(); i++)
      (*Phi)[i]->RestoreDataPointerToLocalVector();
  }

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

Vec3D
SpaceOperator::CalculateGradPhiAtCellInterface(int d/*0,1,2*/, int i, int j, int k, Vec3D*** coords, Vec3D*** dxyz,
                                               int myid, int neighborid, vector<int> *ls_mat_id,
                                               vector<double***> *phi)
{

  int my_ls    = (myid==0) ?       neighborid : myid;
  int neigh_ls = (neighborid==0) ? myid       : neighborid;
  bool found1(false), found2(false);
  for(int s=0; s<ls_mat_id->size(); s++) {
    if(!found1 && (*ls_mat_id)[s] == my_ls) {
      my_ls = s;
      found1 = true;
    }
    if(!found2 && (*ls_mat_id)[s] == neigh_ls) {
      neigh_ls = s;
      found2 = true;
    }
    if(found1 && found2)
      break;
  }
  assert(found1 && found2);

  Vec3D dir;

  if(my_ls == neigh_ls) {// one of the two has matid = 0. only 1 level set function involved

    dir = CalculateGradientAtCellInterface(d,i,j,k,coords,dxyz,(*phi)[my_ls]);

    if(neighborid==0)
      dir *= -1.0;

    if(dir.norm()==0.0) //for whatever reason...
      dir[d] = 1.0; //mesh normal
  } 
  else {// 2 level set functions are involved here

    Vec3D nphi_1 = CalculateGradientAtCellInterface(d,i,j,k,coords,dxyz,(*phi)[my_ls]);
    Vec3D nphi_2 = CalculateGradientAtCellInterface(d,i,j,k,coords,dxyz,(*phi)[neigh_ls]);

    dir = nphi_2 - nphi_1;

    double dir_norm = dir.norm();
    if(dir_norm!=0.0)
      dir /= dir_norm;
    else
      dir[d] = 1.0;
  }

  if(dir[d]<=0.0) {
    fprintf(stderr,"*** Error: (%d,%d,%d)(%d): dir = %e %e %e, myid = %d, neighid = %d.\n", 
            i,j,k, d, dir[0], dir[1], dir[2], myid, neighborid);
    exit(-1);
  }

  return dir;

}

//-----------------------------------------------------

Vec3D
SpaceOperator::CalculateGradientAtCellInterface(int d/*0,1,2*/, int i, int j, int k, Vec3D*** coords,
                                                Vec3D*** dxyz, double*** phi)
{

  Vec3D dir(0.0, 0.0, 0.0); 

  if(d == 0) {
  // calculate the dPhi/dx
    
    dir[0] = (phi[k][j][i] - phi[k][j][i-1])/(coords[k][j][i][0] - coords[k][j][i-1][0]);

    double dy1, dy2;
    if(j<0) {//[j-1] is not available
      dy1 = (phi[k][j+1][i-1] - phi[k][j][i-1])/(coords[k][j+1][i-1][1] - coords[k][j][i-1][1]);
      dy2 = (phi[k][j+1][i]   - phi[k][j][i])  /(coords[k][j+1][i][1]   - coords[k][j][i][1]);
    } else if(j>=NY) {//[j+1] is not available
      dy1 = (phi[k][j][i-1] - phi[k][j-1][i-1])/(coords[k][j][i-1][1] - coords[k][j-1][i-1][1]);
      dy2 = (phi[k][j][i]   - phi[k][j-1][i])  /(coords[k][j][i][1]   - coords[k][j-1][i][1]);
    } else {
      dy1 = CentralDifferenceLocal(   phi[k][j-1][i-1],       phi[k][j][i-1],       phi[k][j+1][i-1],
                                   coords[k][j-1][i-1][1], coords[k][j][i-1][1], coords[k][j+1][i-1][1]);
      dy2 = CentralDifferenceLocal(   phi[k][j-1][i],         phi[k][j][i],         phi[k][j+1][i],
                                   coords[k][j-1][i][1],   coords[k][j][i][1],   coords[k][j+1][i][1]);
    }
    
    double dz1, dz2;
    if(k<0) {//[k-1] is not available
      dz1 = (phi[k+1][j][i-1] - phi[k][j][i-1])/(coords[k+1][j][i-1][2] - coords[k][j][i-1][2]);
      dz2 = (phi[k+1][j][i]   - phi[k][j][i])  /(coords[k+1][j][i][2]   - coords[k][j][i][2]);
    } else if(k>=NZ) {//[k+1] is not available
      dz1 = (phi[k][j][i-1] - phi[k-1][j][i-1])/(coords[k][j][i-1][2] - coords[k-1][j][i-1][2]);
      dz2 = (phi[k][j][i]   - phi[k-1][j][i])  /(coords[k][j][i][2]   - coords[k-1][j][i][2]);
    } else {
      dz1 = CentralDifferenceLocal(   phi[k-1][j][i-1],       phi[k][j][i-1],       phi[k+1][j][i-1],
                                   coords[k-1][j][i-1][2], coords[k][j][i-1][2], coords[k+1][j][i-1][2]);
      dz2 = CentralDifferenceLocal(   phi[k-1][j][i],         phi[k][j][i],         phi[k+1][j][i],
                                   coords[k-1][j][i][2],   coords[k][j][i][2],   coords[k+1][j][i][2]);
    } 

    double x_i_minus_half = coords[k][j][i][0] - 0.5*dxyz[k][j][i][0];
    double dx = coords[k][j][i][0] - coords[k][j][i-1][0];
    double c1 = (coords[k][j][i][0] - x_i_minus_half)/dx;
    double c2 = (x_i_minus_half - coords[k][j][i-1][0])/dx;

    dir[1] = c1*dy1 + c2*dy2;
    dir[2] = c1*dz1 + c2*dz2;

  }
  else if(d == 1) {
  // calculate the dPhi/dy
    
    dir[1] = (phi[k][j][i] - phi[k][j-1][i])/(coords[k][j][i][1] - coords[k][j-1][i][1]);

    double dx1, dx2;
    if(i<0) {//[i-1] is not available
      dx1 = (phi[k][j-1][i+1] - phi[k][j-1][i])/(coords[k][j-1][i+1][0] - coords[k][j-1][i][0]);
      dx2 = (phi[k][j][i+1]   - phi[k][j][i])  /(coords[k][j][i+1][0]   - coords[k][j][i][0]);
    } else if(i>=NX) {//[i+1] is not available
      dx1 = (phi[k][j-1][i] - phi[k][j-1][i-1])/(coords[k][j-1][i][0] - coords[k][j-1][i-1][0]);
      dx2 = (phi[k][j][i]   - phi[k][j][i-1])  /(coords[k][j][i][0]   - coords[k][j][i-1][0]);
    } else {
      dx1 = CentralDifferenceLocal(   phi[k][j-1][i-1],       phi[k][j-1][i],       phi[k][j-1][i+1],
                                   coords[k][j-1][i-1][0], coords[k][j-1][i][0], coords[k][j-1][i+1][0]);
      dx2 = CentralDifferenceLocal(   phi[k][j][i-1],         phi[k][j][i],         phi[k][j][i+1],
                                   coords[k][j][i-1][0],   coords[k][j][i][0],   coords[k][j][i+1][0]);
    }
    
    double dz1, dz2;
    if(k<0) {//[k-1] is not available
      dz1 = (phi[k+1][j-1][i] - phi[k][j-1][i])/(coords[k+1][j-1][i][2] - coords[k][j-1][i][2]);
      dz2 = (phi[k+1][j][i]   - phi[k][j][i])  /(coords[k+1][j][i][2]   - coords[k][j][i][2]);
    } else if(k>=NZ) {//[k+1] is not available
      dz1 = (phi[k][j-1][i] - phi[k-1][j-1][i])/(coords[k][j-1][i][2] - coords[k-1][j-1][i][2]);
      dz2 = (phi[k][j][i]   - phi[k-1][j][i])  /(coords[k][j][i][2]   - coords[k-1][j][i][2]);
    } else {
      dz1 = CentralDifferenceLocal(   phi[k-1][j-1][i],       phi[k][j-1][i],       phi[k+1][j-1][i],
                                   coords[k-1][j-1][i][2], coords[k][j-1][i][2], coords[k+1][j-1][i][2]);
      dz2 = CentralDifferenceLocal(   phi[k-1][j][i],         phi[k][j][i],         phi[k+1][j][i],
                                   coords[k-1][j][i][2],   coords[k][j][i][2],   coords[k+1][j][i][2]);
    } 

    double y_j_minus_half = coords[k][j][i][1] - 0.5*dxyz[k][j][i][1];
    double dy = coords[k][j][i][1] - coords[k][j-1][i][1];
    double c1 = (coords[k][j][i][1] - y_j_minus_half)/dy;
    double c2 = (y_j_minus_half - coords[k][j-1][i][1])/dy;

    dir[0] = c1*dx1 + c2*dx2;
    dir[2] = c1*dz1 + c2*dz2;

  }
  else if(d==2) {
  // calculate the dPhi/dz

    dir[2] = (phi[k][j][i] - phi[k-1][j][i])/(coords[k][j][i][2] - coords[k-1][j][i][2]);

    double dx1, dx2;
    if(i<0) {//[i-1] is not available
      dx1 = (phi[k-1][j][i+1] - phi[k-1][j][i])/(coords[k-1][j][i+1][0] - coords[k-1][j][i][0]);
      dx2 = (phi[k][j][i+1]   - phi[k][j][i])  /(coords[k][j][i+1][0]   - coords[k][j][i][0]);
    } else if(i>=NX) {//[i+1] is not available
      dx1 = (phi[k-1][j][i] - phi[k-1][j][i-1])/(coords[k-1][j][i][0] - coords[k-1][j][i-1][0]);
      dx2 = (phi[k][j][i]   - phi[k][j][i-1])  /(coords[k][j][i][0]   - coords[k][j][i-1][0]);
    } else {
      dx1 = CentralDifferenceLocal(   phi[k-1][j][i-1],       phi[k-1][j][i],       phi[k-1][j][i+1],
                                   coords[k-1][j][i-1][0], coords[k-1][j][i][0], coords[k-1][j][i+1][0]);
      dx2 = CentralDifferenceLocal(   phi[k][j][i-1],         phi[k][j][i],         phi[k][j][i+1],
                                   coords[k][j][i-1][0],   coords[k][j][i][0],   coords[k][j][i+1][0]);
    }
    
    double dy1, dy2;
    if(j<0) {//[j-1] is not available
      dy1 = (phi[k-1][j+1][i] - phi[k-1][j][i])/(coords[k-1][j+1][i][1] - coords[k-1][j][i][1]);
      dy2 = (phi[k][j+1][i]   - phi[k][j][i])  /(coords[k][j+1][i][1]   - coords[k][j][i][1]);
    } else if(j>=NY) {//[j+1] is not available
      dy1 = (phi[k-1][j][i] - phi[k-1][j-1][i])/(coords[k-1][j][i][1] - coords[k-1][j-1][i][1]);
      dy2 = (phi[k][j][i]   - phi[k][j-1][i])  /(coords[k][j][i][1]   - coords[k][j-1][i][1]);
    } else {
      dy1 = CentralDifferenceLocal(   phi[k-1][j-1][i],       phi[k-1][j][i],       phi[k-1][j+1][i],
                                   coords[k-1][j-1][i][1], coords[k-1][j][i][1], coords[k-1][j+1][i][1]);
      dy2 = CentralDifferenceLocal(   phi[k][j-1][i],         phi[k][j][i],         phi[k][j+1][i],
                                   coords[k][j-1][i][1],   coords[k][j][i][1],   coords[k][j+1][i][1]);
    } 

    double z_k_minus_half = coords[k][j][i][2] - 0.5*dxyz[k][j][i][2];
    double dz = coords[k][j][i][2] - coords[k-1][j][i][2];
    double c1 = (coords[k][j][i][2] - z_k_minus_half)/dz;
    double c2 = (z_k_minus_half - coords[k-1][j][i][2])/dz;

    dir[0] = c1*dx1 + c2*dx2;
    dir[1] = c1*dy1 + c2*dy2;

  }

  // normalize the gradient
  double dir_norm = dir.norm();
  if(dir_norm!=0.0)
    dir /= dir_norm;
 
  return dir;

}

//-----------------------------------------------------

void 
SpaceOperator::CheckReconstructedStates(SpaceVariable3D &V,
                                        SpaceVariable3D &Vl, SpaceVariable3D &Vr, SpaceVariable3D &Vb,
                                        SpaceVariable3D &Vt, SpaceVariable3D &Vk, SpaceVariable3D &Vf,
                                        SpaceVariable3D &ID)
{

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
          fprintf(stderr, "\033[0;31m*** Error: Reconstructed state at (%d,%d,%d) violates hyperbolicity. matid = %d.\033[0m\n", i,j,k, myid);
          //fprintf(stderr, "v[%d,%d,%d]  = [%e, %e, %e, %e, %e]\n", i,j,k, v[k][j][i][0], v[k][j][i][1], v[k][j][i][2], v[k][j][i][3], v[k][j][i][4]);
          fprintf(stderr, "vl[%d,%d,%d] = [%e, %e, %e, %e, %e]\n", i,j,k, vl[k][j][i][0], vl[k][j][i][1], vl[k][j][i][2], vl[k][j][i][3], vl[k][j][i][4]);
          fprintf(stderr, "vr[%d,%d,%d] = [%e, %e, %e, %e, %e]\n", i,j,k, vr[k][j][i][0], vr[k][j][i][1], vr[k][j][i][2], vr[k][j][i][3], vr[k][j][i][4]);
          fprintf(stderr, "vb[%d,%d,%d] = [%e, %e, %e, %e, %e]\n", i,j,k, vb[k][j][i][0], vb[k][j][i][1], vb[k][j][i][2], vb[k][j][i][3], vb[k][j][i][4]);
          fprintf(stderr, "vt[%d,%d,%d] = [%e, %e, %e, %e, %e]\n", i,j,k, vt[k][j][i][0], vt[k][j][i][1], vt[k][j][i][2], vt[k][j][i][3], vt[k][j][i][4]);
          fprintf(stderr, "vk[%d,%d,%d] = [%e, %e, %e, %e, %e]\n", i,j,k, vk[k][j][i][0], vk[k][j][i][1], vk[k][j][i][2], vk[k][j][i][3], vk[k][j][i][4]);
          fprintf(stderr, "vf[%d,%d,%d] = [%e, %e, %e, %e, %e]\n", i,j,k, vf[k][j][i][0], vf[k][j][i][1], vf[k][j][i][2], vf[k][j][i][3], vf[k][j][i][4]);
          exit(-1);
        } 
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &nClipped, 1, MPI_INT, MPI_SUM, comm);
  if(nClipped && verbose>0)
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
                                    RiemannSolutions *riemann_solutions, vector<int> *ls_mat_id, vector<SpaceVariable3D*> *Phi)
{

#ifdef LEVELSET_TEST
  return; //testing the level set solver without solving the N-S / Euler equations
#endif

  // -------------------------------------------------
  // calculate fluxes on the left hand side of the equation   
  // -------------------------------------------------
  ComputeAdvectionFluxes(V, ID, R, riemann_solutions, ls_mat_id, Phi);

  if(visco)
    visco->AddDiffusionFluxes(V, ID, R);

  if(symm) //cylindrical or spherical symmetry
    symm->AddSymmetryTerms(V, ID, R); //These terms are placed on the left-hand-side

  // -------------------------------------------------
  // multiply flux by -1, and divide by cell volume (for cells within the actual domain)
  // -------------------------------------------------
  Vec5D***    r = (Vec5D***) R.GetDataPointer();
  double*** vol = (double***)volume.GetDataPointer();

  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++) 
      for(int i=i0; i<imax; i++) {
        r[k][j][i] /= -vol[k][j][i];
      }

  // restore spatial variables
  R.RestoreDataPointerToLocalVector(); //NOTE: although R has been updated, there is no need of 
                                       //      cross-subdomain communications. So, no need to 
                                       //      update the global vec.


  

  volume.RestoreDataPointerToLocalVector();
}


//-----------------------------------------------------


