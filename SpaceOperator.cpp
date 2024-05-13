/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <SpaceOperator.h>
#include <Utils.h>
#include <FluxFcnLLF.h>
#include <Vector3D.h>
#include <Vector5D.h>
#include <GeoTools.h>
#include <DistancePointToSpheroid.h>
#include <DistancePointToParallelepiped.h>
#include <EmbeddedBoundaryDataSet.h>
#include <EmbeddedBoundaryOperator.h>
#include <algorithm> //std::upper_bound
#include <cfloat> //DBL_MAX
#include <KDTree.h>
#include <rbf_interp.hpp>
#include <memory> //unique_ptr
#include <tuple>
using std::cout;
using std::endl;
using std::max;
using std::min;
using std::map;
using std::multimap;
using std::unique_ptr;
using std::tuple;
using std::get;
using namespace GeoTools;

extern int verbose;
extern int INACTIVE_MATERIAL_ID;

//-----------------------------------------------------

SpaceOperator::SpaceOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_,
                             vector<VarFcnBase*> &varFcn_, FluxFcnBase &fluxFcn_,
                             ExactRiemannSolverBase &riemann_,
                             GlobalMeshInfo &global_mesh_,
                             bool screenout) 
  : comm(comm_), dm_all(dm_all_),
    iod(iod_), varFcn(varFcn_), fluxFcn(fluxFcn_), riemann(riemann_),
    interfluxFcn(NULL),
    coordinates(comm_, &(dm_all_.ghosted1_3dof)),
    delta_xyz(comm_, &(dm_all_.ghosted1_3dof)),
    volume(comm_, &(dm_all_.ghosted1_1dof)), global_mesh(global_mesh_),
    rec(comm_, dm_all_, iod_.schemes.ns.rec, coordinates, delta_xyz, &varFcn, &fluxFcn),
    Vl(comm_, &(dm_all_.ghosted1_5dof)),
    Vr(comm_, &(dm_all_.ghosted1_5dof)),
    Vb(comm_, &(dm_all_.ghosted1_5dof)),
    Vt(comm_, &(dm_all_.ghosted1_5dof)),
    Vk(comm_, &(dm_all_.ghosted1_5dof)),
    Vf(comm_, &(dm_all_.ghosted1_5dof)),
    Utmp(comm_, &(dm_all_.ghosted1_5dof)),
    Tag(comm_, &(dm_all_.ghosted1_1dof)),
    symm(NULL), visco(NULL), heat_diffusion(NULL), heo(NULL), smooth(NULL),
    frozen_nodes_ptr(NULL)
{
  
  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);
  coordinates.GetGlobalSize(&NX, &NY, &NZ);

  SetupMesh(global_mesh.x_glob, global_mesh.y_glob, global_mesh.z_glob,
            global_mesh.dx_glob, global_mesh.dy_glob, global_mesh.dz_glob);

  CreateGhostNodeLists(screenout); //create ghost_nodes_inner and ghost_nodes_outer

  rec.Setup(&ghost_nodes_inner, &ghost_nodes_outer); //this function requires mesh info (dxyz)
  
  if(iod.mesh.type == MeshData::SPHERICAL || iod.mesh.type == MeshData::CYLINDRICAL)
    symm = new SymmetryOperator(comm, dm_all, iod.mesh, varFcn, coordinates, delta_xyz, volume);

  if(iod.schemes.ns.smooth.type != SmoothingData::NONE)
    smooth = new SmoothingOperator(comm, dm_all, iod.schemes.ns.smooth, coordinates, delta_xyz, volume);
    
  if(iod.multiphase.flux == MultiPhaseData::LOCAL_LAX_FRIEDRICHS) {
    if(iod.multiphase.phasechange_type == MultiPhaseData::RIEMANN_SOLUTION) {
      print_error("*** Error: MultiPhase::Flux = LocalLaxFriedrichs is incompatible with PhaseChange = RiemannSolution.\n");
      exit_mpi();
    }
    interfluxFcn = new FluxFcnLLF(varFcn, iod);
  }

}

//-----------------------------------------------------

SpaceOperator::~SpaceOperator()
{
  if(symm) delete symm;
  if(visco) delete visco;
  if(heat_diffusion) delete heat_diffusion;
  if(smooth) delete smooth;
  if(interfluxFcn) delete interfluxFcn;
}

//-----------------------------------------------------

void SpaceOperator::Destroy()
{
  rec.Destroy();

  if(symm) symm->Destroy();

  if(visco) visco->Destroy();

  if(heat_diffusion) heat_diffusion->Destroy();

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
  Tag.Destroy();
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
//        fprintf(stdout,"(%d,%d,%d), dx = %e, dy = %e, dz = %e, vol = %e.\n", i,j,k, dxyz[k][j][i][0], dxyz[k][j][i][1], dxyz[k][j][i][2], vol[k][j][i]);
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

          // collect ghost nodes along overset boundaries if any (FACE only)
          if(bcType==MeshData::OVERSET && counter==1)
            ghost_overset.push_back(std::make_pair(Int3(i,j,k), Vec5D(0.0)));

        } 
        else //inside physical domain
          ghost_nodes_inner.push_back(GhostPoint(i,j,k));

      }

  // Find out the owner of each inner ghost
  int mpi_rank(-1);
  MPI_Comm_rank(comm, &mpi_rank);
  Tag.SetConstantValue(mpi_rank, false);
  double *** tag = Tag.GetDataPointer();
  for(auto it = ghost_nodes_inner.begin(); it != ghost_nodes_inner.end(); it++)
    it->owner_proc = (int)tag[it->ijk[2]][it->ijk[1]][it->ijk[0]];

  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++)
        tag[k][j][i] = 0.0; //restore default value (so it can be used for other purposes w/o confusion)
  Tag.RestoreDataPointerToLocalVector();


  int nInner = ghost_nodes_inner.size();
  int nOuter = ghost_nodes_outer.size();
  MPI_Allreduce(MPI_IN_PLACE, &nInner, 1, MPI_INT, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, &nOuter, 1, MPI_INT, MPI_SUM, comm);
  if(screenout) {
    print(comm,"  o Number of ghost nodes inside computational domain (overlapping between subdomains): %d\n",
          nInner);
    print(comm,"  o Number of ghost nodes outside computational domain: %d\n",
          nOuter);
    print(comm,"\n");
  }

/*
  for(int i=0; i<ghost_nodes_outer.size(); i++)
    fprintf(stdout,"Ghost %d: (%d, %d, %d) | Image: (%d, %d, %d) | ProjType = %d | BcType = %d | Proj: (%e, %e, %e), (%e, %e, %e)\n", i, ghost_nodes_outer[i].ijk[0], ghost_nodes_outer[i].ijk[1], ghost_nodes_outer[i].ijk[2], ghost_nodes_outer[i].image_ijk[0], ghost_nodes_outer[i].image_ijk[1], ghost_nodes_outer[i].image_ijk[2], ghost_nodes_outer[i].type_projection, ghost_nodes_outer[i].bcType, ghost_nodes_outer[i].boundary_projection[0], ghost_nodes_outer[i].boundary_projection[1], ghost_nodes_outer[i].boundary_projection[2], ghost_nodes_outer[i].outward_normal[0], ghost_nodes_outer[i].outward_normal[1], ghost_nodes_outer[i].outward_normal[2]);
*/

  // figure out whether the entire domain has overset ghosts...
  int overset_count = ghost_overset.size();
  MPI_Allreduce(MPI_IN_PLACE, &overset_count, 1, MPI_INT, MPI_SUM, comm);
  domain_has_overset = overset_count>0;


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
            fprintf(stdout, "\033[0;31m*** Error: State variables at (%e,%e,%e) violate hyperbolicity." 
                    " matid = %d.\n\033[0m", coords[k][j][i][0],coords[k][j][i][1],coords[k][j][i][2], (int)id[k][j][i]);
            fprintf(stdout, "\033[0;31mv[%d(i),%d(j),%d(k)] = [%e, %e, %e, %e, %e]\n\033[0m", 
                    i,j,k, v[k][j][i][0], v[k][j][i][1], v[k][j][i][2], v[k][j][i][3], v[k][j][i][4]);
            exit(-1);
          }
        }
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &nClipped, 1, MPI_INT, MPI_SUM, comm);
  if(nClipped && verbose>0) {
    print_warning(comm, "Warning: Clipped pressure and/or density in %d cells.\n", nClipped);
    V.RestoreDataPointerAndInsert();
  } else
    V.RestoreDataPointerToLocalVector();

  ID.RestoreDataPointerToLocalVector(); //no changes made

  coordinates.RestoreDataPointerToLocalVector(); //no changes

  return nClipped;
}  

//-----------------------------------------------------
//assign interpolator and gradien calculator (pointers) to the viscosity operator
void SpaceOperator::SetupViscosityOperator(InterpolatorBase *interpolator_, GradientCalculatorBase *grad_,
                                           bool with_embedded_boundary)
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
    visco = new ViscosityOperator(comm, dm_all, iod, varFcn, global_mesh, coordinates, delta_xyz,
                                  volume, *interpolator_, *grad_, with_embedded_boundary);
  }
}

//-----------------------------------------------------

void SpaceOperator::SetupHeatDiffusionOperator(InterpolatorBase *interpolator_, GradientCalculatorBase *grad_)
{
  bool needit = false;
  for(auto it = iod.eqs.materials.dataMap.begin(); it != iod.eqs.materials.dataMap.end(); it++) {
    if(it->second->heat_diffusion.type != HeatDiffusionModelData::NONE) {
      needit = true;
      break; //initialize heat_diffusion if any material has heat diffusion 
    }
  }

  if(needit) {
    assert(interpolator_); //make sure it is not NULL
    assert(grad_);
    heat_diffusion = new HeatDiffusionOperator(comm, dm_all, iod.mesh, iod.eqs, varFcn, coordinates, delta_xyz, volume,
                                               *interpolator_, *grad_);
  }
}

//-----------------------------------------------------
//! Apply boundary conditions by populating the ghost cells.
void SpaceOperator::ApplyBoundaryConditions(SpaceVariable3D &V)
{
  Vec5D*** v = (Vec5D***) V.GetDataPointer();

  for(auto it = ghost_nodes_outer.begin(); it != ghost_nodes_outer.end(); it++) {

    if(it->type_projection != GhostPoint::FACE)
      continue; //corner (edge or vertex) nodes are not populated

    int i(it->ijk[0]), j(it->ijk[1]), k(it->ijk[2]);

    if(it->bcType == (int)MeshData::INLET) {
      v[k][j][i][0] = iod.bc.inlet.density;
      v[k][j][i][1] = iod.bc.inlet.velocity_x;
      v[k][j][i][2] = iod.bc.inlet.velocity_y;
      v[k][j][i][3] = iod.bc.inlet.velocity_z;
      v[k][j][i][4] = iod.bc.inlet.pressure;
    }
    else if(it->bcType == (int)MeshData::INLET2) {
      v[k][j][i][0] = iod.bc.inlet2.density;
      v[k][j][i][1] = iod.bc.inlet2.velocity_x;
      v[k][j][i][2] = iod.bc.inlet2.velocity_y;
      v[k][j][i][3] = iod.bc.inlet2.velocity_z;
      v[k][j][i][4] = iod.bc.inlet2.pressure;
    }
    else if(it->bcType == MeshData::SLIPWALL ||
            it->bcType == MeshData::SYMMETRY) {

      int im_i(it->image_ijk[0]), im_j(it->image_ijk[1]), im_k(it->image_ijk[2]);

      for(int p=0; p<5; p++)
        v[k][j][i][p] = v[im_k][im_j][im_i][p];

      switch (it->side) {
        case GhostPoint::LEFT :
        case GhostPoint::RIGHT :
          v[k][j][i][1] *= -1.0;
          break;
        case GhostPoint::BOTTOM :
        case GhostPoint::TOP :
          v[k][j][i][2] *= -1.0;
          break;
        case GhostPoint::BACK :
        case GhostPoint::FRONT :
          v[k][j][i][3] *= -1.0;
          break;
        default :
          fprintf(stdout,"*** Error: Unknown 'side' tag for GhostPoint (%d).\n", (int)it->side); 
          exit(-1);
      }
    } 
    else if(it->bcType == MeshData::STICKWALL) {
      int im_i(it->image_ijk[0]), im_j(it->image_ijk[1]), im_k(it->image_ijk[2]);
      v[k][j][i][0] =      v[im_k][im_j][im_i][0];
      v[k][j][i][1] = -1.0*v[im_k][im_j][im_i][1];
      v[k][j][i][2] = -1.0*v[im_k][im_j][im_i][2];
      v[k][j][i][3] = -1.0*v[im_k][im_j][im_i][3];
      v[k][j][i][4] =      v[im_k][im_j][im_i][4];
    } 
    else if(it->bcType == MeshData::OUTLET) {
      int im_i(it->image_ijk[0]), im_j(it->image_ijk[1]), im_k(it->image_ijk[2]);
      for(int p=0; p<5; p++)
        v[k][j][i][p] = v[im_k][im_j][im_i][p];
    } 
     else if(it->bcType == MeshData::OVERSET) {

      // nothing to be done here. ghost nodes will be populated below
      
    } else {
      fprintf(stdout,"*** Error: Detected unknown boundary condition type (%d).\n", (int)it->bcType);
      exit(-1);
    }

  }

  // update overset ghosts (if any)
  for(auto&& g : ghost_overset)
    v[g.first[2]][g.first[1]][g.first[0]] = g.second;


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
  
  print(comm, "- Applying a smoothing filter to the N-S state variables. (Iterations: %d)\n", 
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
    if (iod.mesh.bc_x0 == MeshData::INLET    || iod.mesh.bc_x0 == MeshData::INLET2 ||
        iod.mesh.bc_x0 == MeshData::SLIPWALL || iod.mesh.bc_x0 == MeshData::STICKWALL) {

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

        for(int k=kk0; k<kkmax; k++)
          for(int j=jj0; j<jjmax; j++) {

            if(coordinates.OutsidePhysicalDomainAndUnpopulated(ii0,j,k))
              continue; //skip corner nodes

            for(int p=0; p<(int)mydisks.size(); p++) {
              if(IsPointInDisk(coords[k][j][ii0][1], coords[k][j][ii0][2],
                               mydisks[p]->cen_y, mydisks[p]->cen_z, mydisks[p]->radius)){
                v[k][j][ii0][0] = mydisks[p]->state.density;
                v[k][j][ii0][1] = mydisks[p]->state.velocity_x;
                v[k][j][ii0][2] = mydisks[p]->state.velocity_y;
                v[k][j][ii0][3] = mydisks[p]->state.velocity_z;
                v[k][j][ii0][4] = mydisks[p]->state.pressure;
              }
            }

            for(int p=0; p<(int)myrects.size(); p++) {
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
    if (iod.mesh.bc_xmax == MeshData::INLET    || iod.mesh.bc_xmax == MeshData::INLET2 ||
        iod.mesh.bc_xmax == MeshData::SLIPWALL || iod.mesh.bc_xmax == MeshData::STICKWALL) {

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

        for(int k=kk0; k<kkmax; k++)
          for(int j=jj0; j<jjmax; j++) {

            if(coordinates.OutsidePhysicalDomainAndUnpopulated(iimax-1,j,k))
              continue; //skip corner nodes

            for(int p=0; p<(int)mydisks.size(); p++) {
              if(IsPointInDisk(coords[k][j][iimax-1][1], coords[k][j][iimax-1][2],
                               mydisks[p]->cen_y, mydisks[p]->cen_z, mydisks[p]->radius)){
                v[k][j][iimax-1][0] = mydisks[p]->state.density;
                v[k][j][iimax-1][1] = mydisks[p]->state.velocity_x;
                v[k][j][iimax-1][2] = mydisks[p]->state.velocity_y;
                v[k][j][iimax-1][3] = mydisks[p]->state.velocity_z;
                v[k][j][iimax-1][4] = mydisks[p]->state.pressure;
              }
            }

            for(int p=0; p<(int)myrects.size(); p++) {
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
    if (iod.mesh.bc_y0 == MeshData::INLET    || iod.mesh.bc_y0 == MeshData::INLET2 ||
        iod.mesh.bc_y0 == MeshData::SLIPWALL || iod.mesh.bc_y0 == MeshData::STICKWALL) {

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

        for(int k=kk0; k<kkmax; k++)
          for(int i=ii0; i<iimax; i++) {

            if(coordinates.OutsidePhysicalDomainAndUnpopulated(i,jj0,k))
              continue; //skip corner nodes

            for(int p=0; p<(int)mydisks.size(); p++) {
              if(IsPointInDisk(coords[k][jj0][i][2], coords[k][jj0][i][0],
                               mydisks[p]->cen_z, mydisks[p]->cen_x, mydisks[p]->radius)){
                v[k][jj0][i][0] = mydisks[p]->state.density;
                v[k][jj0][i][1] = mydisks[p]->state.velocity_x;
                v[k][jj0][i][2] = mydisks[p]->state.velocity_y;
                v[k][jj0][i][3] = mydisks[p]->state.velocity_z;
                v[k][jj0][i][4] = mydisks[p]->state.pressure;
              }
            }

            for(int p=0; p<(int)myrects.size(); p++) {
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
    if (iod.mesh.bc_ymax == MeshData::INLET    || iod.mesh.bc_ymax == MeshData::INLET2 ||
        iod.mesh.bc_ymax == MeshData::SLIPWALL || iod.mesh.bc_ymax == MeshData::STICKWALL) {

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

        for(int k=kk0; k<kkmax; k++)
          for(int i=ii0; i<iimax; i++) {

            if(coordinates.OutsidePhysicalDomainAndUnpopulated(i,jjmax-1,k))
              continue; //skip corner nodes

            for(int p=0; p<(int)mydisks.size(); p++) {
              if(IsPointInDisk(coords[k][jjmax-1][i][2], coords[k][jjmax-1][i][0],
                               mydisks[p]->cen_z, mydisks[p]->cen_x, mydisks[p]->radius)){
                v[k][jjmax-1][i][0] = mydisks[p]->state.density;
                v[k][jjmax-1][i][1] = mydisks[p]->state.velocity_x;
                v[k][jjmax-1][i][2] = mydisks[p]->state.velocity_y;
                v[k][jjmax-1][i][3] = mydisks[p]->state.velocity_z;
                v[k][jjmax-1][i][4] = mydisks[p]->state.pressure;
              }
            }

            for(int p=0; p<(int)myrects.size(); p++) {
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
    if (iod.mesh.bc_z0 == MeshData::INLET    || iod.mesh.bc_z0 == MeshData::INLET2 ||
        iod.mesh.bc_z0 == MeshData::SLIPWALL || iod.mesh.bc_z0 == MeshData::STICKWALL) {

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

        for(int j=jj0; j<jjmax; j++)
          for(int i=ii0; i<iimax; i++) {

            if(coordinates.OutsidePhysicalDomainAndUnpopulated(i,j,kk0))
              continue; //skip corner nodes

            for(int p=0; p<(int)mydisks.size(); p++) {
              if(IsPointInDisk(coords[kk0][j][i][0], coords[kk0][j][i][1],
                               mydisks[p]->cen_x, mydisks[p]->cen_y, mydisks[p]->radius)){
                v[kk0][j][i][0] = mydisks[p]->state.density;
                v[kk0][j][i][1] = mydisks[p]->state.velocity_x;
                v[kk0][j][i][2] = mydisks[p]->state.velocity_y;
                v[kk0][j][i][3] = mydisks[p]->state.velocity_z;
                v[kk0][j][i][4] = mydisks[p]->state.pressure;
              }
            }

            for(int p=0; p<(int)myrects.size(); p++) {
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
    if (iod.mesh.bc_zmax == MeshData::INLET    || iod.mesh.bc_zmax == MeshData::INLET2 ||
        iod.mesh.bc_zmax == MeshData::SLIPWALL || iod.mesh.bc_zmax == MeshData::STICKWALL) {

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

        for(int j=jj0; j<jjmax; j++)
          for(int i=ii0; i<iimax; i++) {

            if(coordinates.OutsidePhysicalDomainAndUnpopulated(i,j,kkmax-1))
              continue; //skip corner nodes

            for(int p=0; p<(int)mydisks.size(); p++) {
              if(IsPointInDisk(coords[kkmax-1][j][i][0], coords[kkmax-1][j][i][1],
                               mydisks[p]->cen_x, mydisks[p]->cen_y, mydisks[p]->radius)){
                v[kkmax-1][j][i][0] = mydisks[p]->state.density;
                v[kkmax-1][j][i][1] = mydisks[p]->state.velocity_x;
                v[kkmax-1][j][i][2] = mydisks[p]->state.velocity_y;
                v[kkmax-1][j][i][3] = mydisks[p]->state.velocity_z;
                v[kkmax-1][j][i][4] = mydisks[p]->state.pressure;
              }
            }

            for(int p=0; p<(int)myrects.size(); p++) {
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

        myid = id[k][j][i];

        if(myid == INACTIVE_MATERIAL_ID)
          continue;

        for(int p=0; p<5; p++) {
          Vmin[p] = min(Vmin[p], v[k][j][i][p]);
          Vmax[p] = max(Vmax[p], v[k][j][i][p]);
        } 

        c = varFcn[myid]->ComputeSoundSpeedSquare(v[k][j][i][0]/*rho*/, 
                            varFcn[myid]->GetInternalEnergyPerUnitMass(v[k][j][i][0],v[k][j][i][4])/*e*/);

        if(c<0) {
          fprintf(stdout,"*** Error: c^2 (square of sound speed) = %e in SpaceOperator. V = %e, %e, %e, %e, %e, ID = %d.\n",
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

        if(global_mesh.x_glob.size()>1)  dx_over_char_speed_min = min(dx_over_char_speed_min, dxyz[k][j][i][0]/lam_f);
        if(global_mesh.y_glob.size()>1)  dx_over_char_speed_min = min(dx_over_char_speed_min, dxyz[k][j][i][1]/lam_g);
        if(global_mesh.z_glob.size()>1)  dx_over_char_speed_min = min(dx_over_char_speed_min, dxyz[k][j][i][2]/lam_h);
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

void
SpaceOperator::ComputeTimeStepSize(SpaceVariable3D &V, SpaceVariable3D &ID, double &dt, double &cfl,
                                   SpaceVariable3D *LocalDt)
{

  if(LocalDt) { //local time-stepping, dealt with separately.
    assert(iod.ts.timestep<=0.0); //shouldn't have constant time-step size.
    ComputeLocalTimeStepSizes(V,ID,dt,cfl,*LocalDt);
    return;
  }

  double Vmin[5], Vmax[5], cmin, cmax, Machmax, char_speed_max, dx_over_char_speed_min; 
  FindExtremeValuesOfFlowVariables(V, ID, Vmin, Vmax, cmin, cmax, Machmax, char_speed_max, dx_over_char_speed_min);

  if(verbose>1)
    print(comm, "- Maximum values: rho = %e, p = %e, c = %e, Mach = %e, char. speed = %e.\n", 
          Vmax[0], Vmax[4], cmax, Machmax, char_speed_max);

  if(iod.ts.timestep > 0) {
    dt = iod.ts.timestep;
    cfl = dt/dx_over_char_speed_min;
  } else {//apply the CFL number
    cfl = iod.ts.cfl;      
    dt = cfl*dx_over_char_speed_min;
  }

  if(iod.exact_riemann.surface_tension == ExactRiemannSolverData::YES) {
    double dt_surfaceTension = 0.9*ComputeTimeStepSizeSurfaceTension(V, ID); // 0.9 is a safety factor 
    if(dt>dt_surfaceTension) {
      dt = dt_surfaceTension;
      cfl = dt/dx_over_char_speed_min;
    }
  }
}

//-----------------------------------------------------

void
SpaceOperator::ComputeLocalTimeStepSizes(SpaceVariable3D &V, SpaceVariable3D &ID, double &dt, double &cfl,
                                         SpaceVariable3D &LocalDt)
{

  cfl = iod.ts.cfl;

  Vec5D*** v    = (Vec5D***)V.GetDataPointer();
  Vec3D*** dxyz = (Vec3D***)delta_xyz.GetDataPointer();
  double*** id  = ID.GetDataPointer();
  double*** dtl = LocalDt.GetDataPointer(); 

  double cmin(DBL_MAX), dx_over_char_speed_min(DBL_MAX);
  double cmax(-DBL_MAX), Machmax(-DBL_MAX), char_speed_max(-DBL_MAX);
  double Vmin[5], Vmax[5];
  for(int i=0; i<5; i++) {
    Vmin[i] = DBL_MAX; //max. double precision number
    Vmax[i] = -DBL_MAX;
  }

  // Loop through the real domain (excluding the ghost layer)
  double c, mach, lam_f, lam_g, lam_h, dx_over_char_speed_local;
  int myid;
  for(int k=k0; k<kmax; k++) {
    for(int j=j0; j<jmax; j++) {
      for(int i=i0; i<imax; i++) {

        myid = id[k][j][i];

        if(myid == INACTIVE_MATERIAL_ID)
          continue;

        for(int p=0; p<5; p++) {
          Vmin[p] = min(Vmin[p], v[k][j][i][p]);
          Vmax[p] = max(Vmax[p], v[k][j][i][p]);
        }

        c = varFcn[myid]->ComputeSoundSpeedSquare(v[k][j][i][0]/*rho*/,
                            varFcn[myid]->GetInternalEnergyPerUnitMass(v[k][j][i][0],v[k][j][i][4])/*e*/);

        if(c<0) {
          fprintf(stdout,"*** Error: c^2 (square of sound speed) = %e in SpaceOperator. "
                         "V = %e, %e, %e, %e, %e, ID = %d.\n",
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

        dx_over_char_speed_local = DBL_MAX;
        if(global_mesh.x_glob.size()>1)
          dx_over_char_speed_local = min(dx_over_char_speed_local, dxyz[k][j][i][0]/lam_f);
        if(global_mesh.y_glob.size()>1)
          dx_over_char_speed_local = min(dx_over_char_speed_local, dxyz[k][j][i][1]/lam_g);
        if(global_mesh.z_glob.size()>1)
          dx_over_char_speed_local = min(dx_over_char_speed_local, dxyz[k][j][i][2]/lam_h);
        dx_over_char_speed_min = min(dx_over_char_speed_min, dx_over_char_speed_local);

        // calculates local time-step size
        dtl[k][j][i] = cfl*dx_over_char_speed_local;

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

  // global min dt
  dt = cfl*dx_over_char_speed_min; 

  if(verbose>1)
    print(comm, "- Maximum values: rho = %e, p = %e, c = %e, Mach = %e, char. speed = %e.\n",
          Vmax[0], Vmax[4], cmax, Machmax, char_speed_max);

  V.RestoreDataPointerToLocalVector();
  delta_xyz.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();  

  LocalDt.RestoreDataPointerAndInsert();

}

//-----------------------------------------------------

void SpaceOperator::ComputeAdvectionFluxes(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &F,
                                           RiemannSolutions *riemann_solutions, vector<int> *ls_mat_id, 
                                           vector<SpaceVariable3D*> *Phi, vector<SpaceVariable3D*> *KappaPhi,
                                           vector<unique_ptr<EmbeddedBoundaryDataSet> > *EBDS)
{
  //------------------------------------
  // Preparation: Delete previous riemann_solutions
  //------------------------------------
  if(riemann_solutions)
    riemann_solutions->Clear();

  //------------------------------------
  // Reconstruction w/ slope limiters.
  //------------------------------------
  if(TagNodesOutsideConRecDepth(Phi, EBDS, Tag))
    rec.Reconstruct(V, Vl, Vr, Vb, Vt, Vk, Vf, &ID, EBDS, &Tag, false); //false: apply const rec within depth
  else
    rec.Reconstruct(V, Vl, Vr, Vb, Vt, Vk, Vf, &ID, EBDS, NULL); 

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
  vector<double***> phi, kappaPhi;
  if(Phi && ls_mat_id) {
    phi.assign(Phi->size(), NULL);
    for(int i=0; i<(int)phi.size(); i++) {
      phi[i] = (*Phi)[i]->GetDataPointer();
    }
  }
  
  if(iod.exact_riemann.surface_tension == ExactRiemannSolverData::YES) {
    kappaPhi.assign(KappaPhi->size(), NULL);
    for(int i=0; i<(int)kappaPhi.size(); i++) {
      kappaPhi[i] = (*KappaPhi)[i]->GetDataPointer();
    }
  }

  //------------------------------------
  // Extract intersection data
  //------------------------------------
  vector<Vec3D***> xf;  
  vector<Vec3D***> xb;  
  vector<TriangulatedSurface*> surfaces;
  vector<vector<IntersectionPoint>*> intersections;
  if(EBDS) {
    for(auto&& ebds : *EBDS) {
      surfaces.push_back(ebds->surface_ptr);
      intersections.push_back(ebds->intersections_ptr);
      xf.push_back((Vec3D***)ebds->XForward_ptr->GetDataPointer());
      xb.push_back((Vec3D***)ebds->XBackward_ptr->GetDataPointer());
    }
  } 


  //------------------------------------
  // Compute fluxes
  //------------------------------------
  Vec5D localflux1, localflux2;

  // Initialize F to 0
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++) 
      for(int i=ii0; i<iimax; i++) {
        f[k][j][i] = 0.0; //setting f[k][j][i][0] = ... = f[k][j][i][4] = 0.0;
      }


  int myid;
  int neighborid = -1;
  int midid = -1;
  double Vmid[5];
  double Vsm[5], Vsp[5];
  double area = 0.0;
  Int3 ind;

  // Count the riemann solver errors. Note that when this occurs, it could be that
  // (1) the Riemann problem does NOT have a solution, and 
  // (2) the solver returns an "approximate" solution (i.e. code may keep running)
  int riemann_errors = 0; 
  int err;

  bool use_LLF_at_interface = (iod.multiphase.flux == MultiPhaseData::LOCAL_LAX_FRIEDRICHS);
  
  // Loop through the domain interior, and the right, top, and front ghost layers. For each cell, calculate the
  // numerical flux across the left, lower, and back cell boundaries/interfaces
  Vec3D vwallf(0.0), vwallb(0.0), nwallf(0.0), nwallb(0.0);
  for(int k=k0; k<kkmax; k++) {
    for(int j=j0; j<jjmax; j++) {
      for(int i=i0; i<iimax; i++) {

        myid = id[k][j][i];

        //*****************************************
        //calculate flux function F_{i-1/2,j,k}
        //*****************************************
        if(k!=kkmax-1 && j!=jjmax-1) {
 
          neighborid = id[k][j][i-1];

          // See KW's notes for the decision algorithm
          //
          if(EBDS &&
             FindEdgeSurfaceIntersections(0,i,j,k,surfaces,intersections,xf,xb,vwallf,vwallb,nwallf,nwallb)) { 
             //wall

            if(neighborid != INACTIVE_MATERIAL_ID) {
              Vec3D dir = GetNormalForOneSidedRiemann(0, 1, nwallf);
              if(iod.ebm.recon == EmbeddedBoundaryMethodData::CONSTANT) {
                //switch back to constant reconstruction (i.e. v)
                err = riemann.ComputeOneSidedRiemannSolution(dir, v[k][j][i-1], neighborid, vwallf, 
                                                             Vmid, midid, Vsm);
                if(err)
                  riemann_errors++;
                varFcn[neighborid]->ClipDensityAndPressure(Vsm);
                varFcn[neighborid]->CheckState(Vsm);
                fluxFcn.ComputeNumericalFluxAtCellInterface(0/*F*/, v[k][j][i-1]/*Vm*/, Vsm/*Vp*/, neighborid,
                                                            localflux1);
              } 
              else {//linear reconstruction w/ limiter
                err = riemann.ComputeOneSidedRiemannSolution(dir, vr[k][j][i-1], neighborid, vwallf,
                                                             Vmid, midid, Vsm);
                if(err) 
                  riemann_errors++;
                varFcn[neighborid]->ClipDensityAndPressure(Vsm);
                varFcn[neighborid]->CheckState(Vsm);
                fluxFcn.ComputeNumericalFluxAtCellInterface(0/*F*/, vr[k][j][i-1]/*Vm*/, Vsm/*Vp*/, neighborid,
                                                            localflux1);
              }
              //fprintf(stdout,"Vsm = %e %e %e %e %e.\n", Vsm[0], Vsm[1], Vsm[2], Vsm[3], Vsm[4]);
            } else
              localflux1 = 0.0;

            if(myid != INACTIVE_MATERIAL_ID) {
              Vec3D dir = GetNormalForOneSidedRiemann(0,-1, nwallb);
              if(iod.ebm.recon == EmbeddedBoundaryMethodData::CONSTANT) {
                //switch back to constant reconstruction (i.e. v)
                err = riemann.ComputeOneSidedRiemannSolution(dir, v[k][j][i], myid, vwallb, Vmid, midid, Vsp);
                if(err)
                  riemann_errors++;
                varFcn[myid]->ClipDensityAndPressure(Vsp);
                varFcn[myid]->CheckState(Vsp);
                fluxFcn.ComputeNumericalFluxAtCellInterface(0/*F*/, Vsp/*Vm*/, v[k][j][i]/*Vp*/, myid,
                                                            localflux2);
              } 
              else {//linear reconstruction w/ limiter
                err = riemann.ComputeOneSidedRiemannSolution(dir, vl[k][j][i], myid, vwallb, Vmid, midid, Vsp);
                if(err) 
                  riemann_errors++;
                varFcn[myid]->ClipDensityAndPressure(Vsp);
                varFcn[myid]->CheckState(Vsp);
                fluxFcn.ComputeNumericalFluxAtCellInterface(0/*F*/, Vsp/*Vm*/, vl[k][j][i]/*Vp*/, myid,
                                                            localflux2);
              }
              //fprintf(stdout,"Vsp = %e %e %e %e %e.\n", Vsp[0], Vsp[1], Vsp[2], Vsp[3], Vsp[4]);
            } else
              localflux2 = 0.0;
          }

          else if(neighborid!=myid) { //material interface

            if(neighborid!=INACTIVE_MATERIAL_ID && myid!=INACTIVE_MATERIAL_ID) {

              if(use_LLF_at_interface) {
                if(iod.multiphase.recon == MultiPhaseData::CONSTANT)
                  //switch back to constant reconstruction (i.e. v)
                  interfluxFcn->ComputeNumericalFluxAtMaterialInterface(0/*F*/,v[k][j][i-1]/*Vm*/,neighborid,
                                                                        v[k][j][i]/*Vp*/, myid, localflux1);
                else
                  interfluxFcn->ComputeNumericalFluxAtMaterialInterface(0/*F*/,vr[k][j][i-1]/*Vm*/,neighborid,
                                                                        vl[k][j][i]/*Vp*/, myid, localflux1);
                localflux2 = localflux1;
              } 
              else { 
                // determine the axis/direction of the 1D Riemann problem
                Vec3D dir = GetNormalForBimaterialRiemann(0/*i-1/2*/,i,j,k,coords,dxyz,myid,neighborid,
                                                          ls_mat_id,&phi);
                if(iod.exact_riemann.surface_tension == ExactRiemannSolverData::NO) {//no surface tension
                  //Solve 1D Riemann problem
                  if(iod.multiphase.recon == MultiPhaseData::CONSTANT)
                    //switch back to constant reconstruction (i.e. v)
                    err = riemann.ComputeRiemannSolution(dir, v[k][j][i-1], neighborid, v[k][j][i], myid,
                                                         Vmid, midid, Vsm, Vsp);
                  else//linear reconstruction w/ limitor
                    err = riemann.ComputeRiemannSolution(dir, vr[k][j][i-1], neighborid, vl[k][j][i], myid,
                                                         Vmid, midid, Vsm, Vsp);
                }
                else {//with surface tension
		  double curvature = CalculateCurvatureAtCellInterface(0, phi[0], kappaPhi[0], i, j, k);
                  //Solve 1D Riemann problem
                  if(iod.multiphase.recon == MultiPhaseData::CONSTANT)
                    //switch back to constant reconstruction (i.e. v)
                    err = riemann.ComputeRiemannSolution(dir, v[k][j][i-1], neighborid, v[k][j][i], myid,
                                                         Vmid, midid, Vsm, Vsp, curvature);
                  else//linear reconstruction w/ limitor
                    err = riemann.ComputeRiemannSolution(dir, vr[k][j][i-1], neighborid, vl[k][j][i], myid,
                                                         Vmid, midid, Vsm, Vsp, curvature);
                }
                if(err)
                  riemann_errors++;

                //Clip Riemann solution and check it
                varFcn[neighborid]->ClipDensityAndPressure(Vsm);
                varFcn[neighborid]->CheckState(Vsm);
                varFcn[myid]->ClipDensityAndPressure(Vsp);
                varFcn[myid]->CheckState(Vsp);

                if(riemann_solutions && !err) {//store Riemann solution for "phase-change update" 
                  ind[0] = k; ind[1] = j; ind[2] = i;
                  riemann_solutions->left[ind] = std::make_pair((Vec5D)Vsm, neighborid); 
                  ind[2] = i-1;
                  riemann_solutions->right[ind] = std::make_pair((Vec5D)Vsp, myid); 
                }

                if(iod.multiphase.flux == MultiPhaseData::EXACT) { //Godunov-type flux
                  fluxFcn.EvaluateFluxFunction_F(Vmid, midid, localflux1);
                  localflux2 = localflux1;
                } else {//Numerical flux function
                  if(iod.multiphase.recon == MultiPhaseData::CONSTANT) {
                    //switch back to constant reconstruction (i.e. v)
                    fluxFcn.ComputeNumericalFluxAtCellInterface(0/*F*/, v[k][j][i-1]/*Vm*/, Vsm/*Vp*/, 
                                                                neighborid, localflux1);
                    fluxFcn.ComputeNumericalFluxAtCellInterface(0/*F*/, Vsp/*Vm*/, v[k][j][i]/*Vp*/,
                                                                myid, localflux2);
                  } else {//linear reconstruction w/ limiter
                    fluxFcn.ComputeNumericalFluxAtCellInterface(0/*F*/, vr[k][j][i-1]/*Vm*/, Vsm/*Vp*/,
                                                                neighborid, localflux1);
                    fluxFcn.ComputeNumericalFluxAtCellInterface(0/*F*/, Vsp/*Vm*/, vl[k][j][i]/*Vp*/,
                                                                myid, localflux2);
                  }
                }
              }
            } else { // This should only occur when embedded surface has an issue (skip)
              localflux1 = 0.0;
              localflux2 = 0.0;
            }
          }

          else { //neighborid==myid (i.e. same material)
            if(myid!=INACTIVE_MATERIAL_ID) {
              fluxFcn.ComputeNumericalFluxAtCellInterface(0/*F*/, vr[k][j][i-1]/*Vm*/, vl[k][j][i]/*Vp*/,
                                                          myid, localflux1);
              localflux2 = localflux1;
            } else { // This should only occur when embedded surface has an issue (skip)
              localflux1 = 0.0;
              localflux2 = 0.0;
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

          if(EBDS &&
             FindEdgeSurfaceIntersections(1,i,j,k,surfaces,intersections,xf,xb,vwallf,vwallb,nwallf,nwallb)) { 
            //wall

            if(neighborid != INACTIVE_MATERIAL_ID) {
              Vec3D dir = GetNormalForOneSidedRiemann(1, 1, nwallf);
              if(iod.ebm.recon == EmbeddedBoundaryMethodData::CONSTANT) {
                //switch back to constant reconstruction (i.e. v)
                err = riemann.ComputeOneSidedRiemannSolution(dir, v[k][j-1][i], neighborid, vwallf, Vmid,
                                                             midid, Vsm);
                if(err)
                  riemann_errors++;
                varFcn[neighborid]->ClipDensityAndPressure(Vsm);
                varFcn[neighborid]->CheckState(Vsm);
                fluxFcn.ComputeNumericalFluxAtCellInterface(1/*G*/, v[k][j-1][i]/*Vm*/, Vsm/*Vp*/, neighborid,
                                                            localflux1);
              }
              else {//linear reconstruction w/ limiter
                err = riemann.ComputeOneSidedRiemannSolution(dir, vt[k][j-1][i], neighborid, vwallf, Vmid,
                                                             midid, Vsm);
                if(err)
                  riemann_errors++;
                varFcn[neighborid]->ClipDensityAndPressure(Vsm);
                varFcn[neighborid]->CheckState(Vsm);
                fluxFcn.ComputeNumericalFluxAtCellInterface(1/*G*/, vt[k][j-1][i]/*Vm*/, Vsm/*Vp*/, neighborid,
                                                            localflux1);
              }
            } else
              localflux1 = 0.0;

            if(myid != INACTIVE_MATERIAL_ID) {
              Vec3D dir = GetNormalForOneSidedRiemann(1,-1, nwallb);
              if(iod.ebm.recon == EmbeddedBoundaryMethodData::CONSTANT) {
                //switch back to constant reconstruction (i.e. v)
                err = riemann.ComputeOneSidedRiemannSolution(dir, v[k][j][i], myid, vwallb, Vmid, midid, Vsp);
                if(err)
                  riemann_errors++;
                varFcn[myid]->ClipDensityAndPressure(Vsp);
                varFcn[myid]->CheckState(Vsp);
                fluxFcn.ComputeNumericalFluxAtCellInterface(1/*G*/, Vsp/*Vm*/, v[k][j][i]/*Vp*/, myid,
                                                            localflux2);
              }
              else {//linear reconstruction w/ limiter
                err = riemann.ComputeOneSidedRiemannSolution(dir, vb[k][j][i], myid, vwallb, Vmid, midid, Vsp);
                if(err)
                  riemann_errors++;
                varFcn[myid]->ClipDensityAndPressure(Vsp);
                varFcn[myid]->CheckState(Vsp);
                fluxFcn.ComputeNumericalFluxAtCellInterface(1/*G*/, Vsp/*Vm*/, vb[k][j][i]/*Vp*/, myid,
                                                            localflux2);
              }
            } else
              localflux2 = 0.0;
          }

          else if(neighborid!=myid) { //material interface

            if(neighborid!=INACTIVE_MATERIAL_ID && myid!=INACTIVE_MATERIAL_ID) {

              if(use_LLF_at_interface) {
                if(iod.multiphase.recon == MultiPhaseData::CONSTANT)
                  //switch back to constant reconstruction (i.e. v)
                  interfluxFcn->ComputeNumericalFluxAtMaterialInterface(1/*G*/,v[k][j-1][i]/*Vm*/,neighborid,
                                                                        v[k][j][i]/*Vp*/, myid, localflux1);
                else
                  interfluxFcn->ComputeNumericalFluxAtMaterialInterface(1/*G*/,vt[k][j-1][i]/*Vm*/,neighborid,
                                                                        vb[k][j][i]/*Vp*/, myid, localflux1);
                localflux2 = localflux1;
              }
              else {
                // determine the axis/direction of the 1D Riemann problem
                Vec3D dir = GetNormalForBimaterialRiemann(1/*j-1/2*/,i,j,k,coords,dxyz,myid,neighborid,
                                                          ls_mat_id,&phi);
                if(iod.exact_riemann.surface_tension == ExactRiemannSolverData::NO) {//no surface tension
                  //Solve 1D Riemann problem
                  if(iod.multiphase.recon == MultiPhaseData::CONSTANT)
                    //switch back to constant reconstruction (i.e. v)
                    err = riemann.ComputeRiemannSolution(dir, v[k][j-1][i], neighborid, v[k][j][i], myid,
                                                         Vmid, midid, Vsm, Vsp);
                  else
                    err = riemann.ComputeRiemannSolution(dir, vt[k][j-1][i], neighborid, vb[k][j][i], myid,
                                                         Vmid, midid, Vsm, Vsp);
                }
                else { 
		  double curvature = CalculateCurvatureAtCellInterface(1, phi[0], kappaPhi[0], i, j, k);
                  //Solve 1D Riemann problem
                  if(iod.multiphase.recon == MultiPhaseData::CONSTANT)
                    //switch back to constant reconstruction (i.e. v)
                    err = riemann.ComputeRiemannSolution(dir, v[k][j-1][i], neighborid, v[k][j][i], myid,
                                                         Vmid, midid, Vsm, Vsp, curvature);
                  else
                    err = riemann.ComputeRiemannSolution(dir, vt[k][j-1][i], neighborid, vb[k][j][i], myid,
                                                         Vmid, midid, Vsm, Vsp, curvature);
                }

                if(err)
                  riemann_errors++;

                //Clip Riemann solution and check it
                varFcn[neighborid]->ClipDensityAndPressure(Vsm);
                varFcn[neighborid]->CheckState(Vsm);
                varFcn[myid]->ClipDensityAndPressure(Vsp);
                varFcn[myid]->CheckState(Vsp);

                if(riemann_solutions && !err) {//store Riemann solution for "phase-change update"
                  ind[0] = k; ind[1] = j; ind[2] = i;
                  riemann_solutions->bottom[ind] = std::make_pair((Vec5D)Vsm, neighborid); 
                  ind[1] = j-1;
                  riemann_solutions->top[ind] = std::make_pair((Vec5D)Vsp, myid); 
                }

                if(iod.multiphase.flux == MultiPhaseData::EXACT) { //Godunov-type flux
                  fluxFcn.EvaluateFluxFunction_G(Vmid, midid, localflux1);
                  localflux2 = localflux1;
                } else {//Numerical flux function
                  if(iod.multiphase.recon == MultiPhaseData::CONSTANT) {
                    //switch back to constant reconstruction (i.e. v)
                    fluxFcn.ComputeNumericalFluxAtCellInterface(1/*G*/, v[k][j-1][i]/*Vm*/, Vsm/*Vp*/, 
                                                                neighborid, localflux1);
                    fluxFcn.ComputeNumericalFluxAtCellInterface(1/*G*/, Vsp/*Vm*/, v[k][j][i]/*Vp*/,
                                                                myid, localflux2);
                  } else {
                    fluxFcn.ComputeNumericalFluxAtCellInterface(1/*G*/, vt[k][j-1][i]/*Vm*/, Vsm/*Vp*/,
                                                                neighborid, localflux1);
                    fluxFcn.ComputeNumericalFluxAtCellInterface(1/*G*/, Vsp/*Vm*/, vb[k][j][i]/*Vp*/,
                                                                myid, localflux2);
                  }
                }
              }
            } else { // This should only occur when embedded surface has an issue (skip)
              localflux1 = 0.0;
              localflux2 = 0.0;
            }
          }

          else { //neighborid==myid (i.e. same material)
            if(myid!=INACTIVE_MATERIAL_ID) {
              fluxFcn.ComputeNumericalFluxAtCellInterface(1/*G*/, vt[k][j-1][i]/*Vm*/, vb[k][j][i]/*Vp*/,
                                                          myid, localflux1);
              localflux2 = localflux1;
            } else { // This should only occur when embedded surface has an issue (skip)
              localflux1 = 0.0;
              localflux2 = 0.0;
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

          if(EBDS &&
             FindEdgeSurfaceIntersections(2,i,j,k,surfaces,intersections,xf,xb,vwallf,vwallb,nwallf,nwallb)) { 
            //wall

            if(neighborid != INACTIVE_MATERIAL_ID) {
              Vec3D dir = GetNormalForOneSidedRiemann(2, 1, nwallf);
              if(iod.ebm.recon == EmbeddedBoundaryMethodData::CONSTANT) {
                //switch back to constant reconstruction (i.e. v)
                err = riemann.ComputeOneSidedRiemannSolution(dir, v[k-1][j][i], neighborid, vwallf, Vmid,
                                                             midid, Vsm);
                if(err)
                  riemann_errors++;
                varFcn[neighborid]->ClipDensityAndPressure(Vsm);
                varFcn[neighborid]->CheckState(Vsm);
                fluxFcn.ComputeNumericalFluxAtCellInterface(2/*H*/, v[k-1][j][i]/*Vm*/, Vsm/*Vp*/, neighborid,
                                                            localflux1);
              }
              else {//linear reconstruction w/ limiter
                err = riemann.ComputeOneSidedRiemannSolution(dir, vf[k-1][j][i], neighborid, vwallf, Vmid,
                                                             midid, Vsm);
                if(err)
                  riemann_errors++;
                varFcn[neighborid]->ClipDensityAndPressure(Vsm);
                varFcn[neighborid]->CheckState(Vsm);
                fluxFcn.ComputeNumericalFluxAtCellInterface(2/*H*/, vf[k-1][j][i]/*Vm*/, Vsm/*Vp*/, neighborid,
                                                            localflux1);
              }
            } else
              localflux1 = 0.0;

            if(myid != INACTIVE_MATERIAL_ID) {
              Vec3D dir = GetNormalForOneSidedRiemann(2,-1, nwallb);
              if(iod.ebm.recon == EmbeddedBoundaryMethodData::CONSTANT) {
                //switch back to constant reconstruction (i.e. v)
                err = riemann.ComputeOneSidedRiemannSolution(dir, v[k][j][i], myid, vwallb, Vmid, midid, Vsp);
                if(err)
                  riemann_errors++;
                varFcn[myid]->ClipDensityAndPressure(Vsp);
                varFcn[myid]->CheckState(Vsp);
                fluxFcn.ComputeNumericalFluxAtCellInterface(2/*H*/, Vsp/*Vm*/, v[k][j][i]/*Vp*/, myid,
                                                            localflux2);
              }
              else {//linear reconstruction w/ limiter
                err = riemann.ComputeOneSidedRiemannSolution(dir, vk[k][j][i], myid, vwallb, Vmid, midid, Vsp);
                if(err)
                  riemann_errors++;
                varFcn[myid]->ClipDensityAndPressure(Vsp);
                varFcn[myid]->CheckState(Vsp);
                fluxFcn.ComputeNumericalFluxAtCellInterface(2/*H*/, Vsp/*Vm*/, vk[k][j][i]/*Vp*/, myid,
                                                            localflux2);
              }
            } else
              localflux2 = 0.0;
          }

          else if(neighborid!=myid) { //material interface

            if(neighborid!=INACTIVE_MATERIAL_ID && myid!=INACTIVE_MATERIAL_ID) {

              if(use_LLF_at_interface) {
                if(iod.multiphase.recon == MultiPhaseData::CONSTANT)
                  //switch back to constant reconstruction (i.e. v)
                  interfluxFcn->ComputeNumericalFluxAtMaterialInterface(2/*H*/,v[k-1][j][i]/*Vm*/,neighborid,
                                                                        v[k][j][i]/*Vp*/, myid, localflux1);
                else
                  interfluxFcn->ComputeNumericalFluxAtMaterialInterface(2/*H*/,vf[k-1][j][i]/*Vm*/,neighborid,
                                                                        vk[k][j][i]/*Vp*/, myid, localflux1);
                localflux2 = localflux1;
              }
              else {
                // determine the axis/direction of the 1D Riemann problem
                Vec3D dir = GetNormalForBimaterialRiemann(2/*k-1/2*/,i,j,k,coords,dxyz,myid,neighborid,
                                                          ls_mat_id,&phi);
                if(iod.exact_riemann.surface_tension == ExactRiemannSolverData::NO) {//no surface tension
                  //Solve 1D Riemann problem
                  if(iod.multiphase.recon == MultiPhaseData::CONSTANT)
                    //switch back to constant reconstruction (i.e. v)
                    err = riemann.ComputeRiemannSolution(dir, v[k-1][j][i], neighborid, v[k][j][i], myid,
                                                         Vmid, midid, Vsm, Vsp);
                  else
                    err = riemann.ComputeRiemannSolution(dir, vf[k-1][j][i], neighborid, vk[k][j][i], myid,
                                                         Vmid, midid, Vsm, Vsp);
                }
                else {
		  double curvature = CalculateCurvatureAtCellInterface(2, phi[0], kappaPhi[0], i, j, k);
                  //Solve 1D Riemann problem
                  if(iod.multiphase.recon == MultiPhaseData::CONSTANT)
                    //switch back to constant reconstruction (i.e. v)
                    err = riemann.ComputeRiemannSolution(dir, v[k-1][j][i], neighborid, v[k][j][i], myid,
                                                         Vmid, midid, Vsm, Vsp, curvature);
                  else
                    err = riemann.ComputeRiemannSolution(dir, vf[k-1][j][i], neighborid, vk[k][j][i], myid,
                                                         Vmid, midid, Vsm, Vsp, curvature);
                }

                if(err)
                  riemann_errors++;

                //Clip Riemann solution and check it
                varFcn[neighborid]->ClipDensityAndPressure(Vsm);
                varFcn[neighborid]->CheckState(Vsm);
                varFcn[myid]->ClipDensityAndPressure(Vsp);
                varFcn[myid]->CheckState(Vsp);

                if(riemann_solutions && !err) {//store Riemann solution for "phase-change update"
                  ind[0] = k; ind[1] = j; ind[2] = i;
                  riemann_solutions->back[ind] = std::make_pair((Vec5D)Vsm, neighborid); 
                  ind[0] = k-1;
                  riemann_solutions->front[ind] = std::make_pair((Vec5D)Vsp, myid); 
                }

                if(iod.multiphase.flux == MultiPhaseData::EXACT) { //Godunov-type flux
                  fluxFcn.EvaluateFluxFunction_H(Vmid, midid, localflux1);
                  localflux2 = localflux1;
                } else {//Numerical flux function
                  if(iod.multiphase.recon == MultiPhaseData::CONSTANT) {
                    //switch back to constant reconstruction (i.e. v)
                    fluxFcn.ComputeNumericalFluxAtCellInterface(2/*H*/, v[k-1][j][i]/*Vm*/, Vsm/*Vp*/,
                                                                neighborid, localflux1);
                    fluxFcn.ComputeNumericalFluxAtCellInterface(2/*H*/, Vsp/*Vm*/, v[k][j][i]/*Vp*/,
                                                                myid, localflux2);
                  } else {
                    fluxFcn.ComputeNumericalFluxAtCellInterface(2/*H*/, vf[k-1][j][i]/*Vm*/, Vsm/*Vp*/,
                                                                neighborid, localflux1);
                    fluxFcn.ComputeNumericalFluxAtCellInterface(2/*H*/, Vsp/*Vm*/, vk[k][j][i]/*Vp*/,
                                                                myid, localflux2);
                  }
                }
              }
            } else { // This should only occur when embedded surface has an issue (skip)
              localflux1 = 0.0;
              localflux2 = 0.0;
            }
          }

          else { //neighborid==myid (i.e. same material)
            if(myid!=INACTIVE_MATERIAL_ID) {
              fluxFcn.ComputeNumericalFluxAtCellInterface(2/*H*/, vf[k-1][j][i]/*Vm*/, vk[k][j][i]/*Vp*/,
                                                          myid, localflux1);
              localflux2 = localflux1;
            } else { // This should only occur when embedded surface has an issue (skip)
              localflux1 = 0.0;
              localflux2 = 0.0;
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
    print_warning(comm, "Warning: Riemann solver failed to find a bracketing interval or to "
                  "converge on %d edge(s).\n", riemann_errors);


  //------------------------------------
  // Restore Spatial Variables
  //------------------------------------

  if(Phi && ls_mat_id) {
    for(int i=0; i<(int)Phi->size(); i++) {
      (*Phi)[i]->RestoreDataPointerToLocalVector();
    }
  }

  if(iod.exact_riemann.surface_tension == ExactRiemannSolverData::YES) {
    for(int i=0; i<(int)KappaPhi->size(); i++) {
      (*KappaPhi)[i]->RestoreDataPointerToLocalVector();
    }
  }

  if(xf.size()>0) {
    for(auto&& ebds : *EBDS) {
      ebds->XForward_ptr->RestoreDataPointerToLocalVector();
      ebds->XBackward_ptr->RestoreDataPointerToLocalVector();
    }
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

bool
SpaceOperator::TagNodesOutsideConRecDepth(vector<SpaceVariable3D*> *Phi,
                                          vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS,
                                          SpaceVariable3D &Tag0)
{

  // Check whether the user has specified "conRec" depth(s)
  bool phi_tag   = (Phi && Phi->size()>0 && iod.multiphase.conRec_depth>0);
  bool embed_tag = false;
  for(auto it = iod.ebm.embed_surfaces.surfaces.dataMap.begin();
           it != iod.ebm.embed_surfaces.surfaces.dataMap.end(); it++) {
    if(it->second->conRec_depth>0) {
      if(!EBDS || (int)EBDS->size()<=it->first || (*EBDS)[it->first]->Phi_nLayer<1) {
        print_error(comm, "*** Error: Cannot impose constant reconstruction for cells near Embedded Surface %d.\n",
                    it->first);
        exit_mpi();
      }
      embed_tag = true;
      break;
    }
  }

  if(!(phi_tag || embed_tag))
    return false;

  double*** tag = Tag0.GetDataPointer();


  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++)
        tag[k][j][i] = 1;


  // Tag cells specified based on level set(s) 
  if(phi_tag) {

    double depth = iod.multiphase.conRec_depth;  

    int ls_size = Phi->size();

    vector<double***> phi(ls_size, NULL);
    for(int ls=0; ls<ls_size; ls++)
      phi[ls] = (*Phi)[ls]->GetDataPointer(); 

    // loop through subdomain interior
    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++) {
          for(int ls=0; ls<ls_size; ls++) {
            if(fabs(phi[ls][k][j][i])<depth) {
              tag[k][j][i] = 0;
              break;
            }
          }
        }

    for(int ls=0; ls<ls_size; ls++)
      (*Phi)[ls]->RestoreDataPointerToLocalVector();
  }


  // Tag cells specified based on embedded surface(s)
  if(embed_tag) {

    vector<double***> phi;
    vector<double> depth;
    for(auto it = iod.ebm.embed_surfaces.surfaces.dataMap.begin();
             it != iod.ebm.embed_surfaces.surfaces.dataMap.end(); it++) {
      if(it->second->conRec_depth<=0)
        continue;
      phi.push_back((*EBDS)[it->first]->Phi_ptr->GetDataPointer());
      depth.push_back(it->second->conRec_depth);
    }

    // loop through subdomain interior
    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++) {
          if(tag[k][j][i] == 0)
            continue;
          for(int ls=0; ls<(int)phi.size(); ls++) {
            if(fabs(phi[ls][k][j][i])<depth[ls]) {
              tag[k][j][i] = 0;
              break;
            }
          }
        }

    for(auto it = iod.ebm.embed_surfaces.surfaces.dataMap.begin();
             it != iod.ebm.embed_surfaces.surfaces.dataMap.end(); it++) {
      if(it->second->conRec_depth<=0)
        continue;
      (*EBDS)[it->first]->Phi_ptr->RestoreDataPointerToLocalVector();
    }
  } 


  Tag0.RestoreDataPointerToLocalVector(); //modified, but no need to communicate (Reconstructor only checks
                                          //this tag within subdomain interior.)
  return true;
}

//-----------------------------------------------------
// Note that there is also a function in LevelSetOperator (ComputeNormalDirection) that does similar things
Vec3D
SpaceOperator::GetNormalForBimaterialRiemann(int d/*0,1,2*/, int i, int j, int k, Vec3D*** coords,
                                             Vec3D*** dxyz, int myid, int neighborid, vector<int> *ls_mat_id,
                                             vector<double***> *phi)
{
  Vec3D dir(0.0,0.0,0.0);

  if(iod.multiphase.riemann_normal == MultiPhaseData::MESH) //use mesh-based normal
    dir[d] = 1.0;
  else { //LEVEL_SET or AVERAGE
    dir = CalculateGradPhiAtCellInterface(d/*i-1/2,j-1/2,or k-1/2*/,
                                          i,j,k,coords,dxyz,myid,neighborid,ls_mat_id,phi);
    if(iod.multiphase.riemann_normal == MultiPhaseData::AVERAGE) {
      dir[d] += 1.0;
      dir /= dir.norm();
    }
  }
  return dir;
}

//-----------------------------------------------------

Vec3D
SpaceOperator::CalculateGradPhiAtCellInterface(int d/*0,1,2*/, int i, int j, int k, Vec3D*** coords,
                                               Vec3D*** dxyz,
                                               int myid, int neighborid, vector<int> *ls_mat_id,
                                               vector<double***> *phi)
{

  //If any of them is NULL, something is wrong...
  assert(ls_mat_id && ls_mat_id->size()>0);
  assert(phi);

  int my_ls    = (myid==0) ?       neighborid : myid;
  int neigh_ls = (neighborid==0) ? myid       : neighborid;
  bool found1(false), found2(false);
  for(int s=0; s<(int)ls_mat_id->size(); s++) {
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

  if(dir[d]<0.0) {
    if(verbose>=1)  //if dir[d] is very close to zero (but negative), it might be fine 
      fprintf(stdout,"Warning: (%d,%d,%d)(%d): dir = %e %e %e, myid = %d, neighid = %d.\n", 
              i,j,k, d, dir[0], dir[1], dir[2], myid, neighborid);  
    dir[d] *= -1.0; //make sure it is aligned with the positive grid axis
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
  bool clipped;
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

        if(myid == INACTIVE_MATERIAL_ID)
          continue;

        if(boundary==0) {//interior
          clipped = varFcn[myid]->ClipDensityAndPressure(vl[k][j][i]);
          if(clipped) { //go back to constant reconstruction
            vl[k][j][i] = v[k][j][i];
            nClipped++;
          }
          clipped = varFcn[myid]->ClipDensityAndPressure(vr[k][j][i]);
          if(clipped) {
            vr[k][j][i] = v[k][j][i];
            nClipped++;
          }
          clipped = varFcn[myid]->ClipDensityAndPressure(vb[k][j][i]);
          if(clipped) {
            vb[k][j][i] = v[k][j][i];
            nClipped++;
          }
          clipped = varFcn[myid]->ClipDensityAndPressure(vt[k][j][i]);
          if(clipped) {
            vt[k][j][i] = v[k][j][i];
            nClipped++;
          }
          clipped = varFcn[myid]->ClipDensityAndPressure(vk[k][j][i]);
          if(clipped) {
            vk[k][j][i] = v[k][j][i];
            nClipped++;
          }
          clipped = varFcn[myid]->ClipDensityAndPressure(vf[k][j][i]);
          if(clipped) {
            vf[k][j][i] = v[k][j][i];
            nClipped++;
          }

          error = varFcn[myid]->CheckState(vl[k][j][i]) || varFcn[myid]->CheckState(vr[k][j][i]) || 
                  varFcn[myid]->CheckState(vb[k][j][i]) || varFcn[myid]->CheckState(vt[k][j][i]) ||
                  varFcn[myid]->CheckState(vk[k][j][i]) || varFcn[myid]->CheckState(vf[k][j][i]);
        } 
        else {//boundary face
          if(i==ii0) {
            clipped = varFcn[myid]->ClipDensityAndPressure(vr[k][j][i]);
            if(clipped) {
              vr[k][j][i] = v[k][j][i];
              nClipped++;
            }
            error = error || varFcn[myid]->CheckState(vr[k][j][i]);
          } else if (i==iimax-1) {
            clipped = varFcn[myid]->ClipDensityAndPressure(vl[k][j][i]);
            if(clipped) {
              vl[k][j][i] = v[k][j][i];
              nClipped++;
            }
            error = error || varFcn[myid]->CheckState(vl[k][j][i]);
          } else if (j==jj0) {
            clipped = varFcn[myid]->ClipDensityAndPressure(vt[k][j][i]);
            if(clipped) {
              vt[k][j][i] = v[k][j][i];
              nClipped++;
            }
            error = error || varFcn[myid]->CheckState(vt[k][j][i]);
          } else if (j==jjmax-1) {
            clipped = varFcn[myid]->ClipDensityAndPressure(vb[k][j][i]);
            if(clipped) {
              vb[k][j][i] = v[k][j][i];
              nClipped++;
            }
            error = error || varFcn[myid]->CheckState(vb[k][j][i]);
          } else if (k==kk0) {
            clipped = varFcn[myid]->ClipDensityAndPressure(vf[k][j][i]);
            if(clipped) {
              vf[k][j][i] = v[k][j][i];
              nClipped++;
            }
            error = error || varFcn[myid]->CheckState(vf[k][j][i]);
          } else if (k==kkmax-1) {
            clipped = varFcn[myid]->ClipDensityAndPressure(vk[k][j][i]);
            if(clipped) {
              vk[k][j][i] = v[k][j][i];
              nClipped++;
            }
            error = error || varFcn[myid]->CheckState(vk[k][j][i]);
          } 
        }

        if(error) {
          fprintf(stdout, "\033[0;31m*** Error: Reconstructed state at (%d,%d,%d) violates hyperbolicity. matid = %d.\033[0m\n", i,j,k, myid);
          fprintf(stdout, "v[%d,%d,%d]  = [%e, %e, %e, %e, %e]\n", i,j,k, v[k][j][i][0], v[k][j][i][1], v[k][j][i][2], v[k][j][i][3], v[k][j][i][4]);
          fprintf(stdout, "vl[%d,%d,%d] = [%e, %e, %e, %e, %e]\n", i,j,k, vl[k][j][i][0], vl[k][j][i][1], vl[k][j][i][2], vl[k][j][i][3], vl[k][j][i][4]);
          fprintf(stdout, "vr[%d,%d,%d] = [%e, %e, %e, %e, %e]\n", i,j,k, vr[k][j][i][0], vr[k][j][i][1], vr[k][j][i][2], vr[k][j][i][3], vr[k][j][i][4]);
          fprintf(stdout, "vb[%d,%d,%d] = [%e, %e, %e, %e, %e]\n", i,j,k, vb[k][j][i][0], vb[k][j][i][1], vb[k][j][i][2], vb[k][j][i][3], vb[k][j][i][4]);
          fprintf(stdout, "vt[%d,%d,%d] = [%e, %e, %e, %e, %e]\n", i,j,k, vt[k][j][i][0], vt[k][j][i][1], vt[k][j][i][2], vt[k][j][i][3], vt[k][j][i][4]);
          fprintf(stdout, "vk[%d,%d,%d] = [%e, %e, %e, %e, %e]\n", i,j,k, vk[k][j][i][0], vk[k][j][i][1], vk[k][j][i][2], vk[k][j][i][3], vk[k][j][i][4]);
          fprintf(stdout, "vf[%d,%d,%d] = [%e, %e, %e, %e, %e]\n", i,j,k, vf[k][j][i][0], vf[k][j][i][1], vf[k][j][i][2], vf[k][j][i][3], vf[k][j][i][4]);
          exit(-1);
        } 
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &nClipped, 1, MPI_INT, MPI_SUM, comm);
  if(nClipped && verbose>0)
    print_warning(comm, "Warning: Clipped pressure and/or density in %d reconstructed states.\n", nClipped);
 
  V.RestoreDataPointerToLocalVector(); //no changes made
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
                                    [[maybe_unused]] double time,
                                    RiemannSolutions *riemann_solutions, vector<int> *ls_mat_id, 
                                    vector<SpaceVariable3D*> *Phi, vector<SpaceVariable3D*> *KappaPhi,
                                    vector<unique_ptr<EmbeddedBoundaryDataSet> > *EBDS,
                                    SpaceVariable3D *Xi, bool run_heat)
{

#ifdef LEVELSET_TEST
  return; //testing the level set solver without solving the N-S / Euler equations
#endif

  //Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
  //Vec3D*** dxyz = (Vec3D***)delta_xyz.GetDataPointer();

  // -------------------------------------------------
  // calculate fluxes on the left hand side of the equation   
  // -------------------------------------------------
  ComputeAdvectionFluxes(V, ID, R, riemann_solutions, ls_mat_id, Phi, KappaPhi, EBDS);

  if(visco)
    visco->AddDiffusionFluxes(V, ID, EBDS, R); //including extra terms from cylindrical symmetry

  if(heat_diffusion && run_heat)
    heat_diffusion->AddDiffusionFluxes(V, ID, EBDS, R);
  else if(heat_diffusion)
    print("Skipping heat diffusion. \n");

  if(heo) {
    assert(Xi);
#ifdef HYPERELASTICITY_TEST
    heo->PrescribeVelocityForTesting(V, time);
#endif

    heo->AddHyperelasticityFluxes(V, ID, *Xi, EBDS, R);
  }

  if(symm) {//cylindrical or spherical symmetry (symm only handles the sink terms from advective fluxes)
    if(heat_diffusion && run_heat)
      heat_diffusion->AddSymmetryDiffusionTerms(V, ID, R);
    symm->AddSymmetryTerms(V, ID, R); //These terms are placed on the left-hand-side
  }


  Vec5D***    r = (Vec5D***) R.GetDataPointer();
  double*** vol = (double***)volume.GetDataPointer();

  if(frozen_nodes_ptr) {
    for(auto&& fn : *frozen_nodes_ptr)
      r[fn[2]][fn[1]][fn[0]] = 0.0; //re-set residual to 0 for frozen nodes/cells
  }

  // -------------------------------------------------
  // multiply flux by -1, and divide by cell volume (for cells within the actual domain)
  // -------------------------------------------------
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

void
SpaceOperator::UpdateOversetGhostNodes(SpaceVariable3D &V)
{

  if(!domain_has_overset)
    return; //nothing to be done

  assert(V.NumDOF()==5);

  Vec5D*** v = (Vec5D***) V.GetDataPointer();

  for(auto&& g : ghost_overset)
    g.second = v[g.first[2]][g.first[1]][g.first[0]];

  V.RestoreDataPointerToLocalVector();

}


//-----------------------------------------------------

bool
SpaceOperator::FindEdgeSurfaceIntersections(int dir/*0~x,1~y,2~z*/, int i, int j, int k,
                                            vector<TriangulatedSurface*>& surfaces,
                                            vector<vector<IntersectionPoint>*>& intersections,
                                            vector<Vec3D***>& xf, vector<Vec3D***>& xb, 
                                            Vec3D& vwallf, Vec3D& vwallb, Vec3D& nwallf, Vec3D& nwallb)
{

  // Find the correct intersection: <surf. id,  intersection id,  dist-to-intersection>
  tuple<int, int, double> fwd_intersection(std::make_tuple(-1,-1,DBL_MAX));
  tuple<int, int, double> bwd_intersection(std::make_tuple(-1,-1,-DBL_MAX));
  int xid;
  for(int s=0; s<(int)surfaces.size(); s++) {
    xid = xf[s][k][j][i][dir];
    if(xid>=0) {
      double dist = (*intersections[s])[xid].dist;
      if(dist < get<2>(fwd_intersection)) {
        get<0>(fwd_intersection) = s;
        get<1>(fwd_intersection) = xid;
        get<2>(fwd_intersection) = dist;
      }
    }
    xid = xb[s][k][j][i][dir];
    if(xid>=0) {
      double dist = (*intersections[s])[xid].dist;
      if(dist > get<2>(bwd_intersection)) {
        get<0>(bwd_intersection) = s;
        get<1>(bwd_intersection) = xid;
        get<2>(bwd_intersection) = dist;
      }
    }
  }


  // If intersection does not exist, return false.
  if(get<0>(fwd_intersection)<0) {
    assert(get<0>(bwd_intersection)<0);
    return false;
  } else
    assert(get<0>(bwd_intersection)>=0);


  // Intersection exists. 
  int surf;

  // Find normal direction and normal velocity for "forward intersection"
  surf = get<0>(fwd_intersection);
  xid  = get<1>(fwd_intersection);
  IntersectionPoint& xfp((*intersections[surf])[xid]);
  Int3& nf(surfaces[surf]->elems[xfp.tid]);
  nwallf = surfaces[surf]->elemNorm[xfp.tid];
  vwallf = (xfp.xi[0])*surfaces[surf]->Udot[nf[0]]
         + (xfp.xi[1])*surfaces[surf]->Udot[nf[1]]
         + (xfp.xi[2])*surfaces[surf]->Udot[nf[2]];

  // Find normal direction and normal velocity for "backward intersection"
  surf = get<0>(bwd_intersection);
  xid  = get<1>(bwd_intersection);
  IntersectionPoint& xbp((*intersections[surf])[xid]);
  Int3& nb(surfaces[surf]->elems[xbp.tid]);
  nwallb = surfaces[surf]->elemNorm[xbp.tid];
  vwallb = (xbp.xi[0])*surfaces[surf]->Udot[nb[0]]
         + (xbp.xi[1])*surfaces[surf]->Udot[nb[1]]
         + (xbp.xi[2])*surfaces[surf]->Udot[nb[2]];

  return true;

}

//-----------------------------------------------------

Vec3D
SpaceOperator::GetNormalForOneSidedRiemann(int d, int forward_or_backward, Vec3D& nwall)
{
  Vec3D dir(0.0,0.0,0.0);

  switch (iod.ebm.riemann_normal) {
    case EmbeddedBoundaryMethodData::MESH :
      dir[d] = (double)forward_or_backward;
      break;
    case EmbeddedBoundaryMethodData::EMBEDDED_SURFACE :
      dir = nwall;
      if(dir[d]*(double)forward_or_backward<0.0)
        dir[d] *= -1;
      break;
    case EmbeddedBoundaryMethodData::AVERAGE :
      dir = nwall;
      if(dir[d]*(double)forward_or_backward<0.0)
        dir[d] *= -1;
      dir[d] += (double)forward_or_backward;
      dir /= dir.norm(); 
      break;
    default :
      fprintf(stdout, "\033[0;31m*** Error: Undefined Riemann normal code: %d.\n\033[0m", 
              (int)iod.ebm.riemann_normal);
      exit(-1);
  }

  return dir;
}

//-----------------------------------------------------

double
SpaceOperator::CalculateCurvatureAtCellInterface(int d /*0,1,2*/, double*** phi, double*** kappaPhi,
                                                 int i, int j, int k)
{
  double ans = 0.;
  if (d == 0) { // x-direction
    ans = ( kappaPhi[k][j][i-1]*fabs(phi[k][j][i]) + kappaPhi[k][j][i]*fabs(phi[k][j][i-1]) ) /
          ( fabs(phi[k][j][i-1]) + fabs(phi[k][j][i]) ); 
  } else if(d == 1) { // y-direction
    ans = ( kappaPhi[k][j-1][i]*fabs(phi[k][j][i]) + kappaPhi[k][j][i]*fabs(phi[k][j-1][i]) ) /
          ( fabs(phi[k][j-1][i]) + fabs(phi[k][j][i]) ); 
  } else if(d == 2) { // z-direction
    ans = ( kappaPhi[k-1][j][i]*fabs(phi[k][j][i]) + kappaPhi[k][j][i]*fabs(phi[k-1][j][i]) ) /
          ( fabs(phi[k-1][j][i]) + fabs(phi[k][j][i]) ); 
  }
  return ans;
}

//-----------------------------------------------------

double
SpaceOperator::ComputeTimeStepSizeSurfaceTension(SpaceVariable3D &V, SpaceVariable3D &ID)
{

  double sigma = riemann.GetSurfaceTensionCoefficient(); //assume const surface tension coefficient (TODO)
  if(sigma==0.0)
    return DBL_MAX; //no surface tension

  //------------------------------------
  // Extract data
  //------------------------------------
  Vec5D*** v  = (Vec5D***) V.GetDataPointer();
  double*** id = (double***) ID.GetDataPointer();

  int myid;
  int neighborid = -1;
  double my_dx = 0.;  
  double dt_min_surfaceTension = DBL_MAX;  

  // Loop through the domain interior, and the right, top, and front ghost layers. For each pair of cells 
  // occupied by different fluid materials, calculate the time step size due to surface tension.
  for(int k=k0; k<kkmax; k++) {
    for(int j=j0; j<jjmax; j++) {
      for(int i=i0; i<iimax; i++) {

        myid = id[k][j][i];
        if(myid==INACTIVE_MATERIAL_ID)
          continue;

        my_dx = global_mesh.GetMinDXYZ(Int3(i,j,k));

	//**********************************************
	//cell interfaces perpendicular to x-direction
	//**********************************************
        if(k!=kkmax-1 && j!=jjmax-1) {
          neighborid = id[k][j][i-1];
          if(neighborid!=myid && neighborid!=INACTIVE_MATERIAL_ID) { //active material interface
            double neighbor_dx = global_mesh.GetMinDXYZ(Int3(i-1,j,k));
            double dx = min(my_dx, neighbor_dx);
            double dt_surfaceTension = sqrt( (v[k][j][i][0]+v[k][j][i-1][0]) * dx*dx*dx / (4.*M_PI*sigma ) );
            dt_min_surfaceTension = min(dt_min_surfaceTension, dt_surfaceTension);
          }
        }

        //**********************************************
        //cell interfaces perpendicular to y-direction
        //**********************************************
        if(k!=kkmax-1 && i!=iimax-1) {
          neighborid = id[k][j-1][i];
          if(neighborid!=myid && neighborid!=INACTIVE_MATERIAL_ID) { //active material interface
            double neighbor_dx = global_mesh.GetMinDXYZ(Int3(i,j-1,k));
            double dx = min(my_dx, neighbor_dx);
            double dt_surfaceTension = sqrt( (v[k][j][i][0]+v[k][j-1][i][0]) * dx*dx*dx / (4.*M_PI*sigma ) );
            dt_min_surfaceTension = min(dt_min_surfaceTension, dt_surfaceTension); 
          }
        }

        //**********************************************
        //cell interfaces perpendicular to z-direction
        //**********************************************
        if(j!=jjmax-1 && i!=iimax-1) {
          neighborid = id[k-1][j][i];
          if(neighborid!=myid && neighborid!=INACTIVE_MATERIAL_ID) { //active material interface
            double neighbor_dx = global_mesh.GetMinDXYZ(Int3(i,j,k-1));
            double dx = min(my_dx, neighbor_dx);
            double dt_surfaceTension = sqrt( (v[k][j][i][0]+v[k-1][j][i][0]) * dx*dx*dx / (4.*M_PI*sigma ) );
            dt_min_surfaceTension = min(dt_min_surfaceTension, dt_surfaceTension); 
          }
        }
      }
    }
  }
        
  
  MPI_Allreduce(MPI_IN_PLACE, &dt_min_surfaceTension, 1, MPI_DOUBLE, MPI_MIN, comm);

  //------------------------------------
  // Restore Spatial Variables
  //------------------------------------
  ID.RestoreDataPointerToLocalVector(); //no changes
  V.RestoreDataPointerToLocalVector();

  return dt_min_surfaceTension; 

}


  
 

//-----------------------------------------------------

