/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<MultiSurfaceIntersector.h>

using std::vector;
using std::unique_ptr;

//-------------------------------------------------------------------------

MultiSurfaceIntersector::MultiSurfaceIntersector(MPI_Comm &comm_, DataManagers3D &dms_,
                                                 SurfaceIntersectionData &iod_surfX_,
                                                 SpaceVariable3D &coordinates_,
                                                 vector<TriangulatedSurface> &surface_,
                                                 vector<Intersector*> &intersector_,
                                                 vector<GhostPoint> &ghost_nodes_inner_,
                                                 vector<GhostPoint> &ghost_nodes_outer_,
                                                 GlobalMeshInfo &global_mesh_)
                       : comm(comm_), coordinates(coordinates_),
                         ghost_nodes_inner(ghost_nodes_inner_), ghost_nodes_outer(ghost_nodes_outer_),
                         global_mesh(global_mesh_), joint_intersector(NULL), surface_dummy(), iod_surface_dummy()
{

  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);
  coordinates.GetInternalGhostedCornerIndices(&ii0_in, &jj0_in, &kk0_in, &iimax_in, &jjmax_in, &kkmax_in);
  coordinates.GetGlobalSize(&NX, &NY, &NZ);

  //TODO: Currently, limited to two surfaces
  
  int surf1 = iod_surfX_.surface1_id;
  int surf2 = iod_surfX_.surface2_id;
  
  //sanity checks
  if(surf1<0 || surf1>=(int)surface_.size()) {
    print_error("*** Error: Surface %d (in SurfaceIntersection) undefined.\n", surf1);
    exit_mpi();
  }
  if(surf2<0 || surf2>=(int)surface_.size()) {
    print_error("*** Error: Surface %d (in SurfaceIntersection) undefined.\n", surf2);
    exit_mpi();
  }
  if((surf1 == surf2) &&
     iod_surfX_.enclosure_treatment != SurfaceIntersectionData::INACTIVE) {
    print_error("*** Error: Self-intersection (in SurfaceIntersection) only supports Enclosure = Inactive.\n");
    exit_mpi();
  }

  surface.push_back(&surface_[surf1]);
  intersector.push_back(intersector_[surf1]);
  elems_active.push_back(vector<bool>(surface_[surf1].elems.size(), true)); //!< default: active
  if(surf2 != surf1) {
    surface.push_back(&surface_[surf2]);
    intersector.push_back(intersector_[surf2]);
    elems_active.push_back(vector<bool>(surface_[surf2].elems.size(), true)); //!< default: active
  }

  numSurfaces = surface.size();

  ruling_surface_id = -1; //default: inactive within overlapped region
  if(iod_surfX_.enclosure_treatment == SurfaceIntersectionData::IGNORE_SURFACE1)
    ruling_surface_id = 1;
  else if(iod_surfX_.enclosure_treatment == SurfaceIntersectionData::IGNORE_SURFACE2)
    ruling_surface_id = 0;
 
  joint_intersector = new Intersector(comm, dms_, iod_surface_dummy, surface_dummy, coordinates,
                                      ghost_nodes_inner, ghost_nodes_outer, global_mesh);

}

//-------------------------------------------------------------------------

MultiSurfaceIntersector::~MultiSurfaceIntersector()
{
  if(joint_intersector) delete joint_intersector;
}

//-------------------------------------------------------------------------

void
MultiSurfaceIntersector::Destroy()
{ }

//-------------------------------------------------------------------------

bool
MultiSurfaceIntersector::CheckSurfaceIntersections()
{
  if(CheckSurfaceIntersectionsOneWay(true) || (numSurfaces==2 && CheckSurfaceIntersectionsOneWay(false)))
    return true;
  return false;
}

//-------------------------------------------------------------------------

bool
MultiSurfaceIntersector::CheckSurfaceIntersectionsOneWay(bool surf1_surf2)
{
  //Assuming just two surfaces
  int surf1, surf2;
  if(numSurfaces==1) {
    surf1 = surf2 = 0;
  } else if(surf1_surf2) {
    surf1 = 0;
    surf2 = 1;
  } else {
    surf1 = 1;
    surf2 = 0;
  }
    
  // check surf1 edges against surf2 elements
  vector<Vec3D>& Xs(surface[surf1]->X);
  vector<Int3>& Es(surface[surf1]->elems);
  vector<int> surf1_scope_1;
  intersector[surf1]->GetTrianglesInScope1(surf1_scope_1);

  // make sure each processor runs the same number of iterations (for MPI exchange).
  // Note: it may be more efficient to use one-sided comm (RMA), especially MPI_Win_lock/unlock.
  //       But I couldn't make it work as expected.
  int Niter = surf1_scope_1.size();
  MPI_Allreduce(MPI_IN_PLACE, &Niter, 1, MPI_INT, MPI_MAX, comm);
  int commFreq = std::max(10, Niter/10);
  int found = 0;
  for(int i=0; i<Niter; i++) {
    if(i<(int)surf1_scope_1.size()) { //otherwise, do nothing
      Int3 &nod(Es[surf1_scope_1[i]]);
      if(intersector[surf2]->Intersects(Xs[nod[0]], Xs[nod[1]]) ||
         intersector[surf2]->Intersects(Xs[nod[1]], Xs[nod[2]]) ||
         intersector[surf2]->Intersects(Xs[nod[2]], Xs[nod[0]]) ) {
        found = 1;
      }
    }

    if(i%commFreq == commFreq-1) {//collect `found'
      MPI_Allreduce(MPI_IN_PLACE, &found, 1, MPI_INT, MPI_MAX, comm);
      if(found>0) //someone found an intersection
        return true;
    }
  }

  //one last communication
  MPI_Allreduce(MPI_IN_PLACE, &found, 1, MPI_INT, MPI_MAX, comm);

  return found != 0;
}

//-------------------------------------------------------------------------

void
MultiSurfaceIntersector::FloodFillWithMergedIntersections()
{
  assert(numSurfaces==2);
  assert(joint_intersector);

  unique_ptr<EmbeddedBoundaryDataSet> EBDS1    = intersector[0]->GetPointerToResults();
  unique_ptr<EmbeddedBoundaryDataSet> EBDS2    = intersector[1]->GetPointerToResults();
  unique_ptr<EmbeddedBoundaryDataSet> EBDS_jnt = joint_intersector->GetPointerToResults();

  // ------------------------------------------------------------
  // Step 1: Merge intersections and pass them to joint_intersector
  // ------------------------------------------------------------
  Vec3D*** xf1    = (Vec3D***)EBDS1->XForward_ptr->GetDataPointer();
  Vec3D*** xf2    = (Vec3D***)EBDS2->XForward_ptr->GetDataPointer();
  Vec3D*** xf_jnt = (Vec3D***)EBDS_jnt->XForward_ptr->GetDataPointer();

  for(int k=kk0_in; k<kkmax_in; k++)
    for(int j=jj0_in; j<jjmax_in; j++)
      for(int i=ii0_in; i<iimax_in; i++)
        for(int p=0; p<3; p++) {
          if(xf1[k][j][i][p]>=0 || xf2[k][j][i][p]>=0)
            xf_jnt[k][j][i][p] = 0; //xf1 & xf2 contain the intersection index; Not used here. (Just set to 0)
          else
            xf_jnt[k][j][i][p] = -1; //non-blocking
        }

  EBDS_jnt->XForward_ptr->RestoreDataPointerToLocalVector(); //!< no need to sync (already filled ghosts)


  // ------------------------------------------------------------
  // Step 2: Merge occluded nodes (in subD) and pass to joint_intersector
  // ------------------------------------------------------------
  *EBDS_jnt->occluded_ptr = *EBDS1->occluded_ptr; //copy EBDS1 to EBDS_jnt
  EBDS_jnt->occluded_ptr->merge(*EBDS2->occluded_ptr); //merge with EBDS2

  // ------------------------------------------------------------
  // Step 3: Flood fill using joint_intersector
  // ------------------------------------------------------------
  int hasOccluded = joint_intersector->FloodFillColors();



  EBDS1->XForward_ptr->RestoreDataPointerToLocalVector(); //should NOT sync! (See notes in Intersector.cpp)
  EBDS2->XForward_ptr->RestoreDataPointerToLocalVector();
}

//-------------------------------------------------------------------------




//-------------------------------------------------------------------------


//-------------------------------------------------------------------------

