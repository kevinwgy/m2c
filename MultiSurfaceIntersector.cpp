/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<MultiSurfaceIntersector.h>
#include<CommunicationTools.h>

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
                         global_mesh(global_mesh_), joint_intersector(NULL), joint_surface(), iod_surface_dummy()
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
  iod_surface_dummy.surface_thickness = 2.0*intersector_[surf1]->GetSurfaceHalfThickness();

  joint_surface = surface_[surf1]; 

  if(surf2 != surf1) {
    if(intersector_[surf1]->GetSurfaceHalfThickness() != intersector_[surf2]->GetSurfaceHalfThickness()) {
      print_warning("Warning: Surfaces %d and %d (potentiall intersecting) have different thicknesses.\n",
                    surf1, surf2);
      iod_surface_dummy.surface_thickness = std::min(iod_surface_dummy.surface_thickness,
                                                     2.0*intersector_[surf2]->GetSurfaceHalfThickness());
      joint_surface.Append(surface_[surf2]);
    }
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
 


  joint_intersector = new Intersector(comm, dms_, iod_surface_dummy, joint_surface, coordinates,
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
{
  if(joint_intersector)
    joint_intersector->Destroy();
}

//-------------------------------------------------------------------------

void
MultiSurfaceIntersector::UpdateJointSurface()
{
  assert(numSurfaces==1 || numSurfaces==2);
  assert(joint_surface.X0.size() == (numSurfaces==1 ? surface[0]->X0.size() :
                                     surface[0]->X0.size() + surface[1]->X0.size()));
  assert(joint_surface.elems.size() == (numSurfaces==1 ? surface[0]->elems.size() :
                                        surface[0]->elems.size() + surface[1]->elems.size()));

  for(int i=0; i<(int)surface[0]->X0.size(); i++) {
    joint_surface.X0[i]   = surface[0]->X0[i];
    joint_surface.X[i]    = surface[0]->X[i];
    joint_surface.Udot[i] = surface[0]->Udot[i];
  }
  if(numSurfaces == 2) {
    int N1 = surface[0]->X0.size();
    for(int i=0; i<(int)surface[1]->X0.size(); i++) {
      joint_surface.X0[N1+i]   = surface[1]->X0[i];
      joint_surface.X[N1+i]    = surface[1]->X[i];
      joint_surface.Udot[N1+i] = surface[1]->Udot[i];
    }
  }

  for(int i=0; i<(int)surface[0]->elemNorm.size(); i++)
    joint_surface.elemNorm[i] = surface[0]->elemNorm[i];
  if(numSurfaces == 2) {
    int N1 = surface[0]->elemNorm.size();
    for(int i=0; i<(int)surface[1]->elemNorm.size(); i++)
      joint_surface.elemNorm[N1+i] = surface[1]->elemNorm[i];
  }

  for(int i=0; i<(int)surface[0]->elemArea.size(); i++)
    joint_surface.elemArea[i] = surface[0]->elemArea[i];
  if(numSurfaces == 2) {
    int N1 = surface[0]->elemArea.size();
    for(int i=0; i<(int)surface[1]->elemArea.size(); i++)
      joint_surface.elemArea[N1+i] = surface[1]->elemArea[i];
  }

}

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

int
MultiSurfaceIntersector::FindNewEnclosuresByFloodFill()
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
  EBDS1->XForward_ptr->RestoreDataPointerToLocalVector(); //should NOT sync! (See notes in Intersector.cpp)
  EBDS2->XForward_ptr->RestoreDataPointerToLocalVector();


  // ------------------------------------------------------------
  // Step 2: Merge occluded nodes (in subD) and pass to joint_intersector
  // ------------------------------------------------------------
  *EBDS_jnt->occluded_ptr = *EBDS1->occluded_ptr; //copy EBDS1 to EBDS_jnt
  EBDS_jnt->occluded_ptr->merge(*EBDS2->occluded_ptr); //merge with EBDS2

  // ------------------------------------------------------------
  // Step 3: Flood fill using joint_intersector
  // ------------------------------------------------------------
  joint_intersector->FloodFillColors();

  // ------------------------------------------------------------
  // Step 4: Search for new enclosures (topological change)
  // ------------------------------------------------------------
  new_enclosure_color.clear();
  if(EBDS_jnt->nRegions == 0 || 
     EBDS_jnt->nRegions == EBDS1->nRegions + EBDS2->nRegions)
    return 0;
  
  double*** color1    = EBDS1->Color_ptr->GetDataPointer();
  double*** color2    = EBDS2->Color_ptr->GetDataPointer();
  double*** color_jnt = EBDS_jnt->Color_ptr->GetDataPointer();
  DetectNewEnclosures(EBDS1->nPossiblePositiveColors, //a constant (3 at the moment) 
                      EBDS1->nRegions, EBDS2->nRegions,
                      EBDS_jnt->nRegions, color1, color2, color_jnt, new_enclosure_color);
  EBDS1->Color_ptr->RestoreDataPointerToLocalVector();
  EBDS2->Color_ptr->RestoreDataPointerToLocalVector();
  EBDS_jnt->Color_ptr->RestoreDataPointerToLocalVector();
  
  return new_enclosure_color.size();
  //SAVE new_enclosure_color, and move back. Then, write function to detect color boundary, delete intersections (re-order), and update occluded nodes;
  //Re-compute shortest distance

}

//-------------------------------------------------------------------------

int
MultiSurfaceIntersector::DetectNewEnclosures(int nPossiblePositiveColors,
                                             int nRegions1, int nRegions2, int nRegions_jnt,
                                             double*** color1, double*** color2, double*** color_jnt,
                                             std::vector<int> &new_enclosure_color)
{
  if(nRegions_jnt == 0) { //trivial
    new_enclosure_color.clear();
    return 0;
  }
  
  // color1 and color2 must be from two different intersectors/surfaces

  // ---------------------------------------------------------------------
  // Step 1: get color1 --> color_jnt and color2 --> color_jnt maps 
  // ---------------------------------------------------------------------
  vector<vector<int> > color1_new(nRegions1 + 1 + nPossiblePositiveColors);
  vector<vector<int> > color2_new(nRegions2 + 1 + nPossiblePositiveColors);
  //order: 0=>color=0(occluded), 1=>color=-1, ... nRegions1=>color=nRegions1, nRegions1+1=>color=1 (inlet),
  //       nRegions1+2=>color=2, ...

  int c1, c2, c;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        c = color_jnt[k][j][i];

        c1 = color1[k][j][i];
        vector<int> &me1(color1_new[c1<=0 ? -c1 : nRegions1+c1]); //see order above
        if(std::find(me1.begin(), me1.end(), c) == me1.end())
          me1.push_back(c); 

        c2 = color2[k][j][i];
        vector<int> &me2(color2_new[c2<=0 ? -c2 : nRegions2+c2]);
        if(std::find(me2.begin(), me2.end(), c) == me2.end())
          me2.push_back(c); 
      }

  for(auto&& cc : color1_new)
    CommunicationTools::AllGatherVector(comm, cc);
  for(auto&& cc : color2_new)
    CommunicationTools::AllGatherVector(comm, cc);

  // ---------------------------------------------------------------------
  // Step 2: Find new enclosures
  // ---------------------------------------------------------------------
  // Enclosure i (1, 2, ..., nRegions_jnt) obtained from the joint intersector is a new enclosure IFF
  // there exists a map c1->(i, ...) and a map c2->(i, ...), where c1 and c2 are two colors (not necessarily
  // enclosures) from the two intersectors respectively, and (i, ...) means c1 is mapped to "i" AND at least
  // one other color (not necessarily enclosure) from the joint intersector.
  for(int i=1; i<=nRegions_jnt; i++) {
    //check enclosure i
    int counter = 0;
    for(auto&& cc : color1_new) {
      if(std::find(cc.begin(), cc.end(), -i) != cc.end()) {
        counter++;
        if(counter>=2)
          break;
      }
    }
    assert(counter>0); //sanity check
    if(counter==1) 
      continue;

    counter = 0;
    for(auto&& cc : color2_new) {
      if(std::find(cc.begin(), cc.end(), -i) != cc.end()) {
        counter++;
        if(counter>=2)
          break;
      }
    }
    assert(counter>0); //sanity check
      
    if(counter>=2)
      new_enclosure_color.push_back(-i); 
  }

  return new_enclosure_color.size();
}

//-------------------------------------------------------------------------

void
MultiSurfaceIntersector::FindNewEnclosureBoundary(vector<vector<int> > &elem_status)
{
  if(new_enclosure_color.empty())
    return;
  
  elem_status.clear(); //status = 0, 1, 2, or 3 (see Intersector.h)
  for(auto&& cnew : new_enclosure_color) { //repeat the same for each new enclosure (TODO: can be more efficient)
    elem_status.push_back(vector<int>());
    joint_intersector->FindColorBoundary(cnew, elem_status.back());
  }
}

//-------------------------------------------------------------------------

void
MultiSurfaceIntersector::UpdateIntersectionsAndOccludedNodes(vector<vector<int> > &elem_status)
{
  //Collect intersections in joint
}

//-------------------------------------------------------------------------

void
MultiSurfaceIntersector::UpdateShortestDistance()
{
  
}

//-------------------------------------------------------------------------







//-------------------------------------------------------------------------


//-------------------------------------------------------------------------

