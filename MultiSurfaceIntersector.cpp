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
                         global_mesh(global_mesh_), joint_intersector(NULL), joint_surface(), iod_surface_dummy(),
                         found_multi_surf_intersection(false)
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
  if(surf1 == surf2) {
    print_error("*** Error: Self-intersection should be specified under EmbeddedSurface.\n");
    exit_mpi();
  }

  surface_id.push_back(surf1);
  surface.push_back(&surface_[surf1]);
  intersector.push_back(intersector_[surf1]);
  iod_surface_dummy.surface_thickness = 2.0*intersector_[surf1]->GetSurfaceHalfThickness();

  joint_surface = surface_[surf1]; 

  if(intersector_[surf1]->GetSurfaceHalfThickness() != intersector_[surf2]->GetSurfaceHalfThickness()) {
    print_warning("Warning: Surfaces %d and %d (potentiall intersecting) have different thicknesses.\n",
                  surf1, surf2);
    iod_surface_dummy.surface_thickness = std::min(iod_surface_dummy.surface_thickness,
                                                   2.0*intersector_[surf2]->GetSurfaceHalfThickness());
  }

  iod_surface_dummy.allow_self_intersection = EmbeddedSurfaceData::YES;
  joint_surface.Append(surface_[surf2]);

  surface_id.push_back(surf2);
  surface.push_back(&surface_[surf2]);
  intersector.push_back(intersector_[surf2]);

  numSurfaces = surface.size();
  assert(numSurfaces==2);

  intersecting_elems.assign(2, vector<int>());

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
  assert(numSurfaces==2);
  assert(joint_surface.X0.size() == surface[0]->X0.size() + surface[1]->X0.size());
  assert(joint_surface.elems.size() == surface[0]->elems.size() + surface[1]->elems.size());

  for(int i=0; i<(int)surface[0]->X0.size(); i++) {
    joint_surface.X0[i]   = surface[0]->X0[i];
    joint_surface.X[i]    = surface[0]->X[i];
    joint_surface.Udot[i] = surface[0]->Udot[i];
  }
  int N1 = surface[0]->X0.size();
  for(int i=0; i<(int)surface[1]->X0.size(); i++) {
    joint_surface.X0[N1+i]   = surface[1]->X0[i];
    joint_surface.X[N1+i]    = surface[1]->X[i];
    joint_surface.Udot[N1+i] = surface[1]->Udot[i];
  }

  for(int i=0; i<(int)surface[0]->elemNorm.size(); i++)
    joint_surface.elemNorm[i] = surface[0]->elemNorm[i];
  N1 = surface[0]->elemNorm.size();
  for(int i=0; i<(int)surface[1]->elemNorm.size(); i++)
    joint_surface.elemNorm[N1+i] = surface[1]->elemNorm[i];

  for(int i=0; i<(int)surface[0]->elemArea.size(); i++)
    joint_surface.elemArea[i] = surface[0]->elemArea[i];
  N1 = surface[0]->elemArea.size();
  for(int i=0; i<(int)surface[1]->elemArea.size(); i++)
    joint_surface.elemArea[N1+i] = surface[1]->elemArea[i];

}

//-------------------------------------------------------------------------

bool
MultiSurfaceIntersector::CheckSurfaceIntersections()
{
  for(auto&& x : intersecting_elems)
    x.clear();

  bool found1 = CheckSurfaceIntersectionsOneWay(true);
  bool found2 = CheckSurfaceIntersectionsOneWay(false);

  if(found1 || found2) {
    found_multi_surf_intersection = true;
    return true;
  }

  found_multi_surf_intersection = false;
  return false;
}

//-------------------------------------------------------------------------

bool
MultiSurfaceIntersector::CheckSurfaceIntersectionsOneWay(bool surf1_surf2)
{
  //NOTE: This function can be much faster if we just want to know YES or NO intersection
  //      However, we also want to find all the intersecting elements for later use

  //Assuming just two surfaces
  int surf1, surf2;
  if(surf1_surf2) {
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

  int tid;
  for(auto&& eid : surf1_scope_1) {
    Int3 &nod(Es[eid]);
    if(intersector[surf2]->Intersects(Xs[nod[0]], Xs[nod[1]], &tid, true, -1)) {
      intersecting_elems[surf2].push_back(tid);
      intersecting_elems[surf1].push_back(eid);
    }
    if(intersector[surf2]->Intersects(Xs[nod[1]], Xs[nod[2]], &tid, true, -1)) {
      intersecting_elems[surf2].push_back(tid);
      intersecting_elems[surf1].push_back(eid);
    }
    if(intersector[surf2]->Intersects(Xs[nod[2]], Xs[nod[0]], &tid, true, -1)) {
      intersecting_elems[surf2].push_back(tid);
      intersecting_elems[surf1].push_back(eid);
    }
  }

  // Get rid of duplicates and all-gather through MPI
  std::sort(intersecting_elems[surf1].begin(), intersecting_elems[surf1].end());
  intersecting_elems[surf1].erase(std::unique(intersecting_elems[surf1].begin(), intersecting_elems[surf1].end()),
                                  intersecting_elems[surf1].end());
  std::sort(intersecting_elems[surf2].begin(), intersecting_elems[surf2].end());
  intersecting_elems[surf2].erase(std::unique(intersecting_elems[surf2].begin(), intersecting_elems[surf2].end()),
                                  intersecting_elems[surf2].end());

  CommunicationTools::AllGatherVector(comm, intersecting_elems[surf1]);
  CommunicationTools::AllGatherVector(comm, intersecting_elems[surf2]);

  // Remove duplicates again
  std::sort(intersecting_elems[surf1].begin(), intersecting_elems[surf1].end());
  intersecting_elems[surf1].erase(std::unique(intersecting_elems[surf1].begin(), intersecting_elems[surf1].end()),
                                  intersecting_elems[surf1].end());
  std::sort(intersecting_elems[surf2].begin(), intersecting_elems[surf2].end());
  intersecting_elems[surf2].erase(std::unique(intersecting_elems[surf2].begin(), intersecting_elems[surf2].end()),
                                  intersecting_elems[surf2].end());
  
  //intersecting_elems[surf1] and ...[surf2] must be both empty or both non-empty

  if(intersecting_elems[surf1].empty())
    return false;

  return true;
}

//-------------------------------------------------------------------------

int
MultiSurfaceIntersector::FindNewEnclosures()
{
  new_enclosure_color.clear();
  elem_new_status.clear();

  int nNew = FindNewEnclosuresByFloodFill(); // fills new_enclosure_color
  if(nNew>0) {
    UpdateJointSurface();
    FindNewEnclosureBoundary(); // fills elem_new_status
  }
 
  return nNew;
}

//-------------------------------------------------------------------------

int
MultiSurfaceIntersector::UpdateIntersectors()
{
  // Currently, each processor handles the entire surface (can be improved)
  
  if(new_enclosure_color.empty())
    return -1; //nothing to do

  if(ruling_surface_id == -1)
    return -1; //should not drop any intersections, occluded nodes, and change shortest distance
  
  assert(numSurfaces==2);

  int intersector_modified = -1; //!< ID of modified intersector

  if(ruling_surface_id == 0) { //modify intersector[1]
    int N1 = surface[0]->elems.size();
    int N2 = surface[1]->elems.size();
    vector<bool> elem_drop_status(N2, false);
    std::set<int> elem_to_drop;
    for(auto&& status : elem_new_status) {
      assert((int)status.size()-N1 == N2);
      for(int i=N1; i<(int)status.size(); i++) {
        if(status[i]>0 &&
           std::find(intersecting_elems[1].begin(),intersecting_elems[1].end(),i-N1)==intersecting_elems[1].end()) {
          // belong to the new enclosure boundary, and not an intersecting elem
          // (Note: we keep intersecting elems to avoid gaps or 'leaking')
          elem_drop_status[i-N1] = true;
          elem_to_drop.insert(i-N1); 
        }
      }
    }
    if(!elem_to_drop.empty()) {
      intersector_modified = surface_id[1];
      ModifyIntersectionsAndOccludedNodes(1, elem_drop_status, elem_to_drop);
      //tag dropped elements
      std::vector<int>& eTag(surface[1]->elemtag);
      eTag.assign(N2, 0); //default tag = 0
      for(int i=0; i<(int)eTag.size(); i++) {
        if(elem_drop_status[i])
          eTag[i] = 1; //!< means dropped
      }
    } 
  }
  else {
    assert(ruling_surface_id == 1);
    int N1 = surface[0]->elems.size();
    vector<bool> elem_drop_status(N1, false);
    std::set<int> elem_to_drop;
    for(auto&& status : elem_new_status) {
      for(int i=0; i<N1; i++) {
        if(status[i]>0 &&
           std::find(intersecting_elems[0].begin(),intersecting_elems[0].end(),i)==intersecting_elems[0].end()) {
          elem_drop_status[i] = true;
          elem_to_drop.insert(i); 
        }
      }
    }
    if(!elem_to_drop.empty()) {
      intersector_modified = surface_id[0];
      ModifyIntersectionsAndOccludedNodes(0, elem_drop_status, elem_to_drop);
      //tag dropped elements
      std::vector<int>& eTag(surface[0]->elemtag);
      eTag.assign(N1, 0); //default tag = 0
      for(int i=0; i<(int)eTag.size(); i++) {
        if(elem_drop_status[i])
          eTag[i] = 1; //!< means dropped
      }
    }
  }
  
  return intersector_modified;
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
  EBDS_jnt->occluded_ptr->insert(EBDS2->occluded_ptr->begin(), EBDS2->occluded_ptr->end()); //add EBDS2

  // ------------------------------------------------------------
  // Step 3: Merge swept nodes (in subD) and pass to joint_intersector (unnecessary?)
  // ------------------------------------------------------------
  *EBDS_jnt->swept_ptr = *EBDS1->swept_ptr; //copy EBDS1 to EBDS_jnt
  EBDS_jnt->swept_ptr->insert(EBDS2->swept_ptr->begin(), EBDS2->swept_ptr->end()); //add EBDS2

  // ------------------------------------------------------------
  // Step 4: Flood fill using joint_intersector
  // ------------------------------------------------------------
  joint_intersector->FloodFillColors();
  EBDS_jnt = joint_intersector->GetPointerToResults(); //some info (e.g., nRegions) have changed


  // ------------------------------------------------------------
  // Step 5: Search for new enclosures (topological change)
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

  for(auto&& cc : color1_new) {
    CommunicationTools::AllGatherVector(comm, cc);
    //remove duplicates after MPI gathering
    std::sort(cc.begin(), cc.end());
    cc.erase(std::unique(cc.begin(), cc.end()), cc.end());
  }

  for(auto&& cc : color2_new) {
    CommunicationTools::AllGatherVector(comm, cc);
    //remove duplicates after MPI gathering
    std::sort(cc.begin(), cc.end());
    cc.erase(std::unique(cc.begin(), cc.end()), cc.end());
  }


  // ---------------------------------------------------------------------
  // Step 2: Find new enclosures -- Possibility A
  // ---------------------------------------------------------------------
  // Enclosure i (1, 2, ..., nRegions_jnt) obtained from the joint intersector is a new enclosure IF
  // there exists a map c1->(i, ...) and a map c2->(i, ...), where c1 and c2 are two enclosures
  // from the two intersectors respectively, and (i, ...) means c1 is mapped to "i" AND at least
  // one other color (must be enclosure) from the joint intersector.
  for(int i=1; i<=nRegions_jnt; i++) {
    //check enclosure i
    bool found = false;
    for(int j=1; j<=nRegions1; j++) {
      vector<int> &cc(color1_new[j]);
      if(std::find(cc.begin(), cc.end(), -i) != cc.end()) {
        for(auto&& c : cc)
          if(c<0 && c!=-i) {
            found = true;
            break;
          }
      }
      if(found)
        break;
    }
    if(!found)
      continue;

    found = false;
    for(int j=1; j<=nRegions2; j++) {
      vector<int> &cc(color2_new[j]);
      if(std::find(cc.begin(), cc.end(), -i) != cc.end()) {
        for(auto&& c : cc)
          if(c<0 && c!=-i) {
            found = true;
            break;
          }
      }
      if(found)
        break;
    }
      
    if(found)
      new_enclosure_color.push_back(-i); 
  }

  // ---------------------------------------------------------------------
  // Step 3: Find new enclosures -- Possibility B
  // ---------------------------------------------------------------------
  // Enclosure i (1, 2, ..., nRegions_jnt) obtained from the joint intersector is a new enclosure IF
  // there exists a map c1->(i, ...) and a map c2->(i, ...), where c1 and c2 are both colors connected
  // to far-field (i.e., enclosure i belongs to far-field in both intersector 1 & 2
  for(int i=1; i<=nRegions_jnt; i++) {
    //check enclosure i
    bool found = false;
    for(int j=nRegions1+1; j<(int)color1_new.size(); j++) {
      vector<int> &cc(color1_new[j]);
      if(std::find(cc.begin(), cc.end(), -i) != cc.end()) {
        found = true;
        break;
      }
    }
    if(!found)
      continue;

    found = false;
    for(int j=nRegions2+1; j<(int)color2_new.size(); j++) {
      vector<int> &cc(color2_new[j]);
      if(std::find(cc.begin(), cc.end(), -i) != cc.end()) {
        found = true;
        break;
      }
    }
     
    if(found) {
      if(std::find(new_enclosure_color.begin(), new_enclosure_color.end(), -i) 
           == new_enclosure_color.end())
        new_enclosure_color.push_back(-i); //avoid duplicates (won't happen anyway?)
    }
  }


  return new_enclosure_color.size();
}

//-------------------------------------------------------------------------

void
MultiSurfaceIntersector::FindNewEnclosureBoundary()
{
  if(new_enclosure_color.empty())
    return;
  
  elem_new_status.clear(); //status = 0, 1, 2, or 3 (see Intersector.h)
  for(auto&& cnew : new_enclosure_color) { //repeat the same for each new enclosure (TODO: can be more efficient)
    elem_new_status.push_back(vector<int>());
    joint_intersector->FindColorBoundary(cnew, elem_new_status.back());
  }

  for(int i=0; i<(int)elem_new_status[0].size(); i++)
    if(elem_new_status[0][i]>0)
      fprintf(stderr,"status[%d] = %d.\n", i, elem_new_status[0][i]);
    
}

//-------------------------------------------------------------------------

void
MultiSurfaceIntersector::ModifyIntersectionsAndOccludedNodes(int id, vector<bool> elem_drop_status,
                                                             [[maybe_unused]] std::set<int> elem_to_drop)
{
/*
  print("elem_to_drop: ");
  for(auto&& e :elem_to_drop)
    print("%d ", e);
  print("\n");
*/  

  assert(id==0 || id==1); //currently assuming at most two surfaces involved

  unique_ptr<EmbeddedBoundaryDataSet> EBDS = intersector[id]->GetPointerToResults();

  // -------------------------------
  // Step 1: Delete intersections (and LayerTag)
  // -------------------------------
  Vec3D*** xf = (Vec3D***) EBDS->XForward_ptr->GetDataPointer();
  Vec3D*** xb = (Vec3D***) EBDS->XBackward_ptr->GetDataPointer();
  
  vector<IntersectionPoint>& intersections(*EBDS->intersections_ptr);

  int tid_b, tid_f;
  bool double_intersection;
  vector<int> dropped_intersections;
  int special_tag = -314;
  for(int k=kk0_in; k<kkmax_in; k++)
    for(int j=jj0_in; j<jjmax_in; j++)
      for(int i=ii0_in; i<iimax_in; i++) {

        for(int p=0; p<3; p++) {
          // xf and xb must be both -1, or both >=0
          assert((xf[k][j][i][p]>=0) == (xb[k][j][i][p]>=0));
          
          if(xf[k][j][i][p]>=0) {

            double_intersection = (xf[k][j][i][p] != xb[k][j][i][p]);

            //need to get both, xf[k][j][i][p] and xb[k][j][i][p] may be same => 'corrupted' after the first part
            tid_f = intersections[xf[k][j][i][p]].tid;
            tid_b = intersections[xb[k][j][i][p]].tid; 

            if(elem_drop_status[tid_f]) {
              dropped_intersections.push_back(xf[k][j][i][p]);
              intersections[xf[k][j][i][p]].tid = special_tag;
              xf[k][j][i][p] = -1;
            }

            if(elem_drop_status[tid_b]) {
              if(double_intersection) {//otherwise, already inserted
                dropped_intersections.push_back(xb[k][j][i][p]);
                intersections[xb[k][j][i][p]].tid = special_tag;
              }
              xb[k][j][i][p] = -1;
            }
          }
        }

      }

  //re-order intersections (remove the dropped ones)
  vector<int> old2new(intersections.size(),-1);
  vector<int> new2old(intersections.size() - dropped_intersections.size(),-1);
  
  int skipped = 0;
  for(int i=0; i<(int)intersections.size(); i++) {
    if(intersections[i].tid == special_tag) {//to be dropped
      skipped++;
      old2new[i] = special_tag; //should not be used...
    } else {
      old2new[i] = i-skipped;
      new2old[i-skipped] = i;
    }
  }          
  assert(skipped == (int)dropped_intersections.size());
              
  for(int i=0; i<(int)new2old.size(); i++) {
    assert(new2old[i]>=i);
    intersections[i] = intersections[new2old[i]];
  }
  intersections.resize(new2old.size());

  double*** layer = EBDS->LayerTag_ptr->GetDataPointer(); 
  std::set<Int3> &firstLayer(*EBDS->firstLayer_ptr);

  firstLayer.clear(); //reset

  for(int k=kk0_in; k<kkmax_in; k++)
    for(int j=jj0_in; j<jjmax_in; j++)
      for(int i=ii0_in; i<iimax_in; i++) {

        layer[k][j][i] = -1;

        if(xf[k][j][i][0]>=0) {
          xf[k][j][i][0] = old2new[xf[k][j][i][0]];
          xb[k][j][i][0] = old2new[xb[k][j][i][0]];
          layer[k][j][i-1] = layer[k][j][i] = 1;
          if(i-1>=0)
            firstLayer.insert(Int3(i-1,j,k));
          firstLayer.insert(Int3(i,j,k));
        }
        if(xf[k][j][i][1]>=0) {
          xf[k][j][i][1] = old2new[xf[k][j][i][1]];
          xb[k][j][i][1] = old2new[xb[k][j][i][1]];
          layer[k][j-1][i] = layer[k][j][i] = 1;
          if(j-1>=0)
            firstLayer.insert(Int3(i,j-1,k));
          firstLayer.insert(Int3(i,j,k));
        }
        if(xf[k][j][i][2]>=0) {
          xf[k][j][i][2] = old2new[xf[k][j][i][2]];
          xb[k][j][i][2] = old2new[xb[k][j][i][2]];
          layer[k-1][j][i] = layer[k][j][i] = 1;
          if(k-1>=0)
            firstLayer.insert(Int3(i,j,k-1));
          firstLayer.insert(Int3(i,j,k));
        }
      }

  EBDS->XForward_ptr->RestoreDataPointerToLocalVector(); //can NOT sync (because intersections do not sync)
  EBDS->XBackward_ptr->RestoreDataPointerToLocalVector(); 


  // -------------------------------
  // Step 2: Update occluded, firstLayer (and LayerTag)
  // -------------------------------
  double*** occid  = EBDS->OccTriangle_ptr->GetDataPointer(); 
  std::set<Int3>& occluded(*EBDS->occluded_ptr);
  for(auto it = occluded.begin(); it != occluded.end(); ) { //includes internal ghosts (see Intersector.h)
    int tid = occid[(*it)[2]][(*it)[1]][(*it)[0]];
    assert(tid>=0 && tid<(int)elem_drop_status.size());
    if(elem_drop_status[tid]) {//drop this occluded node
      occid[(*it)[2]][(*it)[1]][(*it)[0]] = -1;
      it = occluded.erase(it); //it points to the next element
    } else {
      layer[(*it)[2]][(*it)[1]][(*it)[0]] = 0;
      firstLayer.insert(*it);
      it++;
    }
  }
  EBDS->OccTriangle_ptr->RestoreDataPointerToLocalVector(); //already taken care of internal ghosts

  EBDS->LayerTag_ptr->RestoreDataPointerAndInsert();
}

//-------------------------------------------------------------------------

bool
MultiSurfaceIntersector::FindNewEnclosuresAfterSurfaceUpdate()
{
  elem_new_status.clear();

  int newColor(0);
  if(new_enclosure_color.empty()) {
    int nRegions1(0), nRegions2(0);
    intersector[0]->GetColors(NULL,NULL,NULL,&nRegions1);
    intersector[1]->GetColors(NULL,NULL,NULL,&nRegions2);
    newColor = -(nRegions1 + nRegions2 + 1); //This would be the color for any new enclosure in joint_intersector
  } else
    newColor = new_enclosure_color[0]; //use an existing color

  FindNewEnclosuresByRefill(newColor); // fills new_enclosure_color
  if(!new_enclosure_color.empty()) {
    UpdateJointSurface();
    FindNewEnclosureBoundary(); // fills elem_new_status
  }

  return !new_enclosure_color.empty();
}

//-------------------------------------------------------------------------

void
MultiSurfaceIntersector::FindNewEnclosuresByRefill(int color4new)
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
  EBDS_jnt->occluded_ptr->insert(EBDS2->occluded_ptr->begin(), EBDS2->occluded_ptr->end()); //add EBDS2


  // ------------------------------------------------------------
  // Step 3: Merge swept nodes (in subD) and pass to joint_intersector
  // ------------------------------------------------------------
  *EBDS_jnt->swept_ptr = *EBDS1->swept_ptr; //copy EBDS1 to EBDS_jnt
  EBDS_jnt->swept_ptr->insert(EBDS2->swept_ptr->begin(), EBDS2->swept_ptr->end()); //add EBDS2

  // ------------------------------------------------------------
  // Step 4: Refill using joint_intersector
  // ------------------------------------------------------------
  bool new_enclosure = joint_intersector->RefillAfterSurfaceUpdate(color4new);

  // ------------------------------------------------------------
  // Step 5: Update new_enclosure_color
  // ------------------------------------------------------------
  if(new_enclosure) {
    if(std::find(new_enclosure_color.begin(), new_enclosure_color.end(), color4new) == new_enclosure_color.end())
      new_enclosure_color.push_back(color4new);
  }

}

//-------------------------------------------------------------------------


//-------------------------------------------------------------------------



//-------------------------------------------------------------------------


//-------------------------------------------------------------------------

