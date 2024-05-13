/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<Intersector.h>
#include<GeoTools.h>
#include<Intersections.h>
#include<EmbeddedBoundaryDataSet.h>
using std::pair;
using std::vector;
using std::unique_ptr;

extern int verbose;
extern double domain_diagonal;
//-------------------------------------------------------------------------

Intersector::Intersector(MPI_Comm &comm_, DataManagers3D &dms_, EmbeddedSurfaceData &iod_surface_,
                         TriangulatedSurface &surface_, SpaceVariable3D &coordinates_, 
                         vector<GhostPoint> &ghost_nodes_inner_, vector<GhostPoint> &ghost_nodes_outer_,
                         GlobalMeshInfo &global_mesh_)
           : comm(comm_), iod_surface(iod_surface_), surface(surface_), tree_1(NULL), tree_n(NULL),
             coordinates(coordinates_), 
             ghost_nodes_inner(ghost_nodes_inner_), ghost_nodes_outer(ghost_nodes_outer_),
             global_mesh(global_mesh_),
             BBmin_1(comm_, &(dms_.ghosted1_3dof)),
             BBmax_1(comm_, &(dms_.ghosted1_3dof)),
             BBmin_n(comm_, &(dms_.ghosted1_3dof)),
             BBmax_n(comm_, &(dms_.ghosted1_3dof)),
             TMP(comm_, &(dms_.ghosted1_1dof)),
             TMP2(comm_, &(dms_.ghosted1_1dof)),
             CandidatesIndex_1(comm_, &(dms_.ghosted1_1dof)),
             CandidatesIndex_n(comm_, &(dms_.ghosted1_1dof)),
             ClosestPointIndex(comm_, &(dms_.ghosted1_1dof)),
             XForward(comm_, &(dms_.ghosted1_3dof)),
             XBackward(comm_, &(dms_.ghosted1_3dof)),
             Phi(comm_, &(dms_.ghosted1_1dof)),
             Phi_nLayer(0),
             Color(comm_, &(dms_.ghosted1_1dof)), hasInlet(false), hasOutlet(false), nRegions(0),
             floodfiller(comm_, dms_, ghost_nodes_inner_, ghost_nodes_outer_)
{

  CandidatesIndex_1.SetConstantValue(-1, true);
  CandidatesIndex_n.SetConstantValue(-1, true);
  ClosestPointIndex.SetConstantValue(-1, true);
  XForward.SetConstantValue(-1, true);
  XBackward.SetConstantValue(-1, true);
  Color.SetConstantValue(-1, true);

  half_thickness = 0.5*iod_surface.surface_thickness;
  if(half_thickness==0) {
    print_error("*** Error: Detected an embedded surface with 0 thickness. Even if the physical thickness is 0, \n"
                "           a small numerical tolerance is needed to avoid round-off issues.\n");
    exit_mpi();
  }

  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);
  coordinates.GetInternalGhostedCornerIndices(&ii0_in, &jj0_in, &kk0_in, &iimax_in, &jjmax_in, &kkmax_in);
  coordinates.GetGlobalSize(&NX, &NY, &NZ);

  // Set the capacity of internal vectors, so we don't frequently reallocate memory
  int capacity = (imax-i0)*(jmax-j0)*(kmax-k0)/4; //should be big enough
  intersections.reserve(capacity);
  candidates_1.reserve(capacity*2);
  scope_1.reserve(surface.elems.size());
  candidates_n.reserve(capacity*2);
  scope_n.reserve(surface.elems.size());

  //sanity checks on the triangulated surface
  if(surface.degenerate) {
    print_error("*** Error: Intersector cannot track a degenerate surface.\n");
    exit(-1);
  }
  assert(!(surface.node2node.empty() || 
           surface.node2elem.empty() || 
           surface.elem2elem.empty()));

  closed_surface = surface.CheckSurfaceOrientationAndClosedness();

  //build nodal bounding boxes
  BuildNodalAndSubdomainBoundingBoxes(1, BBmin_1, BBmax_1, subD_bbmin_1, subD_bbmax_1);

  nLayer = -1; //this will be the number of layers in BBmin_n, BBmax_n, etc.

}

//-------------------------------------------------------------------------

Intersector::~Intersector()
{
  if(tree_1)
    delete tree_1;
  if(tree_n)
    delete tree_n;
}

//-------------------------------------------------------------------------

void
Intersector::Destroy()
{
  floodfiller.Destroy();

  BBmin_1.Destroy();
  BBmax_1.Destroy();
  BBmin_n.Destroy();
  BBmax_n.Destroy();
  TMP.Destroy();
  TMP2.Destroy();
  CandidatesIndex_1.Destroy();
  CandidatesIndex_n.Destroy();
  ClosestPointIndex.Destroy();
  XForward.Destroy();
  XBackward.Destroy();
  Phi.Destroy();
  Color.Destroy();
}

//-------------------------------------------------------------------------

unique_ptr<EmbeddedBoundaryDataSet>
Intersector::GetPointerToResults()
{
  unique_ptr<EmbeddedBoundaryDataSet> ebds(new EmbeddedBoundaryDataSet());

  ebds->surface_ptr             = &surface;
  ebds->half_thickness          = half_thickness;
  ebds->XForward_ptr            = &XForward;
  ebds->XBackward_ptr           = &XBackward;
  ebds->Phi_ptr                 = &Phi;
  ebds->Phi_nLayer              = Phi_nLayer;
  ebds->Color_ptr               = &Color;
  ebds->ColorReachesBoundary_ptr= &ColorReachesBoundary;
  ebds->hasInlet                = hasInlet;
  ebds->hasInlet2               = hasInlet2;
  ebds->hasOutlet               = hasOutlet;
  ebds->nRegions                = nRegions;
  ebds->ClosestPointIndex_ptr   = &ClosestPointIndex;
  ebds->closest_points_ptr      = &closest_points;
  ebds->intersections_ptr       = &intersections;
  ebds->occluded_ptr            = &occluded;
  ebds->firstLayer_ptr          = &firstLayer;
  ebds->imposed_occluded_ptr    = &imposed_occluded;
  ebds->swept_ptr               = &swept;

  return ebds;
}

//-------------------------------------------------------------------------

void
Intersector::GetElementsInScope1(std::vector<int> &elems_in_scope)
{
  elems_in_scope.resize(scope_1.size());
  for(int i=0; i<(int)scope_1.size(); i++)
    elems_in_scope[i] = scope_1[i].trId();
}

//-------------------------------------------------------------------------

double
Intersector::TrackSurfaceFullCourse(bool &hasInlet_, bool &hasInlet2_, bool &hasOutlet_,
                                    bool &hasOcc_, int &nRegions_, int phi_layers)
{
  assert(phi_layers>=1);

  // This (smaller) one is for edge-surface intersections
  BuildSubdomainScopeAndKDTree(subD_bbmin_1, subD_bbmax_1, scope_1, &tree_1);

  FindIntersections(); //using scope_1 and tree_1

  bool hasOcc = FloodFillColors();

  // Calculate unsigned distance
  if(nLayer != phi_layers) {
    BuildNodalAndSubdomainBoundingBoxes(phi_layers, BBmin_n, BBmax_n, subD_bbmin_n, subD_bbmax_n); //n layers
    nLayer = phi_layers;
  }
  BuildSubdomainScopeAndKDTree(subD_bbmin_n, subD_bbmax_n, scope_n, &tree_n);
  double dist_max = CalculateUnsignedDistanceNearSurface(phi_layers);

  hasInlet_  = hasInlet;
  hasInlet2_ = hasInlet2;
  hasOutlet_ = hasOutlet;
  hasOcc_    = hasOcc;
  nRegions_  = nRegions;

/*
  Debug
  Vec3D*** xf = (Vec3D***) XForward.GetDataPointer();
  XForward.RestoreDataPointerAndInsert();
  XForward.StoreMeshCoordinates(coordinates);
  XForward.WriteToVTRFile("XForward.vtr", "xf");
  Vec3D*** xb = (Vec3D***) XBackward.GetDataPointer();
  XBackward.RestoreDataPointerAndInsert();
  XBackward.StoreMeshCoordinates(coordinates);
  XBackward.WriteToVTRFile("XBackward.vtr", "xb");
  Color.StoreMeshCoordinates(coordinates);
  Color.WriteToVTRFile("Color.vtr", "color");
  Phi.StoreMeshCoordinates(coordinates);
  Phi.WriteToVTRFile("Phi.vtr", "phi");
  fprintf(stdout,"Got here!\n");
  MPI_Barrier(comm);
  exit_mpi();
*/

  return dist_max;
}

//-------------------------------------------------------------------------

double
Intersector::RecomputeFullCourse(vector<Vec3D> &Xprev, int phi_layers)
{
  assert(phi_layers>=1);

  BuildSubdomainScopeAndKDTree(subD_bbmin_1, subD_bbmax_1, scope_1, &tree_1);
  FindIntersections();
  FindSweptNodes(Xprev);
  RefillAfterSurfaceUpdate();

  if(nLayer != phi_layers) {
    BuildNodalAndSubdomainBoundingBoxes(phi_layers, BBmin_n, BBmax_n, subD_bbmin_n, subD_bbmax_n); //n layers
    nLayer = phi_layers;
  }
  BuildSubdomainScopeAndKDTree(subD_bbmin_n, subD_bbmax_n, scope_n, &tree_n);
  double dist_max = CalculateUnsignedDistanceNearSurface(phi_layers);

  return dist_max;
}

//-------------------------------------------------------------------------

void
Intersector::BuildNodalAndSubdomainBoundingBoxes(int nL, SpaceVariable3D &BBmin, SpaceVariable3D &BBmax,
                                                 Vec3D &subD_bbmin, Vec3D &subD_bbmax)
{
  double tol = 0.1;  //i.e. tolerance = 10% of element size (plus 50% of surface thickness)
  double thicker = (1.0+tol)*half_thickness;

  Vec3D*** bbmin  = (Vec3D***) BBmin.GetDataPointer();
  Vec3D*** bbmax  = (Vec3D***) BBmax.GetDataPointer();

  vector<double> &x_glob(global_mesh.x_glob);
  vector<double> &y_glob(global_mesh.y_glob);
  vector<double> &z_glob(global_mesh.z_glob);
  vector<double> &dx_glob(global_mesh.dx_glob);
  vector<double> &dy_glob(global_mesh.dy_glob);
  vector<double> &dz_glob(global_mesh.dz_glob);

  double delta;
  for(int k=kk0_in; k<kkmax_in; k++)
    for(int j=jj0_in; j<jjmax_in; j++)
      for(int i=ii0_in; i<iimax; i++) {

        delta = tol*dx_glob[i] + thicker;
        bbmin[k][j][i][0] = x_glob[std::max(   0, i-nL)] - delta;
        bbmax[k][j][i][0] = x_glob[std::min(NX-1, i+nL)] + delta;

        delta = tol*dy_glob[j] + thicker;
        bbmin[k][j][i][1] = y_glob[std::max(   0, j-nL)] - delta;
        bbmax[k][j][i][1] = y_glob[std::min(NY-1, j+nL)] + delta;

        delta = tol*dz_glob[k] + thicker;
        bbmin[k][j][i][2] = z_glob[std::max(   0, k-nL)] - delta;
        bbmax[k][j][i][2] = z_glob[std::min(NZ-1, k+nL)] + delta;
      }

  //subD_bb includes the ghost boundary
  for(int p=0; p<3; p++) {
    subD_bbmin[p] = bbmin[kk0_in][jj0_in][ii0_in][p];
    subD_bbmax[p] = bbmax[kkmax_in-1][jjmax_in-1][iimax_in-1][p];
  }

  BBmin.RestoreDataPointerToLocalVector();
  BBmax.RestoreDataPointerToLocalVector();

}

//-------------------------------------------------------------------------

void
Intersector::BuildSubdomainScopeAndKDTree(const Vec3D &subD_bbmin, const Vec3D &subD_bbmax, 
                                          vector<MyTriangle> &scope, KDTree<MyTriangle, 3> **tree)//updating the tree itself
{
  scope.clear();

  vector<Vec3D>& Xs(surface.X);
  vector<Int3>&  Es(surface.elems);
  for(auto it = Es.begin(); it != Es.end(); it++) {
    MyTriangle tri(it - Es.begin(), Xs[(*it)[0]], Xs[(*it)[1]], Xs[(*it)[2]]);
    bool inside = true;
    for(int i=0; i<3; i++) {
      if(tri.val(i) > subD_bbmax[i] || tri.val(i) + tri.width(i) < subD_bbmin[i]) {
        inside = false;
        break;
      }
    }
    if(inside) 
      scope.push_back(tri); //creating a copy of "tri" and store it in scope
  }

//  fprintf(stdout,"scope size: %d.\n", (int)scope.size());

  // build the tree
  if(*tree)
    delete *tree;

  if(scope.size()!=0)
    *tree = new KDTree<MyTriangle,3>(scope.size(), scope.data());
  else
    *tree = NULL;
}

//-------------------------------------------------------------------------

void
Intersector::FindNodalCandidates(SpaceVariable3D &BBmin, SpaceVariable3D &BBmax, KDTree<MyTriangle, 3> *tree,
                                 SpaceVariable3D &CandidatesIndex, 
                                 vector<pair<Int3, vector<MyTriangle> > > &candidates)
{

  candidates.clear();

  int nMaxCand = 1000; //will increase if necessary

  Vec3D*** bbmin   = (Vec3D***) BBmin.GetDataPointer();
  Vec3D*** bbmax   = (Vec3D***) BBmax.GetDataPointer();
  double*** candid = CandidatesIndex.GetDataPointer();

  vector<MyTriangle> tmp(nMaxCand);
  
  // Work on all nodes inside the physical domain, including internal ghost layer
  for(int k=kk0_in; k<kkmax_in; k++)
    for(int j=jj0_in; j<jjmax_in; j++)
      for(int i=ii0_in; i<iimax_in; i++) {

        // find candidates
        int nFound = tree ? FindCandidatesInBox(tree, bbmin[k][j][i], bbmax[k][j][i], tmp, nMaxCand) : 0;

        // update candidates and CandidatesIndex
        if(nFound==0) {
          candid[k][j][i] = -1;
        } else {
          candidates.push_back(std::make_pair(Int3(i,j,k), vector<MyTriangle>())); 
          vector<MyTriangle>& tri(candidates.back().second);
          tri.reserve(nFound);
          for(int i=0; i<nFound; i++)
            tri.push_back(tmp[i]); 
          candid[k][j][i] = candidates.size() - 1; //index of this node in candidates
        }

      }

  BBmin.RestoreDataPointerToLocalVector();
  BBmax.RestoreDataPointerToLocalVector();
  CandidatesIndex.RestoreDataPointerToLocalVector(); //can NOT communicate, because "candidates" do not.

}

//-------------------------------------------------------------------------

void
Intersector::FindIntersections() //also finds occluded and first layer nodes
{

  // Find nodal candidates, layer = 1
  FindNodalCandidates(BBmin_1, BBmax_1, tree_1, CandidatesIndex_1, candidates_1);


  Vec3D*** coords  = (Vec3D***) coordinates.GetDataPointer();
  Vec3D*** xf      = (Vec3D***) XForward.GetDataPointer();
  Vec3D*** xb      = (Vec3D***) XBackward.GetDataPointer();
  double*** candid = CandidatesIndex_1.GetDataPointer();
  double*** color  = Color.GetDataPointer();
  
  double*** occid  = TMP.GetDataPointer(); //occluding triangle id
  double*** layer  = TMP2.GetDataPointer(); //"layer" of each node: 0(occluded), 1, or -1 (unknown)

  //Clear previous values
  intersections.clear();

  previously_occluded_but_not_now.clear();

  //Preparation
  Vec3D tol(half_thickness*1.5, half_thickness*1.5, half_thickness*1.5); //a tolerance

  int max_left   = 1000; //will increase if necessary
  int max_bottom = 1000; 
  int max_back   = 1000; 

  vector<MyTriangle> tmp_left(max_left);
  int found_left;
  vector<MyTriangle> tmp_bottom(max_bottom);
  int found_bottom;
  vector<MyTriangle> tmp_back(max_back);
  int found_back;

  IntersectionPoint xp0, xp1;

  // For completeness (to avoid confusion), we set external ghost cells to -1
  for(auto it = ghost_nodes_outer.begin(); it != ghost_nodes_outer.end(); it++) {
      Int3& ijk(it->ijk);
      xf[ijk[2]][ijk[1]][ijk[0]] = -1; //setting all three components to -1
      xb[ijk[2]][ijk[1]][ijk[0]] = -1;
  }
 

  // ----------------------------------------------------------------------------
  // Find occluded nodes and intersections. Build occluded and firstLayer 
  // ----------------------------------------------------------------------------
  firstLayer.clear();
  // We only deal with edges whose vertices are both in the real domain
  for(int k=kk0_in; k<kkmax_in; k++)
    for(int j=jj0_in; j<jjmax_in; j++)
      for(int i=ii0_in; i<iimax_in; i++) {
  
        // start with a meaningless triangle id
        occid[k][j][i] = -1;

        layer[k][j][i] = -1;

        if(!tree_1) //this subdomain is entirely away from surface, just set color, occid, and layer to default values
          continue; 

        if(k<kmax && j<jmax && i<imax && //candid has a valid value @ i,j,k
           candid[k][j][i] < 0)  //no nodal candidates, intersection impossible, occlusion also impossible
          continue;
 

        //--------------------------------------------
        // Find candidates for left/bottom/back edges
        //--------------------------------------------
        found_left = found_bottom = found_back = 0;

        if(i-1>=ii0_in) { //the edge [k][j][i-1] -> [k][j][i] is inside the physical domain
          found_left = FindCandidatesInBox(tree_1, coords[k][j][i-1] - tol, coords[k][j][i] + tol, tmp_left, max_left);
        }

        if(j-1>=jj0_in) { //the edge [k][j-1][i] -> [k][j][i] is inside the physical domain
          found_bottom = FindCandidatesInBox(tree_1, coords[k][j-1][i] - tol, coords[k][j][i] + tol, tmp_bottom, max_bottom);
        }

        if(k-1>=kk0_in) { //the edge [k-1][j][i] -> [k][j][i] is inside the physical domain
          found_back = FindCandidatesInBox(tree_1, coords[k-1][j][i] - tol, coords[k][j][i] + tol, tmp_back, max_back);
        }

        //--------------------------------------------
        // Check if (i,j,k) is occluded
        //--------------------------------------------
        int tid(-1);
        if ((found_left>0   && IsPointOccludedByTriangles(coords[k][j][i], tmp_left.data(),   found_left,   half_thickness, tid)) ||
            (found_bottom>0 && IsPointOccludedByTriangles(coords[k][j][i], tmp_bottom.data(), found_bottom, half_thickness, tid)) ||
            (found_back>0   && IsPointOccludedByTriangles(coords[k][j][i], tmp_back.data(),   found_back,   half_thickness, tid))) {
          color[k][j][i] = 0;
          occid[k][j][i] = tid;
          layer[k][j][i] = 0;          
        }
        else {
          if(color[k][j][i]==0) {// previously occluded!
            color[k][j][i] = -1; // to be corrected in "Refill"
            previously_occluded_but_not_now.insert(Int3(i,j,k));
          }
        }

        //--------------------------------------------
        // Find intersections (left, bottom, and back edges)
        //--------------------------------------------
        // left
        if(found_left) {
          int count = FindEdgeIntersectionsWithTriangles(coords[k][j][i-1], i-1, j, k, 0, 
                          coords[k][j][i][0] - coords[k][j][i-1][0], tmp_left.data(), found_left, xp0, xp1);
          if(count==0) {
            xf[k][j][i][0] = xb[k][j][i][0] = -1;
          } else if(count==1) {
            intersections.push_back(xp0);
            xf[k][j][i][0] = xb[k][j][i][0] = intersections.size() - 1;
            if(layer[k][j][i-1]==-1)  //it might be 0, meaning occluded. in that case we don't override
              layer[k][j][i-1] = 1;
            if(layer[k][j][i]==-1)
              layer[k][j][i] = 1;
          } else {//more than one intersections
            intersections.push_back(xp0);
            xf[k][j][i][0] = intersections.size() - 1;
            intersections.push_back(xp1);
            xb[k][j][i][0] = intersections.size() - 1;
            if(layer[k][j][i-1]==-1)
              layer[k][j][i-1] =1;
            if(layer[k][j][i]==-1)
              layer[k][j][i] = 1;
          }
        } else {
          xf[k][j][i][0] = xb[k][j][i][0] = -1;
        }

        // bottom 
        if(found_bottom) {
          int count = FindEdgeIntersectionsWithTriangles(coords[k][j-1][i], i, j-1, k, 1, 
                          coords[k][j][i][1] - coords[k][j-1][i][1], tmp_bottom.data(), found_bottom, xp0, xp1);
          if(count==0) {
            xf[k][j][i][1] = xb[k][j][i][1] = -1;
          } else if(count==1) {
            intersections.push_back(xp0);
            xf[k][j][i][1] = xb[k][j][i][1] = intersections.size() - 1;
            if(layer[k][j-1][i]==-1)
              layer[k][j-1][i] = 1;
            if(layer[k][j][i]==-1)
              layer[k][j][i] = 1;
          } else {//more than one intersections
            intersections.push_back(xp0);
            xf[k][j][i][1] = intersections.size() - 1;
            intersections.push_back(xp1);
            xb[k][j][i][1] = intersections.size() - 1;
            if(layer[k][j-1][i]==-1)
              layer[k][j-1][i] = 1;
            if(layer[k][j][i]==-1)
              layer[k][j][i] = 1;
          } 
        } else {
          xf[k][j][i][1] = xb[k][j][i][1] = -1;
        }

        // back 
        if(found_back) {
          int count = FindEdgeIntersectionsWithTriangles(coords[k-1][j][i], i, j, k-1, 2, 
                          coords[k][j][i][2] - coords[k-1][j][i][2], tmp_back.data(), found_back, xp0, xp1);
          if(count==0) {
            xf[k][j][i][2] = xb[k][j][i][2] = -1;
          } else if(count==1) {
            intersections.push_back(xp0);
            xf[k][j][i][2] = xb[k][j][i][2] = intersections.size() - 1;
            if(layer[k-1][j][i]==-1)
              layer[k-1][j][i] = 1;
            if(layer[k][j][i]==-1)
              layer[k][j][i] = 1;
          } else {//more than one intersections
            intersections.push_back(xp0);
            xf[k][j][i][2] = intersections.size() - 1;
            intersections.push_back(xp1);
            xb[k][j][i][2] = intersections.size() - 1;
            if(layer[k-1][j][i]==-1)
              layer[k-1][j][i] = 1;
            if(layer[k][j][i]==-1)
              layer[k][j][i] = 1;
          }
        } else {
          xf[k][j][i][2] = xb[k][j][i][2] = -1;
        }

      }

  // Exchange Color and TMP so internal ghost nodes are accounted for
  CandidatesIndex_1.RestoreDataPointerToLocalVector();
  Color.RestoreDataPointerAndInsert();
  TMP.RestoreDataPointerAndInsert();
  TMP2.RestoreDataPointerAndInsert();

  occid  = TMP.GetDataPointer();

  // ----------------------------------------------------------------------------
  // Make sure all edges connected to occluded nodes have intersections
  // ----------------------------------------------------------------------------
  bool ijk_occluded;
  for(int k=kk0_in; k<kkmax_in; k++)
    for(int j=jj0_in; j<jjmax_in; j++)
      for(int i=ii0_in; i<iimax_in; i++) {
 
        ijk_occluded = (occid[k][j][i]>=0);
        if(i-1>=ii0_in) { //left edge within physical domain
          
          if(occid[k][j][i-1]<0 && !ijk_occluded) { //neither (i-1,j,k) nor (i,j,k) occluded
            //nothing to be done   
          } 
          else { //at least one of the two vertices is occluded

            // insert "imposed" intersection points
            if(occid[k][j][i-1]>=0) { //(i-1,j,k) is occluded
              Vec3D xi;
              int trif = occid[k][j][i-1];
              Int3& nf(surface.elems[trif]);
              // re-run the checker to get xi
              bool is_occluded = GeoTools::IsPointInsideTriangle(coords[k][j][i-1], surface.X[nf[0]],
                                            surface.X[nf[1]], surface.X[nf[2]], half_thickness,
                                            &surface.elemArea[trif], &surface.elemNorm[trif], xi);
              assert(is_occluded);
              intersections.push_back(IntersectionPoint(i-1,j,k, 0, 0.0, trif, xi));
            }
            if(ijk_occluded) { //(i-1,j,k) is occluded
              Vec3D xi;
              int trif = occid[k][j][i];
              Int3& nf(surface.elems[trif]);
              // re-run the checker to get xi
              bool is_occluded = GeoTools::IsPointInsideTriangle(coords[k][j][i], surface.X[nf[0]],
                                            surface.X[nf[1]], surface.X[nf[2]], half_thickness,
                                            &surface.elemArea[trif], &surface.elemNorm[trif], xi);
              assert(is_occluded);
              intersections.push_back(IntersectionPoint(i-1,j,k, 0, coords[k][j][i][0] - coords[k][j][i-1][0], trif, xi));
            }

            // update xb and xf
            if(xf[k][j][i][0]<0) {
              if(occid[k][j][i-1]>=0 && ijk_occluded) { //both vertices are occluded
                xf[k][j][i][0] = intersections.size() - 2;
                xb[k][j][i][0] = intersections.size() - 1;
              } else {// only one of the two vertices is occluded 
                xf[k][j][i][0] = xb[k][j][i][0] = intersections.size() - 1;
              }
            }
            else { //intersection already found. 
              //ensure that the occluded node(s) is the intersection point closest to the occluded node.
              if(occid[k][j][i-1]>=0 && ijk_occluded) {
                xf[k][j][i][0] = intersections.size() - 2;
                xb[k][j][i][0] = intersections.size() - 1;
              } else if(occid[k][j][i-1]>=0) {
                xf[k][j][i][0] = intersections.size() - 1;
                IntersectionPoint &p(intersections[xb[k][j][i][0]]);
                if(p.dist<=half_thickness) {//this is essentially the first vertex
                  xb[k][j][i][0] = intersections.size() - 1;
                }
              } else { //ijk_occluded
                xb[k][j][i][0] = intersections.size() - 1;
                IntersectionPoint &p(intersections[xf[k][j][i][0]]);
                if(p.dist >= coords[k][j][i][0] - coords[k][j][i-1][0] - half_thickness) {//this is essentially the second vertex
                  xf[k][j][i][0] = intersections.size() - 1;
                }
              }
            }
          }
        }


        if(j-1>=jj0_in) { //bottom edge within physical domain
          
          if(occid[k][j-1][i]<0 && !ijk_occluded) { //neither (i,j-1,k) nor (i,j,k) occluded
            //nothing to be done   
          } 
          else { //at least one of the two vertices is occluded

            // insert "imposed" intersection points
            if(occid[k][j-1][i]>=0) { //(i,j-1,k) is occluded
              Vec3D xi;
              int trif = occid[k][j-1][i];
              Int3& nf(surface.elems[trif]);
              // re-run the checker to get xi
              bool is_occluded = GeoTools::IsPointInsideTriangle(coords[k][j-1][i], surface.X[nf[0]],
                                            surface.X[nf[1]], surface.X[nf[2]], half_thickness,
                                            &surface.elemArea[trif], &surface.elemNorm[trif], xi);
              assert(is_occluded);
              intersections.push_back(IntersectionPoint(i,j-1,k, 1, 0.0, trif, xi));
            }
            if(ijk_occluded) { //(i-1,j,k) is occluded
              Vec3D xi;
              int trif = occid[k][j][i];
              Int3& nf(surface.elems[trif]);
              // re-run the checker to get xi
              bool is_occluded = GeoTools::IsPointInsideTriangle(coords[k][j][i], surface.X[nf[0]],
                                            surface.X[nf[1]], surface.X[nf[2]], half_thickness,
                                            &surface.elemArea[trif], &surface.elemNorm[trif], xi);
              assert(is_occluded);
              intersections.push_back(IntersectionPoint(i,j-1,k, 1, coords[k][j][i][1] - coords[k][j-1][i][1], trif, xi));
            }

            // update xb and xf
            if(xf[k][j][i][1]<0) {
              if(occid[k][j-1][i]>=0 && ijk_occluded) { //both vertices are occluded
                xf[k][j][i][1] = intersections.size() - 2;
                xb[k][j][i][1] = intersections.size() - 1;
              } else {// only one of the two vertices is occluded 
                xf[k][j][i][1] = xb[k][j][i][1] = intersections.size() - 1;
              }
            }
            else { //intersection already found. 
              //ensure that the occluded node(s) is the intersection point closest to the occluded node.
              if(occid[k][j-1][i]>=0 && ijk_occluded) {
                xf[k][j][i][1] = intersections.size() - 2;
                xb[k][j][i][1] = intersections.size() - 1;
              } else if(occid[k][j-1][i]>=0) {
                xf[k][j][i][1] = intersections.size() - 1;
                IntersectionPoint &p(intersections[xb[k][j][i][1]]);
                if(p.dist<=half_thickness) {//this is essentially the first vertex
                  xb[k][j][i][1] = intersections.size() - 1;
                }
              } else { //ijk_occluded
                xb[k][j][i][1] = intersections.size() - 1;
                IntersectionPoint &p(intersections[xf[k][j][i][1]]);
                if(p.dist >= coords[k][j][i][1] - coords[k][j-1][i][1] - half_thickness) {//this is essentially the second vertex
                  xf[k][j][i][1] = intersections.size() - 1;
                }
              }
            }
          }
        }


        if(k-1>=kk0_in) { //back edge within physical domain
          
          if(occid[k-1][j][i]<0 && !ijk_occluded) { //neither (i,j,k-1) nor (i,j,k) occluded
            //nothing to be done   
          } 
          else { //at least one of the two vertices is occluded

            // insert "imposed" intersection points
            if(occid[k-1][j][i]>=0) { //(i,j-1,k) is occluded
              Vec3D xi;
              int trif = occid[k-1][j][i];
              Int3& nf(surface.elems[trif]);
              // re-run the checker to get xi
              bool is_occluded = GeoTools::IsPointInsideTriangle(coords[k-1][j][i], surface.X[nf[0]],
                                            surface.X[nf[1]], surface.X[nf[2]], half_thickness,
                                            &surface.elemArea[trif], &surface.elemNorm[trif], xi);
              assert(is_occluded);
              intersections.push_back(IntersectionPoint(i,j,k-1, 2, 0.0, trif, xi));
            }
            if(ijk_occluded) { //(i-1,j,k) is occluded
              Vec3D xi;
              int trif = occid[k][j][i];
              Int3& nf(surface.elems[trif]);
              // re-run the checker to get xi
              bool is_occluded = GeoTools::IsPointInsideTriangle(coords[k][j][i], surface.X[nf[0]],
                                            surface.X[nf[1]], surface.X[nf[2]], half_thickness,
                                            &surface.elemArea[trif], &surface.elemNorm[trif], xi);
              assert(is_occluded);
              intersections.push_back(IntersectionPoint(i,j,k-1, 2, coords[k][j][i][2] - coords[k-1][j][i][2], trif, xi));
            }

            // update xb and xf
            if(xf[k][j][i][2]<0) {
              if(occid[k-1][j][i]>=0 && ijk_occluded) { //both vertices are occluded
                xf[k][j][i][2] = intersections.size() - 2;
                xb[k][j][i][2] = intersections.size() - 1;
              } else {// only one of the two vertices is occluded 
                xf[k][j][i][2] = xb[k][j][i][2] = intersections.size() - 1;
              }
            }
            else { //intersection already found. 
              //ensure that the occluded node(s) is the intersection point closest to the occluded node.
              if(occid[k-1][j][i]>=0 && ijk_occluded) {
                xf[k][j][i][2] = intersections.size() - 2;
                xb[k][j][i][2] = intersections.size() - 1;
              } else if(occid[k-1][j-1][i]>=0) {
                xf[k][j][i][2] = intersections.size() - 1;
                IntersectionPoint &p(intersections[xb[k][j][i][2]]);
                if(p.dist<=half_thickness) {//this is essentially the first vertex
                  xb[k][j][i][2] = intersections.size() - 1;
                }
              } else { //ijk_occluded
                xb[k][j][i][2] = intersections.size() - 1;
                IntersectionPoint &p(intersections[xf[k][j][i][2]]);
                if(p.dist >= coords[k][j][i][2] - coords[k-1][j][i][2] - half_thickness) {//this is essentially the second vertex
                  xf[k][j][i][2] = intersections.size() - 1;
                }
              }
            }
          }
        }

      }



  TMP.RestoreDataPointerToLocalVector();

  XForward.RestoreDataPointerToLocalVector(); //Cannot exchange data, because "intersections" does not communicate
  XBackward.RestoreDataPointerToLocalVector(); //Cannot exchange data, because "intersections" does not communicate

  coordinates.RestoreDataPointerToLocalVector();
/*
  XForward.RestoreDataPointerAndInsert();
  XForward.StoreMeshCoordinates(coordinates);
  XForward.WriteToVTRFile("XForward.vtr", "xf");
  XBackward.RestoreDataPointerAndInsert();
  XBackward.StoreMeshCoordinates(coordinates);
  XBackward.WriteToVTRFile("XBackward.vtr", "xb");
*/
  // ----------------------------------------------------------------------------
  // Build the sets of occluded and firstLayer nodes. Include internal ghost nodes
  // ----------------------------------------------------------------------------
  occluded.clear();
  firstLayer.clear();

  layer  = TMP2.GetDataPointer(); //"layer" of each node: 0(occluded), 1, or -1 (unknown)

  for(int k=kk0_in; k<kkmax_in; k++)
    for(int j=jj0_in; j<jjmax_in; j++)
      for(int i=ii0_in; i<iimax_in; i++) {
        if(layer[k][j][i]==0) {
          occluded.insert(Int3(i,j,k));         
          firstLayer.insert(Int3(i,j,k));
        } else if(layer[k][j][i]==1) 
          firstLayer.insert(Int3(i,j,k));
      }

  TMP2.RestoreDataPointerToLocalVector();

/*
  TMP2.StoreMeshCoordinates(coordinates);
  TMP2.WriteToVTRFile("TMP2.vtr", "layer");
  TMP.StoreMeshCoordinates(coordinates);
  TMP.WriteToVTRFile("TMP.vtr", "occid");
  Color.StoreMeshCoordinates(coordinates);
  Color.WriteToVTRFile("Color.vtr", "color");
*/
}

//-------------------------------------------------------------------------

int
Intersector::FindCandidatesInBox(KDTree<MyTriangle, 3>* mytree, Vec3D bbmin, Vec3D bbmax, 
                                 vector<MyTriangle> &tmp, int& maxCand)
{
  int found = mytree->findCandidatesInBox(bbmin, bbmax, tmp.data(), maxCand);
  if(found>maxCand) {
    maxCand = found; 
    tmp.resize(maxCand);
    found = mytree->findCandidatesInBox(bbmin, bbmax, tmp.data(), maxCand);
  }
  return found;
}

//-------------------------------------------------------------------------

bool
Intersector::IsPointOccludedByTriangles(Vec3D &x0, MyTriangle* tri, int nTri, double my_half_thickness,
                                        int &tid, double* xi) 
{
  vector<Vec3D>&  Xs(surface.X);
  vector<Int3>&   Es(surface.elems);
  vector<Vec3D>&  Ns(surface.elemNorm);
  vector<double>& As(surface.elemArea); 

  for(int i=0; i<nTri; i++) {
    int id = tri[i].trId();
    Int3& nodes(Es[id]);
    if(GeoTools::IsPointInsideTriangle(x0, Xs[nodes[0]], Xs[nodes[1]], Xs[nodes[2]], my_half_thickness,
                                       &As[id], &Ns[id], xi)) {
      tid = id;
      return true;
    }
  }

  return false;
}

//-------------------------------------------------------------------------

bool
Intersector::FloodFillColors()
{
  // ----------------------------------------------------------------
  // Call floodfiller to do the work.
  // ----------------------------------------------------------------
  int nColors = floodfiller.FillBasedOnEdgeObstructions(XForward, -1/*xf==-1 means no intersection*/, occluded, Color);

  // ----------------------------------------------------------------
  // Now we need to convert the color map to what we want: 0~occluded, 1, 2,...: regions connected to Dirichlet
  // boundaries (inlet/inlet2/farfield). 1=inlet, 2=inlet2/farfield, 3=outlet. -1,-2,...: inside enclosures
  // ----------------------------------------------------------------
  double*** color = Color.GetDataPointer();

  // First, we need to find the current colors for inlet, inlet2, outlet
  std::set<int> inlet_color; 
  std::set<int> outlet_color; 
  std::set<int> inlet2_color;  //In most cases, all 3 sets should only contain 0 or 1 member. Otherwise, it means 
                               //one type of boundary (e.g., inlet) is separated by the embedded surface into multiple
                               //disconnected regions. In this case, we have to merge the colors into one.
  for(auto it = ghost_nodes_outer.begin(); it != ghost_nodes_outer.end(); it++) {
    if(it->type_projection != GhostPoint::FACE)
      continue;
    if(it->bcType == MeshData::INLET || it->bcType == MeshData::OVERSET) {
      Int3& ijk(it->image_ijk);
      inlet_color.insert(color[ijk[2]][ijk[1]][ijk[0]]);
    } else if(it->bcType == MeshData::INLET2) {
      Int3& ijk(it->image_ijk);
      inlet2_color.insert(color[ijk[2]][ijk[1]][ijk[0]]);
    } else if(it->bcType == MeshData::OUTLET) {
      Int3& ijk(it->image_ijk);
      outlet_color.insert(color[ijk[2]][ijk[1]][ijk[0]]);
    }
  }

  // Deals with inlet 
  int max_inlet_color(-1);
  for(auto it = inlet_color.begin(); it != inlet_color.end(); it++)
    max_inlet_color = std::max(max_inlet_color, *it);
  MPI_Allreduce(MPI_IN_PLACE, &max_inlet_color, 1, MPI_INT, MPI_MAX, comm);

  vector<int> in_colors;
  if(max_inlet_color!=-1) {
    in_colors.assign(max_inlet_color+1, -1);
    for(auto it = inlet_color.begin(); it != inlet_color.end(); it++)
      in_colors[*it] = 1;
    MPI_Allreduce(MPI_IN_PLACE, in_colors.data(), in_colors.size(), MPI_INT, MPI_MAX, comm);
  }

  if(in_colors.size()==0) 
    print_warning("Warning: Intersector did not find an inlet/farfield boundary of the M2C domain.\n");

  if(in_colors.size()>0 && in_colors[0] == 1 && verbose>1)
    print_warning("Warning: Found occluded node(s) near an inlet or farfield boundary.");


  // Deals with inlet2
  int max_inlet2_color(-1);
  for(auto it = inlet2_color.begin(); it != inlet2_color.end(); it++) 
    max_inlet2_color = std::max(max_inlet2_color, *it);
  MPI_Allreduce(MPI_IN_PLACE, &max_inlet2_color, 1, MPI_INT, MPI_MAX, comm);

  vector<int> in2_colors;
  if(max_inlet2_color!=-1) {
    in2_colors.assign(max_inlet2_color+1, -1);
    for(auto it = inlet2_color.begin(); it != inlet2_color.end(); it++)
      in2_colors[*it] = 1;
    MPI_Allreduce(MPI_IN_PLACE, in2_colors.data(), in2_colors.size(), MPI_INT, MPI_MAX, comm);
  }
 
  if(in2_colors.size()>0 && in2_colors[0] == 1 && verbose>1)
    print_warning("Warning: Found occluded node(s) near an inlet2 or farfield boundary.");


  // Deals with outlet 
  int max_outlet_color(-1);
  for(auto it = outlet_color.begin(); it != outlet_color.end(); it++) 
    max_outlet_color = std::max(max_outlet_color, *it);
  MPI_Allreduce(MPI_IN_PLACE, &max_outlet_color, 1, MPI_INT, MPI_MAX, comm);

  vector<int> out_colors;
  if(max_outlet_color!=-1) {
    out_colors.assign(max_outlet_color+1, -1);
    for(auto it = outlet_color.begin(); it != outlet_color.end(); it++)
      out_colors[*it] = 1;
    MPI_Allreduce(MPI_IN_PLACE, out_colors.data(), out_colors.size(), MPI_INT, MPI_MAX, comm);
  }
 
  if(out_colors.size()>0 && out_colors[0] == 1 && verbose>1)
    print_warning("Warning: Found occluded node(s) near an outlet boundary.");





  // Convert colors 
  std::map<int,int> old2new;
  for(int i=0; i<(int)in_colors.size(); i++)
    if(in_colors[i] == 1)
      old2new[i] = 1; //inlet_color;
  for(int i=0; i<(int)in2_colors.size(); i++) {
    if(in2_colors[i] == 1) {
      if(old2new.find(i) != old2new.end() && verbose>=1) {
        print_warning("Warning: Embedded surface (%d elems) cannot separate Inlet/Farfield and Inlet2/Farfield2.\n"
                      "         Treating both as Inlet/Farfield.\n", surface.elems.size());
        old2new[i] = 1; //inlet_color
      } else
        old2new[i] = 2; //inlet2_color;
    }
  }
  for(int i=0; i<(int)out_colors.size(); i++) {
    if(out_colors[i] == 1) {
      if(old2new.find(i) != old2new.end() && verbose>=1) {
        if(in_colors[i] == 1) {
          print_warning("Warning: Embedded surface (%d elems) cannot separate Inlet/Farfield and Outlet.\n"
                        "         Treating both as Inlet/Farfield.\n", surface.elems.size());
          old2new[i] = 1; //inlet_color
        } else {
          print_warning("Warning: Embedded surface (%d elems) cannot separate Inlet2/Farfield2 and Outlet.\n"
                        "         Treating both as Inlet2/Farfield2.\n", surface.elems.size());
          old2new[i] = 2; //inlet2_color
        }
      } else
        old2new[i] = 3; //outlet_color;
    }
  }



  // give negative colors to enclusures (i.e. regions not connected to farfield)
  int tmp_counter = 0;
  for(int i=1; i<nColors+1; i++)
    if(old2new.find(i) == old2new.end())
      old2new[i] = --tmp_counter;

/*
  for(auto it = old2new.begin(); it != old2new.end(); it++)
    fprintf(stdout,"old2new: %d --> %d.\n", it->first, it->second);
*/

  int total_occluded = 0;

  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        if(color[k][j][i] != 0)
          color[k][j][i] = old2new[color[k][j][i]];
        else
          total_occluded++;
      }

  // statistics
  MPI_Allreduce(MPI_IN_PLACE, &total_occluded, 1, MPI_INT, MPI_SUM, comm);

  bool hasOcc = total_occluded>0; //i.e. one zero color (for occluded)

  hasInlet = hasInlet2 = hasOutlet = false;
  nRegions = 0;
  for(auto it = old2new.begin(); it != old2new.end(); it++)
    if(it->second==1)
      hasInlet = true;
    else if(it->second==2)
      hasInlet2 = true;
    else if(it->second==3)
      hasOutlet = true;
    else if(it->second<0)
      nRegions++;


  // fill ColorReachesBoundary
  ColorReachesBoundary.assign(nRegions, 0);
  for(auto it = ghost_nodes_outer.begin(); it != ghost_nodes_outer.end(); it++) {
    if(it->type_projection != GhostPoint::FACE)
      continue;
    Int3& ijk(it->image_ijk);
    int mycolor = color[ijk[2]][ijk[1]][ijk[0]];
    if(mycolor<0) 
      ColorReachesBoundary[-mycolor] = 1;
  }
  MPI_Allreduce(MPI_IN_PLACE, ColorReachesBoundary.data(), ColorReachesBoundary.size(), MPI_INT, MPI_MAX, comm);
 


  // Finalization

  Color.RestoreDataPointerAndInsert();

  return hasOcc;
}

//-------------------------------------------------------------------------

void
Intersector::RefillAfterSurfaceUpdate()
{
  // "FindSweptNodes" should be called before calling this function. So, in most cases, "nodal_cands_calculated" should be true
  
  // Assuming occluded nodes already have the correct color. We do not deal with them here.
  
  vector<double> &x_glob(global_mesh.x_glob);
  vector<double> &y_glob(global_mesh.y_glob);
  vector<double> &z_glob(global_mesh.z_glob);

  double*** candid = CandidatesIndex_1.GetDataPointer();
  
  //add swept nodes to nodes2fill (including internal ghosts)
  std::set<Int3> nodes2fill = swept;
  std::set<Int3> nodes2fill2;

  // go over swept nodes, correct their colors
  int i0,j0,k0;
  vector<Vec3D>&  Xs(surface.X);
  vector<Int3>&   Es(surface.elems);

  int BAD_SIGN = -999999; //used temporarily.

  int max_it = 100;
  double*** color = NULL;
  int total_remaining_nodes;
  for(int iter = 0; iter < max_it; iter++) {

    total_remaining_nodes = nodes2fill.size();
    MPI_Allreduce(MPI_IN_PLACE, &total_remaining_nodes, 1, MPI_INT, MPI_SUM, comm);
    if(total_remaining_nodes == 0) //Yeah
      break;
  
/*
    if(total_remaining_nodes>0) {
      fprintf(stdout,"Hey total_remaining_nodes = %d.\n", total_remaining_nodes);
      Color.StoreMeshCoordinates(coordinates);
      Color.WriteToVTRFile("Color.vtr", "color");
    }
*/

    if(!color) //first iteration
      color = Color.GetDataPointer();

    nodes2fill2 = nodes2fill;

    for(auto it = nodes2fill.begin(); it != nodes2fill.end(); it++) {
      i0 = (*it)[0];
      j0 = (*it)[1];
      k0 = (*it)[2]; 

      if(!coordinates.IsHere(i0,j0,k0,false))
        continue; //this is an internal ghost. we let its owner fix it.

      if(color[k0][j0][i0] == 0) { //occluded
        nodes2fill2.erase(nodes2fill2.find(*it)); 
        continue;
      }

      color[k0][j0][i0] = BAD_SIGN;

      // if this node is swept, candidates must exist. Otherwise, the surface moved too much!
      assert(candid[k0][j0][i0]>=0);
      vector<MyTriangle> &cands(candidates_1[candid[k0][j0][i0]].second);
      assert(cands.size()>0);

      //go over first layer neighbors, find a nonblocked reliable neighbor
      for(int k=k0-1; k<=k0+1; k++)
        for(int j=j0-1; j<=j0+1; j++)
          for(int i=i0-1; i<=i0+1; i++) {

            if(coordinates.OutsidePhysicalDomain(i,j,k) || (i==i0 && j==j0 && k==k0))
              continue; //this neighbor is out of physical domain, or just the same node

            if(nodes2fill2.find(Int3(i,j,k)) != nodes2fill2.end()) //uses "nodes2fill2", the latest updated list
              continue; //this neighbor is in trouble as well...

            if(color[k][j][i] == 0) //this neighbor is occluded (naturally or "forced"). Either way, it cannot be used.
              continue;

            bool blocked = false;
            for(auto it2 = cands.begin(); it2 != cands.end(); it2++) {
              int id = it2->trId();
              Int3 &nodes(Es[id]);
              blocked = GeoTools::LineSegmentIntersectsTriangle(Vec3D(x_glob[i],y_glob[j],z_glob[k]),
                                                                Vec3D(x_glob[i0],y_glob[j0],z_glob[k0]), 
                                                                Xs[nodes[0]], Xs[nodes[1]], Xs[nodes[2]]);

              if(blocked)
                break;
            }

            if(blocked)
              continue; //this neighbor is blocked from me by the surface

            // If the above checks are all passed, this neighbor is in the same region as me.
            color[k0][j0][i0] = color[k][j][i];
            nodes2fill2.erase(nodes2fill2.find(*it)); //erase from nodes2fill2, NOT from nodes2fill (would mess up pointer!)
            goto DONE_WITH_THIS_NODE; 
         }

      DONE_WITH_THIS_NODE:
      continue; //need this to make "GOTO" work
    }
    nodes2fill = nodes2fill2;

    Color.RestoreDataPointerAndInsert();

    // Get new data
    color = Color.GetDataPointer();

    // remove filled internal ghosts (filled by their owners) from "nodes2fill"
    for(auto it = nodes2fill.begin(); it != nodes2fill.end(); it++) {
      i0 = (*it)[0];
      j0 = (*it)[1];
      k0 = (*it)[2]; 
      if(!coordinates.IsHere(i0,j0,k0,false)) { //this is an internal ghost
        if(color[k0][j0][i0] != BAD_SIGN) //must have been fixed by its owner
          nodes2fill2.erase(nodes2fill2.find(*it)); 
      }
    }
    if(nodes2fill.size() != nodes2fill2.size()) //some node(s) got erased
      nodes2fill = nodes2fill2;

  }

  // fill remaining nodes as FORCED occluded
  imposed_occluded.clear();
  if(total_remaining_nodes>0) {
    print_warning("Warning: Found %d unresolved nodes after performing %d iterations of refill. "
                   "Setting them to be occluded.\n", total_remaining_nodes, max_it);
    imposed_occluded = nodes2fill;
    for(auto it = imposed_occluded.begin(); it != imposed_occluded.end(); it++) 
      color[(*it)[2]][(*it)[1]][(*it)[0]] = 0; //set it to occluded. BUT NO NEW INTERSECTIONS!
      //color must be valid (i.e. not NULL) if total_remaining_nodes>0
  }


  // Update "ColorReachesBoundary"
  if(color) {
    for(auto it = ghost_nodes_outer.begin(); it != ghost_nodes_outer.end(); it++) {
      if(it->type_projection != GhostPoint::FACE)
        continue;
      Int3& ijk(it->image_ijk);
      int mycolor = color[ijk[2]][ijk[1]][ijk[0]];
      if(mycolor<0) 
        ColorReachesBoundary[-mycolor] = 1;
    }
    MPI_Allreduce(MPI_IN_PLACE, ColorReachesBoundary.data(), ColorReachesBoundary.size(), MPI_INT, MPI_MAX, comm);
  }


  if(color)
    Color.RestoreDataPointerToLocalVector();

  CandidatesIndex_1.RestoreDataPointerToLocalVector();

}

//-------------------------------------------------------------------------
 
void
Intersector::FindSweptNodes(std::vector<Vec3D> &X0)
{
  
  vector<double> &x_glob(global_mesh.x_glob);
  vector<double> &y_glob(global_mesh.y_glob);
  vector<double> &z_glob(global_mesh.z_glob);

  swept.clear();

  double*** candid = CandidatesIndex_1.GetDataPointer();

  vector<Vec3D>&  Xs(surface.X);
  vector<Int3>&   Es(surface.elems);
  vector<Vec3D>&  Ns(surface.elemNorm);
  vector<double>& As(surface.elemArea); 

  double collision_time;

  int i,j,k;
  for(auto it = firstLayer.begin(); it != firstLayer.end(); it++) {
    if(occluded.find(*it) != occluded.end())
      continue; //we don't store nodes that are currently occluded

    i = (*it)[0];
    j = (*it)[1];
    k = (*it)[2];

    assert(candid[k][j][i]>=0);
    vector<MyTriangle> &cands(candidates_1[candid[k][j][i]].second);
    assert(cands.size()>0);

    for(auto it2 = cands.begin(); it2 != cands.end(); it2++) {
      int id = it2->trId();
      Int3 &nodes(Es[id]);
      Vec3D coords(x_glob[(*it)[0]], y_glob[(*it)[1]], z_glob[(*it)[2]]); //inside physical domain (safe)
      if(GeoTools::IsPointSweptByTriangle(coords, X0[nodes[0]], X0[nodes[1]], X0[nodes[2]],
                                          Xs[nodes[0]], Xs[nodes[1]], Xs[nodes[2]], &collision_time, half_thickness, 
                                          NULL, NULL, &(As[id]), &(Ns[id]))) {
        swept.insert(*it);
        //fprintf(stdout,"Found swept node: [%d][%d][%d] [%e %e %e].\n", (*it)[2], (*it)[1], (*it)[0],
        //        coords[0], coords[1], coords[2]);
      }
    }
  }
      
  // Verification
  if(verbose>=1) {
    for(auto&& ijk : previously_occluded_but_not_now) {
      if(swept.find(ijk) == swept.end())
        fprintf(stdout,"\033[0;35mWarning: Previous occluded node is not swept (%d %d %d). May be an error.\033[0m\n",
              ijk[0], ijk[1], ijk[2]);
    } 
  }
  
  CandidatesIndex_1.RestoreDataPointerToLocalVector();
}

//-------------------------------------------------------------------------

double
Intersector::CalculateUnsignedDistanceNearSurface(int nL)
{

  assert(nLayer = nL);

  double max_dist = -DBL_MAX;

  FindNodalCandidates(BBmin_n, BBmax_n, tree_n, CandidatesIndex_n, candidates_n);

  vector<double> &x_glob(global_mesh.x_glob);
  vector<double> &y_glob(global_mesh.y_glob);
  vector<double> &z_glob(global_mesh.z_glob);

  Phi_nLayer = nLayer;

  double*** candid = CandidatesIndex_n.GetDataPointer();
  double*** phi    = Phi.GetDataPointer();
  double*** cpi    = ClosestPointIndex.GetDataPointer();
  
  // set const. value to phi, by default
  double default_distance = domain_diagonal;
  for(int k=kk0_in; k<kkmax_in; k++)
    for(int j=jj0_in; j<jjmax_in; j++)
      for(int i=ii0_in; i<iimax_in; i++) {
        phi[k][j][i] = default_distance;
        cpi[k][j][i] = -1;
      }
  closest_points.clear();

  assert(nLayer>=1);

  // loop through layers. 1st layer means nodes connected to intersecting edges (including occluded nodes).
  vector<Vec3D>&  Xs(surface.X);
  vector<Int3>&   Es(surface.elems);
  vector<Vec3D>&  Ns(surface.elemNorm);
  vector<double>& As(surface.elemArea); 
  std::set<Int3> this_layer = firstLayer;

  for(int layer=1; layer<=nLayer; layer++) {

    std::set<Int3> next_layer;

    int i,j,k;
    double bar = 0.99*default_distance;
    double xi[3];
    for(auto it = this_layer.begin(); it != this_layer.end(); it++) {

      i = (*it)[0];
      j = (*it)[1];
      k = (*it)[2];

      assert(candid[k][j][i]>=0);
      vector<MyTriangle> &cands(candidates_n[candid[k][j][i]].second);
      assert(cands.size()>0);

      double dist = DBL_MAX, new_dist;
      ClosestPoint cp(-1,DBL_MAX,xi); //initialize to garbage

      int id;
      for(int tri=0; tri<(int)cands.size(); tri++) {
        id = cands[tri].trId();
        Int3 &nodes(Es[id]);
        Vec3D coords(x_glob[i], y_glob[j], z_glob[k]); //inside physical domain (safe)
        
        new_dist = GeoTools::ProjectPointToTriangle(coords, Xs[nodes[0]], Xs[nodes[1]], Xs[nodes[2]], xi,
                                                    &(As[id]), &(Ns[id]), false);
        if(new_dist<dist) {
          dist = new_dist;
          cp.tid = id;
          cp.dist = new_dist;
          for(int s=0; s<3; s++)
            cp.xi[s] = xi[s];
        }
      }
      phi[k][j][i] = dist;
      closest_points.push_back(std::make_pair(*it, cp));
      cpi[k][j][i] = closest_points.size() - 1;

      if(dist>max_dist)
        max_dist = dist;
/*
      //verification
      double dist_true = fabs(19.4555 - sqrt(x_glob[i]*x_glob[i]+y_glob[j]*y_glob[j]));
      fprintf(stdout,"[%d][%d][%d]: (%e, %e, %e) true: %e, numr: %e.\n", i,j,k, x_glob[i], y_glob[j], z_glob[k],
                     dist_true, dist);
*/

      //insert neighbors to the next layer
      if(layer<nLayer) {
        assert(dist < bar);
        if(i-1>=ii0_in  && phi[k][j][i-1] >= bar) next_layer.insert(Int3(i-1,j,k));
        if(i+1<iimax_in && phi[k][j][i+1] >= bar) next_layer.insert(Int3(i+1,j,k));
        if(j-1>=jj0_in  && phi[k][j-1][i] >= bar) next_layer.insert(Int3(i,j-1,k));
        if(j+1<jjmax_in && phi[k][j+1][i] >= bar) next_layer.insert(Int3(i,j+1,k));
        if(k-1>=kk0_in  && phi[k-1][j][i] >= bar) next_layer.insert(Int3(i,j,k-1));
        if(k+1<kkmax_in && phi[k+1][j][i] >= bar) next_layer.insert(Int3(i,j,k+1));
      }
    }

    this_layer = next_layer;
  }

  MPI_Allreduce(MPI_IN_PLACE, &max_dist, 1, MPI_DOUBLE, MPI_MAX, comm);

  CandidatesIndex_n.RestoreDataPointerToLocalVector();
  ClosestPointIndex.RestoreDataPointerToLocalVector(); //cannot communicate, because "closest_points" do not.
  Phi.RestoreDataPointerAndInsert(); 

  return max_dist;
}

//-------------------------------------------------------------------------

void
Intersector::FindColorBoundary(int this_color, std::vector<int> &status)
{

  // Step 0. Check if the domain (not just this subdomain) has the input color
  if((this_color==1 && !hasInlet) || (this_color==2 && !hasInlet2) || 
     (this_color==3 && !hasOutlet)|| this_color>3 || this_color==0 ||
     this_color<-nRegions) {
    print_error("*** Error: (FindColorBoundary) Unable to find the boundary of an invalid color %d.\n",
                this_color);
    exit_mpi();
  }

  vector<double> &x_glob(global_mesh.x_glob);
  vector<double> &y_glob(global_mesh.y_glob);
  vector<double> &z_glob(global_mesh.z_glob);
  vector<double> &dx_glob(global_mesh.dx_glob);
  vector<double> &dy_glob(global_mesh.dy_glob);
  vector<double> &dz_glob(global_mesh.dz_glob);


  // Step 1. Preparation
  vector<Vec3D>&  Xs(surface.X);
  vector<Int3>&   Es(surface.elems);
  vector<Vec3D>&  Ns(surface.elemNorm);

  vector<int> positive_side(Es.size(), 0); //0 means NOT facing "this_color"
  vector<int> negative_side(Es.size(), 0);


  //Step 2. Build a small (layer=0) local scope. No need to create tree.
  //        Also build a global scope and tree that contain all the triangles
  Vec3D subDmin(0.0), subDmax(0.0);
  subDmin[0] = x_glob[i0] - 0.5*dx_glob[i0] - 2.0*half_thickness;
  subDmin[1] = y_glob[j0] - 0.5*dy_glob[j0] - 2.0*half_thickness;
  subDmin[2] = z_glob[k0] - 0.5*dz_glob[k0] - 2.0*half_thickness;
  subDmax[0] = x_glob[imax-1] + 0.5*dx_glob[imax-1] + 2.0*half_thickness;
  subDmax[1] = y_glob[jmax-1] + 0.5*dy_glob[jmax-1] + 2.0*half_thickness;
  subDmax[2] = z_glob[kmax-1] + 0.5*dz_glob[kmax-1] + 2.0*half_thickness;

  vector<int> local_scope;
  vector<MyTriangle> global_scope;
  global_scope.reserve(Es.size());

  for(int e=0; e<(int)Es.size(); e++) {
    MyTriangle tri(e, Xs[Es[e][0]], Xs[Es[e][1]], Xs[Es[e][2]]);
    global_scope.push_back(tri);

    bool inside = true;
    for(int i=0; i<3; i++) {
      if(tri.val(i) > subDmax[i] || tri.val(i) + tri.width(i) < subDmin[i]) {
        inside = false;
        break;
      }
    }
    if(inside)
      local_scope.push_back(e);
  }

  KDTree<MyTriangle,3> global_tree(global_scope.size(), global_scope.data());
 

  // Step 3. Check both sides 
  Vec3D*** coords  = (Vec3D***) coordinates.GetDataPointer();
  double*** color  = Color.GetDataPointer();

  int bandwidth = 2;
  int nMaxCand = 1000;
  vector<MyTriangle> tmp(nMaxCand);

  for(int side=0; side<2; side++) {

    double disp = 1.5*half_thickness;
    if(side==1) disp *= -1;

    // loop through triangles within the "local scope"
    for(int i=0; i<(int)local_scope.size(); i++) {

      // find point "p" with lofting
      int triangle_id = local_scope[i];
      Int3 &nod(Es[triangle_id]);
      Vec3D p = (Xs[nod[0]]+Xs[nod[1]]+Xs[nod[2]])/3.0 + disp*Ns[triangle_id];

      // locate the node that is closest to p
      Int3 ijk = global_mesh.FindClosestNodeToPoint(p, false);

      // find candidates (nodes with "this_color" for this point)
      vector<Int3> mycands;
      Int3 ijk1;
      for(int dk=-bandwidth; dk<=bandwidth; dk++) {
        ijk1[2] = ijk[2] + dk;
        for(int dj=-bandwidth; dj<=bandwidth; dj++) {
          ijk1[1] = ijk[1] + dj;
          for(int di=-bandwidth; di<=bandwidth; di++) {
            ijk1[0] = ijk[0] + di; 
            if(coordinates.IsHereOrInternalGhost(ijk1[0], ijk1[1], ijk1[2]) &&
               color[ijk1[2]][ijk1[1]][ijk1[0]] == this_color)
              mycands.push_back(ijk1);
          }
        }
      }
      std::sort(mycands.begin(), mycands.end(), 
                [&](Int3 a, Int3 b) {
                   return (coords[a[2]][a[1]][a[0]] - p).norm() < (coords[b[2]][b[1]][b[0]] - p).norm();} );

      // find intersections using the global tree (i.e. all the triangles)
      for(auto it = mycands.begin(); it != mycands.end(); it++) {
       
        Vec3D &q(coords[(*it)[2]][(*it)[1]][(*it)[0]]);

        Vec3D bmin(0.0), bmax(0.0);
        for(int s=0; s<3; s++) {
          bmin[s] = std::min(q[s], p[s]) - half_thickness;
          bmax[s] = std::max(q[s], p[s]) + half_thickness;
        } 
        int nFound = FindCandidatesInBox(&global_tree, bmin, bmax, tmp, nMaxCand);
        bool intersect = false; //if nFound = 0 (unlikely), this should be "false"
        for(int tri=0; tri<nFound; tri++) {
          int id = tmp[tri].trId();
          Int3& nodes(Es[id]);
          if(GeoTools::LineSegmentIntersectsTriangle(p, q, Xs[nodes[0]], Xs[nodes[1]], Xs[nodes[2]])) {
            intersect = true;
            break;  //intersecting one triangle means intersecting the whole surface 
          } 
        }  
        if(!intersect) { //if one <p,q> does not intersect the interface, it means they are on the same "side"
          if(side==0) 
            positive_side[triangle_id] = 1;
          else
            negative_side[triangle_id] = 1;

          break;
        }
      }

    }

    // assemble data
    if(side==0)
      MPI_Allreduce(MPI_IN_PLACE, positive_side.data(), positive_side.size(), MPI_INT, MPI_MAX, comm);
    else
      MPI_Allreduce(MPI_IN_PLACE, negative_side.data(), negative_side.size(), MPI_INT, MPI_MAX, comm);

  }

  
  // Step 4. Finalize output
  status.assign(Es.size(), 0);
  for(int i=0; i<(int)Es.size(); i++) {
    if(positive_side[i]>0) {
      if(negative_side[i]>0)
        status[i] = 3;
      else
        status[i] = 1;
    } else if(negative_side[i]>0)
      status[i] = 2;
  }
  

  // Clean-up
  coordinates.RestoreDataPointerToLocalVector();
  Color.RestoreDataPointerToLocalVector();
 
}

//-------------------------------------------------------------------------

int
Intersector::FindEdgeIntersectionsWithTriangles(Vec3D &x0, int i, int j, int k, int dir, double len, MyTriangle* tri, int nTri,
                                                IntersectionPoint &xf, IntersectionPoint &xb)
{
  vector<Vec3D>&  Xs(surface.X);
  vector<Int3>&   Es(surface.elems);

  vector<pair<double, IntersectionPoint> > X; //pairs distance with intersection point info

  double dist;
  Vec3D xi; //barycentric coords of the projection point
  for(int iTri=0; iTri<nTri; iTri++) {
    int id = tri[iTri].trId();
    Int3& nodes(Es[id]);
    bool found = GeoTools::LineSegmentIntersectsTriangle(x0, dir, len, Xs[nodes[0]], Xs[nodes[1]], Xs[nodes[2]],
                                                         &dist, NULL, &xi);
    if(found) //store the result
      X.push_back(std::make_pair(dist, IntersectionPoint(i,j,k,dir,dist,id,xi)));
  }

  if(X.empty())
    return 0;

  double x0_2_xf = DBL_MAX; 
  double x0_2_xb = -1.0;
  for(auto it = X.begin(); it != X.end(); it++) {
    if(x0_2_xf > it->first) {
      x0_2_xf = it->first;
      xf = it->second;
    }
    if(x0_2_xb < it->first) {
      x0_2_xb = it->first;
      xb = it->second;
    }
  }

  return X.size();
}

//-------------------------------------------------------------------------

bool
Intersector::Intersects(Vec3D &X0, Vec3D &X1)
{

  if(!tree_n) //this subdomain is not even close to the embedded surface.
    return false;

  assert(nLayer>=1);

  // Step 1: Find candidates using the KDTree
  Vec3D bmin(std::min(X0[0],X1[0]), std::min(X0[1],X1[1]), std::min(X0[2],X1[2]));
  Vec3D bmax(std::max(X0[0],X1[0]), std::max(X0[1],X1[1]), std::max(X0[2],X1[2]));
  bmin -= half_thickness;
  bmax += half_thickness;

  int maxCand = 1000;
  vector<MyTriangle> cands(maxCand);
  int found = tree_n->findCandidatesInBox(bmin,bmax,cands.data(),maxCand);
  if(found>maxCand) {
    maxCand = found;
    cands.resize(maxCand);
    found = tree_n->findCandidatesInBox(bmin,bmax,cands.data(),maxCand);
  }

  if(!found)
    return false;

  // Step 2: Check if X0 or X1 is occluded (-->intersection)
  int tid; //not used
  if(IsPointOccludedByTriangles(X0, cands.data(), found, half_thickness, tid))
    return true;
  if(IsPointOccludedByTriangles(X1, cands.data(), found, half_thickness, tid))
    return true;

  // Step 3: Check for intersection
  vector<Vec3D>& Xs(surface.X);
  vector<Int3>&  Es(surface.elems);
  for(int i=0; i<found; i++) {
    Int3 &nodes(Es[cands[i].trId()]);
    if(GeoTools::LineSegmentIntersectsTriangle(X0, X1, Xs[nodes[0]], Xs[nodes[1]], Xs[nodes[2]]))
      return true;
  }

  return false;
}

//-------------------------------------------------------------------------


//-------------------------------------------------------------------------

