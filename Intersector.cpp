#include<Intersector.h>
#include<GeoTools.h>
#include<Intersections.h>
using std::pair;
using std::vector;

extern int verbose;
extern double domain_diagonal;
//-------------------------------------------------------------------------

Intersector::Intersector(MPI_Comm &comm_, DataManagers3D &dms_, EmbeddedSurfaceData &iod_surface_,
                         TriangulatedSurface &surface_, SpaceVariable3D &coordinates_, 
                         vector<GhostPoint> &ghost_nodes_inner_, vector<GhostPoint> &ghost_nodes_outer_,
                         vector<double> &x_, vector<double> &y_, vector<double> &z_,
                         vector<double> &dx_, vector<double> &dy_, vector<double> &dz_)
           : comm(comm_), iod_surface(iod_surface_), surface(surface_), tree(NULL),
             coordinates(coordinates_), 
             ghost_nodes_inner(ghost_nodes_inner_), ghost_nodes_outer(ghost_nodes_outer_),
             x_glob(x_), y_glob(y_), z_glob(z_), dx_glob(dx_), dy_glob(dy_), dz_glob(dz_),
             BBmin(comm_, &(dms_.ghosted1_3dof)),
             BBmax(comm_, &(dms_.ghosted1_3dof)),
             TMP(comm_, &(dms_.ghosted1_1dof)),
             TMP2(comm_, &(dms_.ghosted1_1dof)),
             CandidatesIndex(comm_, &(dms_.ghosted1_1dof)),
             ClosestPointIndex(comm_, &(dms_.ghosted1_1dof)),
             XForward(comm_, &(dms_.ghosted1_3dof)),
             XBackward(comm_, &(dms_.ghosted1_3dof)),
             Phi(comm_, &(dms_.ghosted1_1dof)),
             Sign(comm_, &(dms_.ghosted1_1dof)),
             floodfiller(comm_, dms_, ghost_nodes_inner_, ghost_nodes_outer_)
{

  CandidatesIndex.SetConstantValue(-1, true);
  ClosestPointIndex.SetConstantValue(-1, true);
  XForward.SetConstantValue(-1, true);
  XBackward.SetConstantValue(-1, true);

  half_thickness = 0.5*iod_surface.tracker.surface_thickness;

  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);
  coordinates.GetInternalGhostedCornerIndices(&ii0_in, &jj0_in, &kk0_in, &iimax_in, &jjmax_in, &kkmax_in);
  coordinates.GetGlobalSize(&NX, &NY, &NZ);

  // Set the capacity of internal vectors, so we don't frequently reallocate memory
  int capacity = (imax-i0)*(jmax-j0)*(kmax-k0)/4; //should be big enough
  intersections.reserve(capacity);
  candidates.reserve(capacity*2);
  scope.reserve(surface.elems.size());

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
  BuildNodalAndSubdomainBoundingBoxes(1);

}

//-------------------------------------------------------------------------

Intersector::~Intersector()
{
  if(tree)
    delete tree;
}

//-------------------------------------------------------------------------

void
Intersector::Destroy()
{
  floodfiller.Destroy();

  BBmin.Destroy();
  BBmax.Destroy();
  TMP.Destroy();
  TMP2.Destroy();
  CandidatesIndex.Destroy();
  ClosestPointIndex.Destroy();
  XForward.Destroy();
  XBackward.Destroy();
  Phi.Destroy();
  Sign.Destroy();
}

//-------------------------------------------------------------------------

void
Intersector::TrackSurfaceFullCourse(bool &hasInlet, bool &hasOutlet, bool &hasOcc, int &nRegions, int phi_layers)
{
  assert(phi_layers>=1);

  BuildNodalAndSubdomainBoundingBoxes(1); //1 layer
  BuildSubdomainScopeAndKDTree();
  FindNodalCandidates();
  FindIntersections(true);
  FloodFillColors(hasInlet, hasOutlet, hasOcc, nRegions);
  CalculateUnsignedDistanceNearSurface(phi_layers, phi_layers==1);

  Sign.StoreMeshCoordinates(coordinates);
  Sign.WriteToVTRFile("Sign_filled.vtr", "color");
  Phi.StoreMeshCoordinates(coordinates);
  Phi.WriteToVTRFile("Phi.vtr", "phi");
  fprintf(stderr,"Got here!\n");
  exit_mpi();
}

//-------------------------------------------------------------------------

void
Intersector::BuildNodalAndSubdomainBoundingBoxes(int nLayer)
{
  double tol = 0.1;  //i.e. tolerance = 10% of element size (plus 50% of surface thickness)
  double thicker = (1.0+tol)*half_thickness;

  Vec3D*** bbmin  = (Vec3D***) BBmin.GetDataPointer();
  Vec3D*** bbmax  = (Vec3D***) BBmax.GetDataPointer();

  double delta;
  for(int k=kk0_in; k<kkmax_in; k++)
    for(int j=jj0_in; j<jjmax_in; j++)
      for(int i=ii0_in; i<iimax; i++) {

        delta = tol*dx_glob[i] + thicker;
        bbmin[k][j][i][0] = x_glob[std::max(   0, i-nLayer)] - tol;
        bbmax[k][j][i][0] = x_glob[std::min(NX-1, i+nLayer)] + tol;

        delta = tol*dy_glob[j] + thicker;
        bbmin[k][j][i][1] = y_glob[std::max(   0, j-nLayer)] - tol;
        bbmax[k][j][i][1] = y_glob[std::min(NY-1, j+nLayer)] + tol;

        delta = tol*dz_glob[k] + thicker;
        bbmin[k][j][i][2] = z_glob[std::max(   0, k-nLayer)] - tol;
        bbmax[k][j][i][2] = z_glob[std::min(NZ-1, k+nLayer)] + tol;
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
Intersector::BuildSubdomainScopeAndKDTree()
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

//  fprintf(stderr,"scope size: %d.\n", (int)scope.size());

  // build the tree
  if(tree)
    delete tree;

  if(scope.size()!=0)
    tree = new KDTree<MyTriangle,3>(scope.size(), scope.data());
  else
    tree = NULL;
}

//-------------------------------------------------------------------------

void
Intersector::FindNodalCandidates()
{

  candidates.clear();

  int nMaxCand = 500; //will increase if necessary

  Vec3D*** bbmin   = (Vec3D***) BBmin.GetDataPointer();
  Vec3D*** bbmax   = (Vec3D***) BBmax.GetDataPointer();
  double*** candid = CandidatesIndex.GetDataPointer();

  MyTriangle *tmp = new MyTriangle[nMaxCand];
  
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

  delete [] tmp;
}

//-------------------------------------------------------------------------

void
Intersector::FindIntersections(bool with_nodal_cands) //also finds occluded and first layer nodes
{

  Vec3D*** coords  = (Vec3D***) coordinates.GetDataPointer();
  Vec3D*** xf      = (Vec3D***) XForward.GetDataPointer();
  Vec3D*** xb      = (Vec3D***) XBackward.GetDataPointer();
  double*** candid = with_nodal_cands ? CandidatesIndex.GetDataPointer() : NULL;
  double*** sign   = Sign.GetDataPointer();
  
  double*** occid  = TMP.GetDataPointer(); //occluding triangle id
  double*** layer  = TMP2.GetDataPointer(); //"layer" of each node: 0(occluded), 1, or -1 (unknown)

  //Clear previous values
  intersections.clear();

  //Preparation
  Vec3D tol(half_thickness*1.5, half_thickness*1.5, half_thickness*1.5); //a tolerance

  int max_left   = 500; //will increase if necessary
  int max_bottom = 500; 
  int max_back   = 500; 

  MyTriangle *tmp_left = new MyTriangle[max_left];
  int found_left;
  MyTriangle *tmp_bottom = new MyTriangle[max_bottom];
  int found_bottom;
  MyTriangle *tmp_back = new MyTriangle[max_back];
  int found_back;

  IntersectionPoint xp0, xp1;

  // ----------------------------------------------------------------------------
  // Find occluded nodes and intersections. Build occluded and firstLayer 
  // ----------------------------------------------------------------------------
  firstLayer.clear();
  // We only deal with edges whose vertices are both in the real domain
  for(int k=k0; k<kkmax_in; k++)
    for(int j=j0; j<jjmax_in; j++)
      for(int i=i0; i<iimax_in; i++) {
  
        // start with assuming it is outside
        sign[k][j][i] = 1;

        // start with a meaningless triangle id
        occid[k][j][i] = -1;

        layer[k][j][i] = -1;

        if(!tree) //this subdomain is entirely away from surface, just set sign, occid, and layer to default values
          continue; 

        if(candid && k<kmax && j<jmax && i<imax && //candid has a valid value @ i,j,k
           candid[k][j][i] < 0)  //no nodal candidates, intersection impossible, occlusion also impossible
          continue;
 

        //--------------------------------------------
        // Find candidates for left/bottom/back edges
        //--------------------------------------------
        found_left = found_bottom = found_back = 0;

        if(i-1>=0) { //the edge [k][j][i-1] -> [k][j][i] is inside the physical domain
          found_left = FindCandidatesInBox(tree, coords[k][j][i-1] - tol, coords[k][j][i] + tol, tmp_left, max_left);
        }

        if(j-1>=0) { //the edge [k][j-1][i] -> [k][j][i] is inside the physical domain
          found_bottom = FindCandidatesInBox(tree, coords[k][j-1][i] - tol, coords[k][j][i] + tol, tmp_bottom, max_bottom);
        }

        if(k-1>=0) { //the edge [k-1][j][i] -> [k][j][i] is inside the physical domain
          found_back = FindCandidatesInBox(tree, coords[k-1][j][i] - tol, coords[k][j][i] + tol, tmp_back, max_back);
        }

        if(j==0 && i==6) {
          fprintf(stderr,"found_left = %d, found_bottom = %d, found_back = %d.\n", found_left, found_bottom, found_back);
          fprintf(stderr,"candid = %e.\n", candid[k][j][i]);
        } 

        //--------------------------------------------
        // Check if (i,j,k) is occluded
        //--------------------------------------------
        int tid(-1);
        if ((found_left>0   && IsPointOccludedByTriangles(coords[k][j][i], tmp_left,   found_left,   half_thickness, tid)) ||
            (found_bottom>0 && IsPointOccludedByTriangles(coords[k][j][i], tmp_bottom, found_bottom, half_thickness, tid)) ||
            (found_back>0   && IsPointOccludedByTriangles(coords[k][j][i], tmp_back,   found_back,   half_thickness, tid))) {
          sign[k][j][i] = 0;
          occid[k][j][i] = tid;
          layer[k][j][i] = 0;          
        }

        //--------------------------------------------
        // Find intersections (left, bottom, and back edges)
        //--------------------------------------------
        // left
        if(found_left) {
          int count = FindEdgeIntersectionsWithTriangles(coords[k][j][i-1], i-1, j, k, 0, 
                          coords[k][j][i][0] - coords[k][j][i-1][0], tmp_left, found_left, xp0, xp1);
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
                          coords[k][j][i][1] - coords[k][j-1][i][1], tmp_bottom, found_bottom, xp0, xp1);
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
                          coords[k][j][i][2] - coords[k-1][j][i][2], tmp_back, found_back, xp0, xp1);
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

  // Exchange Sign and TMP so internal ghost nodes are accounted for
  if(candid) CandidatesIndex.RestoreDataPointerToLocalVector();
  Sign.RestoreDataPointerAndInsert();
  TMP.RestoreDataPointerAndInsert();
  TMP2.RestoreDataPointerAndInsert();

  occid  = TMP.GetDataPointer();

  // ----------------------------------------------------------------------------
  // Make sure all edges connected to occluded nodes have intersections
  // ----------------------------------------------------------------------------
  bool ijk_occluded;
  for(int k=k0; k<kkmax_in; k++)
    for(int j=j0; j<jjmax_in; j++)
      for(int i=i0; i<iimax_in; i++) {
 
        ijk_occluded = (occid[k][j][i]>=0);
        if(i-1>=0) { //left edge within physical domain
          
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


        if(j-1>=0) { //bottom edge within physical domain
          
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


        if(k-1>=0) { //back edge within physical domain
          
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
  Sign.StoreMeshCoordinates(coordinates);
  Sign.WriteToVTRFile("Sign.vtr", "sign");
*/
}

//-------------------------------------------------------------------------

int
Intersector::FindCandidatesInBox(KDTree<MyTriangle, 3>* mytree, Vec3D bbmin, Vec3D bbmax, 
                                 MyTriangle* tmp, int& maxCand)
{
  int found = mytree->findCandidatesInBox(bbmin, bbmax, tmp, maxCand);
  if(found>maxCand) {
    maxCand = found; 
    delete [] tmp;
    tmp = new MyTriangle[maxCand];
    found = mytree->findCandidatesInBox(bbmin, bbmax, tmp, maxCand);
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

int
Intersector::FloodFillColors(bool &hasInlet, bool &hasOutlet, bool &hasOcc, int &nRegions)
{
  // ----------------------------------------------------------------
  // Call floodfiller to do the work.
  // ----------------------------------------------------------------
  int nColors = floodfiller.FillBasedOnEdgeObstructions(XForward, -1/*xf==-1 means no intersection*/, occluded, Sign);

  // ----------------------------------------------------------------
  // Now we need to convert the color map to what we want: 0~occluded, 1, 2,...: regions connected to Dirichlet
  // boundaries (inlet/outlet/farfield). 1=inlet, 2=outlet. -1,-2,...: inside enclosures
  // ----------------------------------------------------------------
  double*** sign   = Sign.GetDataPointer();

  // First, we need to find the current colors for inlet, outlet
  std::set<int> inlet_color; 
  std::set<int> outlet_color;  //In most cases, both sets should only contain 0 or 1 member. Otherwise, it means 
                               //one type of boundary (e.g., inlet) is separated by the embedded surface into multiple
                               //disconnected regions. In this case, we have to merge the colors into one.
  for(auto it = ghost_nodes_outer.begin(); it != ghost_nodes_outer.end(); it++) {
    if(it->type_projection != GhostPoint::FACE)
      continue;
    if(it->bcType == MeshData::INLET) {
      Int3& ijk(it->image_ijk);
      inlet_color.insert(sign[ijk[2]][ijk[1]][ijk[0]]);
    } else if(it->bcType == MeshData::OUTLET) {
      Int3& ijk(it->image_ijk);
      outlet_color.insert(sign[ijk[2]][ijk[1]][ijk[0]]);
    }
  }

  int max_inlet_color(-1);
  for(auto it = inlet_color.begin(); it != inlet_color.end(); it++)
    max_inlet_color = std::max(max_inlet_color, *it);
  MPI_Allreduce(MPI_IN_PLACE, &max_inlet_color, 1, MPI_INT, MPI_MAX, comm);

  vector<int> in_colors;
  if(max_inlet_color!=-1) {
    in_colors.resize(max_inlet_color+1, -1);
    for(auto it = inlet_color.begin(); it != inlet_color.end(); it++)
      in_colors[*it] = 1;
    MPI_Allreduce(MPI_IN_PLACE, in_colors.data(), in_colors.size(), MPI_INT, MPI_MAX, comm);
  }

  if(in_colors.size()==0) 
    print_warning("Warning: Intersector did not find an inlet/farfield boundary of the M2C domain.\n");

  if(in_colors.size()>0 && in_colors[0] == 1 && verbose>1)
    print_warning("Warning: Found occluded node(s) near an inlet or farfield boundary.");

  int max_outlet_color(-1);
  for(auto it = outlet_color.begin(); it != outlet_color.end(); it++) 
    max_outlet_color = std::max(max_outlet_color, *it);
  MPI_Allreduce(MPI_IN_PLACE, &max_outlet_color, 1, MPI_INT, MPI_MAX, comm);

  vector<int> out_colors;
  if(max_outlet_color!=-1) {
    out_colors.resize(max_outlet_color+1, -1);
    for(auto it = outlet_color.begin(); it != outlet_color.end(); it++)
      out_colors[*it] = 1;
    MPI_Allreduce(MPI_IN_PLACE, out_colors.data(), out_colors.size(), MPI_INT, MPI_MAX, comm);
  }
 
  if(out_colors.size()>0 && out_colors[0] == 1 && verbose>1)
    print_warning("Warning: Found occluded node(s) near an outlet or farfield boundary.");

  // Convert colors 
  std::map<int,int> old2new;
  for(int i=0; i<in_colors.size(); i++)
    if(in_colors[i] == 1)
     old2new[i] = 1; //inlet_color;
  for(int i=0; i<out_colors.size(); i++)
    if(out_colors[i] == 1) {
     assert(old2new.find(i) == old2new.end());
     old2new[i] = 2; //outlet_color;
    }

  // give negative colors to enclusures (i.e. regions not connected to farfield)
  int tmp_counter = 0;
  for(int i=1; i<nColors+1; i++)
    if(old2new.find(i) == old2new.end())
      old2new[i] = --tmp_counter;

  for(auto it = old2new.begin(); it != old2new.end(); it++)
    fprintf(stderr,"old2new: %d --> %d.\n", it->first, it->second);

  int total_occluded = 0;

  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        if(sign[k][j][i] != 0)
          sign[k][j][i] = old2new[sign[k][j][i]];
        else
          total_occluded++;
      }

  // statistics
  MPI_Allreduce(MPI_IN_PLACE, &total_occluded, 1, MPI_INT, MPI_SUM, comm);

  if(total_occluded>0) 
    hasOcc = true; //i.e. one zero color (for occluded)

  hasInlet = hasOutlet = false;
  nRegions = 0;
  for(auto it = old2new.begin(); it != old2new.end(); it++)
    if(it->second==1)
      hasInlet = true;
    else if(it->second==2)
      hasOutlet = true;
    else if(it->second<0)
      nRegions++;

  fprintf(stderr,"sign[0][0][6] = %e.\n", sign[0][0][6]);

  Sign.RestoreDataPointerAndInsert();

  return int(hasInlet) + int(hasOutlet) + int(hasOcc) + nRegions;
}

//-------------------------------------------------------------------------

void
Intersector::RefillAfterSurfaceUpdate(bool nodal_cands_calculated)
{
  // "FindSweptNodes" should be called before calling this function. So, in most cases, "nodal_cands_calculated" should be true
  
  // Assuming occluded nodes already have the correct sign. We do not deal with them here.
  
  if(!nodal_cands_calculated) {
    BuildNodalAndSubdomainBoundingBoxes(1);
    BuildSubdomainScopeAndKDTree();
    FindNodalCandidates();
  }

  double*** candid = CandidatesIndex.GetDataPointer();
  
  //add swepted nodes to nodes2fill (including internal ghosts)
  std::set<Int3> nodes2fill = swept;
  std::set<Int3> nodes2fill2;

  // go over swept nodes, correct their signs
  int i0,j0,k0;
  vector<Vec3D>&  Xs(surface.X);
  vector<Int3>&   Es(surface.elems);

  int BAD_SIGN = -999999; //used temporarily.

  int max_it = 100;
  double*** sign = NULL;
  int total_remaining_nodes;
  for(int iter = 0; iter < max_it; iter++) {

    total_remaining_nodes = nodes2fill.size();
    MPI_Allreduce(MPI_IN_PLACE, &total_remaining_nodes, 1, MPI_INT, MPI_SUM, comm);
    if(total_remaining_nodes == 0) //Yeah
      break;
  
    if(!sign) //first iteration
      sign = Sign.GetDataPointer();

    nodes2fill2 = nodes2fill;

    for(auto it = nodes2fill.begin(); it != nodes2fill.end(); it++) {
      i0 = (*it)[0];
      j0 = (*it)[1];
      k0 = (*it)[2]; 

      if(!coordinates.IsHere(i0,j0,k0,false))
        continue; //this is an internal ghost. we let its owner fix it.

      sign[k0][j0][i0] = BAD_SIGN;

      // if this node is swept, candidates must exist. Otherwise, the surface moved too much!
      assert(candid[k0][j0][i0]>=0);
      vector<MyTriangle> &cands(candidates[candid[k0][j0][i0]].second);
      assert(cands.size()>0);

      //go over first layer neighbors, find a nonblocked reliable neighbor
      for(int k=k0-1; k<=k0+1; k++)
        for(int j=j0-1; j<=j0+1; j++)
          for(int i=i0-1; i<=i0+1; i++) {

            if(!coordinates.OutsidePhysicalDomain(i,j,k) || (i==i0 && j==j0 && k==k0))
              continue; //this neighbor is out of physical domain, or just the same node

            if(nodes2fill2.find(Int3(i,j,k)) != nodes2fill2.end()) //uses "nodes2fill2", the latest updated list
              continue; //this neighbor is in trouble as well...

            if(sign[k][j][i] == 0) //this neighbor is occluded (naturally or "forced"). Either way, it cannot be used.
              continue;

            bool blocked = false;
            for(auto it2 = cands.begin(); it2 != cands.end(); it2++) {
              int id = it2->trId();
              Int3 &nodes(Es[id]);
              blocked = GeoTools::LineSegmentIntersectsTriangle(Vec3D(x_glob[i],y_glob[j],z_glob[k]),//inside physical domain (safe)
                                                                Vec3D(x_glob[i0],y_glob[j0],z_glob[k0]), 
                                                                Xs[nodes[0]], Xs[nodes[1]], Xs[nodes[2]]);
              if(blocked)
                break;
            }
            if(blocked)
              continue; //this neighbor is blocked from me by the surface

            // If the above checks are all passed, this neighbor is in the same region as me.
            sign[k0][j0][i0] = sign[k][j][i];
            nodes2fill2.erase(nodes2fill2.find(*it)); //erase from nodes2fill2, NOT from nodes2fill (would mess up pointer!)
            goto DONE_WITH_THIS_NODE; 
         }

      DONE_WITH_THIS_NODE:
      continue; //need this to make "GOTO" work
    }
    nodes2fill = nodes2fill2;

    Sign.RestoreDataPointerAndInsert();

    // Get new data
    sign = Sign.GetDataPointer();

    // remove filled internal ghosts (filled by their owners) from "nodes2fill"
    for(auto it = nodes2fill.begin(); it != nodes2fill.end(); it++) {
      i0 = (*it)[0];
      j0 = (*it)[1];
      k0 = (*it)[2]; 
      if(!coordinates.IsHere(i0,j0,k0,false)) { //this is an internal ghost
        if(sign[k0][j0][i0] != BAD_SIGN) //must have been fixed by its owner
          nodes2fill2.erase(nodes2fill2.find(*it)); 
      }
    }
    if(nodes2fill.size() != nodes2fill2.size()) //some node(s) got erased
      nodes2fill = nodes2fill2;

  }

  // fill remaining nodes as FORCED occluded
  imposed_occluded.clear();
  if(total_remaining_nodes>0) {
    print_warning("Warning: Found %d unresolved nodes after performing %d iterations of refill. Setting them to be occluded.\n",
                   total_remaining_nodes, max_it);
    imposed_occluded = nodes2fill;
    for(auto it = imposed_occluded.begin(); it != imposed_occluded.end(); it++) 
      sign[(*it)[2]][(*it)[1]][(*it)[0]] = 0; //set it to occluded. BUT NO NEW INTERSECTIONS!
      //sign must be valid (i.e. not NULL) if total_remaining_nodes>0
  }

  if(sign)
    Sign.RestoreDataPointerToLocalVector();

  CandidatesIndex.RestoreDataPointerToLocalVector();

}

//-------------------------------------------------------------------------
 
void
Intersector::FindSweptNodes(std::vector<Vec3D> &X0, bool nodal_cands_calculated)
{
  if(!nodal_cands_calculated) {
    BuildNodalAndSubdomainBoundingBoxes(1);
    BuildSubdomainScopeAndKDTree();
    FindNodalCandidates();
  }
  
  swept.clear();

  double*** candid = CandidatesIndex.GetDataPointer();

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
    vector<MyTriangle> &cands(candidates[candid[k][j][i]].second);
    assert(cands.size()>0);

    for(auto it2 = cands.begin(); it2 != cands.end(); it2++) {
      int id = it2->trId();
      Int3 &nodes(Es[id]);
      Vec3D coords(x_glob[(*it)[0]], y_glob[(*it)[1]], z_glob[(*it)[2]]); //inside physical domain (safe)
      if(GeoTools::IsPointSweptByTriangle(coords, X0[nodes[0]], X0[nodes[1]], X0[nodes[2]],
                                          Xs[nodes[0]], Xs[nodes[1]], Xs[nodes[2]], &collision_time, half_thickness, 
                                          NULL, NULL, &(As[id]), &(Ns[id])))
        swept.insert(*it);
    }
  }
      
  CandidatesIndex.RestoreDataPointerToLocalVector();
}

//-------------------------------------------------------------------------

void
Intersector::CalculateUnsignedDistanceNearSurface(int nLayer, bool nodal_cands_calculated)
{

  if(!nodal_cands_calculated) {
    BuildNodalAndSubdomainBoundingBoxes(nLayer);
    BuildSubdomainScopeAndKDTree();
    FindNodalCandidates();
  }

  double*** candid = CandidatesIndex.GetDataPointer();
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
      vector<MyTriangle> &cands(candidates[candid[k][j][i]].second);
      assert(cands.size()>0);

      double dist = DBL_MAX, new_dist;
      ClosestPoint cp(-1,DBL_MAX,xi); //initialize to garbage

      int id;
      for(int tri=0; tri<cands.size(); tri++) {
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


  CandidatesIndex.RestoreDataPointerToLocalVector();
  ClosestPointIndex.RestoreDataPointerToLocalVector(); //cannot communicate, because "closest_points" do not.
  Phi.RestoreDataPointerAndInsert(); 

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
  for(int i=0; i<nTri; i++) {
    int id = tri[i].trId();
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

