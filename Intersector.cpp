#include<Intersector.h>
#include<GeoTools.h>
using std::vector;

//-------------------------------------------------------------------------

Intersector::Intersector(MPI_Comm &comm_, DataManagers3D &dms_, EmbeddedSurfaceData &iod_surface_,
                         TriangulatedSurface &surface_, SpaceVariable3D &coordinates_, 
                         SpaceVariable3D &delta_xyz_, SpaceVariable3D &volume_,
                         vector<GhostPoint> &ghost_nodes_inner_, vector<GhostPoint> &ghost_nodes_outer_)
           : comm(comm_), iod_surface(iod_surface_), surface(surface_), tree(NULL),
             coordinates(coordinates_), delta_xyz(delta_xyz_), volume(volume_),
             ghost_nodes_inner(ghost_nodes_inner_), ghost_nodes_outer(ghost_nodes_outer_),
             BBmin(comm_, &(dms_.ghosted1_3dof)),
             BBmax(comm_, &(dms_.ghosted1_3dof)),
             CandidatesIndex(comm_, &(dms_.ghosted1_1dof)),
             XX(comm_, &(dms_.ghosted1_3dof)),
             Phi(comm_, &(dms_.ghosted1_1dof)),
             Sign(comm_, &(dms_.ghosted1_1dof))
{

  CandidatesIndex.SetConstantValue(-1, true);
  XX.SetConstantValue(-1, true);

  half_thickness = 0.5*iod_surface.tracker.surface_thickness;

  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);
  coordinates.GetGlobalSize(&NX, &NY, &NZ);

  // Set the capacity of internal vectors, so we don't frequently reallocate memory
  int capacity = (imax-i0)*(jmax-j0)*(kmax-k0)/4; //should be big enough
  intersections.reserve(capacity);
  occluded.reserve(capacity); 
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
  BuildNodalBoundingBoxes();

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
  BBmin.Destroy();
  BBmax.Destroy();
  XX.Destroy();
  Phi.Destroy();
  Sign.Destroy();
}

//-------------------------------------------------------------------------

void
Intersector::BuildNodalBoundingBoxes()
{
  double tol = 0.01;  //i.e. 1% tolerance

  Vec3D*** coords = (Vec3D***) coordinates.GetDataPointer();
  Vec3D*** bbmin  = (Vec3D***) BBmin.GetDataPointer();
  Vec3D*** bbmax  = (Vec3D***) BBmax.GetDataPointer();

  double delta;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        delta = tol*(coords[k][j][i][0] - coords[k][j][i-1][0]); 
        bbmin[k][j][i][0] = i-1>=0 ? coords[k][j][i-1][0] - tol : coords[k][j][i][0] - tol;
        delta = tol*(coords[k][j][i+1][0] - coords[k][j][i][0]); 
        bbmax[k][j][i][0] = i+1<NX ? coords[k][j][i+1][0] + tol : coords[k][j][i][0] + tol;

        delta = tol*(coords[k][j][i][1] - coords[k][j-1][i][1]); 
        bbmin[k][j][i][1] = j-1>=0 ? coords[k][j-1][i][1] - tol : coords[k][j][i][1] - tol;
        delta = tol*(coords[k][j+1][i][1] - coords[k][j][i][1]); 
        bbmax[k][j][i][1] = j+1<NY ? coords[k][j+1][i][1] + tol : coords[k][j][i][1] + tol;

        delta = tol*(coords[k][j][i][2] - coords[k-1][j][i][2]); 
        bbmin[k][j][i][2] = k-1>=0 ? coords[k-1][j][i][2] - tol : coords[k][j][i][2] - tol;
        delta = tol*(coords[k+1][j][i][2] - coords[k][j][i][2]); 
        bbmax[k][j][i][2] = k+1<NZ ? coords[k+1][j][i][2] + tol : coords[k][j][i][2] + tol;

      }

  subD_bbmin[0] = coords[kk0][jj0][ii0][0]             - tol*(coords[k0][j0][i0][0]           - coords[k0][j0][i0-1][0]);
  subD_bbmax[0] = coords[kkmax-1][jjmax-1][iimax-1][0] + tol*(coords[kmax-1][jmax-1][imax][0] - coords[kmax-1][jmax-1][imax-1][0]);
  subD_bbmin[1] = coords[kk0][jj0][ii0][1]             - tol*(coords[k0][j0][i0][1]           - coords[k0][j0-1][i0][1]);
  subD_bbmax[1] = coords[kkmax-1][jjmax-1][iimax-1][1] + tol*(coords[kmax-1][jmax][imax-1][1] - coords[kmax-1][jmax-1][imax-1][1]);
  subD_bbmin[2] = coords[kk0][jj0][ii0][2]             - tol*(coords[k0][j0][i0][2]           - coords[k0-1][j0][i0][2]);
  subD_bbmax[2] = coords[kkmax-1][jjmax-1][iimax-1][2] + tol*(coords[kmax][jmax-1][imax-1][2] - coords[kmax-1][jmax-1][imax-1][2]);

  coordinates.RestoreDataPointerToLocalVector();
  BBmin.RestoreDataPointerAndInsert();
  BBmax.RestoreDataPointerAndInsert();
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
      if(tri.val[i] > subD_bbmax[i] || tri.val[i] + tri.width[i] < subD_bbmin[i]) {
        inside = false;
        break;
      }
    }
    if(inside) 
      scope.push_back(tri); //creating a copy of "tri" and store it in scope
  }

  // build the tree
  if(tree)
    delete tree;
  tree = new KDTree<MyTriangle,3>(scope.size(), scope.data());
}

//-------------------------------------------------------------------------

void
Intersector::FindNodalCandidates()
{
  assert(tree);

  candidates.clear();

  int nMaxCand = 500; //will increase if necessary

  Vec3D*** coords  = (Vec3D***) coordinates.GetDataPointer();
  Vec3D*** bbmin   = (Vec3D***) BBmin.GetDataPointer();
  Vec3D*** bbmax   = (Vec3D***) BBmax.GetDataPointer();
  double*** candid = CandidatesIndex.GetDataPointer();

  MyTriangle *tmp = new MyTriange[nMaxCand];
  
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        // find candidates
        int nFound = tree.findCandidatesInBox(bbmin[k][j][i], bbmax[k][j][i], tmp, nMaxCand);
        if(nFound>nMaxCand) { //re-do w/ a bigger nMaxCand
          nMaxCand = nFound;
          delete [] tmp;
          tmp = new MyTriangle[nMaxCand];
          tree.findCandidatesInBox(bbmin[k][j][i], bbmax[k][j][i], tmp, nMaxCand);
        }

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

  coordinates.RestoreDataPointerToLocalVector();
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
  Vec3D*** xx      = (Vec3D***) XX.GetDataPointer();
  double*** candid = with_nodal_cands ? CandidatesIndex.GetDataPointer() : NULL;
  double*** sign   = Sign.GetDataPointer();
  

  //Preparation
  Vec3D tol(half_thickness*5, half_thickness*5, half_thickness*5); //a tolerance, more than enough

  int max_left   = 500; //will increase if necessary
  int max_bottom = 500; //will increase if necessary
  int max_back   = 500; //will increase if necessary

  MyTriangle *tmp_left = new MyTriangle[max_left];
  int found_left;
  MyTriangle *tmp_bottom = new MyTriangle[max_bottom];
  int found_bottom;
  MyTriangle *tmp_back= new MyTriangle[max_back];
  int found_back;

  for(int k=k0; k<kkmax; k++)
    for(int j=j0; j<jjmax; j++)
      for(int i=i0; i<iimax; i++) {
  
        if(coordinates.OutsidePhysicalDomain(i,j,k)) //this node is out of physical domain
          continue;

        // start with assuming it is outside
        sign[k][j][i] = 1;

        if(candid && k<kmax && j<jmax && i<imax && //candid has a valid value @ i,j,k
           candid[k][j][i] < 0)  //no nodal candidates
          continue;
 

        found_left = found_bottom = found_back = 0;

        if(i-1>=0) { //the edge [k][j][i-1] -> [k][j][i] is inside the physical domain

          if(candid && k<kmax && j<jmax && //candid has a valid value @ i-1,j,k
             candid[k][j][i-1] < 0) //no nodal candidates
            goto SKIP_LEFT;

          found_left = tree.findCandidatesInBox(coords[k][j][i-1] - tol, coords[k][j][i] + tol, 
                                                tmp_left, max_left);
          if(found_left>max_left) { //re-do w/ a bigger max 
            max_left = found_left;
            delete [] tmp_left;
            tmp_left = new MyTriangle[max_left];
            tree.findCandidatesInBox(coords[k][j][i-1] - tol, coords[k][j][i] + tol, tmp_left, max_left);
          }
        }

        SKIP_LEFT:

        if(j-1>=0) { //the edge [k][j-1][i] -> [k][j][i] is inside the physical domain

          if(candid && k<kmax && i<imax && //candid has a valid value @ i,j-1,k
             candid[k][j-1][i] < 0) //no nodal candidates
            goto SKIP_BOTTOM;

          found_bottom = tree.findCandidatesInBox(coords[k][j-1][i] - tol, coords[k][j][i] + tol, 
                                                  tmp_bottom, max_bottom);
          if(found_bottom>max_bottom) { //re-do w/ a bigger max
            max_bottom = found_bottom;
            delete [] tmp_bottom;
            tmp_bottom = new MyTriangle[max_bottom];
            tree.findCandidatesInBox(coords[k][j-1][i] - tol, coords[k][j][i] + tol, tmp_bottom, max_bottom);
          }
        }

        SKIP_BOTTOM:

        if(k-1>=0) { //the edge [k-1][j][i] -> [k][j][i] is inside the physical domain

          if(candid && j<jmax && i<imax && //candid has a valid value @ i,j,k-1
             candid[k-1][j][i] < 0) //no nodal candidates
            goto SKIP_BACK;

          found_back = tree.findCandidatesInBox(coords[k-1][j][i] - tol, coords[k][j][i] + tol, 
                                                tmp_back, max_back);
          if(found_back>max_back) { //re-do w/ a bigger max
            max_back = found_back;
            delete [] tmp_back;
            tmp_back = new MyTriangle[max_back];
            tree.findCandidatesInBox(coords[k-1][j][i] - tol, coords[k][j][i] + tol, tmp_back, max_back);
          }
        }

        SKIP_BACK:

        // check if (i,j,k) is occluded
        if(found_left>0) {
          GeoTools::IsPointInThickenedTriangle(...) I AM HERE!!!

        }


        // find intersection
      }

}

//-------------------------------------------------------------------------

void
Intersector::FindOcculudedNodes()

void
Intersector::FindIntersections()

void
Intersector::TagFirstLayerNodes();

void
Intersector::FloodFillStatus(); //by default, set sign = 1 (0 for occuluded) only do floodfill for sign = -1

//optional
void
Intersector::FindShortestDistanceForFirstLayer //don't forget to subtract half_distance

void
Intersector::FindShortestDistanceForOtherNodes //call level set reinitializer
