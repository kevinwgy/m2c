#ifndef _INTERSECTOR_H_
#define _INTERSECTOR_H_

#include<IoData.h>
#include<KDTree.h>
#include<TriangulatedSurface.h>
#include<FloodFill.h>

/****************************************************************
 * Class Intersector is responsible for tracking a triangulated
 * surface within a fixed Cartesian mesh. A collision-based
 * algorithm is used, similar to the one presented in
 * Wang et al., IJNMF, 2012. The intersector is able to find
 * (1) nodes covered by the interface (occluded nodes),
 * (2) all the edge-interface intersections,
 * (3) shortest distance from each node to the intersector (w/
 *     the help of level set reinitializer)
 * (4) closest point on the interface to each node (If solution is not
 *     unique, only one solution is found)
 * (5) nodes swept by the (dynamic) surface in one time step
 * (6) (if the interface is a closed surface) the inside/outside
 *     status of each node with respect to the surface.
 * Note: The above results are stored within this class.
 ***************************************************************/

class Intersector {

  //! A nested class that stores information about an intersection point between an edge and a triangle
  struct IntersectionPoint {
    Int3 n0; //!< the first (i.e., left, bottom, or back) node
    int dir; //!< the direction of the edge (0~x, 1~y, 2~z)
    double dist; //!< dist from n0 to the intersection point along dir
    int tid; //!< id of the triangle that it intersects
    double xi[3]; //!< barycentric coords of the intersection point within the triangle it intersects

    /** NOTE: In most cases, the point obtained using tid and xi is IDENTICAL to the point
     *        obtained using n0, dir, and dist. 
     *        However, the two points may not be the same when the intersection is IMPOSED to an edge because
     *        (1) one (or two) of the vertices of the edge is occluded and 
     *        (2) an intersection cannot be identified normally (i.e. w/o imposing thickness).
     *        In this case, the second point (obtained using n0, dir, and dist) would be the occluded vertex,
     *        and the first point is the closest point on the triangle to the occluded vertex.
     *        If both vertices of the edge are occluded. The two vertices are considered as two intersection
     *        points. */

    IntersectionPoint(int i, int j, int k, int dir_, double dist_, int tid_, double* xi_)
      : n0(Int3(i,j,k)), dir(dir_), dist(dist_), tid(tid_) {xi[0] = xi_[0]; xi[1] = xi_[1]; xi[2] = xi_[2];}

    IntersectionPoint &operator=(const IntersectionPoint& p2) {
      n0 = p2.n0;  dir = p2.dir;  dist = p2.dist;  tid = p2.tid;  for(int i=0; i<3; i++) xi[i] = p2.xi[i]; return *this;}
  };


  //! Utility class to find and store bounding boxes for triangles
  class MyTriangle {
    int id;
    double x[3], w[3];
  public:
    MyTriangle() {}
    MyTriangle(int id_, double x_[3], double w_[3]) : id(id_) {
      for(int j=0; j<3; j++) {x[j] = x_[j];  w[j] = w_[j];} 
    }
    MyTriangle(const MyTriangle& t2) : id(t2.id) {
      for(int j=0; j<3; j++) {x[j] = t2.x[j];  w[j] = t2.w[j];}
    }
    MyTriangle(int id_, Vec3D& node1, Vec3D& node2, Vec3D& node3) : id(id_) {
      for(int j=0; j<3; j++) {
        x[j] = std::min(std::min(node1[j], node2[j]), node3[j]);
        w[j] = std::max(std::max(node1[j], node2[j]), node3[j]) - x[j];
      }
    }
    double val(int i) const { return x[i]; }
    double width(int i) const { return w[i]; }
    int trId() const { return id; }
  };


  //! Utility class to store the closest point on the surface to a node (an arbitrary point)
  struct ClosestPoint {
    int tid;
    double dist; //!< dist to the node in question
    double xi[3]; //!< closest point is at xi[0]*node[0] + xi[1]*node[1] + xi[2]*node[2]

    ClosestPoint(int tid_, double dist_, double xi_[3]) : tid(tid_), dist(dist_) {
      for(int i=0; i<3; i++) xi[i] = xi_[i];}

    ClosestPoint &operator=(const ClosestPoint& p2) {
      tid = p2.tid;  dist = p2.dist;
      for(int i=0; i<3; i++) xi[i] = p2.xi[i]; return *this;}
  };

  MPI_Comm& comm;

  EmbeddedSurfaceData &iod_surface;

  //! A tool for flood-fill
  FloodFill floodfiller;

  //! Info about the triangulated surface
  TriangulatedSurface &surface; //!< the surface tracked by the intersector
  bool closed_surface; //!< whether the surface is closed AND normals are consistent
  double half_thickness; //!< half thickness of the surface

  std::vector<MyTriangle> scope; //!< triangles relevant to the current subdomain (no tol for the BB of triangles)
  KDTree<MyTriangle, 3> *tree; //!< a KDTree that organizes the triangles in scope (does not store its own copy)

  //! Mesh info
  SpaceVariable3D& coordinates;

  std::vector<GhostPoint> &ghost_nodes_inner; //!< ghost nodes inside the physical domain (shared with other subd)
  std::vector<GhostPoint> &ghost_nodes_outer; //!< ghost nodes outside the physical domain

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain
  int ii0_in, jj0_in, kk0_in, iimax_in, jjmax_in, kkmax_in; //!< corners of the ghosted subdomain, excluding external ghosts
  int NX, NY, NZ; //!< global size (number of cells in the real domain)

  std::vector<double> &x_glob, &y_glob, &z_glob; //!< x,y,z coords of the entire domain
  std::vector<double> &dx_glob, &dy_glob, &dz_glob; //!< dx, dy, dz of the entire domain

  SpaceVariable3D BBmin, BBmax; /**< The min and max coords of nodal bounding boxes. Only for nodes \n
                                     in the physical domain. For each node, the BB contains the \n
                                     node itself and all the edges (n layers, n TBD) that connect the node \n
                                     with other nodes IN THE PHYSICAL DOMAIN.*/
  Vec3D subD_bbmin, subD_bbmax; //!< bounding box of the subdomain (n layers, n TBD)

  SpaceVariable3D TMP, TMP2; //!< For temporary use.

  /************************
   * Results
   ************************/
  //! "CandidatesIndex" and "candidates" account for internal ghost nodes, but not ghost nodes outside physical domain.
  SpaceVariable3D CandidatesIndex; //!< index in the vector "candidates" (-1 means no candidates)
  std::vector<std::pair<Int3, std::vector<MyTriangle> > > candidates;

  //! XForward/XBackward stores edge-surface intersections where at least one vertex of the edge is inside the subdomain.
  SpaceVariable3D XForward; /**< Edge-surface intersections. X[k][j][i][0]: left-edge, [1]: bottom-edge, [2]: back-edge \n
                                 considers all the edges for which BOTH vertices are within the physical domain \n
                                 -1: no intersection. \n
                                 >=0: intersection. The value is its index in intersections (below) \n */
  SpaceVariable3D XBackward; /**< An edge may intersect multiple triangles. We store two of them, those closest to the \n
                                  two vertices. XForward stores the one that is closest to the left/bottom/back vertex. \n
                                  XBackward stores the one that is closest to the right/top/front vertex. */

  //! Phi and Sign communicate w/ neighbor subdomains. So their values are valid also at internal ghost nodes.
  SpaceVariable3D Phi; //!< unsigned distance from each node to the surface (not thickened). Independent from "occluded".
  SpaceVariable3D Sign; //!< GENERALIZED sign: -N (inside enclosure #N), 0 (occluded), or N (inlet, outlet). N = 1,2,...
                      
  //! Closest point on triangle (for near-field nodes inside subdomain, including internal ghosts nodes)
  SpaceVariable3D ClosestPointIndex; //!< index in the vector closest_points. (-1 means not available)
  std::vector<std::pair<Int3, ClosestPoint> > closest_points;

  //! "intersections" stores edge-surface intersections where at least one vertex of the edge is inside the subdomain
  std::vector<IntersectionPoint> intersections; /**< NOTE: Not all these intersections are registered in XForward \n
                                                     and XBackward. When there are occluded nodes, "intersections" \n
                                                     may contain points that are actually not used/registered! */ 

  //! "occluded" and "firstLayer" account for the internal ghost nodes.
  std::set<Int3> occluded;
  std::set<Int3> firstLayer; //!< nodes that belong to intersecting edges (naturally, including occluded nodes)
  std::set<Int3> imposed_occluded; /**< tracks nodes whose sign cannot be resolved; these nodes are FORCED to have\n
                                        the sign of occluded(0), but intersections from these nodes to neighbors\n
                                        may not exist! Includes internal ghost nodes.*/
                                        

  //! tracks nodes that are swept by the surface during small motion (e.g., in one time step)
  //! Does not include nodes that are occluded at present. Includes internal ghost nodes.
  std::set<Int3> swept;

public:

  Intersector(MPI_Comm &comm_, DataManagers3D &dms_, EmbeddedSurfaceData &iod_surface_,
              TriangulatedSurface &surface_,
              SpaceVariable3D &coordinates_, 
              std::vector<GhostPoint> &ghost_nodes_inner_, std::vector<GhostPoint> &ghost_nodes_outer_,
              std::vector<dobule> &x_, std::vector<double> &y_, std::vector<double> &z_,
              std::vector<double> &dx_, std::vector<double> &dy_, std::vector<double> &dz_);

  ~Intersector();

  void Destroy();

  //! Interface tracking functions
  void TrackSurfaceFullCourse(bool &hasInlet, bool &hasOutlet, bool &hasOcc, int &nRegions);

 

/** Below is like the a la carte menu. Try to use the pre-defined "combos" above as much as you can. 
 *  The functions below are not all independent with each other!*/
public:

  void BuildNodalAndSubdomainBoundingBoxes(int nLayer=1); //!< build BBmin, BBmax, subD_bbmin, subD_bbmax

  void BuildSubdomainScopeAndKDTree(); //!< build "scope" and "tree". Requires subdomain bounding box

  //! Many functions below assume that "BuildNodalBoundingBoxes" and "BuildSubdomainScopeAndKDTree" have been called.

  void FindNodalCandidates(); //!< find nearby triangles for each node based on bounding boxes and KDTree

  void FindIntersections(bool with_nodal_cands = false); //!< find occluded nodes, intersections, and first layer nodes

  int FloodFill(bool &hasInlet, bool &hasOutlet, bool &hasOcc, int &nRegions); /**< determine the generalized sign function ("Sign"). 
                                                                                     Returns the number of "colors" (sum of the four).\n*/
  //! Fill "swept". The inputs are firstLayer nodes and surface nodal coords in the previous time step
  void FindSweptNodes(std::vector<Vec3D> &X0, bool nodal_cands_calculated = false); //!< candidates only need to account for 1 layer

  /** When the structure has moved SLIGHTLY, this "refill" function should be called, not the one above. This function only recomputes
   *  the "Sign" of swept nodes. It is faster, and also maintains the same "signs". Calling the original "FloodFill" function may 
   *  lead to sign(tag) change for the same enclosure.
   *  Note: This function must be called AFTER calling "findSweptNodes"*/
  int RefillAfterSurfaceUpdate(bool &hasInlet, bool &hasOutlet, bool &hasOcc, int &nRegions, bool nodal_cands_calculated = false);

  void CalculateUnsignedDistanceNearSurface(int nLayer, bool nodal_cands_calculated = false); //!< Calculate "Phi" for small "nLayers"

private: 
  //! Utility functions
  //
  //! Use a tree to find candidates. maxCand may change, tmp may be reallocated (if size is insufficient)
  inline int FindCandidatesInBox(KDTree<MyTriangle, 3>* mytree, Vec3D bbmin, Vec3D bbmax, 
                                 MyTriangles* tmp, int& maxCand) {
    int found = mytree->findCandidatesInBox(bbmin, bbmax, tmp, maxCand);
    if(found>maxCand) {
      maxCand = found; delete [] tmp;  tmp = new MyTriangle[maxCand];
      found = mytree->findCandidatesInBox(bbmin, bbmax, tmp, maxCand);
    }
    return found; 
  }

  //! Check if a point is occluded by a set of triangles (thickened)
  bool IsPointOccludedByTriangles(Vec3D &coords, MyTriangle* tri, int nTri, double my_half_thickness
                                  int& tid, double* xi = NULL);

  //! Find the intersections of an edge with a set of triangles. Returns the number of intersections.
  int FindEdgeIntersectionsWithTriangles(Vec3D &x0, int i, int j, int k, int dir/*0~x,1~y,2~z*/, 
                                         double len, MyTriangles* tri, int nTri, 
                                         IntersectionPoint &xf, IntersectionPoint &xb); //!< 2 points, maybe the same

};





#endif