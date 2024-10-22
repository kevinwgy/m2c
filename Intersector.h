/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _INTERSECTOR_H_
#define _INTERSECTOR_H_

#include<IoData.h>
#include<KDTree.h>
#include<TriangulatedSurface.h>
#include<FloodFill.h>
#include<EmbeddedBoundaryDataSet.h>
#include<GlobalMeshInfo.h>
#include<memory> //unique_ptr

/****************************************************************
 * Class Intersector is responsible for tracking a triangulated
 * surface within a fixed Cartesian mesh. A collision-based
 * algorithm is used, similar to the one presented in
 * Wang et al., IJNMF, 2012. The intersector is able to find
 * (1) nodes covered by the interface (occluded nodes),
 * (2) all the edge-interface intersections,
 * (3) shortest distance from each node to the interface (usually
 *     only a few layers of nodes near the surface)
 * (4) closest point on the interface to each node (If solution is not
 *     unique, only one solution is found)
 * (5) nodes swept by the (dynamic) surface in one time step
 * (6) the "color" of each node (i.e. connectivity info)
 * (7) the elements in the embedded surface that form the boundary
 *     of a "color", and their inward-facing side.
 * Note: The above results are stored within this class.
 ***************************************************************/

class Intersector {


  //! Utility class to find and store bounding boxes for triangles
  class MyTriangle {
    int id;
    double x[3], w[3];
  public:
    MyTriangle() {}
    MyTriangle(int id_, double x_[3], double w_[3]) : id(id_) {
      for(int j=0; j<3; j++) {x[j] = x_[j];  w[j] = w_[j];} 
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


  MPI_Comm& comm;

  EmbeddedSurfaceData &iod_surface;

  //! A tool for flood-fill
  FloodFill floodfiller;

  //! Info about the triangulated surface
  TriangulatedSurface &surface; //!< the surface tracked by the intersector
  bool closed_surface; //!< whether the surface is closed AND normals are consistent
  double half_thickness; //!< half thickness of the surface

  //! Mesh info
  SpaceVariable3D& coordinates;

  std::vector<GhostPoint> &ghost_nodes_inner; //!< ghost nodes inside the physical domain (shared with other subd)
  std::vector<GhostPoint> &ghost_nodes_outer; //!< ghost nodes outside the physical domain

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain
  int ii0_in, jj0_in, kk0_in, iimax_in, jjmax_in, kkmax_in; //!< corners of the ghosted subdomain, excluding external ghosts
  int NX, NY, NZ; //!< global size (number of cells in the real domain)

  GlobalMeshInfo &global_mesh;


  //! Infrastructure #1. One layer of neighbors
  SpaceVariable3D BBmin_1, BBmax_1; /**< The min and max coords of nodal bounding boxes. Only for nodes \n
                                     in the physical domain. For each node, the BB contains the \n
                                     node itself and all the edges (n layers, n TBD) that connect the node \n
                                     with other nodes IN THE PHYSICAL DOMAIN.*/
  Vec3D subD_bbmin_1, subD_bbmax_1; //!< bounding box of the subdomain (n layers, n TBD)

  std::vector<MyTriangle> scope_1; //!< triangles relevant to the current subdomain (no tol for the BB of triangles)
  KDTree<MyTriangle, 3> *tree_1; //!< a KDTree that organizes the triangles in scope (does not store its own copy)


  //! Infrastructure #2. N(>1) layer of neighbors
  SpaceVariable3D BBmin_n, BBmax_n; 
  Vec3D subD_bbmin_n, subD_bbmax_n;
  std::vector<MyTriangle> scope_n;
  KDTree<MyTriangle, 3> *tree_n;
  int nLayer; //!< number of layers of neighbors included in the b.b. In most cases, should = Phi_nLayer


  SpaceVariable3D TMP, TMP2; //!< For temporary use.

  /************************
   * Results
   ************************/
  //! "CandidatesIndex" and "candidates" account for internal ghost nodes, but not ghost nodes outside physical domain.
  SpaceVariable3D CandidatesIndex_1; //!< index in the vector "candidates" (-1 means no candidates)
  std::vector<std::pair<Int3, std::vector<MyTriangle> > > candidates_1;

  SpaceVariable3D CandidatesIndex_n;
  std::vector<std::pair<Int3, std::vector<MyTriangle> > > candidates_n;



  //! XForward/XBackward stores edge-surface intersections where both vertices of the edge are inside subdomain or inn. ghost
  SpaceVariable3D XForward; /**< Edge-surface intersections. X[k][j][i][0]: left-edge, [1]: bottom-edge, [2]: back-edge \n
                                 considers all the edges for which BOTH vertices are within the physical domain \n
                                 -1: no intersection. \n
                                 >=0: intersection. The value is its index in intersections (below) \n */
  SpaceVariable3D XBackward; /**< An edge may intersect multiple triangles. We store two of them, those closest to the \n
                                  two vertices. XForward stores the one that is closest to the left/bottom/back vertex. \n
                                  XBackward stores the one that is closest to the right/top/front vertex. */

  //! Phi and Color communicate w/ neighbor subdomains. So their values are valid also at internal ghost nodes.
  SpaceVariable3D Phi; //!< unsigned distance from each node to the surface (not thickened). Independent from "occluded".
  int Phi_nLayer; //!< number of layers of nodes where Phi is calculated.
  SpaceVariable3D Color; //!< GENERALIZED color: -N (inside enclosure #N), 0 (occluded), or N (inlet=1, outlet=2).
  std::vector<int> ColorReachesBoundary; /**< stores whether each negative colored regions (i.e. enclosrures)  touches \n
                                             the domain boundary. Calculated in FloodFillColors and Refill. Note that \n
                                             the size of this vector is "nRegions" calculated in FloodFillColors.*/
  bool hasInlet, hasInlet2, hasOutlet;
  int nRegions;
                      
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
  std::set<Int3> imposed_occluded; /**< tracks nodes whose color cannot be resolved; these nodes are FORCED to have\n
                                        the color of occluded(0), but intersections from these nodes to neighbors\n
                                        may not exist! Includes internal ghost nodes.*/
                                        

  //! tracks nodes that are swept by the surface during small motion (e.g., in one time step)
  //! Does not include nodes that are occluded at present. Includes internal ghost nodes.
  std::set<Int3> swept;

  //! an internally used vector
  std::set<Int3> previously_occluded_but_not_now;

public:

  Intersector(MPI_Comm &comm_, DataManagers3D &dms_, EmbeddedSurfaceData &iod_surface_,
              TriangulatedSurface &surface_,
              SpaceVariable3D &coordinates_, 
              std::vector<GhostPoint> &ghost_nodes_inner_, std::vector<GhostPoint> &ghost_nodes_outer_,
              GlobalMeshInfo &global_mesh_);

  ~Intersector();

  void Destroy();

  //! Get surface half thickness
  inline double GetSurfaceHalfThickness() {return half_thickness;}

  //! Check if a line segment intersects with any triangles inside scope
  bool Intersects(Vec3D &X0, Vec3D &X1);

  //! Interface tracking functions
  double TrackSurfaceFullCourse(bool &hasInlet_, bool &hasInlet2_, bool &hasOutlet_, bool &hasOcc_,
                                int &nRegions_, int phi_layers);

  double RecomputeFullCourse(std::vector<Vec3D> &X0, int phi_layers); 


/** Below is like the a la carte menu. Try to use the pre-defined "combos" above as much as you can. 
 *  The functions below are not all independent with each other!*/
public:

  void BuildNodalAndSubdomainBoundingBoxes(int nL, SpaceVariable3D &BBmin, SpaceVariable3D &BBmax,
                                           Vec3D &subD_bbmin, Vec3D &subD_bbmax); //!< build bounding boxes

  void BuildSubdomainScopeAndKDTree(const Vec3D &subD_bbmin, const Vec3D &subD_bbmax,
                                    std::vector<MyTriangle> &scope, KDTree<MyTriangle, 3> **tree); //!< Requires bounding box

  //! Many functions below assume that bounding boxes, scope, and tree have already been constructed.

  //! find nearby triangles for each node based on bounding boxes and KDTree
  void FindNodalCandidates(SpaceVariable3D &BBmin, SpaceVariable3D &BBmax, KDTree<MyTriangle, 3> *tree,
                           SpaceVariable3D &CandidatesIndex,
                           std::vector<std::pair<Int3, std::vector<MyTriangle> > > &candidates); 

  void FindIntersections(); //!< find occluded nodes, intersections, and first layer nodes

  bool FloodFillColors(); /**< determine the generalized color function ("Color").\n 
                               Returns whether some nodes are occluded.\n"*/
  //! Fill "swept". The inputs are firstLayer nodes and surface nodal coords in the previous time step
  void FindSweptNodes(std::vector<Vec3D> &X0); //!< candidates only need to account for 1 layer

  /** When the structure has moved SLIGHTLY, this "refill" function should be called, not the one above. This function only recomputes
   *  the "Color" of swept nodes. It is faster, and also maintains the same "colors". Calling the original "FloodFill" function may 
   *  lead to color(tag) change for the same enclosure.
   *  Note: This function must be called AFTER calling "findSweptNodes"*/
  void RefillAfterSurfaceUpdate();

  double CalculateUnsignedDistanceNearSurface(int nL); //!< Calculate "Phi" for small "nL"

  //! Find the elements of the embedded surface that constitute the boundary of a "color". For each element in\n
  //! in this set, determine which side(s) of it faces the interior of this color. "status" has the size of\n
  //! surface.elems. For each element, status = 0 means this element does not belong to the boundary of "color",\n
  //! 1 means the positive side (i.e. along "elemNorm") faces the interior, 2 means the negative side faces the\n
  //! interior; 3 means both sides face the interior.
  //! Note: This function assumes intersections and colors "i.e. Color" have been calculated
  void FindColorBoundary(int color, std::vector<int> &status);

  //! Get pointers to all the results
  std::unique_ptr<EmbeddedBoundaryDataSet> GetPointerToResults();

  //! Get "scope_1" (copied to elems_in_scope)
  void GetElementsInScope1(std::vector<int> &elems_in_scope);

  //! Get colors
  void GetColors(bool *hasInlet_ = NULL, bool *hasInlet2_ = NULL, bool *hasOutlet_ = NULL, int *nRegions_ = NULL) {
    if(hasInlet_)   *hasInlet_ = hasInlet;
    if(hasInlet2_) *hasInlet2_ = hasInlet2;
    if(hasOutlet_) *hasOutlet_ = hasOutlet; 
    if(nRegions_)   *nRegions_ = nRegions;
  }

private: 

  //! Utility functions
  //
  //! Use a tree to find candidates. maxCand may change, tmp may be reallocated (if size is insufficient)
  int FindCandidatesInBox(KDTree<MyTriangle, 3>* mytree, Vec3D bbmin, Vec3D bbmax, std::vector<MyTriangle> &tmp, int& maxCand);

  //! Check if a point is occluded by a set of triangles (thickened)
  bool IsPointOccludedByTriangles(Vec3D &coords, MyTriangle* tri, int nTri, double my_half_thickness,
                                  int& tid, double* xi = NULL);

  //! Find the intersections of an edge with a set of triangles. Returns the number of intersections.
  int FindEdgeIntersectionsWithTriangles(Vec3D &x0, int i, int j, int k, int dir/*0~x,1~y,2~z*/, 
                                         double len, MyTriangle* tri, int nTri, 
                                         IntersectionPoint &xf, IntersectionPoint &xb); //!< 2 points, maybe the same


};





#endif
