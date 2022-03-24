#ifndef _INTERSECTOR_H_
#define _INTERSECTOR_H_

#include<IoData.h>
#include<SpaceVariable.h>
#include<KDTree.h>
#include<TriangulatedSurface.h>
#include<GhostPoint.h>

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
 * (5) (if the interface is a closed surface) the inside/outside
 *     status of each node with respect to the surface.
 * Note: The above results are stored within this class.
 ***************************************************************/

class Intersector {

  //! A nested class that stores information about an intersection point
  struct IntersectionPoint {
    Int3 n0; //!< the first (i.e., left, bottom, or back) node
    int dir; //!< the direction of the edge (0~x, 1~y, 2~z)
    int tid; //!< id of the triangle that it intersects
    double xi[3]; //!< barycentric coords of the intersection point within the triangle it intersects
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


  MPI_Comm& comm;

  EmbeddedSurfaceData &iod_surface;

  //! Info about the triangulated surface
  TriangulatedSurface &surface; //!< the surface tracked by the intersector
  bool closed_surface; //!< whether the surface is closed AND normals are consistent
  double half_thickness; //!< half thickness of the surface

  std::vector<MyTriangle> scope; //!< triangles relevant to the current subdomain (no tol for the BB of triangles)
  KDTree<MyTriangle, 3> *tree; //!< a KDTree that organizes the triangles in scope (does not store its own copy)

  //! Mesh info
  SpaceVariable3D& coordinates;
  SpaceVariable3D& delta_xyz;
  SpaceVariable3D& volume;

  std::vector<GhostPoint> &ghost_nodes_inner; //!< ghost nodes inside the physical domain (shared with other subd)
  std::vector<GhostPoint> &ghost_nodes_outer; //!< ghost nodes outside the physical domain

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain
  int NX, NY, NZ; //!< global size (number of cells in the real domain)

  SpaceVariable3D BBmin, BBmax; /**< the min and max coords of nodal bounding boxes. Only for nodes
                                     in the physical domain. For each node, the BB contains the 
                                     node itself and all the edges that connect the node with other
                                     nodes IN THE PHYSICAL DOMAIN.*/
  Vec3D subD_bbmin, subD_bbmax; //!< bounding box of the subdomain

  /************************
   * Results
   ************************/
  SpaceVariable3D CandidatesIndex; //!< index in the vector "candidates" (-1 means no candidates)
  std::vector<std::pair<Int3, std::vector<MyTriangle> > > candidates;

  SpaceVariable3D XX; //!< edge-surface intersections. XX[k][j][i][0]: left-edge, [1]: bottom-edge, [2]: back-edge
                      //!< considers all the edges for which BOTH vertices are within the physical domain 
                      //!< -1: no intersection.
                      //!< >=0: intersection. The value is its index in intersections (below)
  SpaceVariable3D Phi; //!< unsigned distance to the interface
  SpaceVariable3D Sign; //!< -1 (inside), 0 (occluded), or 1 (outside)
                        //!< (-1 is only relevant for closed surfaces with consistent element-orientation)

  std::vector<IntersectionPoint> intersections;
     
  std::vector<Int3> occluded;
  std::vector<Int3> fisrtLayer; //!< includes occluded nodes

public:

  Intersector(MPI_Comm &comm_, DataManagers3D &dms_, EmbeddedSurfaceData &iod_surface_,
              TriangulatedSurface &surface_,
              SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_, SpaceVariable3D &volume_,
              std::vector<GhostPoint> &ghost_nodes_inner_, std::vector<GhostPoint> &ghost_nodes_outer_);

  ~Intersector();

  void Destroy();


private:

  void BuildNodalBoundingBoxes(); //!< build BBmin, BBmax, subD_bbmin, subD_bbmax

  void BuildSubdomainScopeAndKDTree(); //!< build "scope" and "tree"

  void FindNodalCandidates(); //!< find nearby triangles for each node based on bounding boxes and KDTree

  void FindIntersections(bool with_nodal_cands = false); //!< find occluded nodes, intersections, and first layer nodes

};





#endif
