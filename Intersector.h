#ifndef _INTERSECTOR_H_
#define _INTERSECTOR_H_

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
struct IntersectionPoint {
  Int3 n0; //!< the first (i.e., left, bottom, or back) node
  int dir; //!< the direction of the edge (0~x, 1~y, 2~z)
  int tid; //!< id of the triangle that it intersects
  double xi[3]; //!< barycentric coords of the intersection point within the triangle it intersects
};

class Intersector {

  MPI_Comm& comm;

  EmbeddedSurfaceData &iod_surface;

  TriangulatedSurface &surface; //!< the surface tracked by the intersector

  double half_thickness; //!< half thickness of the surface

  //! Mesh info
  SpaceVariable3D& coordinates;
  SpaceVariable3D& delta_xyz;
  SpaceVariable3D& volume;

  std::vector<GhostPoint> &ghost_nodes_inner; //!< ghost nodes inside the physical domain (shared with other subd)
  std::vector<GhostPoint> &ghost_nodes_outer; //!< ghost nodes outside the physical domain

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain

  SpaceVariable3D BBmin, BBMax; //!< the min and max coords of nodal bounding boxes

  /************************
   * Results
   ************************/
  SpaceVariable3D XL, XB, XK; //!< edge-surface intersections. XL: left-edge, XB: bottom-edge, XK: back-edge
                              //!< stores all the edges that involve nodes within the subdomain
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
              SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_, SpaceVariable3D &volume_,
              std::vector<GhostPoint> &ghost_nodes_inner_, std::vector<GhostPoint> &ghost_nodes_outer_);

  ~Intersector();

  void Destroy();



};





#endif
