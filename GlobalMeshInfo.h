#ifndef _GLOBAL_MESH_INFO_H_
#define _GLOBAL_MESH_INFO_H_

#include<vector>
#include<cassert>

struct Int3;
struct Vec3D;

/************************************************
 * class GlobalMeshInfo stores information about
 * the global Cartesian mesh. It also includes
 * functions that calculate various quantities
 * using the global mesh info. This is supposed
 * to be a lightweight, generic class. Special
 * calculations should not be done here.
 *
 * Note: This class assumes the ghost nodes are
 * located with symmetry to the first layer of
 * nodes within the domain interior. This is
 * usually true, but not always. See "ResetGhostLayer"
 * in SpaceOperator for a counterexample.
 ***********************************************/

class GlobalMeshInfo {

public:

  std::vector<double> x_glob, y_glob, z_glob;
  std::vector<double> dx_glob, dy_glob, dz_glob;

public:

  GlobalMeshInfo() {}
  GlobalMeshInfo(std::vector<double> &x_glob_, std::vector<double> &y_glob_, std::vector<double> &z_glob_,
                 std::vector<double> &dx_glob_, std::vector<double> &dy_glob_, std::vector<double> &dz_glob_) {
    x_glob = x_glob_;   y_glob = y_glob_;   z_glob = z_glob_;
    dx_glob = dx_glob_; dy_glob = dy_glob_; dz_glob = dz_glob_;
    assert(x_glob.size() == dx_glob.size());
    assert(y_glob.size() == dy_glob.size());
    assert(z_glob.size() == dz_glob.size());}

  ~GlobalMeshInfo() {}

  //! Get specific info 
  double GetX(int i);
  double GetY(int j);
  double GetZ(int k);
  double GetDx(int i);
  double GetDy(int j);
  double GetDz(int k);

  double GetX(Int3 ijk);
  double GetY(Int3 ijk);
  double GetZ(Int3 ijk);
  double GetDx(Int3 ijk);
  double GetDy(Int3 ijk);
  double GetDz(Int3 ijk);

  Vec3D GetXYZ(Int3 ijk);

  //! Determine is a point is inside the domain (formed by control volumes / cells)
  bool IsPointInDomain(Vec3D &p, bool include_ghost_layer = false);

  //! Determine is a point is inside the "primal" domain determined by nodes
  bool IsPointInNodalMesh(Vec3D &p, bool include_ghost_layer = false);

  //! Find the closest node to an arbitrary point in the 3D space
  Int3 FindClosestNodeToPoint(Vec3D &p, bool include_ghost_layer = false);

  //! Find the control volume / cell that covers a point in 3D (Not the same as the previous function!) 
  bool FindCellCoveringPoint(Vec3D &p, Int3 &ijk, bool include_ghost_layer = false);

  //! Find the "element" (in the primal mesh) that covers a point in 3D.
  bool FindElementCoveringPoint(Vec3D &p, Int3 &ijk0,
                                Vec3D *xi = NULL, //optional output: local coords of "point" within element
                                bool include_ghost_layer = false);


};





#endif
