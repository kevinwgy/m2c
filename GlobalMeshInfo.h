/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _GLOBAL_MESH_INFO_H_
#define _GLOBAL_MESH_INFO_H_

#include<vector>
#include<cassert>
#include<Vector3D.h>
#include<mpi.h>

class DataManagers3D;

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

private:

  /** Neighbors of all processor cores/subdomains, indexed w/ proc id.
      The neighbors of each subdomain follow a certain order. See FindSubdomainInfo(...) in the .cpp file. **/
  std::vector<std::vector<int> > subD_neighbors_27; //!< all 27 neighbors, including self and non-exist ones
  std::vector<std::vector<int> > subD_neighbors_all; //!< all real neighbors, w/o self and non-exist ones
  std::vector<std::vector<int> > subD_neighbors_19; //!< 27 - 8 corners
  std::vector<std::vector<int> > subD_neighbors_face_edge; //!< real neighbors, excluding corners (at most 19)
  std::vector<std::vector<int> > subD_neighbors_7; //!< 19 - 12 edges 
  std::vector<std::vector<int> > subD_neighbors_face; //!< only real neighbors with face-contact (at most 6)

  bool one_dimensional_mesh; //!< set to true if y and z have only one element

  bool two_dimensional_mesh; //!< set to true if (only) z has only one element

public:

  std::vector<double> x_glob, y_glob, z_glob;
  std::vector<double> dx_glob, dy_glob, dz_glob;

  Vec3D xyz_min, xyz_max; //!< boundary of the physical domain (up to cell boundaries)
  int NX, NY, NZ;

  double domain_volume; //!< does not include ghost layer

  std::vector<Vec3D> subD_xyz_min, subD_xyz_max; //!< actual boundaries of subs, up to cell boundaries
  std::vector<Int3> subD_ijk_min, subD_ijk_max; //!< Note: "max" is max index + 1 

public:

  GlobalMeshInfo() {} //needed by M2CTwinMessenger (and maybe others)
  GlobalMeshInfo(std::vector<double> &x_glob_, std::vector<double> &y_glob_, std::vector<double> &z_glob_,
                 std::vector<double> &dx_glob_, std::vector<double> &dy_glob_, std::vector<double> &dz_glob_);

  ~GlobalMeshInfo();

  //! Setup the subD_xxx vectors
  void FindSubdomainInfo(MPI_Comm& comm, DataManagers3D& dms);

  //! Get specific info 
  double GetXmin() {return xyz_min[0];}
  double GetYmin() {return xyz_min[1];}
  double GetZmin() {return xyz_min[2];}
  double GetXmax() {return xyz_max[0];}
  double GetYmax() {return xyz_max[1];}
  double GetZmax() {return xyz_max[2];}

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
  Vec3D GetXYZ(int i, int j, int k);
  Vec3D GetDXYZ(Int3 ijk);
  Vec3D GetDXYZ(int i, int j, int k);

  bool IsMesh1D() {return one_dimensional_mesh;}

  bool IsMesh2D() {return two_dimensional_mesh;}

  //! If mesh is 2D, only consider dx and dy
  double GetMinDXYZ(Int3 ijk);
  double GetMaxDXYZ(Int3 ijk);

  //! Duplicate of function in SpaceVariable (checking nodes/cells)
  inline bool OutsidePhysicalDomain(int i, int j, int k) {
    return (i<0 || i>=NX || j<0 || j>=NY || k<0 || k>=NZ);}

  //! Duplicate of function in SpaceVariable (checking nodes/cells)
  inline bool OutsidePhysicalDomainAndUnpopulated(int i, int j, int k)
  {
    int count = 0;
    if(i<0 || i>=NX) count++;
    if(j<0 || j>=NY) count++;
    if(k<0 || k>=NZ) count++;
    if(count>1)
      return true;
    return false;
  }

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

  //! Find the subdomain/processor core that owns a node (or equiv. cell)
  int GetOwnerOfCell(int i, int j, int k, 
                     bool include_ghost_layer = false); //!< Not optimally fast. Shouldn't call it too frequently.
  int GetOwnerOfNode(int i, int j, int k,
                     bool include_ghost_layer = false) {return GetOwnerOfCell(i,j,k,include_ghost_layer);}

  //! Find the subdomain/processor core that owns a point
  int GetOwnerOfPoint(Vec3D &p, bool include_ghost_layer = false);

  //! Check if a node/cell (i,j,k) is inside a certain subdomain
  bool IsCellInSubdomain(int i, int j, int k, int sub, bool include_ext_ghost_layer = false);
  bool IsNodeInSubdomain(int i, int j, int k, int sub, bool include_ext_ghost_layer = false) 
           {return IsCellInSubdomain(i,j,k,sub,include_ext_ghost_layer);}

  //! Check if a point (x,y,z) is inside a certain subdomain
  bool IsPointInSubdomain(Vec3D &p, int sub, bool include_ext_ghost_layer = false);

  //! Get neighbors
  std::vector<int> &GetAllNeighborsOfSub(int sub); 
  std::vector<int> &GetFaceEdgeNeighborsOfSub(int sub);
  std::vector<int> &GetFaceNeighborsOfSub(int sub);
  std::vector<int> &Get27NeighborhoodOfSub(int sub);
  std::vector<int> &Get19NeighborhoodOfSub(int sub);
  std::vector<int> &Get7NeighborhoodOfSub(int sub);

};





#endif
