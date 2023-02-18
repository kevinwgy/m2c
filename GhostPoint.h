/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _GHOST_POINT_H
#define _GHOST_POINT_H
#include <Vector3D.h>
#include <bits/stdc++.h> //INT_MAX
#include <cfloat> //DBL_MAX

/** This class is responsible for storing information outside each subdomain.
 *  Each 'ghost point' can be a ghost node, or just a location outside the subdomain.
 */
struct GhostPoint {

  bool isNode; //!< whether the ghost point is a node

  int owner_proc; //!< (only used when isNode==true) ID of the processor that owns this ghost node

  Int3  ijk; //!< if not a node, this should be the indices of the lower-left corner of the interpolation cell
  Vec3D xi; //!< local coordinates for trilinear interpolation

  Int3  image_ijk;  //!< if not a node, the indices of the lower-left corner of the interpolation cell
  Vec3D image_xi; //!< local coordinates for trilinear interpolation

  Vec3D boundary_projection;
  Vec3D outward_normal;
  enum ProjectionType {NONE = 0, FACE = 1, EDGE = 2, VERTEX = 3, SIZE = 4} type_projection; 
  
  enum Side {UNDEFINED = 0, LEFT = 1, RIGHT = 2, BOTTOM = 3, TOP = 4, BACK = 5, FRONT = 6} side;

  int bcType; //! type of b.c. of the projection point 

  //! Constructors
  GhostPoint() : 
      isNode(false), ijk(INT_MAX), xi(0.0), image_ijk(INT_MAX), image_xi(0.0),
      boundary_projection(DBL_MAX), outward_normal(DBL_MAX), bcType(0),
      type_projection(NONE), side(UNDEFINED), owner_proc(-1) {}

  GhostPoint(Int3 ijk_) : 
      isNode(true), ijk(ijk_), xi(0.0), image_ijk(INT_MAX), image_xi(0.0),
      boundary_projection(DBL_MAX), outward_normal(DBL_MAX), bcType(0),
      type_projection(NONE), side(UNDEFINED), owner_proc(-1)  {}

  GhostPoint(int i, int j, int k) : 
      isNode(true), ijk(Int3(i,j,k)), xi(0.0), image_ijk(INT_MAX), image_xi(0.0),
      boundary_projection(DBL_MAX), outward_normal(DBL_MAX), bcType(0),
      type_projection(NONE), side(UNDEFINED), owner_proc(-1)  {}

  GhostPoint(Int3 ijk_, Vec3D xi_) : 
      isNode(false), ijk(ijk_), xi(xi_), image_ijk(INT_MAX), image_xi(0.0),
      boundary_projection(DBL_MAX), outward_normal(DBL_MAX), bcType(0),
      type_projection(NONE), side(UNDEFINED), owner_proc(-1)  {}

  GhostPoint(Int3 ijk_, Int3 image_ijk_) : 
      isNode(true), ijk(ijk_), xi(0.0), image_ijk(image_ijk_), image_xi(0.0),
      boundary_projection(DBL_MAX), outward_normal(DBL_MAX), bcType(0),
      type_projection(NONE), side(UNDEFINED), owner_proc(-1)  {}

  GhostPoint(Int3 ijk_, Int3 image_ijk_, ProjectionType projType_, Vec3D proj_, Vec3D normal_, int bcType_) : 
      isNode(true), ijk(ijk_), xi(0.0), image_ijk(image_ijk_), image_xi(0.0),
      boundary_projection(proj_), outward_normal(normal_), bcType(bcType_),
      type_projection(projType_), side(UNDEFINED), owner_proc(-1)  {}

  GhostPoint(Int3 ijk_, Int3 image_ijk_, ProjectionType projType_, Vec3D proj_, Vec3D normal_, int bcType_, Side side_) : 
      isNode(true), ijk(ijk_), xi(0.0), image_ijk(image_ijk_), image_xi(0.0),
      boundary_projection(proj_), outward_normal(normal_), bcType(bcType_),
      type_projection(projType_), side(side_), owner_proc(-1)  {}

  GhostPoint(Int3 ijk_, Vec3D xi_, Int3 image_ijk_, Vec3D image_xi_) : 
      isNode(false), ijk(ijk_), xi(xi_), image_ijk(image_ijk_), image_xi(image_xi_),
      boundary_projection(DBL_MAX), outward_normal(DBL_MAX), bcType(0),
      type_projection(NONE), side(UNDEFINED), owner_proc(-1)  {}


  ~GhostPoint() {}

  bool IsDefined() {return ijk[0]<INT_MAX;}
  bool IsImageDefined() {return image_ijk[0]<INT_MAX;}
  bool IsProjectionDefined() {return boundary_projection[0]<1e20 && outward_normal[0]<1e20 
                                     && type_projection!=NONE;}
  bool IsBcTypeDefined() {return bcType>0;}

  void SetProjection(Vec3D proj_, Vec3D normal_, ProjectionType proj_type_, int bcType_)
  {
    boundary_projection = proj_;  
    outward_normal = normal_;  
    type_projection = proj_type_; 
    bcType = bcType_;
  }

};


#endif
