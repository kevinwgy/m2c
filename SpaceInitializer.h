/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _SPACE_INITIALIZER_H_
#define _SPACE_INITIALIZER_H_

#include <IoData.h>
#include <SpaceVariable.h>
#include <GlobalMeshInfo.h>

class SpaceOperator;
class EmbeddedBoundaryDataSet;
class LevelSetOperator;
class GravityHandler;
struct Vec3D;
struct Vec5D;
struct GhostPoint;

/*****************************************************************************
 * Class SpaceInitializer is responsible for initializing V, ID, and Phi
 ****************************************************************************/

class SpaceInitializer
{
  MPI_Comm&       comm;
  DataManager3D&  dms;
  IoData&         iod;
  GlobalMeshInfo& global_mesh;

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain
  int NX, NY, NZ; //!< global size

public:

  SpaceInitializer(MPI_Comm &comm_, DataManager3D &dms_, IoData &iod_, GlobalMeshInfo &global_mesh_,
                   SpaceVariable3D &coordinates);
  ~SpaceInitialier();

  void Destroy();

  std::multimap<int, std::pair<int,int> >
  SetInitialCondition(SpaceVariable3D &V, SpaceVariable3D &ID, std::vector<SpaceVariable3D*> &Phi,
                      SpaceOperator &spo, std::vector<SpaceVariable3D*> &lso, GravityHandler *ghand = NULL,
                      std::unique_ptr<std::vector<std::unique_ptr<EmbeddedBoundaryDataSet> > > EBDS = nullptr);

private:

  int OrderUserSpecifiedGeometries(std::vector<std::pair<int,int> > &order);

  void AddGeomToVector(int o, int type, int ind, string& name, std::vector<std::pair<int,int> > &order,
                       std::vector<int> &user_specified_order);

  std::multimap<int, std::pair<int,int> >
  InitializeVandID(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceOperator &spo,
                   std::unique_ptr<std::vector<std::unique_ptr<EmbeddedBoundaryDataSet> > > EBDS,
                   std::vector<std::pair<int,int> > &order);

  void ApplyUserSpecifiedInitialConditionFile(Vec3D*** coords, Vec5D*** v, double*** id);

  void InitializeVandIDWithinEnclosure(UserSpecifiedEnclosureData &enclosure, 
                                       SpaceVariable3D &coordinates,
                                       std::vector<GhostPoint> &ghost_nodes_inner,
                                       std::vector<GhostPoint> &ghost_nodes_outer,
                                       Vec5D*** v, double*** id);

  std::pair<int, std::pair<int,int> >
  InitializeVandIDByPoint(PointData& point,
                          std::vector<std::unique_ptr<EmbeddedBoundaryDataSet> > &EBDS,
                          std::vector<double***> &color, Vec5D*** v, double*** id);

  void InitializePhi(SpaceVariable3D &ID, SpaceVariable3D &coordinates, SpaceVariable3D &delta_xyz,
                     std::vector<GhostPoint> &spo_ghost_nodes_inner,
                     std::vector<GhostPoint> &spo_ghost_nodes_outer,
                     std::vector<SpaceVariable3D*> &Phi, std::vector<LevelSetOperator*> &lso,
                     std::unique_ptr<std::vector<std::unique_ptr<EmbeddedBoundaryDataSet> > > EBDS,
                     std::vector<std::pair<int,int> > &order,
                     std::multimap<int, std::pair<int,int> > &id2closure);

  void InitializePhiByDistance(SpaceVariable3D &ID, SpaceVariable3D &coordinates,
                               SpaceVariable3D &delta_xyz,
                               vector<GhostPoint> &spo_ghost_nodes_inner,
                               vector<GhostPoint> &spo_ghost_nodes_outer,
                               SpaceSpaceVariable3D &Phi, LevelSetOperator &lso,
                               unique_ptr<vector<unique_ptr<EmbeddedBoundaryDataSet> > > EBDS = nullptr,
                               vector<pair<int,int> > *surf_and_color = NULL);

  void InitializePhiByReinitialization(SpaceVariable3D &ID,
                                       SpaceSpaceVariable3D &Phi, LevelSetOperator &lso,
                                       vector<pair<int,int> > *surf_and_color = NULL);

  void InitializePhiWithinEnclosure(UserSpecifiedEnclosureData &enclosure,
                                    SpaceVariable3D &coordinates,
                                    SpaceVariable3D &delta_xyz,
                                    std::vector<GhostPoint> &spo_ghost_nodes_inner,
                                    std::vector<GhostPoint> &spo_ghost_nodes_outer,
                                    SpaceVariable3D &Phi, LevelSetOperator &lso);
};

#endif
