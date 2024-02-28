/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _SPACE_INITIALIZER_H_
#define _SPACE_INITIALIZER_H_

#include <IoData.h>
#include <SpaceVariable.h>
#include <GlobalMeshInfo.h>
#include <UserDefinedState.h>
#include <memory> //std::unique_ptr
#include <tuple>

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
  DataManagers3D& dms;
  IoData&         iod;
  GlobalMeshInfo& global_mesh;

  std::tuple<UserDefinedState*, void*, DestroyUDS*> state_calculator; //!< the 1st one is the calculator

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain
  int NX, NY, NZ; //!< global size

public:

  SpaceInitializer(MPI_Comm &comm_, DataManagers3D &dms_, IoData &iod_, GlobalMeshInfo &global_mesh_,
                   SpaceVariable3D &coordinates);
  ~SpaceInitializer();

  void Destroy();

  std::multimap<int, std::pair<int,int> >
  SetInitialCondition(SpaceVariable3D &V, SpaceVariable3D &ID, std::vector<SpaceVariable3D*> &Phi,
                      std::vector<SpaceVariable3D*> &NPhi, std::vector<SpaceVariable3D*> &KappaPhi,
                      SpaceOperator &spo, std::vector<LevelSetOperator*> &lso, GravityHandler *ghand = NULL,
                      std::unique_ptr<std::vector<std::unique_ptr<EmbeddedBoundaryDataSet> > > EBDS = nullptr);

private:

  int OrderUserSpecifiedGeometries(std::vector<std::pair<int,int> > &order);

  void AddGeomToVector(int o, int type, int ind, string name, std::vector<std::pair<int,int> > &order,
                       std::vector<int> &user_specified_order);

  std::multimap<int, std::pair<int,int> >
  InitializeVandID(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceOperator &spo,
                   std::vector<std::unique_ptr<EmbeddedBoundaryDataSet> >* EBDS,
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

  void ApplyUserDefinedVandID(Vec3D*** coords, Vec5D*** v, double*** id);

  void InitializePhi(SpaceVariable3D &ID, SpaceOperator &spo,
                     std::vector<SpaceVariable3D*> &Phi, std::vector<SpaceVariable3D*> &NPhi,
                     std::vector<SpaceVariable3D*> &KappaPhi, std::vector<LevelSetOperator*> &lso,
                     std::vector<std::unique_ptr<EmbeddedBoundaryDataSet> >* EBDS,
                     std::vector<std::pair<int,int> > &order,
                     std::multimap<int, std::pair<int,int> > &id2closure);

  void InitializePhiByDistance(vector<std::pair<int,int> > &order, SpaceVariable3D &coordinates,
                               SpaceVariable3D &delta_xyz,
                               vector<GhostPoint> &spo_ghost_nodes_inner,
                               vector<GhostPoint> &spo_ghost_nodes_outer,
                               SpaceVariable3D &Phi, SpaceVariable3D &NPhi, SpaceVariable3D &KappaPhi,
                               LevelSetOperator &lso,
                               vector<std::unique_ptr<EmbeddedBoundaryDataSet> >* EBDS = NULL,
                               vector<std::pair<int,int> > *surf_and_color = NULL);

  void InitializePhiByReinitialization(SpaceVariable3D &ID,
                                       SpaceVariable3D &Phi, SpaceVariable3D &NPhi,
                                       SpaceVariable3D &KappaPhi, LevelSetOperator &lso,
                                       std::vector<std::pair<int,int> > *surf_and_color = NULL);

  bool InitializePhiWithinEnclosure(UserSpecifiedEnclosureData &enclosure,
                                    SpaceVariable3D &coordinates,
                                    SpaceVariable3D &delta_xyz,
                                    std::vector<GhostPoint> &spo_ghost_nodes_inner,
                                    std::vector<GhostPoint> &spo_ghost_nodes_outer,
                                    SpaceVariable3D &Phi, LevelSetOperator &lso);
};

#endif
