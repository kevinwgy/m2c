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

/*****************************************************************************
 * Class SpaceInitializer is responsible for initializing V, ID, and Phi
 ****************************************************************************/

class SpaceInitializer
{
  MPI_Comm&       comm;
  IoData&         iod;
  GlobalMeshInfo& global_mesh;

public:

  SpaceInitializer(MPI_Comm &comm_, IoData &iod_, GlobalMeshInfo &global_mesh_);
  ~SpaceInitialier();

  void Destroy();

  std::multimap<int, std::pair<int,int> >
  SetInitialCondition(SpaceVariable3D &V, SpaceVariable3D &ID, std::vector<SpaceVariable3D*> &Phi,
                      SpaceOperator &spo, std::vector<SpaceVariable3D*> &lso,
                      std::unique_ptr<std::vector<std::unique_ptr<EmbeddedBoundaryDataSet> > > EBDS = nullptr);

private:


void SortUserSpecifiedGeometries(IoDataGeo &iod_geo,
                                 std::vector<std::pair<int,int> > &order,
                                 std::map<int,GeoTools::DistanceFromPointToPlane*> &plane_cal,
                                 std::map<int,GeoTools::DistanceFromPointToCylinderCone*> &cylindercone_cal,
                                 std::map<int,GeoTools::DistanceFromPointToCylinderSphere*> &cylindersphere_cal,
                                 std::map<int,GeoTools::DistanceFromPointToSphere*> &sphere_cal,
                                 std::map<int,GeoTools::DistanceFromPointToParallelepiped*> &parallelepiped_cal,
                                 std::map<int,GeoTools::DistanceFromPointToSpheroid*> &spheroid_cal);
};

#endif
