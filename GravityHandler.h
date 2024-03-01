/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _GRAVITY_HANDLER_H_
#define _GRAVITY_HANDLER_H_

#include <IoData.h>
#include <SpaceVariable.h>
#include <GlobalMeshInfo.h>
#include <GhostPoint.h>

class EmbeddedBoundaryDataSet;
class LevelSetOperator;

/*****************************************************************************
 * Class GravityHandler handles gravity and free-surface related work. It is
 * generally a tools/functions class. Operate on the data provided, but does not 
 * own them.
 ****************************************************************************/

class GravityHandler
{
  MPI_Comm&                 comm;
  DataManagers3D&           dm_all;
  IoData&                   iod;

  SpaceVariable3D& coordinates;

  vector<GhostPoint>& ghost_nodes_inner; //!< ghost nodes inside the physical domain (shared with other subd)
  vector<GhostPoint>& ghost_nodes_outer; //!< ghost nodes outside the physical domain

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain
  int NX, NY, NZ; //!< global size

  GlobalMeshInfo &global_mesh;


public:

  GravityHandler(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_, SpaceVariable3D &coordinates_,
                 std::vector<GhostPoint> &ghost_nodes_inner_, std::vector<GhostPoint> &ghost_nodes_outer_,
                 GlobalMeshInfo &global_mesh_);

  ~GravityHandler();

  void Destroy();

  void UpdateInitialConditionByFlooding(SpaceVariable3D &V, SpaceVariable3D &ID,
                                        std::vector<LevelSetOperator*> lso, std::vector<SpaceVariable3D*>& Phi,
                                        std::vector<std::unique_ptr<EmbeddedBoundaryDataSet> >* EBDS = NULL);

};

#endif
