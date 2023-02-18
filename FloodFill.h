/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _FLOOD_FILL_H_
#define _FLOOD_FILL_H_

#include<Vector3D.h>
#include<SpaceVariable.h>
#include<GhostPoint.h>

/*******************************************************************
 * Class FloodFill is a utility class that fills continuous regions
 * of the domain with distinct colors (integer flags).
 ******************************************************************/

class FloodFill {

  MPI_Comm& comm;

  SpaceVariable3D TMP; //!< for internal use

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain
  int ii0_in, jj0_in, kk0_in, iimax_in, jjmax_in, kkmax_in; //!< corners of the ghosted subdomain, excluding external ghosts
  int NX, NY, NZ; //!< global size

  std::vector<GhostPoint> ghost_nodes_inner; //!< ghost nodes inside the physical domain (shared with other subd)
  std::vector<GhostPoint> ghost_nodes_outer; //!< ghost nodes outside the physical domain

public:

  FloodFill(MPI_Comm &comm_, DataManagers3D &dms_, std::vector<GhostPoint> &ghost_nodes_inner_, 
            std::vector<GhostPoint> &ghost_nodes_outer_);
  ~FloodFill();

  void Destroy();

  /** Flood-fill nodes in the physical domain (excluding ghost nodes outside the physical domain)
   *  based on edge obstructions. Returns the number of colors for unoccluded nodes. "occluded_nodes" must include
   *  occluded internal ghost nodes. 
   *  Colors are: 0 (occluded), 1, 2, ... */
  int FillBasedOnEdgeObstructions(SpaceVariable3D& Obs, int non_obstruction_flag,
                                  std::set<Int3>& occluded_nodes, SpaceVariable3D& Color);




};



#endif
