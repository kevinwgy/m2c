/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _MULTIPHASE_OPERATOR_H_
#define _MULTIPHASE_OPERATOR_H_
#include<SpaceVariable.h>
#include<PhaseTransition.h>
#include<GhostPoint.h>
#include<GlobalMeshInfo.h>
#include<tuple>
#include<memory>

class Vec5D;
class SpaceOperator;
class LevelSetOperator;
class RiemannSolutions;
class EmbeddedBoundaryDataSet;
class Intersector;

/*******************************************
 * class MultiPhaseOperator contains functions
 * for updating material information and state
 * variables at/around material interfaces
 ******************************************/

class MultiPhaseOperator
{
  MPI_Comm&       comm;
  IoData& iod;

  vector<VarFcnBase*> &varFcn;

  //! Mesh info
  SpaceVariable3D& coordinates;
  SpaceVariable3D& delta_xyz;
  SpaceVariable3D& volume;

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain

  vector<GhostPoint> *ghost_nodes_inner; //!< ghost nodes inside the physical domain (shared with other subd)
  vector<GhostPoint> *ghost_nodes_outer; //!< ghost nodes outside the physical domain

  GlobalMeshInfo &global_mesh;

  //! internal variable for tracking or tagging things.
  SpaceVariable3D Tag;

  //! the material id corresponding to each level set function
  map<int,int> ls2matid;

  //! phase transition 
  vector<vector<PhaseTransitionBase*> > trans;  //!< trans[i] contains all possible destinations of phase i

  //! latent heat reservoir (for modeling phase transition)
  SpaceVariable3D Lambda;

public:
  MultiPhaseOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_,
                     vector<VarFcnBase*> &varFcn_, GlobalMeshInfo &global_mesh_,
                     SpaceOperator &spo, vector<LevelSetOperator*> &lso);
  ~MultiPhaseOperator();

  //! get number of materials
  inline int NumberOfMaterials(bool include_inactive = false) {
    return include_inactive ? varFcn.size() : varFcn.size()-1;}

  //! update material id at (external) ghost nodes (they get the IDs of their images)
  void UpdateMaterialIDAtGhostNodes(SpaceVariable3D &ID);

  //! update material id including the ghost region
  void UpdateMaterialIDByLevelSet(vector<SpaceVariable3D*> &Phi0, vector<SpaceVariable3D*> &Phi,
                                  vector<Intersector*> *intersector, SpaceVariable3D &ID);

  //! update V due to interface motion
  int UpdateStateVariablesAfterInterfaceMotion(SpaceVariable3D &IDn, SpaceVariable3D &ID,
                                               SpaceVariable3D &V, RiemannSolutions &riemann_solutions,
                                               vector<Intersector*> *intersector, vector<Int3> &unresolved);

  int FixUnresolvedNodes(vector<Int3> &unresolved, SpaceVariable3D &IDn, SpaceVariable3D &ID,
                         SpaceVariable3D &V, vector<Intersector*> *intersector, vector<Int3> &still_unresolved,
                         bool apply_failsafe_density); 

  //! detect phase transitions and update Phi, ID, and V
  int UpdatePhaseTransitions(vector<SpaceVariable3D*> &Phi, SpaceVariable3D &ID,
                             SpaceVariable3D &V, vector<int> &phi_updated, vector<Int3> *new_useful_nodes);

  //! add stored latent heat (Lambda) to cells that changed phase due to interface motion
  void AddLambdaToEnthalpyAfterInterfaceMotion(SpaceVariable3D &IDn, SpaceVariable3D &ID, SpaceVariable3D &V);

  //! if the boundaries of multiple material subdomains meet, ensure that the phi's are consistent
  int ResolveConflictsInLevelSets(int time_step, vector<SpaceVariable3D*> &Phi);

  //! update phi values in cells that are failed in "UpdateStateVariablesAfterInterfaceMotion". This often
  //  means the appearance of a narrow gap between subdomain boundaries.
  void UpdateLevelSetsInUnresolvedCells(vector<SpaceVariable3D*> &Phi, vector<Int3> &unresolved);

  //! at each node (including external ghosts), at most one "phi" can be negative.
  int CheckLevelSetOverlapping(vector<SpaceVariable3D*> &Phi);



  //! update ID, V, and Phi due to embedded surface motion. Impose boundary conditions for ID, but not V and Phi
  int UpdateCellsSweptByEmbeddedSurfaces(SpaceVariable3D &V, SpaceVariable3D &ID,
                                         vector<SpaceVariable3D*> &Phi,
                                         std::unique_ptr<vector<std::unique_ptr<EmbeddedBoundaryDataSet> > > EBDS,
                                         vector<Intersector*> *intersector);

  //! update Phi and ID to resolve any conflicts with embedded boundaries
  int ResolveConflictsWithEmbeddedSurfaces(vector<SpaceVariable3D*> &Phi,
          SpaceVariable3D &IDn, SpaceVariable3D &ID,
          vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS, vector<Intersector*> *intersector);

  void Destroy();

protected:
  int UpdateStateVariablesByRiemannSolutions(SpaceVariable3D &IDn, SpaceVariable3D &ID, 
                                             SpaceVariable3D &V, RiemannSolutions &riemann_solutions,
                                             vector<Intersector*> *intersector, vector<Int3> &unresolved);

  int UpdateStateVariablesByExtrapolation(SpaceVariable3D &IDn, SpaceVariable3D &ID, SpaceVariable3D &V,
                                          vector<Intersector*> *intersector, vector<Int3> &unresolved);

  int LocalUpdateByRiemannSolutions(int i, int j, int k, int id, Vec5D &vl, Vec5D &vr, Vec5D &vb, Vec5D &vt,
                                    Vec5D &vk, Vec5D &vf, RiemannSolutions &riemann_solutions, Vec5D &v,
                                    bool upwind = true);

  //! internal function called by UpdatePhaseTransitions
  void UpdatePhiAfterPhaseTransitions(vector<SpaceVariable3D*> &Phi, SpaceVariable3D &ID,
                                      vector<std::tuple<Int3,int,int> > &changed, vector<int> &phi_updated,
                                      vector<Int3> *new_useful_nodes);

  //! internal function called by UpdateCellsSweptByEmbeddedSurfaces
  void FindNeighborsForUpdatingSweptNode(int i, int j, int k, double*** tag, double*** id,
                                         vector<Intersector*> *intersector, std::set<Int3> &occlueded,
                                         std::set<Int3> &imposed_occluded,
                                         vector<std::pair<Int3,bool> > &neighbors);

  //! internal function called by ResolveConflictsWithEmbeddedSurfaces
  bool IsOrphanAcrossEmbeddedSurfaces(int i, int j, int k, double*** idn, double*** id,
                                      vector<Intersector*> *intersector);
};

#endif
