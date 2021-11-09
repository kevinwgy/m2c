#ifndef _MULTIPHASE_OPERATOR_H_
#define _MULTIPHASE_OPERATOR_H_
#include<SpaceVariable.h>
#include<PhaseTransition.h>
#include<tuple>

class Vec5D;
class SpaceOperator;
class LevelSetOperator;
class RiemannSolutions;

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

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain

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
                     vector<VarFcnBase*> &varFcn_, SpaceOperator &spo, vector<LevelSetOperator*> &lso);
  ~MultiPhaseOperator();

  //! update material id including the ghost region
  void UpdateMaterialID(vector<SpaceVariable3D*> &Phi, SpaceVariable3D &ID);

  void UpdateStateVariablesAfterInterfaceMotion(SpaceVariable3D &IDn, SpaceVariable3D &ID,
                                                SpaceVariable3D &V, RiemannSolutions &riemann_solutions);

  //! detect phase transitions and update Phi, ID, and V
  int UpdatePhaseTransitions(vector<SpaceVariable3D*> &Phi, SpaceVariable3D &ID,
                             SpaceVariable3D &V, vector<int> &phi_updated, vector<Int3> *new_useful_nodes);

  void Destroy();

protected:
  void UpdateStateVariablesByRiemannSolutions(SpaceVariable3D &IDn, SpaceVariable3D &ID, 
                                              SpaceVariable3D &V, RiemannSolutions &riemann_solutions);

  void UpdateStateVariablesByExtrapolation(SpaceVariable3D &IDn, SpaceVariable3D &ID, SpaceVariable3D &V);

  int LocalUpdateByRiemannSolutions(int i, int j, int k, int id, Vec5D &vl, Vec5D &vr, Vec5D &vb, Vec5D &vt,
                                    Vec5D &vk, Vec5D &vf, RiemannSolutions &riemann_solutions, Vec5D &v,
                                    bool upwind = true);

  void FixUnresolvedNodes(vector<Int3> &unresolved, SpaceVariable3D &IDn, SpaceVariable3D &ID,
                          SpaceVariable3D &V); 

  //! internal function called by UpdatePhaseTransitions
  void UpdatePhiAfterPhaseTransitions(vector<SpaceVariable3D*> &Phi, SpaceVariable3D &ID,
                                      vector<std::tuple<Int3,int,int> > &changed, vector<int> &phi_updated,
                                      vector<Int3> *new_useful_nodes);

};

#endif
