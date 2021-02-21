#ifndef _MULTIPHASE_OPERATOR_H_
#define _MULTIPHASE_OPERATOR_H_
#include<IoData.h>
#include<SpaceVariable.h>

/*******************************************
 * class MultiPhaseOperator contains functions
 * for updating material information and state
 * variables at/around material interfaces
 ******************************************/
class SpaceOperator;
class LevelSetOperator;

class MultiPhaseOperator
{
  MPI_Comm&       comm;
  IoData &iod;

  //! Mesh info
  SpaceVariable3D& coordinates;
  SpaceVariable3D& delta_xyz;

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain

  //! internal variable for tracking or tagging things.
  SpaceVariable3D Tag;

  //! the material id corresponding to each level set function
  map<int,int> ls2matid;

public:
  MultiPhaseOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_,
                     SpaceOperator &spo, vector<LevelSetOperator*> &lso);
  ~MultiPhaseOperator();

  //update material id including the ghost region
  void UpdateMaterialID(vector<SpaceVariable3D*> &Phi, SpaceVariable3D &ID);

  void UpdateStateVariablesAfterInterfaceMotion(SpaceVariable3D &IDn, SpaceVariable3D &ID,
                                                SpaceVariable3D &V);

  void Destroy();
};

#endif
