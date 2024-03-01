/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _SYMMETRY_OPERATOR_H_
#define _SYMMETRY_OPERATOR_H_
#include <IoData.h>
#include <SpaceVariable.h>
#include <VarFcnBase.h>

/*****************************************************************************
 * Class SymmetryOperator handles the sink terms produced by the advective
 * fluxes in the Navier-Stokes equations associated with spherical or 
 * cylindrical symmetry.
 * Note: Other terms caused by viscosity, heat diffusion etc. are implemented
 *       elsewhere!
 ****************************************************************************/

class SymmetryOperator
{
  MPI_Comm& comm;

  MeshData& iod_mesh;

  vector<VarFcnBase*>& varFcn;

  //! Mesh info.
  SpaceVariable3D &coordinates;
  SpaceVariable3D &delta_xyz;
  SpaceVariable3D &volume;

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain

public:

  SymmetryOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, MeshData &iod_mesh_, 
                   vector<VarFcnBase*> &varFcn_, SpaceVariable3D &coordinates_, 
                   SpaceVariable3D &delta_xyz_, SpaceVariable3D &volume_);

  ~SymmetryOperator();

  void AddSymmetryTerms(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &R);

  void Destroy();

private:

  void AddSphericalSymmetryTerms(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &R);

  void AddCylindricalSymmetryTerms(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &R);

};

#endif
