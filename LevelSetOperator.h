#ifndef _LEVELSETOPERATOR_H_
#define _LEVELSETOPERATOR_H_

#include <IoData.h>
#include <Reconstructor.h>
/*******************************************
 * class LevelSetOperator drives the solution
 * of the level set equation
 ******************************************/
class SpaceOperator;

class LevelSetOperator
{
  MPI_Comm&       comm;
  DataManagers3D& dm_all;
  IoData &iod;
  LevelSetSchemeData& iod_ls; //!< iod may have multiple levelsets, this is the relevant one.

  //! material in the phi<0 region (inside)
  int materialid;
 
  //! Mesh info
  SpaceVariable3D& coordinates;
  SpaceVariable3D& delta_xyz;
  SpaceVariable3D& volume; //!< volume of node-centered control volumes

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain

  //! Class for spatial reconstruction
  Reconstructor rec;

  //! Reconstructed velocity (u, v, w);
  SpaceVariable3D scalar, ul, ur, vb, vt, wk, wf;
  SpaceVariable3D dudx, dvdy, dwdz; //!< derivatives of velocity
  //! Reconstructed signed distance function (Phi)
  SpaceVariable3D Phil, Phir, Phib, Phit, Phik, Phif;

public:
  LevelSetOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_,
                   LevelSetSchemeData &iod_ls_, SpaceOperator &spo);
  ~LevelSetOperator();

  void SetInitialCondition(SpaceVariable3D &Phi);
  void ApplyBoundaryConditions(SpaceVariable3D &Phi);

  void ComputeResidual(SpaceVariable3D &V, SpaceVariable3D &Phi, SpaceVariable3D &R);

  void Reinitialize(SpaceVariable3D &Phi);

  void Destroy();

private:
  // functions for internal use within the class
  void Reconstruct(SpaceVariable3D &V, SpaceVariable3D &Phi);
  void ComputeAdvectionFlux(SpaceVariable3D &R);
  double ComputeLocalAdvectionFlux(double phim, double phip, double um, double up);
  void AddSourceTerm(SpaceVariable3D &Phi, SpaceVariable3D &R);
};

#endif
