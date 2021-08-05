#ifndef _LEVELSETOPERATOR_H_
#define _LEVELSETOPERATOR_H_

#include <IoData.h>
#include <Reconstructor.h>
#include <LevelSetReinitializer.h>
/*******************************************
 * class LevelSetOperator drives the solution
 * of the level set equation
 ******************************************/
class SpaceOperator;

class LevelSetOperator
{
  MPI_Comm&       comm;
  IoData &iod;
  LevelSetSchemeData& iod_ls; //!< iod may have multiple levelsets, this is the relevant one.

  //! material in the phi<0 region (inside)
  int materialid;
 
  //! Mesh info
  SpaceVariable3D& coordinates;
  SpaceVariable3D& delta_xyz;
  SpaceVariable3D& volume; //!< volume of node-centered control volumes

  vector<GhostPoint>& ghost_nodes_inner; //!< ghost nodes inside the physical domain (shared with other subd)
  vector<GhostPoint>& ghost_nodes_outer; //!< ghost nodes outside the physical domain

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain

  //! Class for spatial reconstruction
  Reconstructor rec;

  //! variables related to the narrow-band level set method
  bool narrow_band; //whether the narrow-band lsm is used
  SpaceVariable3D  Level;  //band level of each node (-inf to inf, including ghost boundary)
  SpaceVariable3D  Useful; //all the nodes in the narrow-band (including ghost boundary)
  SpaceVariable3D  Active; //the nodes in the interior of the narrow-band (including ghost boundary)
  vector<Int3> useful_nodes; //vector of useful nodes
  vector<Int3> active_nodes; //vector of active nodes

  //! Class for reinitialization
  LevelSetReinitializer *reinit;

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

  bool Reinitialize(double time, double dt, int time_step,
                    SpaceVariable3D &Phi); //true: reinitialization is done; false: not this time

  int GetMaterialID() {return materialid;}

  void Destroy();

  //! for debugging/testing the level set solver (Euler / N-S solver not activated)
  void PrescribeVelocityFieldForTesting(SpaceVariable3D &V);

private:
  // functions for internal use within the class
  void Reconstruct(SpaceVariable3D &V, SpaceVariable3D &Phi);
  void ReconstructInBand(SpaceVariable3D &V, SpaceVariable3D &Phi); //the narrow-band version
  void ComputeAdvectionFlux(SpaceVariable3D &R);
  void ComputeAdvectionFluxInBand(SpaceVariable3D &R); //the narrow-band version
  double ComputeLocalAdvectionFlux(double phim, double phip, double um, double up);
  void AddSourceTerm(SpaceVariable3D &Phi, SpaceVariable3D &R);
  void AddSourceTermInBand(SpaceVariable3D &Phi, SpaceVariable3D &R); //the narrow-band version

};

#endif
