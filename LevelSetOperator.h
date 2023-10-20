/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _LEVELSETOPERATOR_H_
#define _LEVELSETOPERATOR_H_

#include <IoData.h>
#include <Reconstructor.h>
#include <LevelSetReinitializer.h>
#include <GlobalMeshInfo.h>
#include <memory> //unique_ptr
/*******************************************
 * class LevelSetOperator drives the solution
 * of the level set equation
 ******************************************/
class SpaceOperator;
class EmbeddedBoundaryDataSet;

class LevelSetOperator
{
  MPI_Comm        &comm;
  DataManagers3D  &dms;
  IoData          &iod;

  LevelSetSchemeData& iod_ls; //!< iod may have multiple levelsets, this is the relevant one.


  //! material in the phi<0 region (inside)
  int materialid;
 
  //! Mesh info
  SpaceVariable3D& coordinates;
  SpaceVariable3D& delta_xyz;
  SpaceVariable3D& volume; //!< volume of node-centered control volumes

  //! The ghost node vectors are almost the same as those in SpaceOperator, except bcType
  vector<GhostPoint> ghost_nodes_inner; //!< ghost nodes inside the physical domain (shared with other subd)
  vector<GhostPoint> ghost_nodes_outer; //!< ghost nodes outside the physical domain

  //! These are the same ghost_nodes in SpaceOperator
  vector<GhostPoint>& spo_ghost_nodes_inner; //!< ghost nodes inside the physical domain (shared with other subd)
  vector<GhostPoint>& spo_ghost_nodes_outer; //!< ghost nodes outside the physical domain

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain

  GlobalMeshInfo &global_mesh;

  //! Class for spatial reconstruction (FVM)
  Reconstructor *rec;

  //! Class for calculating spatial gradient (FDM)
  GradientCalculatorBase *grad_minus; //left-biased FD
  GradientCalculatorBase *grad_plus;  //right-biased FD
  
  //! variables related to the narrow-band level set method
  bool narrow_band; //whether the narrow-band lsm is used
  SpaceVariable3D  Level;  //band level of each node (-inf to inf, including ghost boundary)
  SpaceVariable3D  UsefulG2; //all the nodes in the narrow-band (including ghost boundary)
  SpaceVariable3D  Active; //the nodes in the interior of the narrow-band (including ghost boundary)
  vector<Int3> useful_nodes; //vector of useful nodes (i.e. within bandwidth)
  vector<Int3> active_nodes; //vector of active nodes (i.e. within bandwidth - 1)

  //! Class for reinitialization
  LevelSetReinitializer *reinit;

  //! Reconstructed velocity (u, v, w);
  SpaceVariable3D scalar, scalarG2, ul, ur, vb, vt, wk, wf;
  SpaceVariable3D dudx, dvdy, dwdz; //!< derivatives of velocity
  //! Reconstructed signed distance function (Phi)
  SpaceVariable3D Phil, Phir, Phib, Phit, Phik, Phif;

public:
  LevelSetOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_,
                   LevelSetSchemeData &iod_ls_, SpaceOperator &spo);
  ~LevelSetOperator();

  void ApplyBoundaryConditions(SpaceVariable3D &Phi);

  void ApplyBoundaryConditionsNPhi(SpaceVariable3D &NPhi);

  void ApplyBoundaryConditionsKappaPhi(SpaceVariable3D &KappaPhi);

  void ComputeResidual(SpaceVariable3D &V, SpaceVariable3D &Phi, SpaceVariable3D &R, double time);

  void ConstructNarrowBandInReinitializer(SpaceVariable3D &Phi);

  bool Reinitialize(double time, double dt, int time_step,
                    SpaceVariable3D &Phi, int special_maxIts = 0,//!< if >0, will use it instead of iod value
                    bool must_do = false); //!< true: reinitialization is done; false: not this time

  void ReinitializeAfterPhaseTransition(SpaceVariable3D &Phi, vector<Int3> &new_nodes);

  int GetMaterialID() {return materialid;}

  bool NarrowBand() {return narrow_band;}

  bool HasReinitializer() {return reinit!=NULL;}

  LevelSetSchemeData& GetLevelSetSchemeData() {return iod_ls;}
  
  vector<GhostPoint>* GetPointerToInnerGhostNodes() {return &ghost_nodes_inner;}
  vector<GhostPoint>* GetPointerToOuterGhostNodes() {return &ghost_nodes_outer;}

  void AXPlusBY(double a, SpaceVariable3D &X, double b, SpaceVariable3D &Y, bool workOnGhost = false);

  void ComputeNormalDirection(SpaceVariable3D &Phi, SpaceVariable3D &NPhi);

  void ComputeNormal(SpaceVariable3D &Phi, SpaceVariable3D &NPhi);

  void ComputeUnitNormalAndCurvature(SpaceVariable3D &Phi, SpaceVariable3D &NPhi, SpaceVariable3D &KappaPhi);

  void Destroy();

  //! for debugging/testing the level set solver (Euler / N-S solver not activated)
  void PrescribeVelocityFieldForTesting(SpaceVariable3D &V, SpaceVariable3D &Phi, double time);

private:
  // functions for internal use within the class
  
  void CreateGhostNodeLists(); //almost the same as the function in SpaceOperator, except bcType

  //! Apply initial condition of phi within arbitrary enclosure(s) specified using a file
  bool ApplyInitialConditionWithinEnclosure(UserSpecifiedEnclosureData &enclosure, SpaceVariable3D &Phi);

  //! Compute derivatives of phi at nodes
  void ComputeNormalDirectionCentralDifferencing_FullDomain(SpaceVariable3D &Phi, SpaceVariable3D &NPhi);
  void ComputeNormalDirectionCentralDifferencing_NarrowBand(SpaceVariable3D &Phi, SpaceVariable3D &NPhi);

  //! Finite difference method
  void ComputeResidualFDM(SpaceVariable3D &V, SpaceVariable3D &Phi, SpaceVariable3D &R);
  void ComputeResidualFDM_FullDomain(SpaceVariable3D &V, SpaceVariable3D &Phi, SpaceVariable3D &R);
  void ComputeResidualFDM_NarrowBand(SpaceVariable3D &V, SpaceVariable3D &Phi, SpaceVariable3D &R);

  //! Finite volume method
  void ComputeResidualFVM(SpaceVariable3D &V, SpaceVariable3D &Phi, SpaceVariable3D &R);

  void Reconstruct(SpaceVariable3D &V, SpaceVariable3D &Phi);
  void ComputeAdvectionFlux(SpaceVariable3D &R);
  void AddSourceTerm(SpaceVariable3D &Phi, SpaceVariable3D &R);

  void ReconstructInBand(SpaceVariable3D &V, SpaceVariable3D &Phi); //!< the narrow-band version
  void ComputeAdvectionFluxInBand(SpaceVariable3D &R); //!< the narrow-band version
  void AddSourceTermInBand(SpaceVariable3D &Phi, SpaceVariable3D &R); //!< the narrow-band version

  double ComputeLocalAdvectionFlux(double phim, double phip, double um, double up);

  //! Utility
  inline double CentralDifferenceLocal(double phi0, double phi1, double phi2, double x0, double x1, double x2) {
    double c0 = -(x2-x1)/((x1-x0)*(x2-x0));
    double c1 = 1.0/(x1-x0) - 1.0/(x2-x1);
    double c2 = (x1-x0)/((x2-x0)*(x2-x1));
    return c0*phi0 + c1*phi1 + c2*phi2;
  }

  inline double SecondOrderDifference(double phi0, double phi1, double phi2, double x0, double x1, double x2) {
    return 2.0*phi0 / ((x0-x1)*(x0-x2)) + 2.0*phi1 / ((x1-x0)*(x1-x2)) + 2*phi2 / ((x2-x0)*(x2-x1));
  }

};

#endif
