/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _LEVELSET_REINITIALIZER_H_
#define _LEVELSET_REINITIALIZER_H_

#include <GradientCalculatorBase.h>
#include <Interpolator.h>
#include <GhostPoint.h>

/*****************************************************************************
 * Class LevelSetReinitializer handles the reinitialization of the level set
 * function, which means restoring it to a signed distance function (i.e.
 * |grad(phi)| = 1). A pseudo-time integration is performed.
 ****************************************************************************/
class LevelSetReinitializer
{
  MPI_Comm& comm;
  LevelSetSchemeData &iod_ls;

  //! Mesh info
  SpaceVariable3D& coordinates;
  SpaceVariable3D& delta_xyz;
  double min_dxyz;
  double max_dxyz;

  vector<GhostPoint>& ghost_nodes_inner; //!< ghost nodes inside the physical domain (shared with other subd)
  vector<GhostPoint>& ghost_nodes_outer; //!< ghost nodes outside the physical domain

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain

  //! Internal variables
  SpaceVariable3D Tag;
  SpaceVariable3D R;
  SpaceVariable3D Phibk; //!< save input Phi for "failsafe"
  SpaceVariable3D Phi0;
  SpaceVariable3D Phi1;
  SpaceVariable3D Sign;

  SpaceVariable3D NormalDir; //grad(phi)/|grad(phi)|
 
  SpaceVariable3D PhiG2; //this one has 2 ghosted layers

  double phi_max, phi_min; //max pos and neg phi within band
  double phi_out_pos, phi_out_neg; //applied to specify the constant phi outside band

  vector<Int3> useful_nodes_plus1layer; //useful nodes plus one layer outside bandwidth

  //! cfl number for the ficitious time integrator
  double cfl;

  //! a nested structure that stores information of a first-layer node & its neighbors across interface
  struct FirstLayerNode {
    int i,j,k;
    bool s[6]; //neighbor status l,r,b,t,k,f, (in set S or not)
    int ns; //size of set S
    double r[6]; //coefficient for HCR-1
    double r0; //coefficient for HCR-2

    double f; //correction

    Vec3D nphi0; //unit normal calculated using the original phi (grad(phi) / |grad(phi)|)

    FirstLayerNode() : i(0),j(0),k(0),ns(0),r0(0),f(0) {
      for(int ii=0; ii<6; ii++) {
        s[ii] = false;
        r[ii] = 0.0;
      }
    }
    FirstLayerNode(int i_, int j_, int k_) : i(i_),j(j_),k(k_),ns(0),r0(0),f(0) {
      for(int ii=0; ii<6; ii++) {
        s[ii] = false;
        r[ii] = 0.0;
      }
    }
    ~FirstLayerNode() {}
  };


public:

  LevelSetReinitializer(MPI_Comm &comm_, DataManagers3D &dm_all_, LevelSetSchemeData &iod_ls_,
                        SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_,
                        vector<GhostPoint> &ghost_nodes_inner_, vector<GhostPoint> &ghost_nodes_outer_);

  ~LevelSetReinitializer();

  void Destroy();

  void ReinitializeFullDomain(SpaceVariable3D &Phi, int special_maxIts = 0); //reinitialize phi in the entire domain

  void ReinitializeInBand(SpaceVariable3D &Phi, SpaceVariable3D &Level, SpaceVariable3D &UsefulG2,
                          SpaceVariable3D &Active, vector<Int3> &useful_nodes, vector<Int3> &active_nodes,
                          int special_maxIts = 0);

  // For narrow-band level set method
  void ConstructNarrowBand(SpaceVariable3D &Phi,
                           SpaceVariable3D &Level, SpaceVariable3D &UsefulG2, SpaceVariable3D &Active,
                           vector<Int3> &useful_nodes, vector<Int3> &active_nodes);

private:

  // Functions that work for both full-domain and narrow-band level set methods
  //
  void ReinitializeFirstLayerNodes(SpaceVariable3D &Phi0, SpaceVariable3D &Phi, vector<FirstLayerNode> &firstLayer,
                                   SpaceVariable3D *UsefulG2 = NULL, vector<Int3> *useful_nodes = NULL);

  double DifferentiateInFirstLayer(double x0, double x1, double x2,
                                   double tag0, double tag1, double tag2,
                                   double phi0, double phi1, double phi2,
                                   double phi00, double phi3, double eps);

  void ApplyBoundaryConditions(SpaceVariable3D &Phi, SpaceVariable3D *UsefulG2 = NULL);

  void PopulatePhiG2(SpaceVariable3D &Phi, vector<Int3> *useful_nodes = NULL);

  void ComputeNormalDirectionBeyondFirstLayer(SpaceVariable3D &Phi, vector<FirstLayerNode> &firstLayer, 
                              SpaceVariable3D *UsefulG2 = NULL, vector<Int3> *useful_nodes = NULL);

  inline double CentralDifferenceLocal(double phi0, double phi1, double phi2, double x0, double x1, double x2) {
    double c0 = -(x2-x1)/((x1-x0)*(x2-x0));
    double c1 = 1.0/(x1-x0) - 1.0/(x2-x1);
    double c2 = (x1-x0)/((x2-x0)*(x2-x1));
    return c0*phi0 + c1*phi1 + c2*phi2;
  }


  // For full-domain level set method
  bool TagFirstLayerNodes(SpaceVariable3D &Phi, vector<FirstLayerNode> &firstLayer);

  void ApplyCorrectionToFirstLayerNodes(SpaceVariable3D &Phi, vector<FirstLayerNode> &firstLayer, double cfl,
                                        SpaceVariable3D *UsefulG2 = NULL);

  void EvaluateSignFunctionFullDomain(SpaceVariable3D &Phi, double eps/*smoothing factor*/);

  double ComputeResidualFullDomain(SpaceVariable3D &Phi, SpaceVariable3D &R, double cfl);

  double CalculateMaximumRelativeErrorFullDomain(SpaceVariable3D &Phi0, SpaceVariable3D &Phi);


  // For narrow-band level set method
  bool TagFirstLayerNodesInBand(SpaceVariable3D &Phi, vector<Int3> &useful_nodes,
                                vector<FirstLayerNode> &firstLayer, vector<Int3> &firstLayerIncGhost);

  void UpdateNarrowBand(SpaceVariable3D &Phi, vector<Int3> &firstLayerIncGhost,
                        SpaceVariable3D &Level, SpaceVariable3D &UsefulG2, SpaceVariable3D &Active,
                        vector<Int3> &useful_nodes, vector<Int3> &active_nodes);

  void PropagateNarrowBand(SpaceVariable3D &Level, SpaceVariable3D &UsefulG2,
                           SpaceVariable3D &Active, vector<Int3> &useful_nodes,
                           vector<Int3> &active_nodes);

  void CutOffPhiOutsideBand(SpaceVariable3D &Phi, SpaceVariable3D &UsefulG2, vector<Int3> &useful_nodes);

  void EvaluateSignFunctionInBand(SpaceVariable3D &Phi, vector<Int3> &useful_nodes, double eps/*smoothing factor*/);

  double ComputeResidualInBand(SpaceVariable3D &Phi, SpaceVariable3D &UsefulG2,
                               vector<Int3> &useful_nodes, SpaceVariable3D &R, double cfl);

  void UpdatePhiMaxAndPhiMinInBand(SpaceVariable3D &Phi, vector<Int3> &useful_nodes);

  void CreateUsefulNodesPlusOneLayer(vector<Int3> &useful_nodes, SpaceVariable3D &UsefulG2);

  void AXPlusBYInBandPlusOne(double a, SpaceVariable3D &X, double b, SpaceVariable3D &Y, bool workOnGhost = false);

  double CalculateMaximumRelativeErrorInBand(SpaceVariable3D &Phi0, SpaceVariable3D &Phi, vector<Int3> &useful_nodes);

};

#endif
