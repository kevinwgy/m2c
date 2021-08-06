#ifndef _LEVELSET_REINITIALIZER_H_
#define _LEVELSET_REINITIALIZER_H_

#include <GradientCalculator.h>
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

  vector<GhostPoint>& ghost_nodes_inner; //!< ghost nodes inside the physical domain (shared with other subd)
  vector<GhostPoint>& ghost_nodes_outer; //!< ghost nodes outside the physical domain

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain

  //! Default gradient calculator
  InterpolatorBase *interp;
  GradientCalculatorBase *grad;

  //! Internal variables
  SpaceVariable3D Tag;
  SpaceVariable3D R;
  SpaceVariable3D Phi1;
  SpaceVariable3D Sign;
 
  SpaceVariable3D PhiG2; //this one has 2 ghosted layers

  double phi_max, phi_min; //max pos and neg phi within band
  double phi_out_pos, phi_out_neg; //applied to specify the constant phi outside band

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

  void Reinitialize(SpaceVariable3D &Phi); //reinitialize phi in the entire domain

  void ReinitializeInBand(SpaceVariable3D &Phi, SpaceVariable3D &Level, SpaceVariable3D &Useful,
                          SpaceVariable3D &Active, vector<Int3> &useful_nodes, vector<Int3> &active_nodes);

  // For narrow-band level set method
  void ConstructNarrowBand(SpaceVariable3D &Phi,
                           SpaceVariable3D &Level, SpaceVariable3D &Useful, SpaceVariable3D &Active,
                           vector<Int3> &useful_nodes, vector<Int3> &active_nodes);

private:

  void TagFirstLayerNodes(SpaceVariable3D &Phi, vector<FirstLayerNode> &firstLayer);

  void EvaluateSignFunction(SpaceVariable3D &Phi, double eps/*smoothing factor*/);

  void ReinitializeFirstLayerNodes(SpaceVariable3D &Phi0, SpaceVariable3D &Phi, vector<FirstLayerNode> &firstLayer);

  double DifferentiateInFirstLayer(double x0, double x1, double x2,
                                   double tag0, double tag1, double tag2,
                                   double phi0, double phi1, double phi2,
                                   double phi00, double phi3, double eps);

  double ComputeResidual(SpaceVariable3D &Phi, SpaceVariable3D &R, double cfl);

  void ApplyCorrectionToFirstLayerNodes(SpaceVariable3D &Phi, vector<FirstLayerNode> &firstLayer, double cfl);

  void ApplyBoundaryConditions(SpaceVariable3D &Phi, SpaceVariable3D *Useful = NULL);

  void PopulatePhiG2(SpaceVariable3D &Phi);

  // For narrow-band level set method
  void TagFirstLayerNodesInBand(SpaceVariable3D &Phi, vector<Int3> &useful_nodes,
                                vector<FirstLayerNode> &firstLayer, vector<Int3> &firstLayerIncGhost);

  void UpdateNarrowBand(SpaceVariable3D &Phi, vector<Int3> &firstLayerIncGhost,
                        SpaceVariable3D &Level, SpaceVariable3D &Useful, SpaceVariable3D &Active,
                        vector<Int3> &useful_nodes, vector<Int3> &active_nodes);

  void PropagateNarrowBand(SpaceVariable3D &Level, SpaceVariable3D &Useful,
                           SpaceVariable3D &Active, vector<Int3> &useful_nodes,
                           vector<Int3> &active_nodes);

  void CutOffPhiOutsideBand(SpaceVariable3D &Phi, SpaceVariable3D &Useful, vector<Int3> &useful_nodes);

  void EvaluateSignFunctionInBand(SpaceVariable3D &Phi, vector<Int3> &useful_nodes, double eps/*smoothing factor*/);

  double ComputeResidualInBand(SpaceVariable3D &Phi, SpaceVariable3D &Useful,
                               vector<Int3> &useful_nodes, SpaceVariable3D &R, double cfl);

  void UpdatePhiMaxAndPhiMinInBand(SpaceVariable3D &Phi, vector<Int3> &useful_nodes);

};

#endif
