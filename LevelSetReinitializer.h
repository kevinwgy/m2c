#ifndef _LEVELSET_REINITIALIZER_H_
#define _LEVELSET_REINITIALIZER_H_

#include <GradientCalculator.h>

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
  vector<Int3> firstLayer;
 
  SpaceVariable3D PhiG2; //this one has 2 ghosted layers

public:

  LevelSetReinitializer(MPI_Comm &comm_, DataManagers3D &dm_all_, LevelSetSchemeData &iod_ls_,
                        SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_);

  ~LevelSetReinitializer();

  void Destroy();

  void Reinitialize(SpaceVariable3D &Phi);

private:

  void TagFirstLayerNodes(SpaceVariable3D &Phi);

  void EvaluateSignFunction(SpaceVariable3D &Phi);

  void ReinitializeFirstLayerNodes(SpaceVariable3D &Phi0, SpaceVariable3D &Phi);

  double DifferentiateInFirstLayer(double x0, double x1, double x2,
                                   double tag0, double tag1, double tag2,
                                   double phi0, double phi1, double phi2,
                                   double phi00, double phi3, double eps);

  double ComputeResidual(SpaceVariable3D &Phi, SpaceVariable3D &R, double cfl);

  void ApplyBoundaryConditions(SpaceVariable3D &Phi);

  void PopulatePhiG2(SpaceVariable3D &Phi);

};

#endif
