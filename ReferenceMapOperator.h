#ifndef _REFERENCE_MAP_OPERATOR_H_
#define _REFERENCE_MAP_OPERATOR_H_

#include<IoData.h>
#include<GradientCalculatorBase.h>
#include<GhostPoint.h>

/*********************************************************
 * class ReferenceMapOperator is responsible for
 * computing the function that maps current configuration
 * to the original configuration.
 * Similar to LevelSetOperator, except that the reference
 * map is a vector.
 *********************************************************/

class ReferenceMapOperator
{
  MPI_Comm& comm;
  IoData &iod;
  
  //! Mesh info
  SpaceVariable3D& coordinates;
  SpaceVariable3D& delta_xyz;
  SpaceVariable3D& volume; //!< volume of node-centered control volumes

  //! ghost nodes
  vector<GhostPoint> ghost_nodes_inner; //!< ghost nodes inside the physical domain (shared with other subd)
  vector<GhostPoint> ghost_nodes_outer; //!< ghost nodes outside the physical domain

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain

  //! Class for calculating spatial gradient (FDM)
  GradientCalculatorBase *grad_minus; //left-biased FD
  GradientCalculatorBase *grad_plus;  //right-biased FD


public:

  ReferenceMapOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_, SpaceVariable3D &coordinates_,
                       SpaceVariable3D &delta_xyz_, SpaceVariable3D &volume_, 
                       std::vector<GhostPoint> &ghost_nodes_inner_,
                       std::vector<GhostPoint> &ghost_nodes_outer_);  
  ~ReferenceMapOperator();

  void Destroy();

  void SetInitialCondition(SpaceVariable3D &Xi);

  void ApplyBoundaryConditions(SpaceVariable3D &Xi);

  void ComputeResidual(SpaceVariable3D &V, SpaceVariable3D &Xi, SpaceVariable3D &R);

private:

  // Utility
  inline double CentralDifferenceLocal(double phi0, double phi1, double phi2, double x0, double x1, double x2) {
    double c0 = -(x2-x1)/((x1-x0)*(x2-x0));
    double c1 = 1.0/(x1-x0) - 1.0/(x2-x1);
    double c2 = (x1-x0)/((x2-x0)*(x2-x1));
    return c0*phi0 + c1*phi1 + c2*phi2;
  }

};

#endif
