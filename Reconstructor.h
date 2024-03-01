/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _RECONSTRUCTOR_H_
#define _RECONSTRUCTOR_H_
#include <IoData.h>
#include <SpaceVariable.h>
#include <FluxFcnBase.h>
#include <GhostPoint.h>
#include <EmbeddedBoundaryDataSet.h>
using std::min;
using std::max;

/*****************************************************************************
 * Class Reconstructor handles spatial reconstruction and slope limiters. It
 * includes both "general" functions for reconstructing any variables and a
 * special function particularly for reconstructing fluid state variables,
 * including conversion among primitive, conservative, and characteristic
 * variables.
 * Multiple instantiations of Class Reconstructor can be created in the solver
 * for reconstructing different variables.
 ****************************************************************************/
class Reconstructor
{

  MPI_Comm& comm; //!< may not contain all the processors (e.g., see LaserAbsorptionSolver::SetupLoadBalancing)

  ReconstructionData& iod_rec;

  /** Variable functions and flux functions. This is optional, only used for
    * reconstructing the conservative or characteristic fluid state variables.
    * By default, the pointers are set to NULL. */
  vector<VarFcnBase*>* varFcn;
  FluxFcnBase* fluxFcn;

  //! Mesh info
  SpaceVariable3D &coordinates;
  SpaceVariable3D &delta_xyz;
  
  //! Pointer to ghost node lists
  vector<GhostPoint> *ghost_nodes_inner;
  vector<GhostPoint> *ghost_nodes_outer;

  /** Internal variables
   *  CoeffA, CoeffB, and CoeffK are only defined within the real domain
   */
  SpaceVariable3D CoeffA;  //!< constant coefficient for non-uniform meshes, dim=3 for 3D 
  SpaceVariable3D CoeffB;  //!< constant coefficient for non-uniform meshes, dim=3 for 3D 
  SpaceVariable3D CoeffK;  //!< the integer k in Zeng's 2016 paper for Van Albada (stored as real numbers)

  /** Internal variable to tag fixed nodes (only for nodes inside physical domain) */
  SpaceVariable3D* FixedByUser;

  /** Another internal variable for var. conversion*/
  SpaceVariable3D U;

public:
  Reconstructor(MPI_Comm &comm_, DataManagers3D &dm_all_, ReconstructionData &iod_rec_, 
                SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_, 
                vector<VarFcnBase*>* vf_ = NULL, FluxFcnBase* ff_ = NULL);

  ~Reconstructor();

  void Setup(vector<GhostPoint>* inner, vector<GhostPoint>* outer); //!< compute AB and K

  /** Linear reconstruction in 3D 
    * The input and output variables are assumed to be primitive variables. If IoData specifies a
    * different variable to be reconstructed (e.g., primitive or characterstic), a conversion is
    * performed within this function. */
  void Reconstruct(SpaceVariable3D &V, SpaceVariable3D &Vl, SpaceVariable3D &Vr, 
           SpaceVariable3D &Vb, SpaceVariable3D &Vt, SpaceVariable3D &Vk, SpaceVariable3D &Vf,
           SpaceVariable3D *ID = NULL,
           vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS = nullptr,
           SpaceVariable3D *Selected = NULL,
           bool do_nothing_if_not_selected = true); //!< the last input is used (only?) in LevelSetOperator

  /** This function applies reconstruction directly to the input variable U. In other words, no
    * conversions are done inside the function.*/
  void ReconstructIn1D(int dir/*0~x,1~y,2~z*/, SpaceVariable3D &U, SpaceVariable3D &Um, SpaceVariable3D &Up,
                       SpaceVariable3D *Slope = NULL, SpaceVariable3D *Selected = NULL);

  void Destroy(); //!< destroy the SpaceVariable3D variables

private:

  int CalculateSlopeLimiterCoefficientK(double A, double B); //calculate k in Zeng's 2016 paper

  void TagNodesFixedByUser();
  
  /** slope limiter functions
   * (frequently called; need to be as fast as possible)
   */
  inline double GeneralizedMinMod(double A, double B, double alpha, double dq0, double dq1) {
    if(dq0<0.0 && dq1<0.0) 
      return max(alpha*dq0, max(B/(A+1.0)*(dq0+dq1), alpha*dq1));
    else if(dq0>0.0 && dq1>0.0) 
      return min(alpha*dq0, min(B/(A+1.0)*(dq0+dq1), alpha*dq1));
    else
      return 0.0;
  }

  inline double VanAlbada(double A, double B, int k, double dq0, double dq1) {
    if(dq0*dq1<=0) return 0.0;
    if(fabs(dq0*dq1)<1e-16) return 0.0;
    double dq0_power_k = pow(dq0,k);
    double dq1_power_k = pow(dq1,k);
    return B*(dq1*dq0_power_k + dq0*dq1_power_k)/(dq0_power_k + A*dq1_power_k);
  }

};

#endif
