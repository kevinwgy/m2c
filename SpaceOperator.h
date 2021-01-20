#ifndef _SPACEOPERATOR_H_
#define _SPACEOPERATOR_H_
#include <petscdmda.h>
#include <IoData.h>
#include <VarFcnBase.h>
#include <FluxFcnBase.h>
#include <SpaceVariable.h>
#include <Reconstructor.h>

/*******************************************
 * class SpaceOperator drives computations
 * that require domain/mesh information
 ******************************************/
class SpaceOperator
{
  MPI_Comm&       comm;
  DataManagers2D& dm_all;
  IoData&         iod;
  VarFcnBase&     vf; 
  FluxFcnBase&    ff;

  //! Mesh info
  SpaceVariable2D coordinates;
  SpaceVariable2D delta_xy;
  SpaceVariable2D volume; //!< volume of node-centered control volumes (area in 2D)
  
  int i0, j0, imax, jmax; //!< corners of the real subdomain
  int ii0, jj0, iimax, jjmax; //!< corners of the ghosted subdomain

  //! Class for spatial reconstruction
  Reconstructor rec;

  //! Reconstructed primitive state variables at cell boundaries
  SpaceVariable2D Vl, Vr, Vb, Vt;

public:
  SpaceOperator(MPI_Comm &comm_, DataManagers2D &dm_all_, IoData &iod_,
                VarFcnBase &vf_, FluxFcnBase &ff_); 
  ~SpaceOperator();

  void ConservativeToPrimitive(SpaceVariable2D &U, SpaceVariable2D &V, bool workOnGhost = false);
  void PrimitiveToConservative(SpaceVariable2D &V, SpaceVariable2D &U, bool workOnGhost = false);
  int  ClipDensityAndPressure(SpaceVariable2D &V, bool workOnGhost = false, bool checkState = true);

  void SetInitialCondition(SpaceVariable2D &V);
    
  void ApplyBoundaryConditions(SpaceVariable2D &V);

  void ComputeTimeStepSize(SpaceVariable2D &V, double &dt, double &cfl);

  //! Compute the RHS of the ODE system (Only for cells inside the physical domain)
  void ComputeResidual(SpaceVariable2D &V, SpaceVariable2D &R);
  void ComputeAdvectionFluxes(SpaceVariable2D &V, SpaceVariable2D &F);

  SpaceVariable2D& GetMeshCoordinates() {return coordinates;}

  void Destroy();

private:
  void SetupMesh();
  void SetupMeshUniformRectangularDomain();
  void PopulateGhostBoundaryCoordinates();
};


#endif
