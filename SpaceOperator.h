#ifndef _SPACEOPERATOR_H_
#define _SPACEOPERATOR_H_
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
  MPI_Comm&                 comm;
  IoData&                   iod;
  FluxFcnBase&              fluxFcn;

  vector<VarFcnBase*>& varFcn; //!< each material has a varFcn

  //! Mesh info
  SpaceVariable3D coordinates;
  SpaceVariable3D delta_xyz;
  SpaceVariable3D volume; //!< volume of node-centered control volumes
  
  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain

  //! Class for spatial reconstruction
  Reconstructor rec;

  //! Reconstructed primitive state variables at cell boundaries
  SpaceVariable3D Vl, Vr, Vb, Vt, Vk, Vf;

public:
  SpaceOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_,
                vector<VarFcnBase*> &varFcn_, FluxFcnBase &fluxFcn_); 
  ~SpaceOperator();

  void ConservativeToPrimitive(SpaceVariable3D &U, SpaceVariable3D &ID, SpaceVariable3D &V
                               bool workOnGhost = false);
  void PrimitiveToConservative(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &U
                               bool workOnGhost = false);
  int  ClipDensityAndPressure(SpaceVariable3D &V, SpaceVariable3D &ID, 
                              bool workOnGhost = false, bool checkState = true);

  void SetInitialCondition(SpaceVariable3D &V, SpaceVariable3D &ID);
    
  void ApplyBoundaryConditions(SpaceVariable3D &V);

  void FindExtremeValuesOfFlowVariables(SpaceVariable3D &V, SpaceVariable3D &ID,
                                        double *Vmin, double *Vmax, double &cmin, 
                                        double &cmax, double &Machmax, double &char_speed_max,
                                        double &dx_over_char_speed_min);

  void ComputeTimeStepSize(SpaceVariable3D &V, SpaceVariable3D &ID, double &dt, double &cfl);

  //! Compute the RHS of the ODE system (Only for cells inside the physical domain)
  void ComputeResidual(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &R);

  SpaceVariable3D& GetMeshCoordinates() {return coordinates;}
  SpaceVariable3D& GetMeshDeltaXYZ()    {return delta_xyz;}
  SpaceVariable3D& GetMeshCellVolumes() {return volume;}

  void Destroy();

private:
  void SetupMesh();
  void SetupMeshUniformRectangularDomain();
  void PopulateGhostBoundaryCoordinates();

  void ComputeAdvectionFluxes(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &F);
};


#endif
