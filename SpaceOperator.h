#ifndef _SPACEOPERATOR_H_
#define _SPACEOPERATOR_H_
#include <petscdmda.h>
#include <IoData.h>
#include <VarFcnBase.h>
#include <SpaceVariable.h>

class SpaceOperator
{
  MPI_Comm&       comm;
  DataManagers2D& dm_all;
  IoData&         iod;
  VarFcnBase&     vf; 

  // Mesh info
  SpaceVariable2D coordinates;
  SpaceVariable2D delta_xy;
  SpaceVariable2D face_rt; //right and top faces
  SpaceVariable2D volume; //volume of node-centered control volumes (area in 2D)
  
public:
  SpaceOperator(MPI_Comm &comm_, DataManagers2D &dm_all_, IoData &iod_,
                VarFcnBase &vf_);
  ~SpaceOperator();

  void ConservativeToPrimitive(SpaceVariable2D &U, SpaceVariable2D &V, bool workOnGhost = true);
  void PrimitiveToConservative(SpaceVariable2D &V, SpaceVariable2D &U, bool workonGhost = true);

  void SetInitialCondition(SpaceVariable2D &V);
  void applyBoundaryConditions(SpaceVariable2D &V);
    
  void ComputeAdvectionFluxes(&U, &F);

private:
  void SetupMesh();
  void SetupNodalCoordinatesUniformRectangularDomain();
  void PopulateGhostNodalCoordinates();
}


#endif
