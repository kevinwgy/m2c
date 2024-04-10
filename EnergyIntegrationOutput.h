#ifndef _ENERGY_INTEGRATION_OUTPUT_H_
#define _ENERGY_INTEGRATION_OUTPUT_H_

#include <IoData.h>
#include <SpaceVariable.h>
#include <LaserAbsorptionSolver.h>
#include <VarFcnBase.h>
#include<Vector5D.h>

/*********************************************************
 * class EnergyIntegrationOperator integrates the energy 
 * required by user within the specified region
 * ********************************************************/
using namespace std;

class EnergyIntegrationOutput
{
  MPI_Comm &comm;
  OutputData &iod_output;
  std::vector<VarFcnBase*> &vf;

  MeshData& iod_mesh;

  EquationsData& iod_eqs;

  LaserAbsorptionSolver* laser;

  int frequency;
  double frequency_dt;
  double last_snapshot_time;

  FILE *file[EnergyIntegrationData::SIZE]; //!< one file per solution variable

  int numMaterials;

  //! Mesh info
  SpaceVariable3D& coordinates;
  SpaceVariable3D& delta_xyz;
  SpaceVariable3D& cell_volume;


public:

  EnergyIntegrationOutput(MPI_Comm &comm_, IoData &iod, OutputData &iod_output_, 
                          MeshData &iod_mesh_, EquationsData &iod_eqs_, 
                          LaserAbsorptionSolver* laser_,
                          std::vector<VarFcnBase*> &vf_, SpaceVariable3D& coordinates_,
                          SpaceVariable3D& delta_xyz_, SpaceVariable3D& cell_volume_);

  ~EnergyIntegrationOutput();

  void WriteSolutionOfIntegrationEnergy(double time, double dt, int time_step, SpaceVariable3D &V, SpaceVariable3D &ID,
                                       SpaceVariable3D* L, bool force_write);

private:

  void IntegrateVolume(SpaceVariable3D &ID, double* vol);
  void IntegrateMass(SpaceVariable3D &V, SpaceVariable3D &ID, double* mass);
  void IntegrateTotalEnergy(SpaceVariable3D &V, SpaceVariable3D &ID, double* E);
  void IntegrateTotalEnthalpy(SpaceVariable3D &V, SpaceVariable3D &ID, double* H);
  void IntegrateKineticEnergy(SpaceVariable3D &V, SpaceVariable3D &ID, double* kinetic);
  void IntegrateInternalEnergy(SpaceVariable3D &V, SpaceVariable3D &ID, double* internal);
  void IntegratePotentialEnergy(SpaceVariable3D &V, SpaceVariable3D &ID, double* potential);
  void IntegrateLaserRadiation(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D *L, double* radiation);
  

};

#endif
