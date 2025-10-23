/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _INTEGRATION_OUTPUT_H_
#define _INTEGRATION_OUTPUT_H_

#include <IoData.h>
#include <SpaceVariable.h>
#include <LaserAbsorptionSolver.h>
#include <VarFcnBase.h>
#include <Vector5D.h>

/*********************************************************
 * class IntegrationOutput integrates state variables within
 * user-specified regions and output to text files. It is
 * a generalization of EnergyIntegrationOutput and MaterialVolumeOutput.
 * ********************************************************/

class IntegrationOutput
{
  MPI_Comm &comm;
  OutputData &iod_output;
  std::vector<VarFcnBase*> &vf;

  LaserAbsorptionSolver* laser;

  //! User-specified info for all the integration regions
  std::vector<int> frequency;
  std::vector<double> frequency_dt;
  std::vector<double> last_snapshot_time;

  std::vector<std::vector<FILE*> > files; //!< one file per solution variable

  int numMaterials; //!< includes ``ghost/inactive'' id.
  MeshData::Type mesh_type; //!< whether there is spherical or spherical symmetry

  //! Mesh info
  SpaceVariable3D& coordinates;
  SpaceVariable3D& delta_xyz;
  SpaceVariable3D& cell_volume;
  int i0, j0, k0, imax, jmax, kmax;

  //! Tags nodes/cells that need to be integrated
  //The p-th bit (binary digit) of tag[k][j][i] indicates whether [k][j][i] is integrated for
  //the p-th integral.
  SpaceVariable3D Tag;

public:

  IntegrationOutput(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod, LaserAbsorptionSolver* laser_,
                    std::vector<VarFcnBase*> &vf_, SpaceVariable3D& coordinates_,
                    SpaceVariable3D& delta_xyz_, SpaceVariable3D& cell_volume_);

  ~IntegrationOutput();
  void Destroy();

  void WriteIntegrationResults(double time, double dt, int time_step, SpaceVariable3D &V, SpaceVariable3D &ID,
                               SpaceVariable3D* L, bool force_write);

private:

  void SetupIntegrationDomain(double*** tag, int index, IntegrationData &integral); //!< populates one *bit* of tag
  int OrderUserSpecifiedGeometries(IntegrationData &integral, std::vector<std::pair<int,int> > &order);
  void AddGeomToVector(int o, int type, int ind, string name, std::vector<std::pair<int,int> > &order,
                       std::vector<int> &user_specified_order); //!< same as the function in SpaceInitializer


  void IntegrateVolume(int index, double*** tag, Vec3D*** coords, Vec3D*** dxyz, double*** cell,
                       double*** id, double* volume);
  void IntegrateMass(int index, double*** tag, Vec3D*** coords, Vec3D*** dxyz, double*** cell, Vec5D*** v,
                     double*** id, double* mass);
  void IntegrateMomentum(int index, double*** tag, Vec3D*** coords, Vec3D*** dxyz, double*** cell, Vec5D*** v,
                         double*** id, Vec3D* momentum);
  void IntegrateTotalEnergy(int index, double*** tag, Vec3D*** coords, Vec3D*** dxyz, double*** cell, Vec5D*** v,
                            double*** id, double* E);
  void IntegrateTotalEnthalpy(int index, double*** tag, Vec3D*** coords, Vec3D*** dxyz, double*** cell, Vec5D*** v,
                              double*** id, double* H);
  void IntegrateKineticEnergy(int index, double*** tag, Vec3D*** coords, Vec3D*** dxyz, double*** cell, Vec5D*** v,
                              double*** id, double* kinetic);
  void IntegrateInternalEnergy(int index, double*** tag, Vec3D*** coords, Vec3D*** dxyz, double*** cell, Vec5D*** v,
                               double*** id, double* internal);
  void IntegratePotentialEnergy(int index, double*** tag, Vec3D*** coords, Vec3D*** dxyz, double*** cell, Vec5D*** v,
                                double*** id, double* potential);
  void IntegrateLaserRadiation(int index, double*** tag, Vec3D*** coords, Vec3D*** dxyz, double*** cell, Vec5D*** v,
                               double*** id, double*** l, double* radiation);

};

#endif
