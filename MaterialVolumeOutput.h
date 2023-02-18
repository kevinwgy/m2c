/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _MATERIAL_VOLUME_OUTPUT_H_
#define _MATERIAL_VOLUME_OUTPUT_H_
#include <SpaceVariable.h>

class IoData;

/** This class is responsible for calculating and outputing the volume of each
 *  material subdomain.
 */
class MaterialVolumeOutput {

  MPI_Comm &comm;

  SpaceVariable3D& cell_volume;

  int frequency;
  double frequency_dt;
  double last_snapshot_time;

  FILE* file;

  int numMaterials;

public:

  MaterialVolumeOutput(MPI_Comm &comm_, IoData &iod, SpaceVariable3D& cell_volume_);
  ~MaterialVolumeOutput();

  void WriteSolution(double time, double dt, int time_step, SpaceVariable3D& ID, bool force_write);  

private:

  void ComputeMaterialVolumes(SpaceVariable3D& ID, double* vol);

};

#endif
