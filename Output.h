#ifndef _OUTPUT_H_
#define _OUTPUT_H_
#include <IoData.h>
#include <VarFcnBase.h>
#include <SpaceVariable.h>
#include <ProbeOutput.h>
#include <MaterialVolumeOutput.h>
#include <stdio.h>

/** This class is responsible (only) for writing solutions to files. It uses PETSc functionalities
 *  to write VTK files. It is not designed for post-processing results --- we leave this job to
 *  external software (e.g., Paraview). Keep it as simple as possible.
 */
class Output
{
  MPI_Comm& comm;
  IoData& iod;
  vector<VarFcnBase*> &vf;

  //! These variables will temporarily hold solutions before they are printed to file
  SpaceVariable3D scalar;   
  SpaceVariable3D vector3;

  int iFrame; // frame id.

  double last_snapshot_time; //!< latest time when solution snapshot is written to file

  FILE* pvdfile;

  ProbeOutput probe_output;
  std::vector<ProbeOutput*> line_outputs;

  MaterialVolumeOutput matvol_output;

public:
  Output(MPI_Comm &comm_, DataManagers3D &dms, IoData &iod_, vector<VarFcnBase*> &vf_, SpaceVariable3D &cell_volume);
  ~Output();

  void InitializeOutput(SpaceVariable3D &coordinates); //!< attach mesh

  void OutputSolutions(double time, double dt, int time_step, SpaceVariable3D &V,
                       SpaceVariable3D &ID, std::vector<SpaceVariable3D*> &Phi, bool force_write);

  void WriteSolutionSnapshot(double time, int time_step, SpaceVariable3D &V, SpaceVariable3D &ID,
                             vector<SpaceVariable3D*> &Phi); //!< write solution to file

  void FinalizeOutput();

private:
  void OutputMeshInformation(SpaceVariable3D& coordinates);

};

#endif
