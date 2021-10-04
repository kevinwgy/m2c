#ifndef _OUTPUT_H_
#define _OUTPUT_H_
#include <IoData.h>
#include <VarFcnBase.h>
#include <SpaceVariable.h>
#include <ProbeOutput.h>
#include <MaterialVolumeOutput.h>
#include <IonizationOperator.h>
#include <stdio.h>

/** Class Output is responsible  for writing solutions to files. It uses PETSc functionalities
 *  to write VTK files. It should not try to do post-processing work that can be done by some
 *  external software (e.g., Paraview). Keep it as simple as possible.
 */
class Output
{
  MPI_Comm& comm;
  IoData& iod;
  vector<VarFcnBase*> &vf;

  //! Ionization solver (Currently, a post-processer)
  IonizationOperator* ion;

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
  Output(MPI_Comm &comm_, DataManagers3D &dms, IoData &iod_, vector<VarFcnBase*> &vf_, SpaceVariable3D &cell_volume,
         IonizationOperator* ion_ = NULL);
  ~Output();

  void InitializeOutput(SpaceVariable3D &coordinates); //!< attach mesh

  void OutputSolutions(double time, double dt, int time_step, SpaceVariable3D &V,
                       SpaceVariable3D &ID, std::vector<SpaceVariable3D*> &Phi, 
                       SpaceVariable3D *L/*laser radiance*/, bool force_write);

  void WriteSolutionSnapshot(double time, int time_step, SpaceVariable3D &V, SpaceVariable3D &ID,
                             vector<SpaceVariable3D*> &Phi, SpaceVariable3D *L); //!< write solution to file

  void FinalizeOutput();

private:
  void OutputMeshInformation(SpaceVariable3D& coordinates);

};

#endif
