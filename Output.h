#ifndef _OUTPUT_H_
#define _OUTPUT_H_
#include <IoData.h>
#include <VarFcnBase.h>
#include <SpaceVariable.h>
#include <ProbeOutput.h>
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
  //LineOutput  line_output;


public:
  Output(MPI_Comm &comm_, DataManagers3D &dms, IoData &iod_, vector<VarFcnBase*> &vf_);
  ~Output();

  void InitializeOutput(SpaceVariable3D &coordinates); //!< attach mesh

  void WriteSolutionSnapshot(double time, int time_step, SpaceVariable3D &V, SpaceVariable3D &ID,
                             vector<SpaceVariable3D*> &Phi); //!< write solution to file

  void WriteSolutionAtProbes(double time, int time_step, SpaceVariable3D &V, SpaceVariable3D &ID,
                             vector<SpaceVariable3D*> &Phi); //!< write probe solution to file

  bool ToWriteSolutionSnapshot(double time, double dt, int time_step); /**< check whether to write solution 
                                                                        * at this time & time-step */
  inline double GetLastSnapshotTime() {return last_snapshot_time;}

  void FinalizeOutput();

};

#endif
