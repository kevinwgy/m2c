#ifndef _PHASE_TRANSITION_OUTPUT_H_
#define _PHASE_TRANSITION_OUTPUT_H_

#include <MultiPhaseOperator.h>

/********************************************************
 * class PhaseTransitionOutput is responsible for printing
 * phase/material transition statistics to a text file.
 * Computation is mostly done by MultiPhaseOperator.
 *******************************************************/

class PhaseTransitionOutput
{
  MPI_Comm &comm;
  OutputData &iod_output;

  MultiPhaseOperator &mpo;

  FILE* file;

  double last_snapshot_time;

public:

  PhaseTransitionOutput(MPI_Comm &comm_, OutputData &iod_output_, MultiPhaseOperator &mpo_);
  ~PhaseTransitionOutput();

  void WriteStatsToFile(double time, double dt, int time_step, bool force_write);

};

#endif
