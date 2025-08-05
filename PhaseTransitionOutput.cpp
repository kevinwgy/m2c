#include <PhaseTransitionOutput.h>
#include <Utils.h>

// -----------------------------------------------------------

PhaseTransitionOutput::PhaseTransitionOutput(MPI_Comm &comm_, OutputData &iod_output_, MultiPhaseOperator &mpo_)
                     : comm(comm_), iod_output(iod_output_), mpo(mpo_), file(NULL)
{
  assert(strcmp(iod_output.mat_transition.filename, "")); //otherwise, shouldn't construct this object

  last_snapshot_time = -1.0;

  char* filename = new char[strlen(iod_output.prefix)+strlen(iod_output.mat_transition.filename)+1];
  sprintf(filename,"%s%s", iod_output.prefix, iod_output.mat_transition.filename);

  //print header
  file = fopen(filename, "w");
  if(!file) {
    print_error("*** Error: Cannot open file %s for output.\n", filename);
    exit_mpi();
  }
  print(file, "## Phase/Material-Type Transition Statistics\n");
  print(file, "## Time step | Time | Transitioned Energy (total, new) | Dumped Energy (total, new)\n");
  mpi_barrier();
  fflush(file);

  delete [] filename;
}

// -----------------------------------------------------------

PhaseTransitionOutput::~PhaseTransitionOutput()
{
  if(file)
    fclose(file);
}

// -----------------------------------------------------------

void
PhaseTransitionOutput::WriteStatsToFile(double time, double dt, int time_step, bool force_write)
{
  if(!isTimeToWrite(time, dt, time_step, iod_output.mat_transition.frequency_dt,
                    iod_output.mat_transition.frequency, last_snapshot_time, force_write))
    return; //nothing to do

  double energy_transitioned(0.0), energy_transitioned_new(0.0);
  double energy_dumped(0.0), energy_dumped_new(0.0);
  mpo.GetPhaseTransitionStats(energy_transitioned, energy_dumped, energy_transitioned_new,
                              energy_dumped_new);

  if(!force_write &&
     (iod_output.mat_transition.skip_no_transition_timesteps == MaterialTransitionOutputData::YES) &&
     (energy_transitioned_new == 0.0) && (energy_dumped_new == 0.0))
    return; //nothing to do

  //Now, open file and append a line
  assert(file);
  print(file, "%10d  %16.8e  %16.8e  %16.8e  %16.8e  %16.8e\n", time_step, time, energy_transitioned,
        energy_transitioned_new, energy_dumped, energy_dumped_new);
  mpi_barrier();
  fflush(file);

  last_snapshot_time = time;
}

// -----------------------------------------------------------


