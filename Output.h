#ifndef _OUTPUT_H_
#define _OUTPUT_H_
#include <IoData.h>
#include <VarFcnBase.h>
#include <SpaceVariable.h>
#include <stdio.h>

class Output
{
  MPI_Comm& comm;
  IoData& iod;
  VarFcnBase& vf;

  //These variables will temporarily hold solutions before they are printed to file
  SpaceVariable2D scalar;   
  SpaceVariable2D vector3;

  int iFrame;

  FILE* pvdfile;

public:
  Output(MPI_Comm &comm_, DataManagers2D &dms, IoData &iod_, VarFcnBase &vf_);
  ~Output();

  void InitializeOutput(SpaceVariable2D &coordinates); //attach mesh
  void WriteSolutionSnapshot(double time, SpaceVariable2D &V);
  void FinalizeOutput();

};

#endif
