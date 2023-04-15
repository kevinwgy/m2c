/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _SPECIAL_TOOLS_DRIVER_H_
#define _SPECIAL_TOOLS_DRIVER_H_

#include<ConcurrentProgramsHandler.h>
#include<VarFcnBase.h>

/*****************************************
 * class SpecialToolsDriver is the driver
 * that runs special tools implemented
 * in the M2C code.
 ****************************************/

class SpecialToolsDriver 
{

  MPI_Comm& comm;
  IoData& iod;
  std::vector<VarFcnBase*>& vf;
  ConcurrentProgramsHandler& concurrent;

public:
  
  SpecialToolsDriver(IoData &iod_, std::vector<VarFcnBase*>& vf_,
                     MPI_Comm &comm_, ConcurrentProgramsHandler &concurrent_);
  ~SpecialToolsDriver();

  void Run();
 
};

#endif
