/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <SpecialToolsDriver.h>
#include <DynamicLoadCalculator.h>
#include <EOSAnalyzer.h>
#include <cassert>

//------------------------------------------------------------

SpecialToolsDriver::SpecialToolsDriver(IoData &iod_, std::vector<VarFcnBase*>& vf_, MPI_Comm &comm_, 
                                       ConcurrentProgramsHandler &concurrent_)
                  : iod(iod_), vf(vf_), comm(comm_), concurrent(concurrent_)
{
  assert(iod.special_tools.type != SpecialToolsData::NONE);
}

//------------------------------------------------------------

SpecialToolsDriver::~SpecialToolsDriver()
{ }

//------------------------------------------------------------

void SpecialToolsDriver::Run()
{
  if(iod.special_tools.type == SpecialToolsData::DYNAMIC_LOAD_CALCULATION) {
    print("\n");
    print("----------------------------------------------------\n");
    print("- Activated special tool: Dynamic load calculator. -\n");
    print("----------------------------------------------------\n");
    print("\n");

    DynamicLoadCalculator load(iod, comm, concurrent);
    load.Run();
    print("\n");
  } 
  else if(iod.special_tools.type == SpecialToolsData::SpecialToolsData::EOS_TABULATION) {
    print("\n");
    print("----------------------------------------------------\n");
    print("- Activated special tool: EOS Analyzer. -\n");
    print("----------------------------------------------------\n");
    print("\n");

    EOSAnalyzer eos_analyzer(iod.special_tools.eos_tabulationMap, vf);
    eos_analyzer.GenerateAllEOSTables();
    print("\n");
  }
  else {
    print_error("*** Error: Detected unknown type (%d) in SpecialTools.\n", 
                (int)iod.special_tools.type);  
    exit_mpi();
  }
}

//------------------------------------------------------------

