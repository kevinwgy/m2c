#include<SpecialToolsDriver.h>
#include<DynamicLoadCalculator.h>
#include<cassert>

//------------------------------------------------------------

SpecialToolsDriver::SpecialToolsDriver(IoData &iod_, MPI_Comm &comm_, 
                                       ConcurrentProgramsHandler &concurrent_)
                  : iod(iod_), comm(comm_), concurrent(concurrent_)
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
    DynamicLoadCalculator load(iod, comm, concurrent);
    load.Run();
  } else {
    print_error("*** Error: Detected unknown type (%d) in SpecialTools.\n", 
                (int)iod.special_tools.type);  
    exit_mpi();
  }
}

//------------------------------------------------------------

