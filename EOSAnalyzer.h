/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _EOS_ANALYZER_
#define _EOS_ANALYZER_

#include<VarFcnBase.h>
#include<IoData.h>
#include<vector>


/************************************************************************
 * class EOSAnalyzer contains tools for analyzing EOS, including creating
 * 1D or 2D tables.
 ***********************************************************************/

class EOSAnalyzer {

  std::vector<VarFcnBase*>& vf;

  ObjectMap<EOSTabulationData> &eos_tabulationMap; 

public:

  EOSAnalyzer(ObjectMap<EOSTabulationData> &eos_tabulationMap_, std::vector<VarFcnBase*>& vf_);
  ~EOSAnalyzer();

  void GenerateAllEOSTables();

private:

  void GenerateEOSTable(int table_id, EOSTabulationData& tab);

};


#endif
