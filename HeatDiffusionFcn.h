/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<IoData.h>
#include<VarFcnBase.h>
#include<cassert>

/******************************************************************
 * Class HeatDiffusionFcnBase and the derived classes are responsible
 * for calculating heat diffusion fluxes.
 * Each material has a seperate HeatDiffusionFcn
 ******************************************************************/

//-----------------------------------------------------------------------
//
class HeatDiffusionFcnBase {

protected:

  enum Type {NONE = 0, CONSTANT = 1} type;
  double conductivity;

public:

  HeatDiffusionFcnBase() : type(NONE), conductivity(0.0) {}
  virtual ~HeatDiffusionFcnBase() {}   

  inline int GetDiffusionFunctionType() {return (int)type;}
  inline double GetConductivity() {return conductivity;}

};

//--------------------------------------------------------------------------

class HeatDiffusionFcnConstant : public HeatDiffusionFcnBase {

public:
  
  HeatDiffusionFcnConstant(HeatDiffusionModelData &heatdiff) : HeatDiffusionFcnBase() {
    assert(heatdiff.type == HeatDiffusionModelData::CONSTANT);
    type = CONSTANT;
    conductivity = heatdiff.conductivity;
  }

  ~HeatDiffusionFcnConstant() {}

};
