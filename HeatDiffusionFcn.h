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
  double diffusivity;

public:

  HeatDiffusionFcnBase() : type(NONE), diffusivity(0.0) {}
  virtual ~HeatDiffusionFcnBase() {}   

  inline int GetDiffusionFunctionType() {return (int)type;}
  inline double GetDiffusivity() {return diffusivity;}

};

//--------------------------------------------------------------------------

class HeatDiffusionFcnConstant : public HeatDiffusionFcnBase {

public:
  
  HeatDiffusionFcnConstant(HeatDiffusionModelData &heatdiff) : HeatDiffusionFcnBase() {
    assert(heatdiff.type == HeatDiffusionModelData::CONSTANT);
    type = CONSTANT;
    diffusivity = heatdiff.diffusivity;
  }

  ~HeatDiffusionFcnConstant() {}

};
