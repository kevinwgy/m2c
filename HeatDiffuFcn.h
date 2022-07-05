#include<IoData.h>
#include<VarFcnBase.h>
#include<cassert>

/******************************************************************
 * Class HeatDiffuFcnBase and the derived classes are responsible
 * for calculating heat diffusion fluxes.
 * Each material has a seperate HeatDiffuFcn
 ******************************************************************/

//-----------------------------------------------------------------------
//
class HeatDiffuFcnBase{

protected:

   VarFcnBase &vf;

public:

   enum Type {NONE = 0, CONSTANT = 1} type;

   double conduct;

   HeatDiffuFcnBase(VarFcnBase &vf_) : vf(vf_), type(NONE), conduct(0.0) {}
   virtual ~HeatDiffuFcnBase() {}   

};

//--------------------------------------------------------------------------

class HeatDiffuFcnConstant : public HeatDiffuFcnBase{


public:
  
  HeatDiffuFcnConstant(HeatDiffusionModelData &heatdiff, VarFcnBase &vf_) : HeatDiffuFcnBase(vf_) {
    assert(heatdiff.type == HeatDiffusionModelData::CONSTANT);
    type = CONSTANT;
    conduct = heatdiff.diffusivity;
  }
  ~HeatDiffuFcnConstant() {}

};
