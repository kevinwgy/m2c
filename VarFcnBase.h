#ifndef _VAR_FCN_BASE_H_
#define _VAR_FCN_BASE_H_

#include <IoData.h>
#include <Vector3D.h>
#include <cmath>
#include <complex>
typedef std::complex<double> bcomp;
#include <iostream>

//--------------------------------------------------------------------------
// This class is the base class for the VarFcnEOS classes where EOS can be
// a stiffened gas, a Tait EOS or a JWL EOS.
// Only elementary functions are declared and/or defined here.
// All arguments must be pertinent to only a single grid node or a single
// state, since it is assumed that the EOS that must be used at this point 
// is known.
// As a base class, it is also assumed that the Euler equations are 
// considered, meaning there are a priori only 5 variables that
// define each state. Thus for flows with more than 5 variables,
// some transformation operators must be overwritten in
// the appropriate class.
//
// lay-out of the base class is:
//  - 1 -  Transformation Operators
//  - 2 -  General Functions
//  - 3 -  Equations of State Parameters
//  - 4 -  EOS related functions
//--------------------------------------------------------------------------
class VarFcnBase {

public:
  
  const char** pname;
  enum Type{ PERFECTGAS = 0, STIFFENEDGAS = 1, TAIT = 2, JWL = 3} type;

  double rhomin,pmin;
  bool verif_clipping; //KW(05/09/2011): This seems to be a flag for screen output...

  VarFcnBase(FluidModelData &data) {
    rhomin = data.rhomin;
    pmin = data.pmin;
    verif_clipping = false; //true;
  }
  virtual ~VarFcnBase() {}
  
  virtual bool equal(VarFcnBase* oth) { return false; }

  //----- Transformation Operators -----//
  virtual void conservativeToPrimitive(double *U, double *V) = 0; 
  virtual void primitiveToConservative(double *V, double *U) = 0;
  virtual int  conservativeToPrimitiveVerification(int glob, double *U, double *V){
    conservativeToPrimitive(U, V);
    return verification(glob, U, V);
  }
  virtual int  verification(int glob, double *U, double *V) { return 0; }
  int verification_negative_rho_p(double*V, double rho_min, double p_min) {return (V[0] <= rho_min || V[4] <= p_min); }
  virtual void computedVdU(double *V, double *dVdU) = 0;
  virtual void computedUdV(double *V, double *dUdV) = 0;

  virtual void computeFofU(double n[3], double *U, double *F) {}
  virtual void computeFofV(double n[3], double *V, double *F) {}
  virtual void computedFdU(double n[3], double *V, double *dFdU) {}
  virtual void computedFdV(double n[3], double *V, double *dFdV) {}
  
  virtual void multiplyBydVdU(double *V, double *vec, double *res);
  virtual void multiplyBydVdU(double *V, bcomp *vec, bcomp *res);
  virtual void multiplyBydVdUT(double *V, double *vec, double *res);
  virtual void multiplyBydVdUT(double *V, bcomp *vec, bcomp *res);
  virtual void postMultiplyBydUdV(double *V, double *mat, double *res);
  virtual void postMultiplyBydVdU(double *V, double *mat, double *res);
  virtual void postMultiplyBydUdV(double *V, bcomp *mat, bcomp *res);
  virtual void preMultiplyBydUdV(double *V, double *mat, double *res);


  virtual void extrapolatePrimitive(double un, double c, double *Vb, double *Vinter, double *V) {
    fprintf(stderr, "*** Error:  extrapolatePrimitive Function not defined\n");
    exit(1); }
  virtual void extrapolateCharacteristic(double n[3], double un, double c, double *Vb, double *dV) {
    fprintf(stderr, "*** Error:  extrapolateCharacteristic Function not defined\n");
    exit(1); }
  virtual void primitiveToCharacteristicVariations(double n[3], double *V, double *dV, double *dW){
    fprintf(stderr, "*** Error:  primitiveToCharacteristic Function not defined\n");
    exit(1); }
  virtual void characteristicToPrimitiveVariations(double n[3], double *V, double *dW, double *dV){
    fprintf(stderr, "*** Error:  characteristicToPrimitive Function not defined\n");
    exit(-1);}


  virtual void primitiveToConservativeDerivative(double *V, double *dV, double *U, double *dU) {
    fprintf(stderr, "*** Error:  primitiveToConservativeDerivative Function not defined\n");
    exit(1); }
  virtual void conservativeToPrimitiveDerivative(double *U, double *dU, double *V, double *dV) {
    fprintf(stderr, "*** Error:  conservativeToPrimitiveDerivative Function not defined\n");
    exit(1); }
  virtual void computeConservativeToPrimitiveDerivativeOperators(double*, double*,
                                                                 double dVdU[5][5], double dVdPstiff[5]) {
    fprintf(stderr, "*** Error:  computeComputeConservativeToPrimitiveDerivativeOperators5 Function not defined\n");
    exit(1); }


  //----- General Functions -----//
  virtual int getType() const{ return type; }
  virtual bool doVerification() const{ return (rhomin!=-1.e9||pmin!=-1.e9); }

  virtual double getDensity(double *V)   const{ return V[0]; }
  virtual Vec3D  getVelocity(double *V)  const{ return Vec3D(V[1], V[2], V[3]); }
  virtual double getVelocityX(double *V) const{ return V[1]; }
  virtual double getVelocityY(double *V) const{ return V[2]; }
  virtual double getVelocityZ(double *V) const{ return V[3]; }
  virtual double getPressure(double *V)  const{ return V[4]; }
  virtual void computedPdV(double *dPdV)  const{ 
    for(int i=0; i<5; ++i) dPdV[i] = 0.0;
    dPdV[4] = 1.0;    
  }

  virtual void setDensity(const double density, double *V) { V[0] = density; }
  virtual void setPressure(const double p, double *V)  { V[4] = p;   }
  virtual void setDensity(double *V, double *Vorig) { V[0] = Vorig[0]; }
  virtual void setPressure(double *V, double *Vorig){ V[4] = Vorig[4]; }

  //checks that the Euler equations are still hyperbolic
  virtual double checkPressure(double *V) const{
    fprintf(stderr, "*** Error:  checkPressure Function not defined\n");
    exit(1); }
  virtual bool checkReconstructedValues(double *V, int nodeNum, int otherNodeNum, int phi, int otherPhi, int failsafe) const{
    fprintf(stderr, "*** Error:  VarFcnBase::checkReconstructedValues not implemented for this equation of state (integer tag = %d)\n", type);
    return false;
  }

  virtual double computeTemperature(double *V) const{
    fprintf(stderr, "*** Error:  computeTemperature Function not defined\n");
    exit(1); }
  virtual void computeTemperatureGradient(double *V,double* Tg) const{
    fprintf(stderr, "*** Error:  computeTemperatureGradient Function not defined\n");
    exit(1); }
  virtual void computeTemperatureHessian(double *V,double& Trr, double& Trp, 
                                 double& Tpp) const { 
    fprintf(stderr, "*** Error:  computeTemperatureHessian Function not defined\n");
    exit(1); }
  virtual void getV4FromTemperature(double *V, double T) const{
    fprintf(stderr, "*** Error:  getV4FromTemperature Function not defined\n");
    exit(1); }
  virtual void getdV4FromdTemperature(double *V, double T) const{
    fprintf(stderr, "*** Error:  getdV4FromdTemperature Function not defined\n");
    exit(1); }
  virtual double computeRhoEnergy(double *V)   const{
    fprintf(stderr, "*** Error:  computeRhoEnergy Function not defined\n");
    exit(1); }
  virtual double computeRhoEpsilon(double *V)  const{ //this function computes the internal energy (=rho*e-0.5*rho*u^2)
    fprintf(stderr, "*** Error:  computeRhoEpsilon Function not defined\n");
    exit(1); }
  virtual double getVelocitySquare(double *V)  const{ return V[1]*V[1]+V[2]*V[2]+V[3]*V[3]; }
  virtual double getVelocityNorm(double *V)    const{ return sqrt(getVelocitySquare(V)); }

  virtual double computeMachNumber(double *V)  const{ return getVelocityNorm(V)/computeSoundSpeed(V); }
  virtual double computeSoundSpeed(double *V)  const{ 
    fprintf(stderr, "*** Error:  computeSoundSpeed Function not defined\n");
    exit(1); }
  virtual double computeSoundSpeed(const double density, const double entropy) const{
    fprintf(stderr, "*** Error:  computeSoundSpeed(entropy) Function not defined\n");
    exit(1); }
  virtual double computeEntropy(const double density, const double pressure) const{
    fprintf(stderr, "*** Error:  computeEntropy Function not defined\n");
    exit(1); }
  virtual double computeIsentropicPressure(const double entropy, const double density) const{
    fprintf(stderr, "*** Error:  computeIsentropicEntropy Function not defined\n");
    exit(1); }
  virtual double computeTotalPressure(double machr, double* V) const{
    fprintf(stderr, "*** Error:  computeTotalPressure Function not defined\n");
    exit(1); }

  //virtual void rstVar(IoData &iod) {}
  //virtual void rV(IoData &iod) { pmin = iod.eqs.fluidModel.pmin; }
  virtual Vec3D getDerivativeOfVelocity(double *dV) const { return Vec3D(dV[1], dV[2], dV[3]); }
  virtual double computeDerivativeOfTemperature(double *V, double *dV) const { return 0.0; }
  virtual void computeDerivativeOperatorsOfTemperature(double *V, double dTdV[5]) const {
    fprintf(stderr, "*** Error:  computeDerivativeOperatorsOfTemperature Function not defined\n");
    exit(-1);  }
  virtual double computeDerivativeOfMachNumber(double *V, double *dV, double dMach) const { return 0.0; }
  virtual double computeDerivativeOfSoundSpeed(double *V, double *dV, double dMach) const { return 0.0; }
  virtual double computeDerivativeOfTotalPressure(double machr, double dmachr, double* V, double* dV, double dmach) const { return 0.0; }
  virtual double getDerivativeOfPressureConstant() const { return 0.0; }
  virtual double getDerivativeOfVelocityNorm(double *V, double *dV) const { 
    double vnorm=sqrt(V[1]*V[1]+V[2]*V[2]+V[3]*V[3]);
    if (vnorm < 100*std::numeric_limits<float>::min())//TODO this will be incorrect on stationary points and for shape sensitivity
      return 0.0;
    else
      return (V[1]*dV[1]+V[2]*dV[2]+V[3]*dV[3])/vnorm;
  }

  virtual double specificHeatCstPressure() const{ 
    fprintf(stderr, "*** Error: specificHeatCstPressure not defined\n");
    exit(1);
  }
  virtual double computePressureCoefficient(double *V, double pinfty, 
                                            double mach, bool dimFlag) const{ 
    fprintf(stderr, "*** Error: computePressureCoefficient not defined\n");
    exit(1);
  }


  //----- Equation of State Parameters -----//
  /* PERFECT/STIFFENED GAS EOS */
  virtual double getGamma() const{
    fprintf(stderr, "*** Error:  getGamma Function not defined\n");
    exit(1); }
  virtual double getGamma1() const{
    fprintf(stderr, "*** Error:  getGamma1 Function not defined\n");
    exit(1); }
  virtual double getPressureConstant() const{
    fprintf(stderr, "*** Error:  getPressureConstant Function not defined\n");
    exit(1); }
  /* TAIT EOS */
  virtual double getCv() const{
    fprintf(stderr, "*** Error:  getCv Function not defined\n");
    exit(1); }
  virtual double getAlphaWater() const{
    fprintf(stderr, "*** Error: getAlphaWater Function not defined\n");
    exit(1);}
  virtual double getBetaWater() const{
    fprintf(stderr, "*** Error: getBetaWater Function not defined\n");
    exit(1);}
  virtual double getPrefWater() const{
    fprintf(stderr, "*** Error: getPrefWater Function not defined\n");
    exit(1);}
  /* JWL EOS */
  virtual double getOmega() const{
    fprintf(stderr, "*** Error: getOmega Function not defined\n");
    exit(1);}
  virtual double getOmegap1() const{
    fprintf(stderr, "*** Error: getOmegap1 Function not defined\n");
    exit(1);}
  virtual double getInvOmega() const{
    fprintf(stderr, "*** Error: getInvOmega Function not defined\n");
    exit(1);}
  virtual double getInvOmegap1() const{
    fprintf(stderr, "*** Error: getInvOmegap1 Function not defined\n");
    exit(1);}
  virtual double getA1() const{
    fprintf(stderr, "*** Error: getA1 Function not defined\n");
    exit(1);}
  virtual double getA2() const{
    fprintf(stderr, "*** Error: getA2 Function not defined\n");
    exit(1);}
  virtual double getR1() const{
    fprintf(stderr, "*** Error: getR1 Function not defined\n");
    exit(1);}
  virtual double getR2() const{
    fprintf(stderr, "*** Error: getR2 Function not defined\n");
    exit(1);}
  virtual double getRhoref() const{
    fprintf(stderr, "*** Error: getRhoref Function not defined\n");
    exit(1);}
  virtual double getR1r() const{
    fprintf(stderr, "*** Error: getR1r Function not defined\n");
    exit(1);}
  virtual double getR2r() const{
    fprintf(stderr, "*** Error: getR2r Function not defined\n");
    exit(1);}

  //----- EOS related functions -----//
  virtual double computeExponentials(const double density) const{
    fprintf(stderr, "*** Error: computeExponentials Function not defined\n");
    exit(1);}
  virtual double computeDerivativeOfExponentials(const double density) const{
    fprintf(stderr, "*** Error: computeDerivativeOfExponentials Function not defined\n");
    exit(1);}
  virtual double computeExponentials2(const double density) const{
    fprintf(stderr, "*** Error: computeExponentials2 Function not defined\n");
    exit(1);}
  virtual double computeDerivativeOfExponentials2(const double density) const{
    fprintf(stderr, "*** Error: computeDerivativeOfExponentials2 Function not defined\n");
    exit(1);}
  virtual double computeFrho(double *V) const{
    fprintf(stderr, "*** Error: computeFrho Function not defined\n");
    exit(1);}
  virtual double computeFrho(const double density) const{
    fprintf(stderr, "*** Error: computeFrho Function not defined\n");
    exit(1);}
  virtual double computeFrhop(double *V) const{
    fprintf(stderr, "*** Error: computeFrhop Function not defined\n");
    exit(1);}
  virtual double computeFrhop(const double density) const{
    fprintf(stderr, "*** Error: computeFrhop Function not defined\n");
    exit(1);}

};
//------------------------------------------------------------------------------
inline
void VarFcnBase::multiplyBydVdU(double *V, double *vec, double *res) {

  double dVdU[25];
  computedVdU(V, dVdU);

  res[0] = dVdU[0]*vec[0];
  res[1] = dVdU[5]*vec[0]+dVdU[6]*vec[1];
  res[2] = dVdU[10]*vec[0]+dVdU[12]*vec[2];
  res[3] = dVdU[15]*vec[0]+dVdU[18]*vec[3];
  res[4] = dVdU[20]*vec[0]+dVdU[21]*vec[1]+dVdU[22]*vec[2]+dVdU[23]*vec[3]+dVdU[24]*vec[4];

}
//------------------------------------------------------------------------------
inline
void VarFcnBase::multiplyBydVdU(double *V, bcomp *vec, bcomp *res) {

  double dVdU[25];
  computedVdU(V, dVdU);

  res[0] = dVdU[0]*vec[0];
  res[1] = dVdU[5]*vec[0]+dVdU[6]*vec[1];
  res[2] = dVdU[10]*vec[0]+dVdU[12]*vec[2];
  res[3] = dVdU[15]*vec[0]+dVdU[18]*vec[3];
  res[4] = dVdU[20]*vec[0]+dVdU[21]*vec[1]+dVdU[22]*vec[2]+dVdU[23]*vec[3]+dVdU[24]*vec[4];

}

//------------------------------------------------------------------------------
inline
void VarFcnBase::multiplyBydVdUT(double *V, double *vec, double *res) {

  double dVdU[25];
  computedVdU(V, dVdU);

  res[0] = dVdU[0]*vec[0] + dVdU[5]*vec[1] + dVdU[10]*vec[2] + dVdU[15]*vec[3] + dVdU[20]*vec[4];
  res[1] = dVdU[6]*vec[1] + dVdU[21]*vec[4];
  res[2] = dVdU[12]*vec[2] + dVdU[22]*vec[4];
  res[3] = dVdU[18]*vec[3] + dVdU[23]*vec[4];
  res[4] = dVdU[24]*vec[4];

}
//------------------------------------------------------------------------------
inline
void VarFcnBase::multiplyBydVdUT(double *V, bcomp *vec, bcomp *res) {

  double dVdU[25];
  computedVdU(V, dVdU);

  res[0] = dVdU[0]*vec[0] + dVdU[5]*vec[1] + dVdU[10]*vec[2] + dVdU[15]*vec[3] + dVdU[20]*vec[4];
  res[1] = dVdU[6]*vec[1] + dVdU[21]*vec[4];
  res[2] = dVdU[12]*vec[2] + dVdU[22]*vec[4];
  res[3] = dVdU[18]*vec[3] + dVdU[23]*vec[4];
  res[4] = dVdU[24]*vec[4];

}
//------------------------------------------------------------------------------
inline
void VarFcnBase::postMultiplyBydUdV(double *V, double *mat, double *res) {

  double dUdV[25];
  computedUdV(V, dUdV);

  res[0] = mat[0]*dUdV[0]+mat[1]*dUdV[5]+mat[2]*dUdV[10]+mat[3]*dUdV[15]+mat[4]*dUdV[20];
  res[1] = mat[1]*dUdV[6]+mat[4]*dUdV[21];
  res[2] = mat[2]*dUdV[12]+mat[4]*dUdV[22];
  res[3] = mat[3]*dUdV[18]+mat[4]*dUdV[23];
  res[4] = mat[4]*dUdV[24];
  res[5] = mat[5]*dUdV[0]+mat[6]*dUdV[5]+mat[7]*dUdV[10]+mat[8]*dUdV[15]+mat[9]*dUdV[20];
  res[6] = mat[6]*dUdV[6]+mat[9]*dUdV[21];
  res[7] = mat[7]*dUdV[12]+mat[9]*dUdV[22];
  res[8] = mat[8]*dUdV[18]+mat[9]*dUdV[23];
  res[9] = mat[9]*dUdV[24];
  res[10] = mat[10]*dUdV[0]+mat[11]*dUdV[5]+mat[12]*dUdV[10]+mat[13]*dUdV[15]+mat[14]*dUdV[20];  res[11] = mat[11]*dUdV[6]+mat[14]*dUdV[21];
  res[12] = mat[12]*dUdV[12]+mat[14]*dUdV[22];
  res[13] = mat[13]*dUdV[18]+mat[14]*dUdV[23];
  res[14] = mat[14]*dUdV[24];
  res[15] = mat[15]*dUdV[0]+mat[16]*dUdV[5]+mat[17]*dUdV[10]+mat[18]*dUdV[15]+mat[19]*dUdV[20];  res[16] = mat[16]*dUdV[6]+mat[19]*dUdV[21];
  res[17] = mat[17]*dUdV[12]+mat[19]*dUdV[22];
  res[18] = mat[18]*dUdV[18]+mat[19]*dUdV[23];
  res[19] = mat[19]*dUdV[24];
  res[20] = mat[20]*dUdV[0]+mat[21]*dUdV[5]+mat[22]*dUdV[10]+mat[23]*dUdV[15]+mat[24]*dUdV[20];  res[21] = mat[21]*dUdV[6]+mat[24]*dUdV[21];
  res[22] = mat[22]*dUdV[12]+mat[24]*dUdV[22];
  res[23] = mat[23]*dUdV[18]+mat[24]*dUdV[23];
  res[24] = mat[24]*dUdV[24];

}
//------------------------------------------------------------------------------
inline
void VarFcnBase::postMultiplyBydVdU(double *V, double *mat, double *res) {

  double dVdU[25];
  computedVdU(V, dVdU);

  res[0] = mat[0]*dVdU[0]+mat[1]*dVdU[5]+mat[2]*dVdU[10]+mat[3]*dVdU[15]+mat[4]*dVdU[20];
  res[1] = mat[1]*dVdU[6]+mat[4]*dVdU[21];
  res[2] = mat[2]*dVdU[12]+mat[4]*dVdU[22];
  res[3] = mat[3]*dVdU[18]+mat[4]*dVdU[23];
  res[4] = mat[4]*dVdU[24];
  res[5] = mat[5]*dVdU[0]+mat[6]*dVdU[5]+mat[7]*dVdU[10]+mat[8]*dVdU[15]+mat[9]*dVdU[20];
  res[6] = mat[6]*dVdU[6]+mat[9]*dVdU[21];
  res[7] = mat[7]*dVdU[12]+mat[9]*dVdU[22];
  res[8] = mat[8]*dVdU[18]+mat[9]*dVdU[23];
  res[9] = mat[9]*dVdU[24];
  res[10] = mat[10]*dVdU[0]+mat[11]*dVdU[5]+mat[12]*dVdU[10]+mat[13]*dVdU[15]+mat[14]*dVdU[20];  res[11] = mat[11]*dVdU[6]+mat[14]*dVdU[21];
  res[12] = mat[12]*dVdU[12]+mat[14]*dVdU[22];
  res[13] = mat[13]*dVdU[18]+mat[14]*dVdU[23];
  res[14] = mat[14]*dVdU[24];
  res[15] = mat[15]*dVdU[0]+mat[16]*dVdU[5]+mat[17]*dVdU[10]+mat[18]*dVdU[15]+mat[19]*dVdU[20];  res[16] = mat[16]*dVdU[6]+mat[19]*dVdU[21];
  res[17] = mat[17]*dVdU[12]+mat[19]*dVdU[22];
  res[18] = mat[18]*dVdU[18]+mat[19]*dVdU[23];
  res[19] = mat[19]*dVdU[24];
  res[20] = mat[20]*dVdU[0]+mat[21]*dVdU[5]+mat[22]*dVdU[10]+mat[23]*dVdU[15]+mat[24]*dVdU[20];  res[21] = mat[21]*dVdU[6]+mat[24]*dVdU[21];
  res[22] = mat[22]*dVdU[12]+mat[24]*dVdU[22];
  res[23] = mat[23]*dVdU[18]+mat[24]*dVdU[23];
  res[24] = mat[24]*dVdU[24];

}
//------------------------------------------------------------------------------
inline
void VarFcnBase::postMultiplyBydUdV(double *V, bcomp *mat, bcomp *res) {

  double dUdV[25];
  computedUdV(V, dUdV);

  res[0] = mat[0]*dUdV[0]+mat[1]*dUdV[5]+mat[2]*dUdV[10]+mat[3]*dUdV[15]+mat[4]*dUdV[20];
  res[1] = mat[1]*dUdV[6]+mat[4]*dUdV[21];
  res[2] = mat[2]*dUdV[12]+mat[4]*dUdV[22];
  res[3] = mat[3]*dUdV[18]+mat[4]*dUdV[23];
  res[4] = mat[4]*dUdV[24];
  res[5] = mat[5]*dUdV[0]+mat[6]*dUdV[5]+mat[7]*dUdV[10]+mat[8]*dUdV[15]+mat[9]*dUdV[20];
  res[6] = mat[6]*dUdV[6]+mat[9]*dUdV[21];
  res[7] = mat[7]*dUdV[12]+mat[9]*dUdV[22];
  res[8] = mat[8]*dUdV[18]+mat[9]*dUdV[23];
  res[9] = mat[9]*dUdV[24];
  res[10] = mat[10]*dUdV[0]+mat[11]*dUdV[5]+mat[12]*dUdV[10]+mat[13]*dUdV[15]+mat[14]*dUdV[20];  res[11] = mat[11]*dUdV[6]+mat[14]*dUdV[21];
  res[12] = mat[12]*dUdV[12]+mat[14]*dUdV[22];
  res[13] = mat[13]*dUdV[18]+mat[14]*dUdV[23];
  res[14] = mat[14]*dUdV[24];
  res[15] = mat[15]*dUdV[0]+mat[16]*dUdV[5]+mat[17]*dUdV[10]+mat[18]*dUdV[15]+mat[19]*dUdV[20];  res[16] = mat[16]*dUdV[6]+mat[19]*dUdV[21];
  res[17] = mat[17]*dUdV[12]+mat[19]*dUdV[22];
  res[18] = mat[18]*dUdV[18]+mat[19]*dUdV[23];
  res[19] = mat[19]*dUdV[24];
  res[20] = mat[20]*dUdV[0]+mat[21]*dUdV[5]+mat[22]*dUdV[10]+mat[23]*dUdV[15]+mat[24]*dUdV[20];  res[21] = mat[21]*dUdV[6]+mat[24]*dUdV[21];
  res[22] = mat[22]*dUdV[12]+mat[24]*dUdV[22];
  res[23] = mat[23]*dUdV[18]+mat[24]*dUdV[23];
  res[24] = mat[24]*dUdV[24];

}
//------------------------------------------------------------------------------
inline
void VarFcnBase::preMultiplyBydUdV(double *V, double *mat, double *res) {

  double dUdV[25];
  computedUdV(V, dUdV);

  res[0] = dUdV[0]*mat[0];
  res[1] = dUdV[0]*mat[1];
  res[2] = dUdV[0]*mat[2];
  res[3] = dUdV[0]*mat[3];
  res[4] = dUdV[0]*mat[4];
  res[5] = dUdV[5]*mat[0]+dUdV[6]*mat[5];
  res[6] = dUdV[5]*mat[1]+dUdV[6]*mat[6];
  res[7] = dUdV[5]*mat[2]+dUdV[6]*mat[7];
  res[8] = dUdV[5]*mat[3]+dUdV[6]*mat[8];
  res[9] = dUdV[5]*mat[4]+dUdV[6]*mat[9];
  res[10] = dUdV[10]*mat[0]+dUdV[12]*mat[10];
  res[11] = dUdV[10]*mat[1]+dUdV[12]*mat[11];
  res[12] = dUdV[10]*mat[2]+dUdV[12]*mat[12];
  res[13] = dUdV[10]*mat[3]+dUdV[12]*mat[13];
  res[14] = dUdV[10]*mat[4]+dUdV[12]*mat[14];
  res[15] = dUdV[15]*mat[0]+dUdV[18]*mat[15];
  res[16] = dUdV[15]*mat[1]+dUdV[18]*mat[16];
  res[17] = dUdV[15]*mat[2]+dUdV[18]*mat[17];
  res[18] = dUdV[15]*mat[3]+dUdV[18]*mat[18];
  res[19] = dUdV[15]*mat[4]+dUdV[18]*mat[19];
  res[20] = dUdV[20]*mat[0]+dUdV[21]*mat[5]+dUdV[22]*mat[10]+dUdV[23]*mat[15]+dUdV[24]*mat[20];
  res[21] = dUdV[20]*mat[1]+dUdV[21]*mat[6]+dUdV[22]*mat[11]+dUdV[23]*mat[16]+dUdV[24]*mat[21];
  res[22] = dUdV[20]*mat[2]+dUdV[21]*mat[7]+dUdV[22]*mat[12]+dUdV[23]*mat[17]+dUdV[24]*mat[22];
  res[23] = dUdV[20]*mat[3]+dUdV[21]*mat[8]+dUdV[22]*mat[13]+dUdV[23]*mat[18]+dUdV[24]*mat[23];
  res[24] = dUdV[20]*mat[4]+dUdV[21]*mat[9]+dUdV[22]*mat[14]+dUdV[23]*mat[19]+dUdV[24]*mat[24];

}
//------------------------------------------------------------------------------

#endif
