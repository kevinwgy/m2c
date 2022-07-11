#ifndef _VAR_FCN_ANEOS_BASE_H_
#define _VAR_FCN_ANEOS_BASE_H_

#include <VarFcnBase.h>

/****************************************************************************
 * This class is the base class for a class of (complete) EOS, known as ANEOS,
 * where "AN" represents "Analytic". See KW's notes. This type of EOS was 
 * developed by Sandia National Lab, primarily for solid materials. The main
 * idea is to model the Helmholtz free energy, more specifically, as a summation
 * of multiple contributions, then model each term "analytically".
 * This type of EOS is based on defining specific Helmholtz free energy as a function
 * of density (rho) and temperature (T), i.e. F(rho,T). (Note: This is specific
 * Helmholtz free energy, denoted as \bar{\Psi} in KW's notes.)
 ***************************************************************************/

class VarFcnANEOSBase : public VarFcnBase {

public:
  
 
  //----- EOS-Specific Functions -----//
  //! get pressure from density (rho) and internal energy per unit mass (e)
  double GetPressure(double rho, double e) const{
    double T = GetTemperature(rho, e);
    return GetPressure(rho,T);
  }

  //! get e (internal energy per unit mass) from density (rho) and pressure (p)
  double GetInternalEnergyPerUnitMass(double rho, double p) const{
    fprintf(stderr,"\033[0;31m*** Error:  GetInternalEnergyPerUnitMass Function not defined\n\033[0m");
    exit(-1); return 0.0;}

  //! get e - e0 from density (rho) and pressure (p)
  virtual double GetReferenceInternalEnergyPerUnitMass() const{
    fprintf(stderr,"\033[0;31m*** Error:  GetReferenceInternalEnergyPerUnitMass Function not defined\n\033[0m");
    exit(-1); return 0.0;}

  //! get rho (density) from p (pressure) and p (internal energy per unit mass)
  virtual double GetDensity(double p, double e) const{
    fprintf(stderr,"\033[0;31m*** Error:  GetDensity Function not defined\n\033[0m");
    exit(-1); return 0.0;}

  //! dpdrho = \frac{\partial p(\rho,e)}{\partial \rho}
  virtual double GetDpdrho(double rho, double e) const{
    fprintf(stderr,"\033[0;31m*** Error:  GetDpdrho Function not defined\n\033[0m");
    exit(-1); return 0.0;}

  //! BigGamma = 1/rho*(\frac{\partial p(\rho,e)}{\partial e})
  //  It is called "BigGamma" to distinguish it from the small "gamma" in perfect and stiffened EOS.
  virtual double GetBigGamma(double rho, double e) const{
    fprintf(stderr,"\033[0;31m*** Error:  GetBigGamma Function not defined\n\033[0m");
    exit(-1); return 0.0;}

  //! temperature law, defined separately for each EOS
  virtual double GetTemperature(double rho, double e) const{
    fprintf(stderr,"\033[0;31m*** Error:  GetTemperature Function not defined\n\033[0m");
    exit(-1); return 0.0;}

  //! temperature law, defined separately for each EOS
  virtual double GetReferenceTemperature() const{
    fprintf(stderr,"\033[0;31m*** Error:  GetReferenceTemperature Function not defined\n\033[0m");
    exit(-1); return 0.0;}

  //! temperature law, defined separately for each EOS
  virtual double GetInternalEnergyPerUnitMassFromTemperature(double rho, double T) const{
    fprintf(stderr,"\033[0;31m*** Error:  GetInternalEnergyPerUnitMassFromTemperature Function not defined\n\033[0m");
    exit(-1); return 0.0;}

  //! calculate e from rho and h
  virtual double GetInternalEnergyPerUnitMassFromEnthalpy(double rho, double h) const{
    fprintf(stderr,"\033[0;31m*** Error:  GetInternalEnergyPerUnitMassFromEnthalpy Function not defined\n\033[0m");
    exit(-1); return 0.0;}

  //checks that the Euler equations are still hyperbolic
  virtual bool CheckState(double rho, double p, bool silence = false) const{
    if(!std::isfinite(rho) || !std::isfinite(p)) {
      if(!silence)
        fprintf(stderr, "*** Error: CheckState failed. rho = %e, p = %e.\n\033[0m", rho, p);
      return true;
    }
    if(rho <= 0.0) {
      if(!silence && verbose>1)
        fprintf(stdout, "Warning: Negative density or violation of hyperbolicity. rho = %e, p = %e.\n", rho, p);
      return true;
    }
    double e = GetInternalEnergyPerUnitMass(rho,p);
    double c2 = GetDpdrho(rho, e) + p/rho*GetBigGamma(rho, e);
    if(c2<=0){
      if(!silence && verbose>1)
        fprintf(stdout, "Warning: Negative density or violation of hyperbolicity. rho = %e, p = %e.\n", rho, p);
      return true;
    }
    return false;
  }

  //checks that the Euler equations are still hyperbolic
  virtual bool CheckState(double *V, bool silence = false) const{
    if(!std::isfinite(V[0]) || !std::isfinite(V[1]) || !std::isfinite(V[2]) || !std::isfinite(V[3]) || !std::isfinite(V[4])){
      if(!silence)
        fprintf(stderr, "\033[0;31m*** Error: CheckState failed. V = %e %e %e %e %e\n\033[0m", V[0], V[1], V[2], V[3], V[4]);
      return true;
    }
    return CheckState(V[0], V[4]); 
  }
 
  //check for phase transitions
  virtual bool CheckPhaseTransition(int id/*id of the other phase*/) const{
    return false; //by default, phase transition is not allowed/considered
  }

  //----- Transformation Operators -----//
  virtual void ConservativeToPrimitive(double *U, double *V); 
  virtual void PrimitiveToConservative(double *V, double *U);

  //----- General Functions -----//
  inline int GetType() const{ return type; }

  virtual double ComputeSoundSpeed(double rho, double e);
  virtual double ComputeSoundSpeedSquare(double rho, double e); //!< this one does not crash on negative c^2
  virtual double ComputeMachNumber(double *V);
  virtual double ComputeEnthalpyPerUnitMass(double rho, double p); //!< h = e + p/rho
  virtual double ComputeTotalEnthalpyPerUnitMass(double *V); //!< H = 1/rho*(E + p)

  // Clipping
  virtual bool ClipDensityAndPressure(double *V, double *U = 0);
};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

inline
void VarFcnBase::ConservativeToPrimitive(double *U, double *V)
{
  V[0] = U[0];

  double invRho = 1.0 / U[0];

  V[1] = U[1] * invRho;
  V[2] = U[2] * invRho;
  V[3] = U[3] * invRho;

  double e = (U[4] - 0.5*V[0]*(V[1]*V[1]+V[2]*V[2]+V[3]*V[3])) * invRho;
  V[4] = GetPressure(V[0], e);
}

//------------------------------------------------------------------------------

inline
void VarFcnBase::PrimitiveToConservative(double *V, double *U)
{
  U[0] = V[0];

  U[1] = V[0] * V[1];
  U[2] = V[0] * V[2];
  U[3] = V[0] * V[3];

  double e = GetInternalEnergyPerUnitMass(V[0],V[4]); //pass in rho and p
  U[4] = V[0]*(e + 0.5*(V[1]*V[1]+V[2]*V[2]+V[3]*V[3]));
}

//------------------------------------------------------------------------------

inline
double VarFcnBase::ComputeSoundSpeed(double rho, double e)
{
  double c2 = GetDpdrho(rho, e) + GetPressure(rho,e)/rho*GetBigGamma(rho, e);
  if(c2<=0) {
    fprintf(stderr,"\033[0;31m*** Error: Cannot calculate speed of sound (Square-root of a negative number): rho = %e, e = %e.\n\033[0m",
            rho, e);
    exit(-1);
  }
  return sqrt(c2);
}

//------------------------------------------------------------------------------

inline
double VarFcnBase::ComputeSoundSpeedSquare(double rho, double e)
{
  return GetDpdrho(rho, e) + GetPressure(rho,e)/rho*GetBigGamma(rho, e);
}

//------------------------------------------------------------------------------

inline
double VarFcnBase::ComputeMachNumber(double *V)
{
  double e = GetInternalEnergyPerUnitMass(V[0],V[4]); 
  double c = ComputeSoundSpeedSquare(V[0], e);

  if(c<0) {
    fprintf(stderr,"\033[0;31m*** Error: c^2 (square of sound speed) = %e in ComputeMachNumber. V = %e, %e, %e, %e, %e.\n\033[0m",
            c, V[0], V[1], V[2], V[3], V[4]);
    exit(-1);
  } else
    c = sqrt(c);

  return sqrt(V[1]*V[1]+V[2]*V[2]+V[3]*V[3])/c;
}

//------------------------------------------------------------------------------

inline
double VarFcnBase::ComputeEnthalpyPerUnitMass(double rho, double p)
{
  return GetInternalEnergyPerUnitMass(rho,p) + p/rho;
}

//------------------------------------------------------------------------------

inline
double VarFcnBase::ComputeTotalEnthalpyPerUnitMass(double *V) //!< H = 1/rho*(E + p)
{
  double e = GetInternalEnergyPerUnitMass(V[0],V[4]); 
  return e + 0.5*(V[1]*V[1]+V[2]*V[2]+V[3]*V[3]) + V[4]/V[0];
}

//------------------------------------------------------------------------------

inline
bool VarFcnBase::ClipDensityAndPressure(double *V, double *U)
{
//verification of density and pressure value
//if pressure/density < pmin/rhomin, set pressure/density to pmin/rhomin
//and rewrite V and U!!
  bool clip = false;

  if(V[0]<rhomin){
//    if(verbose)
//      fprintf(stderr,"clip density from %e to %e.\n", V[0], rhomin);
    V[0] = rhomin;
    clip = true;
  }

  if(V[4]<pmin){
//    if(verbose)
//      fprintf(stdout, "clip pressure from %e to %e\n", V[4], pmin);
    V[4] = pmin;
    clip = true;
  }

  if(V[0]>rhomax) {
    V[0] = rhomax;
    clip = true;
  }

  if(V[4]>pmax) {
    V[4] = pmax;
    clip = true;
  }

  if(clip && U) //also modify U
    PrimitiveToConservative(V,U);

  return clip;
}

//------------------------------------------------------------------------------



#endif
