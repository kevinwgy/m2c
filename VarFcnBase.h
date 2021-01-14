#ifndef _VAR_FCN_BASE_H_
#define _VAR_FCN_BASE_H_

#include <IoData.h>
#include <Vector3D.h>
#include <cmath>
#include <iostream>

/****************************************************************************
 * This class is the base class for the VarFcnEOS classes where EOS can be
 * a stiffened gas, a Tait EOS or a JWL EOS.
 * Only elementary functions are declared and/or defined here.
 * All arguments must be pertinent to only a single grid node or a single
 * state, since it is assumed that the EOS that must be used at this point 
 * is known.
 * As a base class, it is also assumed that the Euler equations are 
 * considered, meaning there are a priori only 5 variables that
 * define each state. Thus for flows with more than 5 variables,
 * some transformation operators must be overwritten in
 * the appropriate class.
 *
 * lay-out of the base class is:
 *  - 1 -  Transformation Operators
 *  - 2 -  General Functions
 *  - 3 -  Equations of State Parameters
 *  - 4 -  EOS related functions
 ***************************************************************************/
class VarFcnBase {

public:
  
  enum Type{ STIFFENEDGAS = 0, MIEGRUNEISEN = 1} type;

  double rhomin,pmin;

  VarFcnBase(FluidModelData &data) {
    rhomin = data.rhomin;
    pmin = data.pmin;
  }
  virtual ~VarFcnBase() {}
 
  //----- EOS-Specific Functions -----//
  //! get pressure from density (rho) and internal energy per unit mass (e)
  virtual double GetPressure(double rho, double e) const{
    fprintf(stderr, "*** Error:  GetPressure Function not defined\n");
    exit(1); }

  //! get e (internal energy per unit mass) from density (rho) and pressure (p)
  virtual double GetInternalEnergyPerUnitMass(double rho, double p) const{
    fprintf(stderr, "*** Error:  GetInternalEnergyPerUnitMass Function not defined\n");
    exit(1); }

  //! get rho (density) from p (pressure) and p (internal energy per unit mass)
  virtual double GetDensity(double p, double e) const{
    fprintf(stderr, "*** Error:  GetDensity Function not defined\n");
    exit(1); }

  //! dpdrho = \frac{\partial p(\rho,e)}{\partial \rho}
  virtual double GetDpdrho(double rho, double e) const{
    fprintf(stderr, "*** Error:  GetDpdrho Function not defined\n");
    exit(1); }

  //! BigGamma = 1/rho*(\frac{\partial p(\rho,e)}{\partial e})
  //  It is called "BigGamma" to distinguish it from the small "gamma" in perfect and stiffened EOS.
  virtual double GetBigGamma(double rho, double e) const{
    fprintf(stderr, "*** Error:  GetBigGamma Function not defined\n");
    exit(1); }

  //checks that the Euler equations are still hyperbolic
  virtual bool CheckState(double *V) const{
    fprintf(stderr, "*** Error:  VarFcnBase::checkState not implemented for this equation of state (type = %d)\n", type);
    return false;
  }

 
  //----- Transformation Operators -----//
  inline void ConservativeToPrimitive(double *U, double *V); 
  inline void PrimitiveToConservative(double *V, double *U);

  //----- General Functions -----//
  inline int GetType() const{ return type; }

  inline double ComputeSoundSpeed(double rho, double e);
  inline double ComputeMachNumber(double *V);

  // Clipping
  inline bool ClipDensityAndPressure(double *V, double *U = 0);
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
    fprintf(stderr,"Error: Cannot calculate speed of sound (Square-root of a negative number): rho = %e, e = %e.\n",
            rho, e);
    exit(-1);
  }
  return sqrt(c2);
}

//------------------------------------------------------------------------------

inline 
double VarFcnBase::ComputeMachNumber(double *V)
{
  double e = GetInternalEnergyPerUnitMass(V[0],V[4]); 
  double c = ComputeSoundSpeed(V[0], e);
  return sqrt(V[1]*V[1]+V[2]*V[2]+V[3]*V[3])/c;
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
    fprintf(stderr,"clip density from %e to %e.\n", V[0], rhomin);
    V[0] = rhomin;
    clip = true;
  }

  if(V[4]<pmin){
    fprintf(stdout, "clip pressure from %e to %e\n", V[4], pmin);
    V[4] = pmin;
    clip = true;
  }

  if(clip && U) //also modify U
    PrimitiveToConservative(V,U);

  return clip;
}

//------------------------------------------------------------------------------



#endif
