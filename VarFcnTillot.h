/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _VAR_FCN_TILLOT_H_
#define _VAR_FCN_TILLOT_H_

#include <VarFcnBase.h>

/********************************************************************************
 * This class is the VarFcn class for the Tillotson equation of state (EOS)
 * Only elementary functions are declared and/or defined here.
 * All arguments must be pertinent to only a single grid node or a single
 * state.
 *
 * The Tillotson EOS presented in Brundage, 2013 is implemented here. For details,
 * see KW's notes.
 *
 *   p      : pressure
 *   e      : internal energy per unit mass.
 *   rho    : density
 *   eta    : rho/rho0
 *   mu     : eta - 1 
 *   omega  : rho0/rho - 1
 *   chi    : 1/(e/(e0*eta*eta)+1)
 *
 *   rho0   : density in the ambient state
 *   e0     : internal energy (per unit mass) in the ambient state
 *   a, b   : non-D model parameters. (a+b/4: Gruneisen param in ambient state; a+b: ... at low T (e=0))
 *   A, B   : model parameters. dim: [force]/[length]^2. (A: bulk modulus at p = 0 and e = 0)
 *   alpha  : non-D model parameters (in the formula for "hot expanded states")
 *   beta   : non-D model parameters (in the formula for "hot expanded states")
 *   rhoIV  : "incipient vaporization" density (constant model param.)
 *   eIV    : "incipient vaporization" internal energy (constant model param.)
 *   eCV    : "complete vaporization" internal energy (constant model param.) (eCV-eIV: latent heat per unit mass)
 *
 *   cv     : specific heat at constant volume
 *   T0     : temperature at rho0 and e0
 *   temperature_depends_on_density: whether T depends on both rho and e, or just e.
 *
 * References: KW's notes 
 ********************************************************************************/

class VarFcnTillot : public VarFcnTillot {

private:

  double rho0, e0, a, b, A, B, alpha, beta;
  double rhoIV, eIV, eCV;
  double cv, T0;
  bool temperature_depends_on_density;


public:

  VarFcnTillot(MaterialModelData &data)
  ~VarFcnTillot() {}

  inline double GetPressure(double rho, double e) {return GetPressureCase[GetCase(rho,e)](rho,e);}

  inline double GetInternalEnergyPerUnitMass(double rho, double p);

  inline double GetDensity(double p, double e);

  inline double GetDpdrho(double rho, double e) {return GetDpdrhoCase[GetCase(rho,e)](rho,e);}

  inline double GetBigGamma(double rho, double e) {return GetGammaCase[GetCase(rho,e)](rho,e);}

  I AM HERE

private:

  inline double GetChiWithEta(double eta, double e) {return 1.0/(e/(e0*eta*eta)+1.0);}
  inline double GetChiWithOmega(double omega, double e) {return 1.0/(e/e0*(omega+1.0)*(omega+1)+1.0);}

  //! Determine the case id. Returns 0, 1, 2, or 3 for Case 1, 2, 3, and 1|2
  inline int GetCase(double rho, double e) {
    if(rho<=0.0 || e<0.0) {
      fprintf(stdout,"\033[0;31m*** Error: VarFcnTillot::GetCase detected negative rho (%e) or e (%e).\033[0m\n",
              rho, e);
      exit(-1);
    }
    if(rho>=rho0) return 0; // Case 1 
    if(e>=eCV)    return 1; // Case 2
    if(rho<rhoIV) return 2; // Case 3
    if(e<=eIV)    return 0; // Case 1
    return 3; // Case 1|2
  }
   
  
  /********************************
   *            Case 1
   *******************************/
  inline double GetPressure1(double rho, double e) {
    double eta = rho/rho0;
    double mu  = eta - 1.0;
    return (a + b*GetChiWithEta(eta,e))*rho*e + (A + B*mu)*mu;
  }

  inline double GetGamma1(double rho, double e) {
    double chi = GetChiWithEta(rho/rho0, e);
    return a + b*chi*chi;
  }

  inline double GetDpdrho1(double rho, double e) {
    double eta = rho/rho0;
    double mu  = eta - 1.0;
    double chi = GetChiWithEta(eta, e);
    return a*e + b*e*chi*chi*(1.0 + 3.0*e/(e0*eta*eta)) + (A + 2.0*B*mu)/rho0;
  }
         

  /********************************
   *            Case 2
   *******************************/
  inline double GetPressure2(double rho, double e) {
    double mu    = rho/rho0 - 1.0;
    double omega = rho0/rho - 1.0;
    double rho_e = rho*e;
    return a*rho_e + (b*rho_e*GetChiWithOmega(omega,e) + A*mu*exp(-beta*mu))*exp(-alpha*omega*omega);
  }
 
  inline double GetGamma2(double rho, double e) {
    double omega = rho0/rho - 1.0;
    double chi   = GetChiWithOmega(omega, e);
    return a + b*chi*chi*exp(-alpha*omega*omega);  
  }

  inline double GetDpdrho2(double rho, double e) {
    double omega = rho0/rho - 1.0;
    double chi   = GetChiWithOmega(omega, e);
    double omom  = omega*(1.0+omega);
    return a*e + A/rho0*(1.0 - omom*(beta+2.0*alpha*omega))*exp(-omega*(beta+alpha*omega))
               + b*e*chi*chi*(1.0 + 2.0*alpha*omom + e/e0*(1.0+omega)*(1.0+omega)*(3.0+2.0*alpha*omom))*exp(-alpha*omega*omega);
  }


  /********************************
   *            Case 3 
   *******************************/
  inline double GetPressure3(double rho, double e) {
    double eta = rho/rho0;
    return (a + b*GetChiWithEta(eta,e))*rho*e + A*(eta - 1.0);
  }

  inline double GetGamma3(double rho, double e) { //same as GetGamma1
    double chi = GetChiWithEta(rho/rho0, e);
    return a + b*chi*chi;
  }

  inline double GetDpdrho3(double rho, double e) {
    double eta = rho/rho0;
    double mu  = eta - 1.0;
    double chi = GetChiWithEta(eta, e);
    return a*e + b*e*chi*chi*(1.0 + 3.0*e/(e0*eta*eta)) + A/rho0;
  }


  /********************************
   *           Case 1|2 
   *******************************/
  inline double GetPressure12(double rho, double e) {
    return ((eCV-e)*GetPressure1(rho,e) + (e-eIV)*GetPressure2(rho,e))/(eCV-eIV);
  }

  inline double GetGamma12(double rho, double e) {
    return (  (GetPressure2(rho,e)-GetPressure1(rho,e))/rho 
            + (eCV-e)*GetGamma1(rho,e) + (e-eIV)*GetGamma2(rho,e) )/(eCV-eIV);
  }

  inline double GetDpdrho12(double rho, double e) {
    return ((eCV-e)*GetDpdrho1(rho,e) + (e-eIV)*GetDpdrho2(rho,e))/(eCV-eIV);
  }


  /********************************
   * Function Arrays (All Cases)      
   *******************************/
  typedef double (*DoubleFunction) (double rho, double e);
  DoubleFunction GetPressureCase[] = {GetPressure1, GetPressure2, GetPressure3, GetPressure12};
  DoubleFunction GetGammaCase[]    = {GetGamma1, GetGamma2, GetGamma3, GetGamma12};
  DoubleFunction GetDpdrhoCase[]   = {GetDpdrho1, GetDpdrho2, GetDpdrho3, GetDpdrho12};

};

#endif
