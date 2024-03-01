/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _VAR_FCN_ANEOS_BASE_H_
#define _VAR_FCN_ANEOS_BASE_H_

#include <VarFcnBase.h>

/****************************************************************************
 * This class is the base class for a class of (complete) EOS, known as ANEOS,
 * where "AN" represents "Analytic". See KW's notes. This is also referred to
 * as "Physically-based" EOS. The main idea is to model the Helmholtz free energy, 
 * more specifically, as a summation of multiple contributions, then model each 
 * term "analytically".
 * This type of EOS is based on defining specific Helmholtz free energy as a function
 * of density (rho) and temperature (T), i.e. F(rho,T). 
 ***************************************************************************/

class VarFcnANEOSBase : public VarFcnBase {

protected:

  //! Numerical parameters
  double finite_difference_step; //!< non-dimensional step size for finite difference

public:
  
  //! Constructor (Note: Parameter values must be set in derived classes!)
  VarFcnANEOSBase(MaterialModelData &data) : VarFcnBase(data) {finite_difference_step = 1.0e-6;}
  ~VarFcnANEOSBase() {}

  // ------------------------------------------------------------------------------
  //! Compute derivatives of p(rho,e) by finite difference (central difference)
  // ------------------------------------------------------------------------------

  //! dpdrho = \frac{\partial p(\rho,e)}{\partial \rho}
  virtual double GetDpdrho(double rho, double e) {
    assert(rho>0);
    double drho = finite_difference_step*rho;
    double pminus = GetPressure(rho - drho, e);
    double pplus  = GetPressure(rho + drho, e);
    return (pplus - pminus)/(2.0*drho);
  }

  //! BigGamma = 1/rho*(\frac{\partial p(\rho,e)}{\partial e})   (Gruneisen parameter)
  virtual double GetBigGamma(double rho, double e) {
    double de = finite_difference_step*(e==0.0 ? 1.0 : fabs(e));
    double pminus = GetPressure(rho, e - de);
    double pplus  = GetPressure(rho, e + de);
    return (pplus - pminus)/(2.0*de*rho);
  }

};

#endif
