/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _VISCO_FCN_H_
#define _VISCO_FCN_H_

#include<IoData.h>
#include<VarFcnBase.h>
#include<cassert>

/****************************************************
 * Class ViscoFcnBase and the derived classes are
 * responsible for calculating viscous stresses.
 * Each material has a separate ViscoFcn
 ****************************************************/

//---------------------------------------------------------------------------------
//
class ViscoFcnBase {

protected:

  VarFcnBase &vf;

public:

  enum Type {NONE = 0, CONSTANT = 1, SUTHERLAND = 2, ARTIFICIAL_RODIONOV = 3} type; 

  ViscoFcnBase(VarFcnBase &vf_) : vf(vf_), type(NONE) {}
  virtual ~ViscoFcnBase() {}

  //! Get viscosity coefficients (may depend on state variable and even local mesh resolution)
  virtual double GetMu([[maybe_unused]] double* V3 = NULL, [[maybe_unused]] double rho = 0.0, [[maybe_unused]] double p = 0.0,
                       [[maybe_unused]] double div = 0.0, [[maybe_unused]] double h = 0.0) {return 0.0;}
  virtual double GetLambda([[maybe_unused]] double* V3 = NULL, [[maybe_unused]] double rho = 0.0, [[maybe_unused]] double p = 0.0,
                           [[maybe_unused]] double div = 0.0, [[maybe_unused]] double h = 0.0) {return 0.0;}

  //! Compute the stress tensor components

  virtual double GetTauXX([[maybe_unused]] double* dVdx, [[maybe_unused]] double* dVdy, [[maybe_unused]] double* dVdz,
                          [[maybe_unused]] double* V3 = NULL, [[maybe_unused]] double rho = 0.0, [[maybe_unused]] double p = 0.0,
                          [[maybe_unused]] double* h = NULL) {
    return 0.0;}

  virtual double GetTauXY([[maybe_unused]] double* dVdx, [[maybe_unused]] double* dVdy, [[maybe_unused]] double* dVdz,
                          [[maybe_unused]] double* V3 = NULL, [[maybe_unused]] double rho = 0.0, [[maybe_unused]] double p = 0.0,
                          [[maybe_unused]] double* h = NULL) {
    return 0.0;}

  virtual double GetTauXZ([[maybe_unused]] double* dVdx, [[maybe_unused]] double* dVdy, [[maybe_unused]] double* dVdz, 
                          [[maybe_unused]] double* V3= NULL, [[maybe_unused]] double rho = 0.0,  [[maybe_unused]] double p = 0.0,
                          [[maybe_unused]] double* h = NULL) {
    return 0.0;}

  virtual double GetTauYY([[maybe_unused]] double* dVdx, [[maybe_unused]] double* dVdy, [[maybe_unused]] double* dVdz,
                          [[maybe_unused]] double* V3 = NULL, [[maybe_unused]] double rho = 0.0, [[maybe_unused]] double p = 0.0,
                          [[maybe_unused]] double* h = NULL) {
    return 0.0;}

  virtual double GetTauYZ([[maybe_unused]] double* dVdx, [[maybe_unused]] double* dVdy, [[maybe_unused]] double* dVdz,
                          [[maybe_unused]] double* V3 = NULL, [[maybe_unused]] double rho = 0.0, [[maybe_unused]] double p = 0.0,
                          [[maybe_unused]] double* h = NULL) {
    return 0.0;}

  virtual double GetTauZZ([[maybe_unused]] double* dVdx, [[maybe_unused]] double* dVdy, [[maybe_unused]] double* dVdz,
                          [[maybe_unused]] double* V3 = NULL, [[maybe_unused]] double rho = 0.0, [[maybe_unused]] double p = 0.0,
                          [[maybe_unused]] double* h = NULL) {
    return 0.0;}


  //! Compute the viscous flux function
  void EvaluateViscousFluxFunction_F(double* flux/*output*/, double* dVdx, double* dVdy, double* dVdz, 
                                     double* V3, double rho = 0.0, double p = 0.0, double* h = NULL) {
    if(type==NONE) {
      flux[0] = flux[1] = flux[2] = flux[3] = flux[4] = 0.0;
      return;
    }

    double tauxx = GetTauXX(dVdx,dVdy,dVdz,V3,rho,p,h);
    double tauxy = GetTauXY(dVdx,dVdy,dVdz,V3,rho,p,h);
    double tauxz = GetTauXZ(dVdx,dVdy,dVdz,V3,rho,p,h);
    flux[0] = 0.0;
    flux[1] = tauxx;
    flux[2] = tauxy;
    flux[3] = tauxz;
    flux[4] = V3[0]*tauxx + V3[1]*tauxy + V3[2]*tauxz;
  }

  //! Compute the viscous flux function
  void EvaluateViscousFluxFunction_G(double* flux/*output*/, double* dVdx, double* dVdy, double* dVdz, 
                                     double* V3, double rho = 0.0, double p = 0.0, double* h = NULL) {
    if(type==NONE) {
      flux[0] = flux[1] = flux[2] = flux[3] = flux[4] = 0.0;
      return;
    }

    double tauxy = GetTauXY(dVdx,dVdy,dVdz,V3,rho,p,h);
    double tauyy = GetTauYY(dVdx,dVdy,dVdz,V3,rho,p,h);
    double tauyz = GetTauYZ(dVdx,dVdy,dVdz,V3,rho,p,h);
    flux[0] = 0.0;
    flux[1] = tauxy;
    flux[2] = tauyy;
    flux[3] = tauyz;
    flux[4] = V3[0]*tauxy + V3[1]*tauyy + V3[2]*tauyz;
  }

  //! Compute the viscous flux function
  void EvaluateViscousFluxFunction_H(double* flux/*output*/, double* dVdx, double* dVdy, double* dVdz, 
                                     double* V3, double rho = 0.0, double p = 0.0, double* h = NULL) {
    if(type==NONE) {
      flux[0] = flux[1] = flux[2] = flux[3] = flux[4] = 0.0;
      return;
    }

    double tauxz = GetTauXZ(dVdx,dVdy,dVdz,V3,rho,p,h);
    double tauyz = GetTauYZ(dVdx,dVdy,dVdz,V3,rho,p,h);
    double tauzz = GetTauZZ(dVdx,dVdy,dVdz,V3,rho,p,h);
    flux[0] = 0.0;
    flux[1] = tauxz;
    flux[2] = tauyz;
    flux[3] = tauzz;
    flux[4] = V3[0]*tauxz + V3[1]*tauyz + V3[2]*tauzz;
  }

};

//---------------------------------------------------------------------------------

class ViscoFcnConstant : public ViscoFcnBase {

  double dyna, bulk; //!< See Kevin's notes

public:

  ViscoFcnConstant(ViscosityModelData &vis, VarFcnBase &vf_) : ViscoFcnBase(vf_) { 
    assert(vis.type == ViscosityModelData::CONSTANT);
    type = CONSTANT;
    dyna = vis.dynamicViscosity;
    bulk = vis.bulkViscosity;
  }
  ~ViscoFcnConstant() {}

  inline double GetMu([[maybe_unused]] double* V3 = NULL, [[maybe_unused]] double rho = 0.0, [[maybe_unused]] double p = 0.0,
                      [[maybe_unused]] double div = 0.0, [[maybe_unused]] double h = 0.0) {return dyna;}

  inline double GetLambda([[maybe_unused]] double* V3 = NULL, [[maybe_unused]] double rho = 0.0, [[maybe_unused]] double p = 0.0,
                          [[maybe_unused]] double div = 0.0, [[maybe_unused]] double h = 0.0) {return bulk-2.0/3.0*dyna;}

  //! Compute the stress tensor components

  inline double GetTauXX(double* dVdx, double* dVdy, double* dVdz, [[maybe_unused]] double* V3 = NULL, 
                         [[maybe_unused]] double rho = 0.0, [[maybe_unused]] double p = 0.0, [[maybe_unused]] double* h = NULL) {
    double div = dVdx[0]+dVdy[1]+dVdz[2];
    return bulk*div + dyna*(2.0*dVdx[0] - 2.0/3.0*div);
  }

  inline double GetTauXY(double* dVdx, double* dVdy, [[maybe_unused]] double* dVdz, [[maybe_unused]] double* V3 = NULL,
                         [[maybe_unused]] double rho = 0.0, [[maybe_unused]] double p = 0.0, [[maybe_unused]] double* h = NULL) {
    return dyna*(dVdx[1]+dVdy[0]);
  }

  inline double GetTauXZ(double* dVdx, [[maybe_unused]] double* dVdy, double* dVdz, [[maybe_unused]] double* V3 = NULL,
                         [[maybe_unused]] double rho = 0.0, [[maybe_unused]] double p = 0.0, [[maybe_unused]] double* h = NULL) {
    return dyna*(dVdx[2]+dVdz[0]);
  }

  inline double GetTauYY(double* dVdx, double* dVdy, double* dVdz, [[maybe_unused]] double* V3 = NULL,
                         [[maybe_unused]] double rho = 0.0, [[maybe_unused]] double p = 0.0, [[maybe_unused]] double* h = NULL) {
    double div = dVdx[0]+dVdy[1]+dVdz[2];
    return bulk*div + dyna*(2.0*dVdy[1] - 2.0/3.0*div);
  }
 
  inline double GetTauYZ([[maybe_unused]] double* dVdx, double* dVdy, double* dVdz, [[maybe_unused]] double* V3 = NULL,
                         [[maybe_unused]] double rho = 0.0, [[maybe_unused]] double p = 0.0, [[maybe_unused]] double* h = NULL) {
    return dyna*(dVdy[2]+dVdz[1]);
  }

  inline double GetTauZZ(double* dVdx, double* dVdy, double* dVdz, [[maybe_unused]] double* V3 = NULL,
                         [[maybe_unused]] double rho = 0.0, [[maybe_unused]] double p = 0.0, [[maybe_unused]] double* h = NULL) {
    double div = dVdx[0]+dVdy[1]+dVdz[2];
    return bulk*div + dyna*(2.0*dVdz[2] - 2.0/3.0*div);
  }

};

//---------------------------------------------------------------------------------

class ViscoFcnSutherland : public ViscoFcnBase {

  double mu0, T0, Smu; //!< See Kevin's notes
  double bulk;

public:

  ViscoFcnSutherland(ViscosityModelData &vis, VarFcnBase &vf_) : ViscoFcnBase(vf_) { 
    assert(vis.type == ViscosityModelData::SUTHERLAND);
    type = SUTHERLAND;
    mu0 = vis.sutherlandMu0;
    T0 = vis.sutherlandT0;
    Smu= vis.sutherlandConstant;
    bulk = vis.bulkViscosity;
  }
  ~ViscoFcnSutherland() {}

  //! Compute the stress tensor components

  inline double GetTauXX(double* dVdx, double* dVdy, double* dVdz, double* V3,
                         double rho, double p, [[maybe_unused]] double* h = NULL) { //V3 is not used
    double dyna = GetMu(V3,rho,p);
    double div = dVdx[0]+dVdy[1]+dVdz[2];
    return bulk*div + dyna*(2.0*dVdx[0] - 2.0/3.0*div);
  }

  inline double GetTauXY(double* dVdx, double* dVdy, [[maybe_unused]] double* dVdz, double* V3,
                         double rho, double p, [[maybe_unused]] double* h = NULL) { //V3 is not used
    double dyna = GetMu(V3,rho,p);
    return dyna*(dVdx[1]+dVdy[0]);
  }

  inline double GetTauXZ(double* dVdx, [[maybe_unused]] double* dVdy, double* dVdz, double* V3,
                         double rho, double p, [[maybe_unused]] double* h = NULL) { //V3 is not used
    double dyna = GetMu(V3,rho,p);
    return dyna*(dVdx[2]+dVdz[0]);
  }

  inline double GetTauYY(double* dVdx, double* dVdy, double* dVdz, double* V3,
                         double rho, double p, [[maybe_unused]] double* h = NULL) { //V3 is not used
    double dyna = GetMu(V3,rho,p);
    double div = dVdx[0]+dVdy[1]+dVdz[2];
    return bulk*div + dyna*(2.0*dVdy[1] - 2.0/3.0*div);
  }
 
  inline double GetTauYZ([[maybe_unused]] double* dVdx, double* dVdy, double* dVdz, double* V3,
                         double rho, double p, [[maybe_unused]] double* h = NULL) { //V3 is not used
    double dyna = GetMu(V3,rho,p);
    return dyna*(dVdy[2]+dVdz[1]);
  }

  inline double GetTauZZ(double* dVdx, double* dVdy, double* dVdz, double* V3,
                         double rho, double p, [[maybe_unused]] double* h = NULL) { //V3 is not used
    double dyna = GetMu(V3,rho,p);
    double div = dVdx[0]+dVdy[1]+dVdz[2];
    return bulk*div + dyna*(2.0*dVdz[2] - 2.0/3.0*div);
  }



  inline double GetMu([[maybe_unused]] double* V3, double rho, double p, [[maybe_unused]] double div = 0.0, [[maybe_unused]] double h = 0.0) {
    assert(rho>0);
    double T = vf.GetTemperature(rho,p);
    return mu0*pow(T/T0,1.5)*(T0+Smu)/(T+Smu);
  }
  
  inline double GetLambda(double* V3, double rho, double p, [[maybe_unused]] double div = 0.0, [[maybe_unused]] double h = 0.0) {
    return bulk - 2.0/3.0*GetMu(V3,rho,p);
  }

};

//---------------------------------------------------------------------------------
//This is an artificial viscosity model by Alex Rodionov (2017,2018), an extension of the
//classical model by von Neumann and Richtmyer (1950)
class ViscoFcnRodionov : public ViscoFcnBase {

  double Cav, Cth; //!< See Kevin's notes
  double bulk; 

public:

  ViscoFcnRodionov(ViscosityModelData &vis, VarFcnBase &vf_) : ViscoFcnBase(vf_) { 
    assert(vis.type == ViscosityModelData::ARTIFICIAL_RODIONOV);
    type = ARTIFICIAL_RODIONOV;
    Cav = vis.Cav;
    Cth = vis.Cth;
    bulk = vis.bulkViscosity; //!< generally set to 0
  }
  ~ViscoFcnRodionov() {}

  //! Compute the stress tensor components

  inline double GetTauXX(double* dVdx, double* dVdy, double* dVdz, double* V3,
                         double rho, double p, double* h) { //V3 not used
    assert(h);
    double div = dVdx[0]+dVdy[1]+dVdz[2];
    double dyna = GetMu(V3, rho, p, div, *h);
    return bulk*div + dyna*(2.0*dVdx[0] - 2.0/3.0*div);
  }

  inline double GetTauXY(double* dVdx, double* dVdy, double* dVdz, double* V3,
                         double rho, double p, double* h) { //V3 not used
    assert(h);
    double div = dVdx[0]+dVdy[1]+dVdz[2];
    double dyna = GetMu(V3, rho, p, div, *h);
    return dyna*(dVdx[1]+dVdy[0]);
  }

  inline double GetTauXZ(double* dVdx, double* dVdy, double* dVdz, double* V3,
                         double rho, double p, double* h) { //V3 not used
    assert(h);
    double div = dVdx[0]+dVdy[1]+dVdz[2];
    double dyna = GetMu(V3, rho, p, div, *h);
    return dyna*(dVdx[2]+dVdz[0]);
  }

  inline double GetTauYY(double* dVdx, double* dVdy, double* dVdz, double* V3,
                         double rho, double p, double* h) { //V3 not used
    assert(h);
    double div = dVdx[0]+dVdy[1]+dVdz[2];
    double dyna = GetMu(V3, rho, p, div, *h);
    return bulk*div + dyna*(2.0*dVdy[1] - 2.0/3.0*div);
  }
 
  inline double GetTauYZ(double* dVdx, double* dVdy, double* dVdz, double* V3,
                         double rho, double p, double* h) { //V3 not used
    assert(h);
    double div = dVdx[0]+dVdy[1]+dVdz[2];
    double dyna = GetMu(V3, rho, p, div, *h);
    return dyna*(dVdy[2]+dVdz[1]);
  }

  inline double GetTauZZ(double* dVdx, double* dVdy, double* dVdz, double* V3,
                         double rho, double p, double* h) { //V3 not used
    assert(h);
    double div = dVdx[0]+dVdy[1]+dVdz[2];
    double dyna = GetMu(V3, rho, p, div, *h);
    return bulk*div + dyna*(2.0*dVdz[2] - 2.0/3.0*div);
  }




  inline double GetMu([[maybe_unused]] double* V3, double rho, double p, double div, double h) { //V3 not used
    assert(h>0.0 && rho>0.0);
    double e = vf.GetInternalEnergyPerUnitMass(rho,p);
    double c = vf.ComputeSoundSpeed(rho, e);
    double pp = Cth*c/h;
    if(div < -pp)
      return Cav*rho*h*h*sqrt(div*div - pp*pp);
    else
      return 0.0;
  }
 
  inline double GetLambda(double* V3, double rho, double p, double div, double h) {
    return bulk - 2.0/3.0*GetMu(V3,rho,p,div,h);
  }
 
};


#endif
