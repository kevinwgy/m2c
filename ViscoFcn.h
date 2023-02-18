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
  virtual double GetMu(double* V = NULL, double div = 0.0, double h = 0.0) {return 0.0;}
  virtual double GetLambda(double* V = NULL, double div = 0.0, double h = 0.0) {return 0.0;}

  //! Compute the stress tensor components

  virtual double GetTauXX(double* dVdx, double* dVdy, double* dVdz, double* V = NULL, double* h = NULL) {
    return 0.0;}

  virtual double GetTauXY(double* dVdx, double* dVdy, double* dVdz, double* V = NULL, double* h = NULL) {
    return 0.0;}

  virtual double GetTauXZ(double* dVdx, double* dVdy, double* dVdz, double* V = NULL, double* h = NULL) {
    return 0.0;}

  virtual double GetTauYY(double* dVdx, double* dVdy, double* dVdz, double* V = NULL, double* h = NULL) {
    return 0.0;}

  virtual double GetTauYZ(double* dVdx, double* dVdy, double* dVdz, double* V = NULL, double* h = NULL) {
    return 0.0;}

  virtual double GetTauZZ(double* dVdx, double* dVdy, double* dVdz, double* V = NULL, double* h = NULL) {
    return 0.0;}


  //! Compute the viscous flux function
  void EvaluateViscousFluxFunction_F(double* flux/*output*/, double* dVdx, double* dVdy, double* dVdz, 
                                     double* V, double* h = NULL) {
    if(type==NONE) {
      flux[0] = flux[1] = flux[2] = flux[3] = flux[4] = 0.0;
      return;
    }

    double tauxx = GetTauXX(dVdx,dVdy,dVdz,V,h);
    double tauxy = GetTauXY(dVdx,dVdy,dVdz,V,h);
    double tauxz = GetTauXZ(dVdx,dVdy,dVdz,V,h);
    flux[0] = 0.0;
    flux[1] = tauxx;
    flux[2] = tauxy;
    flux[3] = tauxz;
    flux[4] = V[0]*tauxx + V[1]*tauxy + V[2]*tauxz;
  }

  //! Compute the viscous flux function
  void EvaluateViscousFluxFunction_G(double* flux/*output*/, double* dVdx, double* dVdy, double* dVdz, 
                                     double* V, double* h = NULL) {
    if(type==NONE) {
      flux[0] = flux[1] = flux[2] = flux[3] = flux[4] = 0.0;
      return;
    }

    double tauxy = GetTauXY(dVdx,dVdy,dVdz,V,h);
    double tauyy = GetTauYY(dVdx,dVdy,dVdz,V,h);
    double tauyz = GetTauYZ(dVdx,dVdy,dVdz,V,h);
    flux[0] = 0.0;
    flux[1] = tauxy;
    flux[2] = tauyy;
    flux[3] = tauyz;
    flux[4] = V[0]*tauxy + V[1]*tauyy + V[2]*tauyz;
  }

  //! Compute the viscous flux function
  void EvaluateViscousFluxFunction_H(double* flux/*output*/, double* dVdx, double* dVdy, double* dVdz, 
                                     double* V, double* h = NULL) {
    if(type==NONE) {
      flux[0] = flux[1] = flux[2] = flux[3] = flux[4] = 0.0;
      return;
    }

    double tauxz = GetTauXZ(dVdx,dVdy,dVdz,V,h);
    double tauyz = GetTauYZ(dVdx,dVdy,dVdz,V,h);
    double tauzz = GetTauZZ(dVdx,dVdy,dVdz,V,h);
    flux[0] = 0.0;
    flux[1] = tauxz;
    flux[2] = tauyz;
    flux[3] = tauzz;
    flux[4] = V[0]*tauxz + V[1]*tauyz + V[2]*tauzz;
  }

};

//---------------------------------------------------------------------------------

class ViscoFcnConstant : public ViscoFcnBase {

  double dyna, bulk; //See Kevin's notes

public:

  ViscoFcnConstant(ViscosityModelData &vis, VarFcnBase &vf_) : ViscoFcnBase(vf_) { 
    assert(vis.type == ViscosityModelData::CONSTANT);
    type = CONSTANT;
    dyna = vis.dynamicViscosity;
    bulk = vis.bulkViscosity;
  }
  ~ViscoFcnConstant() {}

  inline double GetMu(double* V = NULL, double div = 0.0, double h = 0.0) {return dyna;}

  inline double GetLambda(double* V = NULL, double div = 0.0, double h = 0.0) {return bulk-2.0/3.0*dyna;}

  //! Compute the stress tensor components

  inline double GetTauXX(double* dVdx, double* dVdy, double* dVdz, double* V = NULL, double* h = NULL) {
    double div = dVdx[0]+dVdy[1]+dVdz[2];
    return bulk*div + dyna*(2.0*dVdx[0] - 2.0/3.0*div);
  }

  inline double GetTauXY(double* dVdx, double* dVdy, double* dVdz, double* V = NULL, double* h = NULL) {
    return dyna*(dVdx[1]+dVdy[0]);
  }

  inline double GetTauXZ(double* dVdx, double* dVdy, double* dVdz, double* V = NULL, double* h = NULL) {
    return dyna*(dVdx[2]+dVdz[0]);
  }

  inline double GetTauYY(double* dVdx, double* dVdy, double* dVdz, double* V = NULL, double* h = NULL) {
    double div = dVdx[0]+dVdy[1]+dVdz[2];
    return bulk*div + dyna*(2.0*dVdy[1] - 2.0/3.0*div);
  }
 
  inline double GetTauYZ(double* dVdx, double* dVdy, double* dVdz, double* V = NULL, double* h = NULL) {
    return dyna*(dVdy[2]+dVdz[1]);
  }

  inline double GetTauZZ(double* dVdx, double* dVdy, double* dVdz, double* V = NULL, double* h = NULL) {
    double div = dVdx[0]+dVdy[1]+dVdz[2];
    return bulk*div + dyna*(2.0*dVdz[2] - 2.0/3.0*div);
  }

};

//---------------------------------------------------------------------------------

class ViscoFcnSutherland : public ViscoFcnBase {

  double mu0, T0, Smu; //See Kevin's notes
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

  inline double GetTauXX(double* dVdx, double* dVdy, double* dVdz, double* V, double* h = NULL) {
    double dyna = GetMu(V);
    double div = dVdx[0]+dVdy[1]+dVdz[2];
    return bulk*div + dyna*(2.0*dVdx[0] - 2.0/3.0*div);
  }

  inline double GetTauXY(double* dVdx, double* dVdy, double* dVdz, double* V, double* h = NULL) {
    double dyna = GetMu(V);
    return dyna*(dVdx[1]+dVdy[0]);
  }

  inline double GetTauXZ(double* dVdx, double* dVdy, double* dVdz, double* V, double* h = NULL) {
    double dyna = GetMu(V);
    return dyna*(dVdx[2]+dVdz[0]);
  }

  inline double GetTauYY(double* dVdx, double* dVdy, double* dVdz, double* V, double* h = NULL) {
    double dyna = GetMu(V);
    double div = dVdx[0]+dVdy[1]+dVdz[2];
    return bulk*div + dyna*(2.0*dVdy[1] - 2.0/3.0*div);
  }
 
  inline double GetTauYZ(double* dVdx, double* dVdy, double* dVdz, double* V, double* h = NULL) {
    double dyna = GetMu(V);
    return dyna*(dVdy[2]+dVdz[1]);
  }

  inline double GetTauZZ(double* dVdx, double* dVdy, double* dVdz, double* V, double* h = NULL) {
    double dyna = GetMu(V);
    double div = dVdx[0]+dVdy[1]+dVdz[2];
    return bulk*div + dyna*(2.0*dVdz[2] - 2.0/3.0*div);
  }



  inline double GetMu(double* V, double div = 0.0, double h = 0.0) {
    assert(V);
    double T = vf.GetTemperature(V[0]/*rho*/,V[4]/*p*/);
    return mu0*pow(T/T0,1.5)*(T0+Smu)/(T+Smu);
  }
  
  inline double GetLambda(double* V, double div = 0.0, double h = 0.0) {
    return bulk - 2.0/3.0*GetMu(V,div,h);
  }

};

//---------------------------------------------------------------------------------
//This is an artificial viscosity model by Alex Rodionov (2017,2018), an extension of the
//classical model by von Neumann and Richtmyer (1950)
class ViscoFcnRodionov : public ViscoFcnBase {

  double Cav, Cth; //See Kevin's notes
  double bulk; 

public:

  ViscoFcnRodionov(ViscosityModelData &vis, VarFcnBase &vf_) : ViscoFcnBase(vf_) { 
    assert(vis.type == ViscosityModelData::ARTIFICIAL_RODIONOV);
    type = ARTIFICIAL_RODIONOV;
    Cav = vis.Cav;
    Cth = vis.Cth;
    bulk = vis.bulkViscosity; //generally set to 0
  }
  ~ViscoFcnRodionov() {}

  //! Compute the stress tensor components

  inline double GetTauXX(double* dVdx, double* dVdy, double* dVdz, double* V, double* h) {
    double div = dVdx[0]+dVdy[1]+dVdz[2];
    double dyna = GetMu(V, div, *h);
    return bulk*div + dyna*(2.0*dVdx[0] - 2.0/3.0*div);
  }

  inline double GetTauXY(double* dVdx, double* dVdy, double* dVdz, double* V, double* h) {
    double div = dVdx[0]+dVdy[1]+dVdz[2];
    double dyna = GetMu(V, div, *h);
    return dyna*(dVdx[1]+dVdy[0]);
  }

  inline double GetTauXZ(double* dVdx, double* dVdy, double* dVdz, double* V, double* h) {
    double div = dVdx[0]+dVdy[1]+dVdz[2];
    double dyna = GetMu(V, div, *h);
    return dyna*(dVdx[2]+dVdz[0]);
  }

  inline double GetTauYY(double* dVdx, double* dVdy, double* dVdz, double* V, double* h) {
    double div = dVdx[0]+dVdy[1]+dVdz[2];
    double dyna = GetMu(V, div, *h);
    return bulk*div + dyna*(2.0*dVdy[1] - 2.0/3.0*div);
  }
 
  inline double GetTauYZ(double* dVdx, double* dVdy, double* dVdz, double* V, double* h) {
    double div = dVdx[0]+dVdy[1]+dVdz[2];
    double dyna = GetMu(V, div, *h);
    return dyna*(dVdy[2]+dVdz[1]);
  }

  inline double GetTauZZ(double* dVdx, double* dVdy, double* dVdz, double* V, double* h) {
    double div = dVdx[0]+dVdy[1]+dVdz[2];
    double dyna = GetMu(V, div, *h);
    return bulk*div + dyna*(2.0*dVdz[2] - 2.0/3.0*div);
  }




  inline double GetMu(double* V, double div, double h) {
    assert(V);
    assert(h>0.0);
    double e = vf.GetInternalEnergyPerUnitMass(V[0]/*rho*/, V[4]/*p*/);
    double c = vf.ComputeSoundSpeed(V[0], e);
    double p = Cth*c/h;
    if(div < -p)
      return Cav*V[0]*h*h*sqrt(div*div - p*p);
    else
      return 0.0;
  }
 
  inline double GetLambda(double* V, double div, double h) {
    return bulk - 2.0/3.0*GetMu(V,div,h);
  }
 
};


#endif
