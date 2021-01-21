#ifndef _VAR_FCN_MG_H
#define _VAR_FCN_MG_H

#include <VarFcnBase.h>
#include <fstream>
#include <Utils.h>

/********************************************************************************
 * This class is the VarFcn class for the Mie-Gruneisen EOS in Euler
 * Equations. Only elementary functions are declared and/or defined here.
 * All arguments must be pertinent to only a single grid node or a single
 * state.
 *
 * EOS: Pressure = rho0*C0*C0*(1-Gamma0/2*Xi)/(1-s*Xi)^2 + Gamm0*rho0*e
 * where
 *   e      : internal energy per unit mass.
 *   rho    : density
 *
 *   rho0   : ref. density
 *   C0     : bulk speed of sound
 *   Gamma0 : Gruneisen coefficient at ref. state
 *   s      : slope of shock Hugoniot
 *   Xi     : 1 - rho0/rho
 *
 *   Note: A temperature law can be defined as e = cv*(T-T0), where cv is a constant
 *         specific heat (at constant volume), and T0 a ref. temperature
 *
 ********************************************************************************/
class VarFcnMG : public VarFcnBase {

private:
  double rho0;
  double C0;
  double s;
  double Gamma0;
  double cv;

  double rho0_C0_C0; //rho0*C0*C0
  double Gamma0_over_2; //Gamma0/2
  double Gamma0_rho0; //Gamma0*rho0

public:
  VarFcnMG(MaterialModelData &data, bool verbose_ = true);
  ~VarFcnMG() {}

  //! ----- EOS-Specific Functions -----
  inline double GetPressure(double rho, double e) const {
    double Xi = 1.0 - rho0/rho;
    return (rho0_C0_C0*Xi*(1.0 - Gamma0_over_2*Xi))/((1.0-s*Xi)*(1.0-s*Xi)) + Gamma0_rho0*e;}

  inline double GetInternalEnergyPerUnitMass(double rho, double p) const {
    print_error("To be implemented!\n"); exit_mpi(); return 0.0;}

  inline double GetDensity(double p, double e) const {
    print_error("To be implemented!\n"); exit_mpi(); return 0.0;}

  inline double GetDpdrho(double rho, double e) const{
    print_error("To be implemented!\n"); exit_mpi(); return 0.0;}

  inline double GetBigGamma(double rho, double e) const {return Gamma0_rho0/rho;}

  //! Verify hyperbolicity (i.e. c^2 > 0): Report error if rho < 0 or p + Pstiff < 0 (Not p + gamma*Pstiff). 
  inline bool CheckState(double *V) const{
    double e = GetInternalEnergyPerUnitMass(V[0],V[4]);
    double c2 = GetDpdrho(V[0], e) + GetPressure(V[0], e)/V[0]*GetBigGamma(V[0], e);
    if(V[0] <= 0.0 || c2 <= 0.0){
      if(verbose)
        fprintf(stdout, "Warning:  found negative density (%e) or loss of hyperbolicity (p = %e).\n", V[0], V[4]);
      return true;
    }
    return false;
  }

};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
inline
VarFcnMG::VarFcnMG(MaterialModelData &data, bool verbose_) : VarFcnBase(data,verbose_) {

  if(data.eos != MaterialModelData::MIE_GRUNEISEN){
    fprintf(stderr, "*** Error: MaterialModelData is not of type Mie-Gruneisen.\n");
    exit(1);
  }

  type   = MIE_GRUNEISEN;

  rho0   = data.mgModel.rho0;
  C0     = data.mgModel.C0;
  s      = data.mgModel.s;
  Gamma0 = data.mgModel.Gamma0;

  cv     = data.mgModel.cv;

  rho0_C0_C0 = rho0*C0*C0;
  Gamma0_over_2 = 0.5*Gamma0;
  Gamma0_rho0 = Gamma0*rho0;
}

//------------------------------------------------------------------------------

#endif
