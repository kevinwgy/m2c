/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _VAR_FCN_HOMO_INCOMP_H
#define _VAR_FCN_HOMO_INCOMP_H

#include <VarFcnBase.h>
#include <fstream>

/********************************************************************************
 * This class is the VarFcn class for homogeneous incompressible materials
 * Temperature law: e - e0 = c*(T - T0)
 ********************************************************************************/
class VarFcnHomoIncomp : public VarFcnBase {

private:

  double rho0, p0, c, T0, e0, invc;

public:

  VarFcnHomoIncomp(MaterialModelData &data) : VarFcnBase(data),
                                              rho0(data.incompModel.rho0), p0(data.incompModel.p0),
                                              c(data.incompModel.c), T0(data.incompModel.T0),
                                              e0(data.incompModel.e0) {
    type = HOMOGENEOUS_INCOMPRESSIBLE;
    invc = c==0.0 ? 0.0 : 1.0/c;
  } 
  ~VarFcnHomoIncomp() {}

  //! ----- EOS-Specific Functions -----
  inline double GetPressure([[maybe_unused]] double rho, [[maybe_unused]] double e) {
    fprintf(stdout,"\033[0;31m*** Error: GetPressure undefined for incompressible materials.\033[0m\n"); 
    exit(-1); return 0.0;}

  inline double GetReferencePressure() {return p0;}

  inline double GetInternalEnergyPerUnitMass([[maybe_unused]] double rho, [[maybe_unused]] double p) {
    fprintf(stdout,"\033[0;31m*** Error: GetInternalEnergyPerUnitMass undefined for incompressible "
            "materials.\033[0m\n"); 
    exit(-1); return 0.0;}

  inline double GetReferenceInternalEnergyPerUnitMass() {return e0;}

  inline double GetDensity([[maybe_unused]] double p, [[maybe_unused]] double e) {return rho0;}

  inline double GetDpdrho([[maybe_unused]] double rho, [[maybe_unused]] double e) {
    fprintf(stdout,"\033[0;31m*** Error: GetDpdrho undefined for incompressible materials.\033[0m\n"); 
    exit(-1); return 0.0;}

  inline double GetBigGamma([[maybe_unused]] double rho, [[maybe_unused]] double e) {
    fprintf(stdout,"\033[0;31m*** Error: GetBigGamma undefined for incompressible materials.\033[0m\n"); 
    exit(-1); return 0.0;}

  inline double GetTemperature([[maybe_unused]] double rho, double e) {return T0 + invc*(e-e0);}

  inline double GetReferenceTemperature() {return T0;}

  inline double GetInternalEnergyPerUnitMassFromTemperature([[maybe_unused]] double rho, double T) {
    return e0 + c*(T-T0);}

  inline double GetInternalEnergyPerUnitMassFromEnthalpy([[maybe_unused]] double rho,
                                                         [[maybe_unused]] double h) {
    fprintf(stdout,"\033[0;31m*** Error: GetInternalEnergyPerUnitMassFromEnthalpy undefined for"
            " incompressible materials.\033[0m\n"); 
    exit(-1); return 0.0;}

  bool CheckState(double rho, double p, bool silence = false) {
    if(!std::isfinite(rho) || !std::isfinite(p)) {
      if(!silence)
        fprintf(stdout, "*** Error: CheckState failed. rho = %e, p = %e.\033[0m\n", rho, p);
      return true;
    }
    if(rho <= 0.0) {
      if(!silence && verbose>1)
        fprintf(stdout, "Warning: Negative density. rho = %e, p = %e.\n", rho, p);
      return true;
    }
    return false;
  }

  bool CheckState(double *V, bool silence = false) {
    if(!std::isfinite(V[0]) || !std::isfinite(V[1]) || !std::isfinite(V[2]) || !std::isfinite(V[3]) || 
       !std::isfinite(V[4])) {
      if(!silence)
        fprintf(stdout, "\033[0;31m*** Error: CheckState failed. V = %e %e %e %e %e.\033[0m\n",
                V[0], V[1], V[2], V[3], V[4]);
      return true;
    }
    return CheckState(V[0], V[4], silence);
  }

  inline bool CheckPhaseTransition([[maybe_unused]] int id/*id of the other phase*/) {return false;}

  //! Overwrite the calculations done in the base class
  inline void ConservativeToPrimitive([[maybe_unused]] double *U, [[maybe_unused]] double *V) {
    fprintf(stdout,"\033[0;31m*** Error: ConservativeToPrimitive undefined for"
            " incompressible materials.\033[0m\n"); 
    exit(-1);;}

  inline void PrimitiveToConservative([[maybe_unused]] double *V, [[maybe_unused]] double *U) {
    fprintf(stdout,"\033[0;31m*** Error: PrimitiveToConservative undefined for"
            " incompressible materials.\033[0m\n"); 
    exit(-1);}

  inline double ComputeSoundSpeed([[maybe_unused]] double rho, [[maybe_unused]] double e) {
    fprintf(stdout,"\033[0;31m*** Error: ComputeSoundSpeed undefined for incompressible materials.\033[0m\n"); 
    exit(-1); return 0.0;}

  inline double ComputeSoundSpeedSquare([[maybe_unused]] double rho, [[maybe_unused]] double e) {
    fprintf(stdout,"\033[0;31m*** Error: ComputeSoundSpeedSquare undefined for"
            " incompressible materials.\033[0m\n"); 
    exit(-1); return 0.0;}

  inline double ComputeMachNumber([[maybe_unused]] double *V) {
    fprintf(stdout,"\033[0;31m*** Error: ComputeMachNumber undefined for incompressible materials.\033[0m\n"); 
    exit(-1); return 0.0;}

  inline double ComputeEnthalpyPerUnitMass([[maybe_unused]] double rho, [[maybe_unused]] double p) {
    fprintf(stdout,"\033[0;31m*** Error: ComputeEnthalpyPerUnitMass undefined for "
                   "incompressible materials.\033[0m\n"); 
    exit(-1); return 0.0;}

  inline double ComputeTotalEnthalpyPerUnitMass([[maybe_unused]] double *V) {
    fprintf(stdout,"\033[0;31m*** Error: ComputeTotalEnthalpyPerUnitMass undefined for "
                   "incompressible materials.\033[0m\n"); 
    exit(-1); return 0.0;}

  bool ClipDensityAndPressure(double *V, double *U) {
    if(U) {fprintf(stdout,"\033[0;31m*** Error: ClipDensityAndPressure cannot handle 'U' for "
                   "incompressible material.\033[0m\n"); exit(-1);}
    bool clip = false;
    if(V[4]<pmin) {V[4] = pmin; clip = true;}
    if(V[4]>pmax) {V[4] = pmax; clip = true;}
    return clip;
  }

};

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

#endif
