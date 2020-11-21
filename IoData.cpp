#include <IoData.h>
#include <parser/Assigner.h>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cstring>
#include <cmath>
#include <unistd.h>
//#include <dlfcn.h>

//------------------------------------------------------------------------------

MeshData::MeshData()
{
  type = CYLINDRICAL;
  x0 = 0.0;
  xmax = 1.0;
  y0 = 0.0;
  ymax = 1.0;
  Nx = 50;
  Ny = 50;
}

//------------------------------------------------------------------------------

void MeshData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 7, father);

  new ClassToken<MeshData>(ca, "Type", this,
                               reinterpret_cast<int MeshData::*>(&MeshData::type), 3,
                               "TwoDimensional", 0, "Cylindrical", 1, "ThreeDimensional");
  new ClassDouble<MeshData>(ca, "X0", this, &MeshData::x0);
  new ClassDouble<MeshData>(ca, "Xmax", this, &MeshData::xmax);
  new ClassDouble<MeshData>(ca, "Y0", this, &MeshData::y0);
  new ClassDouble<MeshData>(ca, "Ymax", this, &MeshData::ymax);
  new ClassInt<MeshData>(ca, "NumberOfElementsX", this, &MeshData::Nx);
  new ClassInt<MeshData>(ca, "NumberOfElementsY", this, &MeshData::Ny);
 } 

//------------------------------------------------------------------------------

GasModelData::GasModelData()
{

  type = IDEAL;
  specificHeatRatio = 1.4;
  idealGasConstant = 287.1;
  specificHeatPressure = -1.0;
  pressureConstant = 0.0;

}

//------------------------------------------------------------------------------

void GasModelData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 5, father);

  new ClassToken<GasModelData>(ca, "Type", this,
                               reinterpret_cast<int GasModelData::*>(&GasModelData::type), 2,
                               "Ideal", 0, "Stiffened", 1);
  new ClassDouble<GasModelData>(ca, "SpecificHeatRatio", this,
                                &GasModelData::specificHeatRatio);
  new ClassDouble<GasModelData>(ca, "IdealGasConstant", this,
                                &GasModelData::idealGasConstant);
  new ClassDouble<GasModelData>(ca, "SpecificHeatAtConstantPressure", this,
                                &GasModelData::specificHeatPressure);
  new ClassDouble<GasModelData>(ca, "PressureConstant", this,
                                &GasModelData::pressureConstant);

}

//------------------------------------------------------------------------------

FluidModelData::FluidModelData()
{

  fluid = STIFFENED_GAS;
  rhomin = -1.e9; // note: if these defaults are changed then doVerification()
  pmin = -1.e9;   //       in VarFcnBase.h must also be changed.

}

//------------------------------------------------------------------------------

void FluidModelData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 3, father);

  new ClassToken<FluidModelData>(ca, "Fluid", this,
                                 reinterpret_cast<int FluidModelData::*>(&FluidModelData::fluid), 1,
                                 "StiffenedGas", FluidModelData::STIFFENED_GAS);
  new ClassDouble<FluidModelData>(ca, "DensityCutOff", this, &FluidModelData::rhomin);
  new ClassDouble<FluidModelData>(ca, "PressureCutOff", this, &FluidModelData::pmin);

  gasModel.setup("GasModel", ca);

};

//------------------------------------------------------------------------------

EquationsData::EquationsData()
{
  numPhase = 1;
}

//------------------------------------------------------------------------------

void EquationsData::setup(const char *name, ClassAssigner *father)
{
  fluidModelMap.setup("FluidModel", 0);
}

//------------------------------------------------------------------------------

SchemeData::SchemeData(int af) : allowsFlux(af)
{

  flux = ROE;
  reconstruction = LINEAR;
  limiter = NONE;
  dissipation = SECOND_ORDER;

  beta = 1.0/3.0;
  gamma = 1.0;
  xiu  = -2.0/15.0;
  xic = -1.0/30.0;
  eps = 0.1;

  xirho = 1.0;
  xip = 1.0;

  vel_fac = sqrt(5.0);
}

//------------------------------------------------------------------------------

void SchemeData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner* ca;
  if (allowsFlux)
    ca = new ClassAssigner(name, 12, father);
  else
    ca = new ClassAssigner(name, 11, father);

  if (allowsFlux) {
    new ClassToken<SchemeData>
      (ca, "Flux", this,
       reinterpret_cast<int SchemeData::*>(&SchemeData::flux), 4,
       "Roe", 0, "VanLeer", 1, "HLLE", 2, "HLLC", 3);
  }

  new ClassToken<SchemeData>
    (ca, "Reconstruction", this,
     reinterpret_cast<int SchemeData::*>(&SchemeData::reconstruction), 2,
     "Constant", 0, "Linear", 1);

  new ClassToken<SchemeData>
    (ca, "Limiter", this,
     reinterpret_cast<int SchemeData::*>(&SchemeData::limiter), 6,
     "None", 0, "VanAlbada", 1, "Barth", 2, "Venkatakrishnan", 3, "PressureSensor", 4,
     "ExtendedVanAlbada",5);

  new ClassToken<SchemeData>
    (ca, "Dissipation", this,
     reinterpret_cast<int SchemeData::*>(&SchemeData::dissipation), 2,
     "SecondOrder", 0, "SixthOrder", 1);

  new ClassDouble<SchemeData>(ca, "Beta", this, &SchemeData::beta);
  new ClassDouble<SchemeData>(ca, "Gamma", this, &SchemeData::gamma);
  new ClassDouble<SchemeData>(ca, "XiU", this, &SchemeData::xiu);
  new ClassDouble<SchemeData>(ca, "XiC", this, &SchemeData::xic);
  new ClassDouble<SchemeData>(ca, "Eps", this, &SchemeData::eps);

  new ClassDouble<SchemeData>(ca, "XiRho", this, &SchemeData::xirho);
  new ClassDouble<SchemeData>(ca, "XiP", this, &SchemeData::xip);
  new ClassDouble<SchemeData>(ca, "VelFac", this, &SchemeData::vel_fac);

}

//------------------------------------------------------------------------------

BoundarySchemeData::BoundarySchemeData()
{
  type = STEGER_WARMING;
}

//------------------------------------------------------------------------------

void BoundarySchemeData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 1, father);

  new ClassToken<BoundarySchemeData>(ca, "Type", this,
         reinterpret_cast<int BoundarySchemeData::*>(&BoundarySchemeData::type), 5,
         "StegerWarming", 0, "ConstantExtrapolation", 1, "LinearExtrapolation", 2,
         "Ghidaglia", 3, "ModifiedGhidaglia", 4);

}

//------------------------------------------------------------------------------

SchemesData::SchemesData() : ls(0), tm(0)
{

}

//------------------------------------------------------------------------------

void SchemesData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 3, father);

  ns.setup("NavierStokes", ca);
  ls.setup("LevelSet",ca);
  bc.setup("Boundaries", ca);

  return ca;
}

//------------------------------------------------------------------------------

ExplicitData::ExplicitData()
{
  type = RUNGE_KUTTA_4;
}

//------------------------------------------------------------------------------

void ExplicitData::setup(const char *name, ClassAssigner *father)
{

 ClassAssigner *ca = new ClassAssigner(name, 1, father);

  new ClassToken<ExplicitData>
    (ca, "Type", this,
     reinterpret_cast<int ExplicitData::*>(&ExplicitData::type), 3,
     "RungeKutta4", 0, "RungeKutta2", 1, "ForwardEuler", 2);

}

//------------------------------------------------------------------------------

TsData::TsData()
{

  type = EXPLICIT;
  maxIts = 100;
  timestep = -1.0;
  maxTime = 1.e99;

}

//------------------------------------------------------------------------------

void TsData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 4, father);

  new ClassToken<TsData>(ca, "Type", this,
                         reinterpret_cast<int TsData::*>(&TsData::type), 2,
                         "Explicit", 0, "Implicit", 1);
  new ClassInt<TsData>(ca, "MaxIts", this, &TsData::maxIts);
  new ClassDouble<TsData>(ca, "TimeStep", this, &TsData::timestep);
  new ClassDouble<TsData>(ca, "MaxTime", this, &TsData::maxTime);

  expl.setup("Explicit", ca);
}

//------------------------------------------------------------------------------

void Input::readCmdLine(int argc, char** argv)
{
  if(argc==1) {
    fprintf(stderr,"ERROR: Input file not provided!\n");
    exit(-1);
  }
  cmdFileName = argv[1];
}

//------------------------------------------------------------------------------

void Input::readCmdFile()
{
  extern FILE *yyCmdfin;
  extern int yyCmdfparse();

  setupCmdFileVariables();
//  cmdFilePtr = freopen(cmdFileName, "r", stdin);
  yyCmdfin = cmdFilePtr = fopen(cmdFileName, "r");

  if (!cmdFilePtr) {
    fprintf(stderr,"*** Error: could not open \'%s\'\n", cmdFileName);
    exit(-1);
  }

  int error = yyCmdfparse();
  if (error) {
    fprintf(stderr,"*** Error: command file contained parsing errors\n");
    exit(error);
  }
  fclose(cmdFilePtr);
}

//------------------------------------------------------------------------------

void IoData::setupCmdFileVariables()
{
  domain.setup("Domain");
  eqs.setup("Equations");
  schemes.setup("Space");
  ts.setup("Time");
}

//------------------------------------------------------------------------------
