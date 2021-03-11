#include <Utils.h>
#include <IoData.h>
#include <parser/Assigner.h>
#include <parser/Dictionary.h>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cfloat>
#include <climits>
#include <cmath>
#include <unistd.h>
//#include <dlfcn.h>
using namespace std;

RootClassAssigner *nullAssigner = new RootClassAssigner;

//------------------------------------------------------------------------------

StateVariable::StateVariable()
{

  materialid  = 0;
  density     = -1.0;
  velocity_x  = 0.0;
  velocity_y  = 0.0;
  velocity_z  = 0.0;
  pressure    = -1.0;
  temperature = -1.0;

}

//------------------------------------------------------------------------------

void StateVariable::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 7, father);

  new ClassInt<StateVariable>(ca, "MaterialID", this, &StateVariable::materialid);
  new ClassDouble<StateVariable>(ca, "Density", this, &StateVariable::density);
  new ClassDouble<StateVariable>(ca, "VelocityX", this, &StateVariable::velocity_x);
  new ClassDouble<StateVariable>(ca, "VelocityY", this, &StateVariable::velocity_y);
  new ClassDouble<StateVariable>(ca, "VelocityZ", this, &StateVariable::velocity_z);
  new ClassDouble<StateVariable>(ca, "Pressure", this, &StateVariable::pressure);
  new ClassDouble<StateVariable>(ca, "Temperature", this, &StateVariable::temperature);

}

//------------------------------------------------------------------------------

PointData::PointData()
{
  x  = 0.0;
  y  = 0.0;
  z  = 0.0;
}

//------------------------------------------------------------------------------

Assigner *PointData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 4, nullAssigner);

  new ClassDouble<PointData>
    (ca, "X", this, &PointData::x);
  new ClassDouble<PointData>
    (ca, "Y", this, &PointData::y);
  new ClassDouble<PointData>
    (ca, "Z", this, &PointData::z);

  initialConditions.setup("InitialState", ca);

  return ca;
}

//------------------------------------------------------------------------------


PlaneData::PlaneData()
{

  cen_x  = 0.0;
  cen_y  = 0.0;
  cen_z  = 0.0;
  nx     = 0.0;
  ny     = 0.0;
  nz     = 0.0;

}

//------------------------------------------------------------------------------

Assigner *PlaneData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 7, nullAssigner);

  new ClassDouble<PlaneData> (ca, "Point_x", this, &PlaneData::cen_x);
  new ClassDouble<PlaneData> (ca, "Point_y", this, &PlaneData::cen_y);
  new ClassDouble<PlaneData> (ca, "Point_z", this, &PlaneData::cen_z);
  new ClassDouble<PlaneData> (ca, "Normal_x", this, &PlaneData::nx);
  new ClassDouble<PlaneData> (ca, "Normal_y", this, &PlaneData::ny);
  new ClassDouble<PlaneData> (ca, "Normal_z", this, &PlaneData::nz);

  initialConditions.setup("InitialState", ca);

  return ca;
}

//------------------------------------------------------------------------------

SphereData::SphereData()
{
  
  cen_x  = 0.0;
  cen_y  = 0.0;
  cen_z  = 0.0;
  radius = -1.0;

}

//------------------------------------------------------------------------------

Assigner *SphereData::getAssigner()
{
  
  ClassAssigner *ca = new ClassAssigner("normal", 5, nullAssigner);
  
  new ClassDouble<SphereData> (ca, "Center_x", this, &SphereData::cen_x);
  new ClassDouble<SphereData> (ca, "Center_y", this, &SphereData::cen_y);
  new ClassDouble<SphereData> (ca, "Center_z", this, &SphereData::cen_z);
  new ClassDouble<SphereData> (ca, "Radius", this, &SphereData::radius);
  
  initialConditions.setup("InitialState", ca);
  
  return ca;
}

//------------------------------------------------------------------------------

CylinderConeData::CylinderConeData() {

  cen_x  = 0.0;
  cen_y  = 0.0;
  cen_z  = 0.0;
  nx     = 0.0;
  ny     = 0.0;
  nz     = 1.0;

  r = 1.0;
  L = 0.0;

  cone_height = 0.0;
  opening_angle_degrees = 45.0;
}

//------------------------------------------------------------------------------

Assigner *CylinderConeData::getAssigner()
{
  ClassAssigner *ca = new ClassAssigner("normal", 11, nullAssigner);

  new ClassDouble<CylinderConeData> (ca, "Axis_x", this, &CylinderConeData::nx);
  new ClassDouble<CylinderConeData> (ca, "Axis_y", this, &CylinderConeData::ny);
  new ClassDouble<CylinderConeData> (ca, "Axis_z", this, &CylinderConeData::nz);
  new ClassDouble<CylinderConeData> (ca, "BaseCenter_x", this, &CylinderConeData::cen_x);
  new ClassDouble<CylinderConeData> (ca, "BaseCenter_y", this, &CylinderConeData::cen_y);
  new ClassDouble<CylinderConeData> (ca, "BaseCenter_z", this, &CylinderConeData::cen_z);
  new ClassDouble<CylinderConeData> (ca, "CylinderRadius", this, &CylinderConeData::r);
  new ClassDouble<CylinderConeData> (ca, "CylinderHeight", this, &CylinderConeData::L);

  new ClassDouble<CylinderConeData> (ca, "ConeOpeningAngleInDegrees", this, &CylinderConeData::opening_angle_degrees);
  new ClassDouble<CylinderConeData> (ca, "ConeHeight", this, &CylinderConeData::cone_height);

  initialConditions.setup("InitialState", ca);

  return ca;
}

//------------------------------------------------------------------------------

void MultiInitialConditionsData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 4, father);
  pointMap.setup("Point", ca);
  planeMap.setup("Plane", ca);
  sphereMap.setup("Sphere", ca);
  cylinderconeMap.setup("CylinderAndCone", ca);
}

//------------------------------------------------------------------------------

MeshData::MeshData()
{
  type = THREEDIMENSIONAL;
  x0 = 0.0;
  xmax = 1.0;
  y0 = 0.0;
  ymax = 1.0;
  z0 = 0.0;
  zmax = 1.0;
  Nx = 50;
  Ny = 50;
  Nz = 50;

  bc_x0   = NONE;
  bc_xmax = NONE;
  bc_y0   = NONE;
  bc_ymax = NONE;
  bc_z0   = NONE;
  bc_zmax = NONE;
}

//------------------------------------------------------------------------------

void MeshData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 16, father);

  new ClassToken<MeshData>(ca, "Type", this,
                               reinterpret_cast<int MeshData::*>(&MeshData::type), 2,
                               "ThreeDimensional", 0, "Cylindrical", 1);
  new ClassDouble<MeshData>(ca, "X0", this, &MeshData::x0);
  new ClassDouble<MeshData>(ca, "Xmax", this, &MeshData::xmax);
  new ClassDouble<MeshData>(ca, "Y0", this, &MeshData::y0);
  new ClassDouble<MeshData>(ca, "Ymax", this, &MeshData::ymax);
  new ClassDouble<MeshData>(ca, "Z0", this, &MeshData::z0);
  new ClassDouble<MeshData>(ca, "Zmax", this, &MeshData::zmax);
  new ClassInt<MeshData>(ca, "NumberOfCellsX", this, &MeshData::Nx);
  new ClassInt<MeshData>(ca, "NumberOfCellsY", this, &MeshData::Ny);
  new ClassInt<MeshData>(ca, "NumberOfCellsZ", this, &MeshData::Nz);

  new ClassToken<MeshData>(ca, "BoundaryConditionX0", this,
                               reinterpret_cast<int MeshData::*>(&MeshData::bc_x0), 5,
                               "None", 0, "Inlet", 1, "Outlet", 2, "Wall", 3, "Symmetry", 4);
  new ClassToken<MeshData>(ca, "BoundaryConditionXmax", this,
                               reinterpret_cast<int MeshData::*>(&MeshData::bc_xmax), 5,
                               "None", 0, "Inlet", 1, "Outlet", 2, "Wall", 3, "Symmetry", 4);
  new ClassToken<MeshData>(ca, "BoundaryConditionY0", this,
                               reinterpret_cast<int MeshData::*>(&MeshData::bc_y0), 5,
                               "None", 0, "Inlet", 1, "Outlet", 2, "Wall", 3, "Symmetry", 4);
  new ClassToken<MeshData>(ca, "BoundaryConditionYmax", this,
                               reinterpret_cast<int MeshData::*>(&MeshData::bc_ymax), 5,
                               "None", 0, "Inlet", 1, "Outlet", 2, "Wall", 3, "Symmetry", 4);
  new ClassToken<MeshData>(ca, "BoundaryConditionZ0", this,
                               reinterpret_cast<int MeshData::*>(&MeshData::bc_z0), 5,
                               "None", 0, "Inlet", 1, "Outlet", 2, "Wall", 3, "Symmetry", 4);
  new ClassToken<MeshData>(ca, "BoundaryConditionZmax", this,
                               reinterpret_cast<int MeshData::*>(&MeshData::bc_zmax), 5,
                               "None", 0, "Inlet", 1, "Outlet", 2, "Wall", 3, "Symmetry", 4);
 } 

//------------------------------------------------------------------------------

StiffenedGasModelData::StiffenedGasModelData()
{

  specificHeatRatio = 1.4;
  idealGasConstant = 287.1;
  specificHeatPressure = -1.0;
  pressureConstant = 0.0;

}

//------------------------------------------------------------------------------

void StiffenedGasModelData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 4, father);

  new ClassDouble<StiffenedGasModelData>(ca, "SpecificHeatRatio", this,
                                &StiffenedGasModelData::specificHeatRatio);
  new ClassDouble<StiffenedGasModelData>(ca, "IdealGasConstant", this,
                                &StiffenedGasModelData::idealGasConstant);
  new ClassDouble<StiffenedGasModelData>(ca, "SpecificHeatAtConstantPressure", this,
                                &StiffenedGasModelData::specificHeatPressure);
  new ClassDouble<StiffenedGasModelData>(ca, "PressureConstant", this,
                                &StiffenedGasModelData::pressureConstant);

}

//------------------------------------------------------------------------------

MieGruneisenModelData::MieGruneisenModelData() 
{
  //default values are for copper
  rho0 = 8.96e-3;       // unit: g/mm3
  cv = 3.90e8;          // unit: mm2/(s2.K)
  C0 = 3.933e6;         // unit: mm/s
  s = 1.5;              // non-dimensional
  Gamma0 = 1.99;        // non-dimensional

}

//------------------------------------------------------------------------------

void MieGruneisenModelData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 5, father);

  new ClassDouble<MieGruneisenModelData>(ca, "ReferenceDensity", this, 
                                         &MieGruneisenModelData::rho0);
  new ClassDouble<MieGruneisenModelData>(ca, "SpecificHeatAtConstantVolume", this, 
                                         &MieGruneisenModelData::cv);
  new ClassDouble<MieGruneisenModelData>(ca, "BulkSpeedOfSound", this, 
                                         &MieGruneisenModelData::C0);
  new ClassDouble<MieGruneisenModelData>(ca, "HugoniotSlope", this, 
                                         &MieGruneisenModelData::s);
  new ClassDouble<MieGruneisenModelData>(ca, "Gamma0", this, 
                                         &MieGruneisenModelData::Gamma0);

}

//------------------------------------------------------------------------------

MaterialModelData::MaterialModelData()
{

  eos = STIFFENED_GAS;
  rhomin = -DBL_MAX; // By default, no clipping
  pmin = -DBL_MAX;   // By default, no clipping

}

//------------------------------------------------------------------------------

Assigner *MaterialModelData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 5, nullAssigner);

  new ClassToken<MaterialModelData>(ca, "EquationOfState", this,
                                 reinterpret_cast<int MaterialModelData::*>(&MaterialModelData::eos), 2,
                                 "StiffenedGas", MaterialModelData::STIFFENED_GAS, 
                                 "MieGruneisen", MaterialModelData::MIE_GRUNEISEN);
  new ClassDouble<MaterialModelData>(ca, "DensityCutOff", this, &MaterialModelData::rhomin);
  new ClassDouble<MaterialModelData>(ca, "PressureCutOff", this, &MaterialModelData::pmin);

  sgModel.setup("StiffenedGasModel", ca);
  mgModel.setup("MieGruneisenModel", ca);

  return ca;
};

//------------------------------------------------------------------------------

EquationsData::EquationsData()
{

}

//------------------------------------------------------------------------------

void EquationsData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 1, father); 

  materials.setup("Material", ca);
}

//------------------------------------------------------------------------------

ReconstructionData::ReconstructionData() 
{
  type = LINEAR;
  limiter = GENERALIZED_MINMOD;
  generalized_minmod_coeff = 2.0; //The MC Limiter
}

//------------------------------------------------------------------------------

void ReconstructionData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 3, father); 

  new ClassToken<ReconstructionData>
    (ca, "Type", this,
     reinterpret_cast<int ReconstructionData::*>(&ReconstructionData::type), 2,
     "Constant", 0, "Linear", 1);

  new ClassToken<ReconstructionData>
    (ca, "Limiter", this,
     reinterpret_cast<int ReconstructionData::*>(&ReconstructionData::limiter), 3,
     "None", 0, "GeneralizedMinMod", 1, "VanAlbada", 2);

  new ClassDouble<ReconstructionData>(ca, "GeneralizedMinModCoefficient", this, 
    &ReconstructionData::generalized_minmod_coeff);
}

//------------------------------------------------------------------------------

SchemeData::SchemeData() 
{
  flux = ROE;
  delta = 0.2; //the coefficient in Harten's entropy fix.
}

//------------------------------------------------------------------------------

void SchemeData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner* ca;
  ca = new ClassAssigner(name, 3, father);

  new ClassToken<SchemeData>
    (ca, "Flux", this,
     reinterpret_cast<int SchemeData::*>(&SchemeData::flux), 3,
     "Roe", 0, "LocalLaxFriedrichs", 1, "HLLC", 2);

  new ClassDouble<SchemeData>(ca, "EntropyFixCoefficient", this, &SchemeData::delta);

  rec.setup("Reconstruction", ca);

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

LevelSetSchemeData::LevelSetSchemeData() 
{
  materialid = -1;

  flux = ROE;

  delta = 0.2; //the coefficient in Harten's entropy fix.

  reinitialization_freq = -1;
}

//------------------------------------------------------------------------------

Assigner *LevelSetSchemeData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 5, nullAssigner);

  new ClassInt<LevelSetSchemeData>(ca, "MaterialID", this, 
    &LevelSetSchemeData::materialid);

  new ClassToken<LevelSetSchemeData>
    (ca, "Flux", this,
     reinterpret_cast<int LevelSetSchemeData::*>(&LevelSetSchemeData::flux), 3,
     "Roe", 0, "LocalLaxFriedrichs", 1, "Upwind", 2);

  new ClassDouble<LevelSetSchemeData>(ca, "EntropyFixCoefficient", this, &LevelSetSchemeData::delta);

  new ClassInt<LevelSetSchemeData>(ca, "ReinitializationFrequency", this, 
    &LevelSetSchemeData::reinitialization_freq);

  rec.setup("Reconstruction", ca);

  return ca;
}

//------------------------------------------------------------------------------



//------------------------------------------------------------------------------

SchemesData::SchemesData() 
{

}

//------------------------------------------------------------------------------

void SchemesData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 3, father);

  ns.setup("NavierStokes", ca);

  bc.setup("Boundaries", ca);

  ls.setup("LevelSet", ca);
}

//------------------------------------------------------------------------------

ExactRiemannSolverData::ExactRiemannSolverData()
{
  maxIts_main = 200;
  maxIts_shock = 200;
  numSteps_rarefaction = 100;
  tol_main = 1.0e-5;
  tol_shock = 1.0e-6;
  tol_rarefaction = 1.0e-5;
}

//------------------------------------------------------------------------------

void ExactRiemannSolverData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 6, father);

  new ClassInt<ExactRiemannSolverData>(ca, "MaxIts", this, 
                                       &ExactRiemannSolverData::maxIts_main);

  new ClassInt<ExactRiemannSolverData>(ca, "MaxItsShock", this, 
                                       &ExactRiemannSolverData::maxIts_shock);

  new ClassInt<ExactRiemannSolverData>(ca, "IntegrationSteps", this, 
                                       &ExactRiemannSolverData::numSteps_rarefaction);

  new ClassDouble<ExactRiemannSolverData>(ca, "Tolerance", this, 
                                          &ExactRiemannSolverData::tol_main);

  new ClassDouble<ExactRiemannSolverData>(ca, "ToleranceShock", this, 
                                          &ExactRiemannSolverData::tol_shock);

  new ClassDouble<ExactRiemannSolverData>(ca, "ToleranceRarefaction", this, 
                                          &ExactRiemannSolverData::tol_rarefaction);
}

//------------------------------------------------------------------------------

MultiPhaseData::MultiPhaseData()
{
  flux = EXACT;
}

//------------------------------------------------------------------------------

void MultiPhaseData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 1, father);

  new ClassToken<MultiPhaseData>
    (ca, "Flux", this,
     reinterpret_cast<int MultiPhaseData::*>(&MultiPhaseData::flux), 2,
     "Exact", 0, "Numerical", 1);

}

//------------------------------------------------------------------------------

ExplicitData::ExplicitData()
{
  type = RUNGE_KUTTA_2;
}

//------------------------------------------------------------------------------

void ExplicitData::setup(const char *name, ClassAssigner *father)
{

 ClassAssigner *ca = new ClassAssigner(name, 1, father);

  new ClassToken<ExplicitData>
    (ca, "Type", this,
     reinterpret_cast<int ExplicitData::*>(&ExplicitData::type), 3,
     "ForwardEuler", 0, "RungeKutta2", 1, "RungeKutta3", 2);

}

//------------------------------------------------------------------------------

TsData::TsData()
{

  type = EXPLICIT;
  maxIts = INT_MAX;
  timestep = -1.0;
  cfl = 0.5;
  maxTime = 1e6;

}

//------------------------------------------------------------------------------

void TsData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 6, father);

  new ClassToken<TsData>(ca, "Type", this,
                         reinterpret_cast<int TsData::*>(&TsData::type), 2,
                         "Explicit", 0, "Implicit", 1);
  new ClassInt<TsData>(ca, "MaxIts", this, &TsData::maxIts);
  new ClassDouble<TsData>(ca, "TimeStep", this, &TsData::timestep);
  new ClassDouble<TsData>(ca, "CFL", this, &TsData::cfl);
  new ClassDouble<TsData>(ca, "MaxTime", this, &TsData::maxTime);

  expl.setup("Explicit", ca);
}

//------------------------------------------------------------------------------

BcsData::BcsData()
{

}

//------------------------------------------------------------------------------

void BcsData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 3, father);

  inlet.setup("Inlet", ca);
  outlet.setup("Outlet", ca);
  wall.setup("Wall", ca);

}

//------------------------------------------------------------------------------

BcsWallData::BcsWallData()
{

  type = ADIABATIC;
  temperature = -1.0;

}

//------------------------------------------------------------------------------

void BcsWallData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 2, father);

  new ClassToken<BcsWallData>(ca, "Type", this,
                              reinterpret_cast<int BcsWallData::*>(&BcsWallData::type), 2,
                              "Isothermal", 0, "Adiabatic", 1);
  new ClassDouble<BcsWallData>(ca, "Temperature", this, &BcsWallData::temperature);

}

//------------------------------------------------------------------------------

IcData::IcData()
{
  user_specified_ic = "";

  type = NONE;

  for(int i=0; i<SIZE; i++)
    specified[i] = 0;
}

//------------------------------------------------------------------------------

void IcData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 2, father);

  new ClassStr<IcData>(ca, "UserDataFile", this, &IcData::user_specified_ic);

  multiInitialConditions.setup("GeometricEntities");
}

//------------------------------------------------------------------------------

void IcData::readUserSpecifiedIC()
{
  if(!strcmp(user_specified_ic, "")) //no user_specified_ic
    return;

  std::fstream input;
  input.open(user_specified_ic, std::fstream::in);
  if (!input.is_open()) {
    print_error("Error: could not open user-specified initial condition file %s.\n", user_specified_ic);
    exit_mpi();
  } else
    print("- Reading user-specified initial condition file: %s.\n", user_specified_ic);

  std::string word, line;

  // Read the first line of user-specified file
  input.ignore(2,' '); //This line should start with ## (so the file can be plotted
                       //directly using gnuplot). It must then specify the type of initial cond.

  input >> word;
  // just comparing the first 4 letters is sufficient
  if(!(word.compare(0,4,"Planar",0,4) && 
       word.compare(0,4,"PLANAR",0,4) && 
       word.compare(0,4,"planar",0,4))) {
    type = PLANAR;
    input.ignore(256,'\n'); //done with line 1
    readUserSpecifiedIC_Planar(input);
  } 
  else if(!(word.compare(0,4,"Cylindrical",0,4) && 
            word.compare(0,4,"CYLINDRICAL",0,4) && 
            word.compare(0,4,"cylindrical",0,4))) {
    type = CYLINDRICAL;
    input.ignore(256,'\n'); //done with line 1
    readUserSpecifiedIC_Cylindrical(input);
  } 
  else if(!(word.compare(0,4,"Spherical",0,4) && 
            word.compare(0,4,"SPHERICAL",0,4) && 
            word.compare(0,4,"spherical",0,4))) {
    type = SPHERICAL;
    input.ignore(256,'\n'); //done with line 1
    readUserSpecifiedIC_Spherical(input);
  } 
  else {
    print_error("Error: Unknown initial condition type %s.\n", word);
    exit_mpi();
  }

  input.close();
  //print("Read user-specified initial condition.\n");
}

//------------------------------------------------------------------------------

// Read planar initial condition from file
void IcData::readUserSpecifiedIC_Planar(std::fstream &input)
{
  std::string word, line;

  // Read the second line of user-specified file
  input.ignore(2,' '); //This line should start with ## 
                       //It must then contain 3 real numbers corresponding to the (x,y,z) coordinates
                       //of the "0" in this data file within the actual mesh
  input >> x0[0] >> x0[1] >> x0[2];
  input.ignore(256,'\n'); //done with line 2

  // Read the next line
  input.ignore(2,' '); //This line should start with ##
                       //It must contain 3 real numbers corresponding to the direction of "+x" in the 1D data
  input >> dir[0] >> dir[1] >> dir[2];
  dir /= dir.norm();
  input.ignore(256,'\n'); //done with this line

  // Read the next line, which specifies the variable of each column
  input.ignore(2,' '); //This line should start with ##
                       //It must contain between 2 and SIZE words indicating the field
  for(int i=0; i<SIZE; i++)
    specified[i] = 0;
  std::map<int,int> column2var; //maps the column number to "Vars" index 
  getline(input, line);

  std::istringstream iss(line);
  int column = 0; 
  for(;iss>>word;) {
    // just check the first four letters
    if(!(word.compare(0,4,"Coordinate",0,4) && 
         word.compare(0,4,"COORDINATE",0,4) && 
         word.compare(0,4,"coordinate",0,4))) {
      column2var[column] = COORDINATE;
      specified[COORDINATE] = 1;
    } 
    else if(!(word.compare(0,4,"Density",0,4) && 
              word.compare(0,4,"DENSITY",0,4) && 
              word.compare(0,4,"density",0,4))) {
      column2var[column] = DENSITY;
      specified[DENSITY] = 1;
    } 
    else if(!(word.compare(0,4,"Velocity",0,4) && 
              word.compare(0,4,"VELOCITY",0,4) && 
              word.compare(0,4,"velocity",0,4))) {
      column2var[column] = VELOCITY;
      specified[VELOCITY] = 1;
    } 
    else if(!(word.compare(0,4,"Pressure",0,4) && 
              word.compare(0,4,"PRESSURE",0,4) && 
              word.compare(0,4,"pressure",0,4))) {
      column2var[column] = PRESSURE;
      specified[PRESSURE] = 1;
    } 
    else if(!(word.compare(0,4,"LevelSet",0,4) && 
              word.compare(0,4,"LEVELSET",0,4) && 
              word.compare(0,4,"levelset",0,4))) {
      column2var[column] = LEVELSET;
      specified[LEVELSET] = 1;
    } 
    else if(!(word.compare(0,4,"MaterialID",0,4) && 
              word.compare(0,4,"MATERIALID",0,4) && 
              word.compare(0,4,"materialid",0,4))) {
      column2var[column] = MATERIALID;
      specified[MATERIALID] = 1;
    } 
    else if(!(word.compare(0,4,"Temperature",0,4) && 
              word.compare(0,4,"TEMPERATURE",0,4) && 
              word.compare(0,4,"temperature",0,4))) {
      column2var[column] = TEMPERATURE;
      specified[TEMPERATURE] = 1;
    } else {
      print_error("Error: I do not understand the word '%s' in the user-specified initial condition file.\n", word.c_str());
      exit_mpi();
    }
    column++;
  }

  if(column<2 || !specified[COORDINATE]) {
    print_error("Error: Need additional data in the initial condition file.\n");
    exit_mpi();
  }

  // Now start reading the actual data (until end of file)
  int MaxRows = 100000; //should be a lot more than enough
  double data;
  for(int r=0; r<MaxRows; r++) {
    getline(input, line);
    std::istringstream is(line);
    for(int i=0; i<column; i++) {
      is >> data; 
      if(is.fail())
        break;
      //cout << "row " << r << ", column " << i << ": " << data << endl;
      user_data[column2var[i]].push_back(data);
    }
    if(input.eof())
      break;
  }

}

//------------------------------------------------------------------------------

// Read cylindrical initial condition from file
void IcData::readUserSpecifiedIC_Cylindrical(std::fstream &input)
{
  std::string word, line;

  // Read the second line of user-specified file
  input.ignore(2,' '); //This line should start with ## 
                       //It must then contain 3 real numbers corresponding to the (x,y,z) coordinates
                       //of the "0" in this data file within the actual mesh
  input >> x0[0] >> x0[1] >> x0[2];
  input.ignore(256,'\n'); //done with line 2

  // Read the next line
  input.ignore(2,' '); //This line should start with ##
                       //It must contain 3 real numbers corresponding to the positive direction of axis of symmetry
  input >> dir[0] >> dir[1] >> dir[2];
  dir /= dir.norm();
  input.ignore(256,'\n'); //done with this line

  // Read the next line, which specifies the variable of each column
  input.ignore(2,' '); //This line should start with ##
                       //It must contain between 2 and SIZE words indicating the field
  for(int i=0; i<SIZE; i++)
    specified[i] = 0;
  std::map<int,int> column2var; //maps the column number to "Vars" index 
  getline(input, line);

  std::istringstream iss(line);
  int column = 0; 
  for(;iss>>word;) {
    // just check the first four letters
    if(!(word.compare(0,4,"Coordinate",0,4) && 
         word.compare(0,4,"COORDINATE",0,4) && 
         word.compare(0,4,"coordinate",0,4))) {
      column2var[column] = COORDINATE;
      specified[COORDINATE] = 1;
    } 
    else if(!(word.compare(0,4,"Density",0,4) && 
              word.compare(0,4,"DENSITY",0,4) && 
              word.compare(0,4,"density",0,4))) {
      column2var[column] = DENSITY;
      specified[DENSITY] = 1;
    } 
    else if(!(word.compare(0,4,"Velocity",0,4) && 
              word.compare(0,4,"VELOCITY",0,4) && 
              word.compare(0,4,"velocity",0,4))) {
      column2var[column] = VELOCITY;
      specified[VELOCITY] = 1;
    } 
    else if(!(word.compare(0,4,"Pressure",0,4) && 
              word.compare(0,4,"PRESSURE",0,4) && 
              word.compare(0,4,"pressure",0,4))) {
      column2var[column] = PRESSURE;
      specified[PRESSURE] = 1;
    } 
    else if(!(word.compare(0,4,"LevelSet",0,4) && 
              word.compare(0,4,"LEVELSET",0,4) && 
              word.compare(0,4,"levelset",0,4))) {
      column2var[column] = LEVELSET;
      specified[LEVELSET] = 1;
    } 
    else if(!(word.compare(0,4,"MaterialID",0,4) && 
              word.compare(0,4,"MATERIALID",0,4) && 
              word.compare(0,4,"materialid",0,4))) {
      column2var[column] = MATERIALID;
      specified[MATERIALID] = 1;
    } 
    else if(!(word.compare(0,4,"Temperature",0,4) && 
              word.compare(0,4,"TEMPERATURE",0,4) && 
              word.compare(0,4,"temperature",0,4))) {
      column2var[column] = TEMPERATURE;
      specified[TEMPERATURE] = 1;
    } else {
      print_error("Error: I do not understand the word '%s' in the user-specified initial condition file.\n", word.c_str());
      exit_mpi();
    }
    column++;
  }

  if(column<2 || !specified[COORDINATE]) {
    print_error("Error: Need additional data in the initial condition file.\n");
    exit_mpi();
  }


  // Read the next line: should start with ##, then say "AXIAL"
  input.ignore(2,' '); 
  input >> word;
  if(word.compare(0,4,"Axial",0,4) && 
     word.compare(0,4,"AXIAL",0,4) && 
     word.compare(0,4,"axial",0,4)) {
    print_error("Error: Expect keyword 'Axial' in user-specified initial condition.\n");
    exit_mpi();
  } 
  input.ignore(256,'\n'); //done with this line


  // Now start reading the data in the axial direction
  int MaxRows = 100000; //should be a lot more than enough
  double data;
  bool found_radial = false;
  for(int r=0; r<MaxRows; r++) {

    getline(input, line);

    // check if we got the line that says "## Radial"
    std::istringstream check(line);
    check >> word; check >> word;
    if(!(word.compare(0,4,"Radial",0,4) && 
         word.compare(0,4,"RADIAL",0,4) && 
         word.compare(0,4,"radial",0,4))) {
      found_radial = true;
      break; 
    }

    std::istringstream is(line);
    for(int i=0; i<column; i++) {
      is >> data; 
      if(is.fail())
        break;
      //cout << "row " << r << ", column " << i << ": " << data << endl;
      user_data[column2var[i]].push_back(data);
    }

    if(input.eof())
      break;
  }

  if(found_radial) { //read radial variation
    for(int r=0; r<MaxRows; r++) {
      getline(input, line);

      std::istringstream is(line);
      for(int i=0; i<column; i++) {
        is >> data; 
        if(is.fail())
          break;
        //cout << "row " << r << ", column " << i << ": " << data << endl;
        user_data2[column2var[i]].push_back(data);
      }

      if(input.eof())
        break;
    }
  }

}

//------------------------------------------------------------------------------

// Read spherical initial condition from file
// Note: This is almost identical to reading planar data, except the variable "dir" does not need be read
void IcData::readUserSpecifiedIC_Spherical(std::fstream &input)
{
  std::string word, line;

  // Read the second line of user-specified file
  input.ignore(2,' '); //This line should start with ## 
                       //It must then contain 3 real numbers corresponding to the (x,y,z) coordinates
                       //of the "0" in this data file within the actual mesh
  input >> x0[0] >> x0[1] >> x0[2];
  input.ignore(256,'\n'); //done with line 2

  // Read the next line, which specifies the variable of each column
  input.ignore(2,' '); //This line should start with ##
                       //It must contain between 2 and SIZE words indicating the field
  for(int i=0; i<SIZE; i++)
    specified[i] = 0;
  std::map<int,int> column2var; //maps the column number to "Vars" index 
  getline(input, line);

  std::istringstream iss(line);
  int column = 0; 
  for(;iss>>word;) {
    // just check the first four letters
    if(!(word.compare(0,4,"Coordinate",0,4) && 
         word.compare(0,4,"COORDINATE",0,4) && 
         word.compare(0,4,"coordinate",0,4))) {
      column2var[column] = COORDINATE;
      specified[COORDINATE] = 1;
    } 
    else if(!(word.compare(0,4,"Density",0,4) && 
              word.compare(0,4,"DENSITY",0,4) && 
              word.compare(0,4,"density",0,4))) {
      column2var[column] = DENSITY;
      specified[DENSITY] = 1;
    } 
    else if(!(word.compare(0,4,"Velocity",0,4) && 
              word.compare(0,4,"VELOCITY",0,4) && 
              word.compare(0,4,"velocity",0,4))) {
      column2var[column] = VELOCITY;
      specified[VELOCITY] = 1;
    } 
    else if(!(word.compare(0,4,"Pressure",0,4) && 
              word.compare(0,4,"PRESSURE",0,4) && 
              word.compare(0,4,"pressure",0,4))) {
      column2var[column] = PRESSURE;
      specified[PRESSURE] = 1;
    } 
    else if(!(word.compare(0,4,"LevelSet",0,4) && 
              word.compare(0,4,"LEVELSET",0,4) && 
              word.compare(0,4,"levelset",0,4))) {
      column2var[column] = LEVELSET;
      specified[LEVELSET] = 1;
    } 
    else if(!(word.compare(0,4,"MaterialID",0,4) && 
              word.compare(0,4,"MATERIALID",0,4) && 
              word.compare(0,4,"materialid",0,4))) {
      column2var[column] = MATERIALID;
      specified[MATERIALID] = 1;
    } 
    else if(!(word.compare(0,4,"Temperature",0,4) && 
              word.compare(0,4,"TEMPERATURE",0,4) && 
              word.compare(0,4,"temperature",0,4))) {
      column2var[column] = TEMPERATURE;
      specified[TEMPERATURE] = 1;
    } else {
      print_error("Error: I do not understand the word '%s' in the user-specified initial condition file.\n", word.c_str());
      exit_mpi();
    }
    column++;
  }

  if(column<2 || !specified[COORDINATE]) {
    print_error("Error: Need additional data in the initial condition file.\n");
    exit_mpi();
  }

  // Now start reading the actual data (until end of file)
  int MaxRows = 100000; //should be a lot more than enough
  double data;
  for(int r=0; r<MaxRows; r++) {
    getline(input, line);
    std::istringstream is(line);
    for(int i=0; i<column; i++) {
      is >> data; 
      if(is.fail())
        break;
      //cout << "row " << r << ", column " << i << ": " << data << endl;
      user_data[column2var[i]].push_back(data);
    }
    if(input.eof())
      break;
  }

}
//------------------------------------------------------------------------------

OutputData::OutputData()
{
  prefix = "";
  solution_filename_base = "solution";

  frequency = 0;
  frequency_dt = -1.0;

  density = OFF;
  velocity = OFF;
  pressure = OFF;
  materialid = OFF;
  temperature = OFF;
  levelset0 = OFF;
  levelset1 = OFF;
  levelset2 = OFF;
  levelset3 = OFF;
  levelset4 = OFF;

  for(int i=0; i<MAXLS; i++)
    levelset[i] = OFF;

  verbose = OFF;
}

//------------------------------------------------------------------------------

void OutputData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 10+MAXLS, father);

  new ClassStr<OutputData>(ca, "Prefix", this, &OutputData::prefix);
  new ClassStr<OutputData>(ca, "Solution", this, &OutputData::solution_filename_base);

  new ClassInt<OutputData>(ca, "Frequency", this, &OutputData::frequency);
  new ClassDouble<OutputData>(ca, "TimeInterval", this, &OutputData::frequency_dt);

  new ClassToken<OutputData>(ca, "Density", this,
                               reinterpret_cast<int OutputData::*>(&OutputData::density), 2,
                               "Off", 0, "On", 1);
  new ClassToken<OutputData>(ca, "Velocity", this,
                               reinterpret_cast<int OutputData::*>(&OutputData::velocity), 2,
                               "Off", 0, "On", 1);
  new ClassToken<OutputData>(ca, "Pressure", this,
                               reinterpret_cast<int OutputData::*>(&OutputData::pressure), 2,
                               "Off", 0, "On", 1);
  new ClassToken<OutputData>(ca, "MaterialID", this,
                               reinterpret_cast<int OutputData::*>(&OutputData::materialid), 2,
                               "Off", 0, "On", 1);
  new ClassToken<OutputData>(ca, "Temperature", this,
                               reinterpret_cast<int OutputData::*>(&OutputData::temperature), 2,
                               "Off", 0, "On", 1);

  new ClassToken<OutputData>(ca, "LevelSet0", this,
                               reinterpret_cast<int OutputData::*>(&OutputData::levelset0), 2,
                               "Off", 0, "On", 1);
  new ClassToken<OutputData>(ca, "LevelSet1", this,
                               reinterpret_cast<int OutputData::*>(&OutputData::levelset1), 2,
                               "Off", 0, "On", 1);
  new ClassToken<OutputData>(ca, "LevelSet2", this,
                               reinterpret_cast<int OutputData::*>(&OutputData::levelset2), 2,
                               "Off", 0, "On", 1);
  new ClassToken<OutputData>(ca, "LevelSet3", this,
                               reinterpret_cast<int OutputData::*>(&OutputData::levelset3), 2,
                               "Off", 0, "On", 1);
  new ClassToken<OutputData>(ca, "LevelSet4", this,
                               reinterpret_cast<int OutputData::*>(&OutputData::levelset4), 2,
                               "Off", 0, "On", 1);


  new ClassToken<OutputData>(ca, "VerboseScreenOutput", this,
                               reinterpret_cast<int OutputData::*>(&OutputData::verbose), 2,
                               "Off", 0, "On", 1);
}

//------------------------------------------------------------------------------

IoData::IoData(int argc, char** argv)
{
  readCmdLine(argc, argv);
  readCmdFile();
}

//------------------------------------------------------------------------------

void IoData::readCmdLine(int argc, char** argv)
{
  if(argc==1) {
    print_error("Error: Input file not provided!\n");
    exit_mpi();
  }
  cmdFileName = argv[1];
}

//------------------------------------------------------------------------------

void IoData::readCmdFile()
{
  extern FILE *yyCmdfin;
  extern int yyCmdfparse();

  setupCmdFileVariables();
//  cmdFilePtr = freopen(cmdFileName, "r", stdin);
  yyCmdfin = cmdFilePtr = fopen(cmdFileName, "r");

  if (!cmdFilePtr) {
    print_error("*** Error: could not open \'%s\'\n", cmdFileName);
    exit_mpi();
  }

  int error = yyCmdfparse();
  if (error) {
    print_error("*** Error: command file contained parsing errors\n");
    exit(error);
  }
  fclose(cmdFilePtr);

  //READ ADDITIONAL FILES
  if(strcmp(ic.user_specified_ic, ""))
    ic.readUserSpecifiedIC(); 

  //FIX Levelset output (TODO: need a better way...)
  output.levelset[0] = output.levelset0;
  output.levelset[1] = output.levelset1;
  output.levelset[2] = output.levelset2;
  output.levelset[3] = output.levelset3;
  output.levelset[4] = output.levelset4;
}

//------------------------------------------------------------------------------

void IoData::setupCmdFileVariables()
{
  eqs.setup("Equations");
  ic.setup("InitialCondition");
  bc.setup("BoundaryConditions");

  mesh.setup("Mesh");

  schemes.setup("Space");

  exact_riemann.setup("ExactRiemannSolution");

  ts.setup("Time");

  multiphase.setup("MultiPhase");

  output.setup("Output");
}

//------------------------------------------------------------------------------
