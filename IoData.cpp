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
#include <cmath>
#include <unistd.h>
//#include <dlfcn.h>
using namespace std;

//------------------------------------------------------------------------------

MeshData::MeshData()
{
  type = CYLINDRICAL;
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
                               reinterpret_cast<int MeshData::*>(&MeshData::type), 3,
                               "TwoDimensional", 0, "Cylindrical", 1, "ThreeDimensional", 2);
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
  ClassAssigner *ca = new ClassAssigner(name, 2, father); //not used.

  fluid1.setup("FluidModel1");
}

//------------------------------------------------------------------------------

SchemeData::SchemeData(int af) : allowsFlux(af)
{

  flux = ROE;
  reconstruction = LINEAR;
  limiter = GENERALIZED_MINMOD;

  generalized_minmod_coeff = 2.0; //The MC Limiter
}

//------------------------------------------------------------------------------

void SchemeData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner* ca;
  if (allowsFlux)
    ca = new ClassAssigner(name, 4, father);
  else
    ca = new ClassAssigner(name, 3, father);

  if (allowsFlux) {
    new ClassToken<SchemeData>
      (ca, "Flux", this,
       reinterpret_cast<int SchemeData::*>(&SchemeData::flux), 4,
       "Roe", 0, "HLLE", 1, "HLLC", 2, "KurganovTadmor", 3);
  }

  new ClassToken<SchemeData>
    (ca, "Reconstruction", this,
     reinterpret_cast<int SchemeData::*>(&SchemeData::reconstruction), 2,
     "Constant", 0, "Linear", 1);

  new ClassToken<SchemeData>
    (ca, "Limiter", this,
     reinterpret_cast<int SchemeData::*>(&SchemeData::limiter), 4,
     "None", 0, "GeneralizedMinMod", 1, "VanAlbada", 2, "ModifiedVanAlbada", 3);

  new ClassDouble<SchemeData>(ca, "GeneralizedMinModCoefficient", this, 
    &SchemeData::generalized_minmod_coeff);

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

SchemesData::SchemesData() : ls(0), ns(1)
{

}

//------------------------------------------------------------------------------

void SchemesData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 3, father);

  ns.setup("NavierStokes", ca);
  ls.setup("LevelSet",ca);
  bc.setup("Boundaries", ca);

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

BcsFreeStreamData::BcsFreeStreamData()
{

  density = -1.0;
  velocity_x = 0.0;
  velocity_y = 0.0;
  velocity_z = 0.0;
  pressure = -1.0;
  temperature = -1.0;

}

//------------------------------------------------------------------------------

void BcsFreeStreamData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 6, father);

  new ClassDouble<BcsFreeStreamData>(ca, "Density", this, &BcsFreeStreamData::density);
  new ClassDouble<BcsFreeStreamData>(ca, "VelocityX", this, &BcsFreeStreamData::velocity_x);
  new ClassDouble<BcsFreeStreamData>(ca, "VelocityY", this, &BcsFreeStreamData::velocity_y);
  new ClassDouble<BcsFreeStreamData>(ca, "VelocityZ", this, &BcsFreeStreamData::velocity_z);
  new ClassDouble<BcsFreeStreamData>(ca, "Pressure", this, &BcsFreeStreamData::pressure);
  new ClassDouble<BcsFreeStreamData>(ca, "Temperature", this, &BcsFreeStreamData::temperature);

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
  user_specified_ic2 = "";

  type = NONE;

  for(int i=0; i<SIZE; i++)
    specified[i] = 0;
}

//------------------------------------------------------------------------------

void IcData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 2, father);

  new ClassStr<IcData>(ca, "UserDataFile", this, &IcData::user_specified_ic);
  new ClassStr<IcData>(ca, "UserDataFile2", this, &IcData::user_specified_ic2);
}

//------------------------------------------------------------------------------

void IcData::readUserSpecifiedIC()
{
  if(!strcmp(user_specified_ic, "")) //no user_specified_ic
    return;

  std::fstream input;
  input.open(user_specified_ic, std::fstream::in);
  if (!input.is_open()) {
    print_error("ERROR: could not open user-specified initial condition file %s.\n", user_specified_ic);
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
       word.compare(0,4,"planar",0,4)))
    type = PLANAR;
  else if(!(word.compare(0,4,"Cylindrical",0,4) && 
            word.compare(0,4,"CYLINDRICAL",0,4) && 
            word.compare(0,4,"cylindrical",0,4)))
    type = CYLINDRICAL;
  else if(!(word.compare(0,4,"Spherical",0,4) && 
            word.compare(0,4,"SPHERICAL",0,4) && 
            word.compare(0,4,"spherical",0,4)))
    type = SPHERICAL;
  else {
    print_error("ERROR: Unknown initial condition type %s.\n", word);
    exit_mpi();
  }
  input.ignore(256,'\n'); //done with line 1

  // Read the second line of user-specified file
  input.ignore(2,' '); //This line should start with ## 
                       //It must then contain 3 real numbers corresponding to the (x,y,z) coordinates
                       //of the "0" in this data file within the actual mesh
  input >> x0[0] >> x0[1] >> x0[2];
  input.ignore(256,'\n'); //done with line 2

  // Read the next line, which should be provided only if type is "planar" or "cylindrical"
  input.ignore(2,' '); //This line should start with ##
                       //It must contain 3 real numbers corresponding to the direction of "+x"
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
      print_error("ERROR: I do not understand the word '%s' in the user-specified initial condition file.\n", word.c_str());
      exit_mpi();
    }
    column++;
  }

  if(column<2 || !specified[COORDINATE]) {
    print_error("ERROR: Need additional data in the initial condition file.\n");
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

  input.close();
  //print("Read user-specified initial condition.\n");
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
  levelset = OFF;
  materialid = OFF;
  temperature = OFF;
}

//------------------------------------------------------------------------------

void OutputData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 10, father);

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
  new ClassToken<OutputData>(ca, "LevelSet", this,
                               reinterpret_cast<int OutputData::*>(&OutputData::levelset), 2,
                               "Off", 0, "On", 1);
  new ClassToken<OutputData>(ca, "MaterialID", this,
                               reinterpret_cast<int OutputData::*>(&OutputData::materialid), 2,
                               "Off", 0, "On", 1);
  new ClassToken<OutputData>(ca, "Temperature", this,
                               reinterpret_cast<int OutputData::*>(&OutputData::temperature), 2,
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
    print_error("ERROR: Input file not provided!\n");
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
}

//------------------------------------------------------------------------------

void IoData::setupCmdFileVariables()
{
  eqs.setup("Equations");
  ic.setup("InitialCondition");
  bc.setup("BoundaryConditions");

  mesh.setup("Mesh");

  schemes.setup("Space");
  ts.setup("Time");

  output.setup("Output");
}

//------------------------------------------------------------------------------
