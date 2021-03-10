#ifndef _IO_DATA_H_
#define _IO_DATA_H_

#include <cstdio>
#include <map>
#include "parser/Assigner.h"
#include "parser/Dictionary.h"
#include <Vector3D.h>
#include <Utils.h>

using std::map;
using std::pair;
using std::vector;

/*********************************************************************
 * class IoData reads and processes the input data provided by the user
 *********************************************************************
*/
//------------------------------------------------------------------------------

template<class DataType>
class ObjectMap {

public:
  map<int, DataType *> dataMap;

  void setup(const char *name, ClassAssigner *p) {
    SysMapObj<DataType> *smo = new SysMapObj<DataType>(name, &dataMap);
    if (p) p->addSmb(name, smo);
    else addSysSymbol(name, smo);
  }

  ~ObjectMap() {
    for(typename map<int, DataType *>::iterator it=dataMap.begin();it!=dataMap.end();++it)
      delete it->second;
  }
};

//------------------------------------------------------------------------------

struct StateVariable {

  int    materialid;
  double density;
  double velocity_x;
  double velocity_y;
  double velocity_z;
  double pressure;
  double temperature;

  StateVariable();
  ~StateVariable() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct PointData {

  double x,y,z;
  StateVariable initialConditions;

  PointData();
  ~PointData() {}
  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct PlaneData {

  double cen_x, cen_y, cen_z, nx, ny, nz;
  StateVariable initialConditions;

  PlaneData();
  ~PlaneData() {}
  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct SphereData {

  double cen_x, cen_y, cen_z, radius;
  StateVariable initialConditions;

  SphereData();
  ~SphereData() {}
  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct CylinderConeData {

  //! info about the cylinder
  double cen_x, cen_y, cen_z, nx, ny, nz, r, L;
  //! info about the cone (connected to the top of the cylinder)
  double opening_angle_degrees, cone_height;
 
  StateVariable initialConditions;

  CylinderConeData();
  ~CylinderConeData() {}
  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct MultiInitialConditionsData {

  ObjectMap<PointData>    pointMap;
  ObjectMap<PlaneData>    planeMap;
  ObjectMap<SphereData>   sphereMap;
  ObjectMap<CylinderConeData> cylinderconeMap;

  void setup(const char *, ClassAssigner * = 0);
};

//------------------------------------------------------------------------------

struct MeshData {

  enum Type {THREEDIMENSIONAL = 0, CYLINDRICAL = 1} type;
  double x0, xmax, y0, ymax, z0, zmax;
  int Nx, Ny, Nz;

  enum BcType {NONE = 0, INLET = 1, OUTLET = 2, WALL = 3, SYMMETRY = 4, SIZE = 5};
  BcType bc_x0, bc_xmax, bc_y0, bc_ymax, bc_z0, bc_zmax;

  MeshData();
  ~MeshData() {} 

  void setup(const char *, ClassAssigner * = 0);
};

//------------------------------------------------------------------------------

struct StiffenedGasModelData {

  double specificHeatRatio;
  double idealGasConstant;
  double pressureConstant;
  double specificHeatPressure;

  StiffenedGasModelData();
  ~StiffenedGasModelData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct MieGruneisenModelData {

  double rho0;
  double cv; 
  double C0;
  double s;
  double Gamma0;

  MieGruneisenModelData();
  ~MieGruneisenModelData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct MaterialModelData {

  int id;
  enum EOS {STIFFENED_GAS = 0, MIE_GRUNEISEN = 1} eos;
  double rhomin;
  double pmin;

  StiffenedGasModelData sgModel;
  MieGruneisenModelData mgModel;

  MaterialModelData();
  ~MaterialModelData() {}
  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct EquationsData {

  ObjectMap<MaterialModelData> materials;

  EquationsData();
  ~EquationsData() {}

  void setup(const char *, ClassAssigner * = 0);
};

//------------------------------------------------------------------------------

struct ReconstructionData {

  enum Reconstruction {CONSTANT = 0, LINEAR = 1} type;
  enum Limiter {NONE = 0, GENERALIZED_MINMOD = 1, VANALBADA = 2} limiter; 
  double generalized_minmod_coeff;

  ReconstructionData();
  ~ReconstructionData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct SchemeData {

  enum Flux {ROE = 0, LOCAL_LAX_FRIEDRICHS = 1, HLLC = 2} flux;
 
  double delta; //! The coeffient in Harten's entropy fix.

  ReconstructionData rec;

  SchemeData();
  ~SchemeData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct BoundarySchemeData {

  enum Type { STEGER_WARMING = 0,
              CONSTANT_EXTRAPOLATION = 1,
              LINEAR_EXTRAPOLATION = 2,
              GHIDAGLIA = 3, MODIFIED_GHIDAGLIA = 4} type;

  BoundarySchemeData();
  ~BoundarySchemeData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct LevelSetSchemeData {

  int materialid; //! The material in the phi<0 region ("inside")

  enum Flux {ROE = 0, LOCAL_LAX_FRIEDRICHS = 1} flux;

  ReconstructionData rec;
  
  double delta; //! The coeffient in Harten's entropy fix.

  int reinitialization_freq; 

  LevelSetSchemeData();
  ~LevelSetSchemeData() {}
  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct SchemesData {

  SchemeData ns;

  BoundarySchemeData bc;

  ObjectMap<LevelSetSchemeData> ls;

  SchemesData();
  ~SchemesData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct ExactRiemannSolverData {

  int maxIts_main;
  int maxIts_shock;
  int numSteps_rarefaction;
  double tol_main;
  double tol_shock;
  double tol_rarefaction;

  ExactRiemannSolverData();
  ~ExactRiemannSolverData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct MultiPhaseData {

  enum Flux {EXACT = 0, NUMERICAL = 1} flux;

  MultiPhaseData();
  ~MultiPhaseData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct ExplicitData {

//!time-integration scheme used
  enum Type {FORWARD_EULER = 0, RUNGE_KUTTA_2 = 1, RUNGE_KUTTA_3 = 2} type;

  ExplicitData();
  ~ExplicitData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct TsData {

  enum Type {EXPLICIT = 0, IMPLICIT = 1} type;
  int maxIts;
  double timestep;
  double cfl;
  double maxTime;
  ExplicitData expl;

  TsData();
  ~TsData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct BcsWallData {

  enum Type {ISOTHERMAL = 0, ADIABATIC = 1} type;
  double temperature;

  BcsWallData();
  ~BcsWallData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct BcsData {

  StateVariable inlet;
  StateVariable outlet;
  BcsWallData wall;

  BcsData();
  ~BcsData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct IcData {

  //! initial condition specified using simple geometric entities (e.g., point, plane,)
  MultiInitialConditionsData multiInitialConditions;

  //! user-specified file
  const char *user_specified_ic;

  enum Type {NONE = 0, PLANAR = 1, CYLINDRICAL = 2, SPHERICAL = 3} type;
  Vec3D x0;
  Vec3D dir; //!< relevant for PLANAR and CYLINDRICAL
  enum Vars {COORDINATE = 0, DENSITY = 1, VELOCITY = 2, PRESSURE = 3,
             LEVELSET = 4, MATERIALID = 5, TEMPERATURE = 6, SIZE = 7};

  int specified[SIZE];  //!< 0~unspecified, 1~specified

  vector<double> user_data[SIZE];

  vector<double> user_data2[SIZE]; //!< for radial variation 

  IcData();
  ~IcData() {}

  void setup(const char *, ClassAssigner * = 0);

  void readUserSpecifiedIC();
  void readUserSpecifiedIC_Planar(std::fstream &input);
  void readUserSpecifiedIC_Cylindrical(std::fstream &input);
  void readUserSpecifiedIC_Spherical(std::fstream &input);
};

//------------------------------------------------------------------------------

struct OutputData {

  const char *prefix; //!< path
  const char *solution_filename_base; //!< filename without path

  enum Options {OFF = 0, ON = 1};
  Options density, velocity, pressure, materialid, temperature;
  Options verbose;

  const static int MAXLS = 5;
  Options levelset[MAXLS];

  Options levelset0;
  Options levelset1;
  Options levelset2;
  Options levelset3;
  Options levelset4;

  int frequency;
  double frequency_dt; //!< -1 by default. To activate it, set it to a positive number

  OutputData();
  ~OutputData() {}

  void setup(const char *, ClassAssigner * = 0);
};

//------------------------------------------------------------------------------

class IoData {

  char *cmdFileName;
  FILE *cmdFilePtr;

public:

  MeshData mesh;

  EquationsData eqs;
  IcData ic;  //!< initial condition
  BcsData bc;  //!< boundary condition

  SchemesData schemes;

  ExactRiemannSolverData exact_riemann;

  MultiPhaseData multiphase;

  TsData ts;

  OutputData output;

public:

  IoData() {}
  IoData(int, char**);
  ~IoData() {}

  void readCmdLine(int, char**);
  void setupCmdFileVariables();
  void readCmdFile();

};
#endif
