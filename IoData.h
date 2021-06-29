#ifndef _IO_DATA_H_
#define _IO_DATA_H_

#include <cstdio>
#include <map>
#include "parser/Assigner.h"
#include "parser/Dictionary.h"
#include <Vector2D.h>
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

struct SpheroidData {

  double cen_x, cen_y, cen_z;
  double axis_x, axis_y, axis_z;
  double length, diameter;

  StateVariable initialConditions;

  SpheroidData();
  ~SpheroidData() {}
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
  ObjectMap<SpheroidData> spheroidMap;
  ObjectMap<CylinderConeData> cylinderconeMap;

  void setup(const char *, ClassAssigner * = 0);
};

//------------------------------------------------------------------------------

struct MeshResolution1DPointData {

  double coord;
  double h; 

  MeshResolution1DPointData();
  ~MeshResolution1DPointData() {}
  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct MeshData {

  enum Type {THREEDIMENSIONAL = 0, CYLINDRICAL = 1} type;
  double x0, xmax, y0, ymax, z0, zmax;
  int Nx, Ny, Nz;

  // mesh resolution info
  ObjectMap<MeshResolution1DPointData>  xpoints_map;
  ObjectMap<MeshResolution1DPointData>  ypoints_map;
  ObjectMap<MeshResolution1DPointData>  zpoints_map;

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
  double c0;
  double Gamma0;
  double s;
  double e0;

  double cv; 

  MieGruneisenModelData();
  ~MieGruneisenModelData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct JonesWilkinsLeeModelData {

  double omega;
  double A1;
  double A2;
  double R1;
  double R2;
  double rho0;

  JonesWilkinsLeeModelData();
  ~JonesWilkinsLeeModelData() {}

  void setup(const char *, ClassAssigner * = 0);  

};

//------------------------------------------------------------------------------

struct MaterialModelData {

  int id;
  enum EOS {STIFFENED_GAS = 0, MIE_GRUNEISEN = 1, JWL = 2} eos;
  double rhomin;
  double pmin;

  StiffenedGasModelData    sgModel;
  MieGruneisenModelData    mgModel;
  JonesWilkinsLeeModelData jwlModel;

  MaterialModelData();
  ~MaterialModelData() {}
  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct ViscosityModelData {

  enum Type {NONE = 0, CONSTANT = 0, SUTHERLAND = 1, ARTIFICIAL_RODIONOV = 2} type;

  // constant
  double dynamicViscosity;
  double bulkViscosity;

  // Sutherland
  double sutherlandConstant;
  double sutherlandReferenceTemperature;

  // Artificial viscosity (Rodionov)
  double Cav, Cth; 


  ViscosityModelData();

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct EquationsData {

  ObjectMap<MaterialModelData> materials;

  ViscosityModelData viscosity;

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

  enum Flux {ROE = 0, LOCAL_LAX_FRIEDRICHS = 1, UPWIND = 2} flux;

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

  enum ReconstructionAtInterface {CONSTANT = 0, LINEAR = 1} recon;

  enum PhaseChangeType {RIEMANN_SOLUTION = 0, EXTRAPOLATION = 1} phasechange_type;

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

struct DiskData {

  double cen_x, cen_y, cen_z, radius;
  double normal_x, normal_y, normal_z;
  StateVariable state;

  DiskData();
  ~DiskData() {}
  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct RectangleData {

  double cen_x, cen_y, cen_z, a, b;
  double normal_x, normal_y, normal_z;
  StateVariable state;

  RectangleData();
  ~RectangleData() {}
  Assigner *getAssigner();

};

//------------------------------------------------------------------------------


struct MultiBoundaryConditionsData {

  ObjectMap<DiskData>      diskMap;
  ObjectMap<RectangleData> rectangleMap;

  void setup(const char *, ClassAssigner * = 0);
};

//------------------------------------------------------------------------------

struct BcsData {

  StateVariable inlet;
  StateVariable outlet;
  BcsWallData wall;

  MultiBoundaryConditionsData multiBoundaryConditions;

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

  enum Type {NONE = 0, PLANAR = 1, CYLINDRICAL = 2, SPHERICAL = 3, 
             GENERALCYLINDRICAL = 4} type;
  Vec3D x0;
  Vec3D dir; //!< relevant for PLANAR, CYLINDRICAL and GENERALCYLINDRICAL

  Vec2D xmin, xmax; //!< bounding box for interpolation, only for GENERALCYLINDRICAL

  enum Vars {COORDINATE = 0, RADIALCOORDINATE = 1, DENSITY = 2, VELOCITY = 3, 
             RADIALVELOCITY = 4, PRESSURE = 5, LEVELSET = 6, MATERIALID = 7, 
             TEMPERATURE = 8, SIZE = 9};

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
  void readUserSpecifiedIC_GeneralCylindrical(std::fstream &input);
};

//------------------------------------------------------------------------------

struct ProbeNode {

  double locationX;
  double locationY;
  double locationZ;

  ProbeNode();
  ~ProbeNode() {}

  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct Probes {

  ObjectMap<ProbeNode> myNodes;

  int frequency;

  enum Vars  {DENSITY = 0, VELOCITY_X = 1, VELOCITY_Y = 2, VELOCITY_Z = 3, PRESSURE = 4, TEMPERATURE = 5, 
              MATERIALID = 6, LEVELSET0 = 7, LEVELSET1 = 8, LEVELSET2 = 9, LEVELSET3 = 10, LEVELSET4 = 11, SIZE = 12};

  const char *density;
  const char *velocity_x;
  const char *velocity_y;
  const char *velocity_z;
  const char *pressure;
  const char *temperature;
  const char *materialid;
  const char *levelset0;
  const char *levelset1;
  const char *levelset2;
  const char *levelset3;
  const char *levelset4;

  Probes();
  ~Probes() {}

  void setup(const char *, ClassAssigner * = 0);
};

//------------------------------------------------------------------------------

struct LinePlot {

  LinePlot();
  ~LinePlot() {}

  Assigner *getAssigner();

  const char *filename_base; //!< filename without path

  double x0,y0,z0;
  double x1,y1,z1;

  int numPoints;
  int frequency;

};

//------------------------------------------------------------------------------

struct MaterialVolumes {

  const char *filename;
  int frequency;

  MaterialVolumes();
  ~MaterialVolumes() {}

  void setup(const char *, ClassAssigner * = 0);
};

//------------------------------------------------------------------------------

struct OutputData {

  const char *prefix; //!< path
  const char *solution_filename_base; //!< filename without path

  enum Options {OFF = 0, ON = 1};
  Options density, velocity, pressure, materialid, internal_energy, temperature;
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


  Probes probes;

  ObjectMap<LinePlot> linePlots;

  MaterialVolumes materialVolumes;

  const char *mesh_filename; //!< file for nodal coordinates

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
