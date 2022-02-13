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

struct CylinderSphereData {

  double cen_x, cen_y, cen_z, nx, ny, nz, r, L;
 
  enum OnOff {Off = 0, On = 1};
  OnOff front_cap;
  OnOff back_cap;

  StateVariable initialConditions;

  CylinderSphereData();
  ~CylinderSphereData() {}
  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct MultiInitialConditionsData {

  ObjectMap<PointData>    pointMap;
  ObjectMap<PlaneData>    planeMap;
  ObjectMap<SphereData>   sphereMap;
  ObjectMap<SpheroidData> spheroidMap;
  ObjectMap<CylinderConeData> cylinderconeMap;
  ObjectMap<CylinderSphereData> cylindersphereMap;

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

  enum Type {THREEDIMENSIONAL = 0, SPHERICAL = 1, CYLINDRICAL = 2} type;
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

  void check(); //!< check input parameters (for spherical & cylindrical domains)
};

//------------------------------------------------------------------------------

struct StiffenedGasModelData {

  double specificHeatRatio;
  double pressureConstant;

  //! parameters related to temperature
  double cv; //!< specific heat at constant volume
  double T0;  //!< temperature is T0 when internal energy (per mass) is e0
  double e0;  //!< internal energy per specific mass at T0

  double cp; //!< specific heat at constant pressure
  double h0; //!< enthalpy per specific mass at T0
  //NOTE: For stiffened (non-perfect) gas, calculating temperature using (cv, T0, e0) is NOT
  //      equivalent to using (cp, T0, h0). See KW's note. By default, cv is used. But if cv
  //      is 0 while cp>0, cp will be used.
 

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

  //! parameters related to temperature
  double cv; //!< specific heat at constant volume
  double T0;  //!< temperature is T0 when internal energy (per mass) is e0

  double cp; //!< specific heat at constant pressure
  double h0; //!< enthalpy per specific mass at T0

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

struct ViscosityModelData {

  enum Type {NONE = 0, CONSTANT = 1, SUTHERLAND = 2, ARTIFICIAL_RODIONOV = 3} type;

  // constant
  double dynamicViscosity;
  double bulkViscosity;

  // Sutherland
  double sutherlandConstant;
  double sutherlandT0;
  double sutherlandMu0;

  // Artificial viscosity (Rodionov)
  double Cav, Cth; 

  ViscosityModelData();

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct MaterialModelData {

  int id;
  enum EOS {STIFFENED_GAS = 0, MIE_GRUNEISEN = 1, JWL = 2} eos;
  double rhomin;
  double pmin;
  double rhomax;
  double pmax;

  double failsafe_density; //for updating phase change -- last resort

  StiffenedGasModelData    sgModel;
  MieGruneisenModelData    mgModel;
  JonesWilkinsLeeModelData jwlModel;

  ViscosityModelData viscosity;

  MaterialModelData();
  ~MaterialModelData() {}
  Assigner *getAssigner();

};


//------------------------------------------------------------------------------

struct MaterialTransitionData {

  int from_id, to_id;

  double temperature_lowerbound; //!< transition occurs if temperature is lower than this value
  double temperature_upperbound; //!< transition occurs if temperature is higher than this value
  double pressure_lowerbound;
  double pressure_upperbound;

  double latent_heat;

  MaterialTransitionData();
  ~MaterialTransitionData() {}

  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct EquationsData {

  ObjectMap<MaterialModelData> materials;

  ObjectMap<MaterialTransitionData> transitions;

  EquationsData();
  ~EquationsData() {}

  void setup(const char *, ClassAssigner * = 0);
};

//------------------------------------------------------------------------------

struct FixData {

  ObjectMap<SphereData>   sphereMap;
  ObjectMap<SpheroidData> spheroidMap;
  ObjectMap<CylinderConeData> cylinderconeMap;
  ObjectMap<CylinderSphereData> cylindersphereMap;

  FixData();
  ~FixData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct ReconstructionData {

  enum Reconstruction {CONSTANT = 0, LINEAR = 1} type;
  enum Limiter {NONE = 0, GENERALIZED_MINMOD = 1, VANALBADA = 2} limiter; 
  enum SlopeNearInterface {ZERO = 0, NONZERO = 1} slopeNearInterface;

  double generalized_minmod_coeff;

  // for nonlinear conservation laws, reconstruction can be done for variables
  // of different forms (primitive, conservative, characteristic), and the effect
  // can be different.
  enum VariableType {PRIMITIVE = 0, CONSERVATIVE = 1, PRIMITIVE_CHARACTERISTIC = 2,
                     CONSERVATIVE_CHARACTERISTIC = 3} varType;

  //User-specified regions in which type is reset to CONSTANT
  FixData fixes;

  ReconstructionData();
  ~ReconstructionData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct SmoothingData {

  enum Type {NONE = 0, BOX = 1, GAUSSIAN = 2} type;

  int iteration; //number of smoothing iterations applied each time

  double sigma_factor; //coefficient for Gaussian smoothing filter (multipled by dx)

  int frequency;
  double frequency_dt;

  enum ONOFF {OFF = 0, ON = 1} conservation;
  double conservation_tol;

  SmoothingData();
  ~SmoothingData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct SchemeData {

  enum Flux {ROE = 0, LOCAL_LAX_FRIEDRICHS = 1, HLLC = 2, GODUNOV = 3} flux;
 
  double delta; //! The coeffient in Harten's entropy fix.

  ReconstructionData rec;

  SmoothingData smooth;

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

struct LevelSetReinitializationData {

  int frequency; 
  double frequency_dt;

  int maxIts;

  double cfl;

  double convergence_tolerance;
  
  enum FirstLayerTreatment {FIXED = 0, 
                            CONSTRAINED1 = 1, CONSTRAINED2 = 2, //CR-1 and CR-2,Hartmann (2008)
                            ITERATIVE_CONSTRAINED1 = 3, ITERATIVE_CONSTRAINED2 = 4} //HCR-1&2,Hartmann 2010
           firstLayerTreatment;

  LevelSetReinitializationData();
  ~LevelSetReinitializationData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct LevelSetSchemeData {

  int materialid; //! The material in the phi<0 region ("inside")

  enum Solver {FINITE_VOLUME = 0, FINITE_DIFFERENCE = 1} solver;

  enum FiniteDifferenceMethod {UPWIND_CENTRAL_3 = 0} fd;

  enum Flux {ROE = 0, LOCAL_LAX_FRIEDRICHS = 1, UPWIND = 2} flux;
  ReconstructionData rec;

  enum BcType {NONE = 0, ZERO_NEUMANN = 1, LINEAR_EXTRAPOLATION = 2, NON_NEGATIVE = 3, SIZE = 4};
  BcType bc_x0, bc_xmax, bc_y0, bc_ymax, bc_z0, bc_zmax;

  
  double delta; //! The coeffient in Harten's entropy fix.

  int bandwidth; //number of layers of nodes on each side of interface

  LevelSetReinitializationData reinit;

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
  int maxIts_bracket;
  int maxIts_shock;
  int numSteps_rarefaction;
  double tol_main;
  double tol_shock;
  double tol_rarefaction;

  double min_pressure; //this is to guide the finding of bracketing interval (set it to
                       //be a low (maybe negative) pressure that is *clearly* out of bound

  double failure_threshold; //when to apply a fixed pressure (at failure)
  double pressure_at_failure; //this is a fixed pressure to be specified as ps when the solver fails to
                              //find a bracketing interval and the best approximation obtained is poor.
                              //this is the last resort. Usually it can be set to a very low but physical pressure

  ExactRiemannSolverData();
  ~ExactRiemannSolverData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct MultiPhaseData {

  enum Flux {EXACT = 0, NUMERICAL = 1} flux;

  enum ReconstructionAtInterface {CONSTANT = 0, LINEAR = 1} recon;

  double conRec_depth; //!< depth (fabs(phi)) where constant reconstruction is applied (default: 0)

  enum PhaseChangeType {RIEMANN_SOLUTION = 0, EXTRAPOLATION = 1} phasechange_type;

  enum RiemannNormal {LEVEL_SET = 0, MESH = 1, AVERAGE = 2} riemann_normal;

  enum OnOff {Off = 0, On = 1};
  OnOff latent_heat_transfer; //!< whether stored latent heat would be added to the enthalpy
                              //Note: In the case of a "physical" phase transition, the
                              //      latent heat is always added to the enthalpy. The option here
                              //      is about whether this operation will be done if a phase
                              //      change occurs due to the motion of material interface(s).

  int levelset_correction_frequency; //!< frequency of eliminating small gaps or inconsistencies
                                     //   between multiple level set functions

  OnOff apply_failsafe_density;

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

  //-----------------------------------------------------------------------
  //! initial condition specified using simple geometric entities (e.g., point, plane,)
  MultiInitialConditionsData multiInitialConditions;
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  //! user-specified file
  const char *user_specified_ic;

  enum RadialBasisFunction {MULTIQUADRIC = 0, INVERSE_MULTIQUADRIC = 1, 
                            THIN_PLATE_SPLINE = 2, GAUSSIAN = 3} rbf; //radial basis function for interpolation

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
  //-----------------------------------------------------------------------

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

struct LaserAbsorptionCoefficient {

  int materialid;
  double slope;
  double T0;
  double alpha0;

  LaserAbsorptionCoefficient();
  ~LaserAbsorptionCoefficient() {}

  Assigner *getAssigner();
};

//------------------------------------------------------------------------------

struct LaserData {

  // physical parameters
  double source_intensity; //!< radiance
  enum SourceDistribution {CONSTANT = 0, GAUSSIAN = 1} source_distribution;
  double source_power;
  const char *source_power_timehistory_file;
  double source_center_x, source_center_y, source_center_z;
  double source_dir_x, source_dir_y, source_dir_z;
  double source_radius;
  double source_beam_waist;
  double focusing_angle_degrees; //!< divering if <0
  double range;
  double lmin; //!< inside the laser domain (and ghosts), L>=lmin. (should be a tiny pos number.)
  ObjectMap<LaserAbsorptionCoefficient> abs; //!< absorption coefficients

  // parallel solution approach
  enum Parallelization {ORIGINAL = 0, BALANCED = 1} parallel;
  int min_cells_per_core;

  // numerical parameters
  double source_depth;
  double alpha;
  double convergence_tol;
  double max_iter;
  double relax_coeff;

  LaserData();
  ~LaserData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct AtomicIonizationModel {

  double molar_fraction;
  int    atomic_number;
  double molar_mass;

  int max_charge;

  const char* ionization_energy_filename; //!< file that stores ionization energies
  
  const char* excitation_energy_files_prefix; //!< prefix of files that contain excitation energies (for different excited states)
  const char* excitation_energy_files_suffix;

  const char* degeneracy_files_prefix; //!< prefix of files that contain degeneracy (for different excited states)
  const char* degeneracy_files_suffix;

  AtomicIonizationModel();
  ~AtomicIonizationModel() {}

  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct MaterialIonizationModel{

  enum Type {NONE = 0, SAHA_IDEAL = 1, SAHA_NONIDEAL = 2, SIZE = 3} type;

  int maxIts;
  double convergence_tol;

  enum PartitionFunctionEvaluation {ON_THE_FLY = 0, CUBIC_SPLINE_INTERPOLATION = 1,
                                    LINEAR_INTERPOLATION = 2} partition_evaluation;
  
  //! Tmin is still used as a threshold, below which
  //! calculation will not be performed (i.e. ionization will not happen)
  double ionization_Tmin;

  // numerical parameters for sampling & interpolating the partition function
  double sample_Tmin, sample_Tmax;
  int sample_size;

  ObjectMap<AtomicIonizationModel> elementMap;
  
  MaterialIonizationModel();
  ~MaterialIonizationModel() {}

  Assigner *getAssigner();
};

//------------------------------------------------------------------------------

struct IonizationData {

  // global constants
  double planck_constant;
  double electron_charge; //needed?
  double electron_mass;
  double boltzmann_constant;
  
  ObjectMap<MaterialIonizationModel> materialMap;
  
  IonizationData();
  ~IonizationData() {}

  void setup(const char *, ClassAssigner * = 0);

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
  double frequency_dt;

  enum Vars  {DENSITY = 0, VELOCITY_X = 1, VELOCITY_Y = 2, VELOCITY_Z = 3, PRESSURE = 4, TEMPERATURE = 5, 
              DELTA_TEMPERATURE = 6, MATERIALID = 7, LASERRADIANCE = 8, LEVELSET0 = 9, LEVELSET1 = 10, 
              LEVELSET2 = 11, LEVELSET3 = 12, LEVELSET4 = 13, IONIZATION = 14, SIZE = 15};

  const char *density;
  const char *velocity_x;
  const char *velocity_y;
  const char *velocity_z;
  const char *pressure;
  const char *temperature;
  const char *delta_temperature;
  const char *materialid;
  const char *laser_radiance;
  const char *levelset0;
  const char *levelset1;
  const char *levelset2;
  const char *levelset3;
  const char *levelset4;
  const char *ionization_result;

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
  double frequency_dt;

};

//------------------------------------------------------------------------------

struct MaterialVolumes {

  const char *filename;
  int frequency;
  double frequency_dt;

  MaterialVolumes();
  ~MaterialVolumes() {}

  void setup(const char *, ClassAssigner * = 0);
};

//------------------------------------------------------------------------------

struct OutputData {

  const char *prefix; //!< path
  const char *solution_filename_base; //!< filename without path

  enum Options {OFF = 0, ON = 1};
  Options density, velocity, pressure, materialid, internal_energy, temperature, delta_temperature, laser_radiance;

  enum VerbosityLevel {LOW = 0, MEDIUM = 1, HIGH = 2} verbose;

  const static int MAXLS = 5;
  Options levelset[MAXLS];
  Options levelset0;
  Options levelset1;
  Options levelset2;
  Options levelset3;
  Options levelset4;

  //! ionization-related outputs
  Options mean_charge; //!< Zav
  Options heavy_particles_density; //!< number density of heavy particles
  Options electron_density; //!< number density of electrons
  int max_charge_number; //!< all the ions with higher charge number will be lumped together in the output
  const static int MAXSPECIES = 5; //!< This only limits outputing. MaterialIonizationModel can have more species.
  Options molar_fractions[MAXSPECIES];
  Options molar_fractions0;
  Options molar_fractions1;
  Options molar_fractions2;
  Options molar_fractions3;
  Options molar_fractions4;
  inline bool ionization_output_requested() {
    return mean_charge==ON || heavy_particles_density==ON || electron_density==ON || molar_fractions0==ON ||
           molar_fractions1==ON || molar_fractions2==ON || molar_fractions3==ON || molar_fractions4==ON;}

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

struct EmbeddedSurfaceData {

  const char *filename; //!< file for nodal coordinates and elements
  
  enum Type {None = 0, Wall = 1, Symmetry = 2, DirectState = 3, MassFlow = 4, PorousWall = 5,
             Size = 6} type;
             
  enum ThermalCondition {Adiabatic = 0, Isothermal = 1, Source = 2} thermal;

  double heat_source;

  EmbeddedSurfaceData();
  ~EmbeddedSurfaceData() {}

  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct EmbeddedSurfacesData {

  ObjectMap<EmbeddedSurfaceData> surfaces;

  EmbeddedSurfacesData();
  ~EmbeddedSurfacesData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct AerosCouplingData {

  enum FSICouplingAlgorithm {NONE = 0, BY_AEROS = 1, C0 = 2, A6 = 3} fsi_algo;
  enum Fracture {OFF = 0, ON = 1} fracture;

  AerosCouplingData();
  ~AerosCouplingData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct ConcurrentProgramsData {

  AerosCouplingData aeros;

  ConcurrentProgramsData();
  ~ConcurrentProgramsData() {} 

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

class IoData {

  char *cmdFileName;
  FILE *cmdFilePtr;

public:

  ConcurrentProgramsData concurrent;

  EmbeddedSurfacesData embed_surfaces;

  MeshData mesh;

  EquationsData eqs;
  IcData ic;  //!< initial condition
  BcsData bc;  //!< boundary condition

  SchemesData schemes;

  ExactRiemannSolverData exact_riemann;

  MultiPhaseData multiphase;

  LaserData laser;

  IonizationData ion;

  TsData ts;

  OutputData output;

public:

  IoData() {}
  IoData(int, char**);
  ~IoData() {}

  void readCmdLine(int, char**);
  void setupCmdFileVariables();
  void readCmdFile();
  void finalize();

};
#endif
