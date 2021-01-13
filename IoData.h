#ifndef _IO_DATA_H_
#define _IO_DATA_H_

#include <cstdio>
#include <map>
#include "parser/Assigner.h"
#include "parser/Dictionary.h"
#include <Vector3D.h>

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
  void setup(const char *name, ClassAssigner *);
  ~ObjectMap()
    {
      for(typename map<int, DataType *>::iterator it=dataMap.begin();it!=dataMap.end();++it)
      {
        delete it->second;
      }
    }
};

//------------------------------------------------------------------------------

struct MeshData {

  enum Type {TWODIMENSIONAL = 0, CYLINDRICAL = 1, THREEDIMENSIONAL = 2} type;
  double x0, xmax, y0, ymax, z0, zmax;
  int Nx, Ny, Nz;

  enum BcType {NONE = 0, INLET = 1, OUTLET = 2, WALL = 3, SYMMETRY = 4, SIZE = 5};
  BcType bc_x0, bc_xmax, bc_y0, bc_ymax, bc_z0, bc_zmax;

  MeshData();
  ~MeshData() {} 

  void setup(const char *, ClassAssigner * = 0);
};

//------------------------------------------------------------------------------

struct GasModelData {

  enum Type {IDEAL = 0, STIFFENED = 1} type;

  double specificHeatRatio;
  double idealGasConstant;
  double pressureConstant;
  double specificHeatPressure;

  GasModelData();
  ~GasModelData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct FluidModelData {

  enum Fluid { STIFFENED_GAS = 0} fluid;
  double rhomin;
  double pmin;

  GasModelData gasModel;

  FluidModelData();
  ~FluidModelData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct EquationsData {

  int numPhase;
  FluidModelData fluid1;

  EquationsData();
  ~EquationsData() {}

  void setup(const char *, ClassAssigner * = 0);
};

//------------------------------------------------------------------------------

struct SchemeData {

  enum Flux {ROE = 0, HLLE = 2, HLLC = 3, KURGANOV_TADMOR = 4} flux;

  enum Reconstruction {CONSTANT = 0, LINEAR = 1} reconstruction;

  enum Limiter {NONE = 0, GENERALIZED_MINMOD = 1, VANALBADA = 2, MODIFIED_VANALBADA = 3} limiter; 

  double generalized_minmod_coeff;
  

  int allowsFlux;

  //! allowsFlux = 0 for levelset equation (the choice of flux for the levelset is hardcoded)
  SchemeData(int allowsFlux = 1);
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

struct SchemesData {

  SchemeData ns;
  SchemeData ls;
  BoundarySchemeData bc;

  SchemesData();
  ~SchemesData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct ExplicitData {

//!time-integration scheme used
  enum Type {RUNGE_KUTTA_4 = 0, RUNGE_KUTTA_2 = 1, FORWARD_EULER = 2} type;

  ExplicitData();
  ~ExplicitData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct TsData {

  enum Type {EXPLICIT = 0, IMPLICIT = 1} type;
  int maxIts;
  double timestep;
  double maxTime;
  ExplicitData expl;

  TsData();
  ~TsData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct BcsFreeStreamData {

  double density;
  double velocity_x;
  double velocity_y;
  double velocity_z;
  double pressure;
  double temperature;

  BcsFreeStreamData();
  ~BcsFreeStreamData() {}

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

  BcsFreeStreamData inlet;
  BcsFreeStreamData outlet;
  BcsWallData wall;

  BcsData();
  ~BcsData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct IcData {

  //! user-specified file
  const char *user_specified_ic;

  enum Type {NONE = 0, PLANAR = 1, CYLINDRICAL = 2, SPHERICAL = 3} type;
  Vec3D x0;
  Vec3D dir; //!< relevant for PLANAR and CYLINDRICAL
  enum Vars {COORDINATE = 0, DENSITY = 1, VELOCITY = 2, PRESSURE = 3,
             LEVELSET = 4, MATERIALID = 5, TEMPERATURE = 6, SIZE = 7};

  int specified[SIZE];  //!< 0~unspecified, 1~specified

  vector<double> user_data[SIZE];

  const char *user_specified_ic2; //!< for radial variation in the case of cylindrical
  vector<double> user_data2[SIZE]; //!< for radial variation 

  IcData();
  ~IcData() {}

  void setup(const char *, ClassAssigner * = 0);
  void readUserSpecifiedIC();
};

//------------------------------------------------------------------------------

struct OutputData {

  const char *prefix; //!< path
  const char *solution_filename_base; //!< filename without path

  enum Options {OFF = 0, ON = 1};
  Options density, velocity, pressure, levelset, materialid, temperature;
 
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
