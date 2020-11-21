#ifndef _IO_DATA_H_
#define _IO_DATA_H_

#include <cstdio>
#include <map>
#include "parser/ParseTree.h"
#include "parser/Dictionary.h"

using std::map;

class Assigner;
class ClassAssigner;
class Communicator;

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
  double x0, xmax, y0, ymax;
  int Nx, Ny;

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
  ObjectMap<FluidModelData> fluidModelMap;

  EquationsData();
  ~EquationsData() {}

  void setup(const char *, ClassAssigner * = 0);
};

//------------------------------------------------------------------------------

struct SchemeData {

  enum Flux {ROE = 0, VANLEER = 1, HLLE = 2, HLLC = 3} flux;

  enum Reconstruction {CONSTANT = 0, LINEAR = 1} reconstruction;

  enum Limiter {NONE = 0, VANALBADA = 1, BARTH = 2, VENKAT = 3, P_SENSOR = 4,
                EXTENDEDVANALBADA = 5} limiter;
  enum Dissipation {SECOND_ORDER = 0, SIXTH_ORDER = 1} dissipation;

  double beta;
  double gamma;
  double xiu;
  double xic;
  double eps;

  double xirho;
  double xip;
  double vel_fac;

  int allowsFlux;

  // allowsFlux = 0 for levelset equation (the choice of flux for the levelset is hardcoded)
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

//time-integration scheme used
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

class IoData {

  char *cmdFileName;
  FILE *cmdFilePtr;

public:

  MeshData mesh;
  EquationsData eqs;
  SchemesData schemes;
  TsData ts;

public:

  IoData() {}
  ~IoData() {}

  void readCmdLine(int, char**);
  void setupCmdFileVariables();
  void readCmdFile();

}
#endif
