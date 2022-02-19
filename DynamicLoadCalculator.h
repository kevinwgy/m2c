#ifndef _DYNAMIC_LOAD_CALCULATOR_H_
#define _DYNAMIC_LOAD_CALCULATOR_H_

#include<ConcurrentProgramsHandler.h>

/**********************************************
 * class DynamicLoadCalculator is a special
 * tool that reads solution time-history from
 * user-specified files, calculates (dynamic)
 * loads on an embedded structure, and sends
 * them to a structural dynamics solver (e.g.,
 * Aero-S).
 *********************************************/

class DynamicLoadCalculator
{

  MPI_Comm& comm;
  IoData& iod;
  ConcurrentProgramsHandler& concurrent;

  //! Information about transient input data
  std::string prefix, suffix;
  double tmin, tmax;
  enum Var {NODE_NUMBER = 0, COORDINATES = 1, DENSITY = 2, VELOCITY = 3, PRESSURE = 4, 
            FORCE = 5, SIZE = 6};
  std::map<int,Var> column2var; 
  std::map<Var,int> var2column;

  //! Pair of time and the file label between prefix and suffix
  std::vector<std::pair<double,std::string> > stamp;


public:

  DynamicLoadCalculator(IoData &iod_, MPI_Comm &comm_, ConcurrentProgramsHandler &concurrent_);
  ~DynamicLoadCalculator();

  void Run();

private:

  void RunForAeroS();

  void ReadMetaFile(std::string filename);
  void ReadSnapshot(std::string filename);

  void BuildKDTree(vector<vector<double> >& S, KDTree<PointIn3D,3> *tree);
  void InterpolateInSpace(vector<vector<double> >& S, KDTree<PointIn3D,3>* tree,
                          vector<Vec3D>& X, Var var, int var_dim, double* output);
  void InterpolateInTime(double t1, double* input1, double t2, double* input2,
                         double t, double* output, int size);
};

#endif
