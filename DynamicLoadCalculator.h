#ifndef _DYNAMIC_LOAD_CALCULATOR_H_
#define _DYNAMIC_LOAD_CALCULATOR_H_

#include<ConcurrentProgramsHandler.h>
#include<KDTree.h>
#include<memory> //shared_ptr

struct TriangulatedSurface;

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

  //! Internal variables storing the snapshots currently stored in memory
  int id0, id1;
  std::shared_ptr<std::vector<std::vector<double> > > S0, S1; //!< use smart pointers (automatically deleted)
  std::shared_ptr<KDTree<PointIn3D,3> > tree0, tree1; 
  std::vector<Vec3D> F0, F1; //!< interpolated forces (using S0 and S1)

public:

  DynamicLoadCalculator(IoData &iod_, MPI_Comm &comm_, ConcurrentProgramsHandler &concurrent_);
  ~DynamicLoadCalculator();

  void Run();

private:

  void RunForAeroS();

  void ComputeForces(TriangulatedSurface *surface, std::vector<Vec3D> *force, double t);

  void ReadMetaFile(std::string filename);
  void ReadSnapshot(std::string filename, std::vector<std::vector<double> >& S);

  void BuildKDTree(std::vector<std::vector<double> >& S, KDTree<PointIn3D,3>* &tree);
  void InterpolateInSpace(std::vector<std::vector<double> >& S, KDTree<PointIn3D,3>* tree,
                          std::vector<Vec3D>& X, int active_nodes, Var var, int var_dim, double* output);
  void InterpolateInTime(double t1, double* input1, double t2, double* input2,
                         double t, double* output, int size);
};

#endif
