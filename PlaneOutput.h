/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _PLANE_OUTPUT_H_
#define _PLANE_OUTPUT_H_

#include<IonizationOperator.h>
#include<GlobalMeshInfo.h>

/** Class PlaneOutput is responsible for interpolating solutions on a plane,
 *  and printing the solution to a file. Upon initialization, the mesh of
 *  the plane is outputted in "top" format. A solution file is also initialized.
 *  Each solution snapshot is appended to the end of the existing solution file.
 *  To visualize the result, the user may use "xp2exo" to combine the mesh and
 *  the solutions into an exo (Exodus II) file.
 *  Note: The mesh of the plane have duplicated nodes, when simulation is performed
 *  on more than one processor core. Nonetheless, the solutions at the duplicates
 *  should be identical.
 */
class PlaneOutput {

  MPI_Comm &comm;
  PlanePlot &iod_pplot;

  std::vector<VarFcnBase*> &vf;

  GlobalMeshInfo &global_mesh;

  IonizationOperator *ion; //!< for post-processing

  int frequency;
  double frequency_dt;
  int iFrame;
  double last_snapshot_time;

  std::string filename[PlanePlot::SIZE]; //!< one file per solution variable
  std::string filename_mesh;
 
  //! points and interpolated solutions within subdomain
  std::vector<Vec3D> points;
  std::vector<double> sol1; 
  std::vector<Vec3D> sol3; 
  std::vector<Int3> ijk;
  std::vector<std::pair<int, std::array<bool,8> > > ijk_valid;
  std::vector<Vec3D> trilinear_coords;

  //! interpolated solutions for the entire mesh
  std::vector<double> sol1_global;
  std::vector<Vec3D> sol3_global;

  //! for MPI communication (only defined on proc 0)
  std::vector<int> package_disp1;
  std::vector<int> package_size1;
  std::vector<int> package_disp3; 
  std::vector<int> package_size3;

public:

  PlaneOutput(MPI_Comm &comm_, OutputData &iod_output_, PlanePlot &iod_pplot_, 
              std::vector<VarFcnBase*> &vf_, GlobalMeshInfo &global_mesh_, 
              IonizationOperator *ion_);
  
  ~PlaneOutput();

  void InitializeOutput(SpaceVariable3D &coordinates); //!< Setup mesh & interpolation, print mesh

  void WriteSolutionOnPlane(double time, double dt, int time_step, SpaceVariable3D &V, 
                            SpaceVariable3D &ID, std::vector<SpaceVariable3D*> &Phi, 
                            SpaceVariable3D* L, bool force_write);

private:

  void GetScalarSolutionAndWrite(double time, double*** v, int dim, int p, std::string& fname);

  void GetVec3SolutionAndWrite(double time, double*** v, int dim, int p, std::string& fname);

  void GetTemperatureAndWrite(double time, double*** v, double*** id, std::string& fname);

  void GetDeltaTemperatureAndWrite(double time, double*** v, double*** id, std::string& fname);

  void GetIonizationAndWrite(double time, double*** v, double*** id, std::string& fname);


  void WriteScalarSolutionToFile(double time, std::vector<double>& S, std::string& fname);

  void WriteVec3SolutionToFile(double time, std::vector<Vec3D>& S, std::string& fname);

};

#endif
