#ifndef _PLANE_OUTPUT_H_
#define _PLANE_OUTPUT_H_

#include<ProbeOutput.h>

/** Class PlaneOutput is responsible for interpolating solutions on a plane,
 *  and printing the solution to a file. Upon initialization, the mesh of
 *  the plane is outputted in "top" format. A solution file is also initialized.
 *  Each solution snapshot is appended to the end of the existing solution file.
 *  To visualize the result, the user may use "xp2exo" to combine the mesh and
 *  the solutions into an exo (Exodus II) file.
 */
class PlaneOutput {

  MPI_Comm &comm;
  PlanePlot &iod_pplot;

  std::vector<VarFcnBase*> &vf;

  GlobalMeshInfo &global_mesh;

  //! using "ProbeOutput" to post-process and interpolate solutions
  ProbeOutput probe_util; 

  int frequency;
  double frequency_dt;
  int iFrame;
  double last_snapshot_time;

  FILE *file[PlanePlot::SIZE]; //!< one file per solution variable
  
  //! points and interpolated solutions within subdomain
  std::vector<Vec3D> points;
  std::vector<double> sol1; 
  std::vector<Vec3D> sol3; 
  std::set<Int3> ijk;
  std::vector<std::pair<int, std::array<bool,8> > > ijk_valid;
  std::vector<Vec3D> trilinear_coords;
  std::vector<int> loc2glob; //!< point map: local to global

  //! interpolated solutions for the entire mesh
  std::vector<double> sol1_global;
  std::vector<Vec3D> sol3_global;

public:

  PlaneOutput(MPI_Comm &comm_, OutputData &iod_output_, PlanePlot &iod_pplot_, 
              std::vector<VarFcnBase*> &vf_, GlobalMeshInfo &global_mesh_, 
              IonizationOperator *ion_);
  
  ~PlaneOutput();

  void InitializeOutput(SpaceVariable3D &coordinates); //!< Setup mesh and interpolation

  void WriteSolutionOnPlane(double time, double dt, int time_step, SpaceVariable3D &V, 
                            SpaceVariable3D &ID, std::vector<SpaceVariable3D*> &Phi, 
                            SpaceVariable3D* L, bool force_write);


};

#endif
