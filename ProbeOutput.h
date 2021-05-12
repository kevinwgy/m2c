#ifndef _PROBE_OUTPUT_H_
#define _PROBE_OUTPUT_H_
#include <IoData.h>
#include <SpaceVariable.h>

/** This class is responsible for interpolating solutions at probe locations and outputing
 *  the interpolated solutions to files. It is owned by class Output
 *  Each solution variable (e.g., density, pressure) is written to a separate file
 */
class ProbeOutput {

  MPI_Comm &comm;

  int numNodes;
  int frequency;

  std::vector<Vec3D> locations;
  FILE *file[Probes::SIZE]; //one file per solution variable

  //! For each probe node, ijk are the lower nodal indices of the element that contains the node
  std::vector<Int3> ijk;
  //! For each probe node, trilinear_coords contains the local x,y,z coordinates w/i the element
  std::vector<Vec3D> trilinear_coords;

public:
  ProbeOutput(MPI_Comm &comm_, OutputData &iod_output);
  ~ProbeOutput();

  void SetupInterpolation(SpaceVariable3D &coordinates);
  void WriteSolutionAtProbes(double time, int time_step, SpaceVariable3D &V, SpaceVariable3D &ID,
                             std::vector<SpaceVariable3D*> &Phi); //!< write probe solution to file

private:
  double InterpolateSolutionAtProbe(Int3& ijk, Vec3D &trilinear_coords, double ***v, int dim, int p);

};




#endif
