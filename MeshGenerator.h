#ifndef _MESH_GENENERATOR_H_
#define _MESH_GENENERATOR_H_

#include <IoData.h>

/*********************************************************************************
 * class MeshGenerator is responsible for calculating the nodal positions 
 * and cell sizes of a uniform or non-uniform Cartesian grid. (EXCLUDING
 * the ghost boundary layer outside of the real domain)
 *
 * Rule: (using the x-direction as example)
 * Case | Nx specified? | res. points specified? | Output
 *  1.        YES                 NO               Uniform mesh
 *  2.        YES                 YES              THROW ERROR AND QUIT
 *  3.        NO                  NO               THROW ERROR AND QUIT
 *  4.        NO                  YES              Non-uniform, driven by res. points
 *********************************************************************************/

class MeshGenerator
{

public:
  MeshGenerator() {}
  ~MeshGenerator() {}

  //! Function that calculate nodal positions and cell sizes (dx,dy,dz)
  void ComputeMeshCoordinatesAndDeltas(MeshData &iod_mesh,
                                       vector<double> &x, vector<double> &y,
                                       vector<double> &z, vector<double> &dx,
                                       vector<double> &dy, vector<double> &dz);

protected:
  //! Internal functions
  void ComputeMesh1DUniform(double x0, double xmax, double Nx, vector<dobule> &x, vector<double> &dy);
  void ComputeMesh1DNonUniform(double x0, double xmax, vector<pair<double,double> > &xpoints,
                               vector<double> &x, vector<double> &dx);
};




#endif
