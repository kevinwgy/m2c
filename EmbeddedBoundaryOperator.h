#ifndef _EMBEDDED_BOUNDARY_OPERATOR_H_
#define _EMBEDDED_BOUNDARY_OPERATOR_H_

#include<TriangulatedSurface.h>
#include<IoData.h>
#include<SpaceVariable.h>
#include<cassert>
#include<map>

/******************************************************************
 * Class EmbeddedBoundaryOperator stores data related to one or
 * multiple embedded surfaces, which can be provided either by a
 * structural dynamics solver (e.g., AERO-S) or directly by the
 * user in a file. The class also contains functions for tracking
 * the embedded surfaces and enforcing interface/boundary conditions
 *****************************************************************/

class EmbeddedBoundaryOperator {

  IoData &iod;

  bool hasSurfFromOtherSolver;

  vector<TriangulatedSurface> surfaces; //embedded surfaces
  vector<vector<Vec3D> > F; //forces

//  SurfaceTracker tracker;
 
public:
   
  EmbeddedBoundaryOperator(IoData &iod_, bool surface_from_other_solver = false);
  ~EmbeddedBoundaryOperator();

  void Destroy();

  int  NumberOfSurfaces() {return surfaces.size();}
  bool HasSurfaceFromOtherSolver() {return hasSurfFromOtherSolver;} //!< index of surface from other solver is 0

  vector<TriangulatedSurface> *GetPointerToSurfaces() {return &surfaces;}
  TriangulatedSurface         *GetPointerToSurface(int i) {assert(i>=0 && i<surfaces.size()); return &surfaces[i];}
  vector<vector<Vec3D> >      *GetPointerToForces() {return &F;}
  vector<Vec3D>               *GetPointerToForcesOnSurface(int i) {assert(i>=0 && i<F.size()); return &F[i];}

  void ComputeForces(SpaceVariable3D &V, SpaceVariable3D &ID);

  void TrackSurfaces();
  void TrackUpdatedSurfaceFromOtherSolver();

private:

  void ReadMeshFile(const char *filename, std::string& nodeSetName,
                    std::map<int, EmbeddedSurfaceData::Type>& boundaryConditionsMap, 
                    std::vector<int>& faceIDs, int& numStNodes, int& numStElems, 
                    std::vector<Vec3D>& Xs, int *&surfaceID, vector<Int3> &stElem, int *&faceID, int& maxFaceID);

};




#endif
