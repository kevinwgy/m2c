#ifndef _EMBEDDED_SURFACE_H_
#define _EMBEDDED_SURFACE_H_
#include<Vector3D.h>
#include<vector>

//----------------------------------------------------------------
// Class EmbeddedSurface stores information about a dynamic 
// surface embedded in the computational domain.
//----------------------------------------------------------------

class EmbeddedSurface {

  TriangulatedSurface S0; //original mesh
  TriangulatedSurface S; //deformed/current mesh

  double t0, t; //time corresponding to S0 and S

public:
  EmbeddedSurface();
  ~EmbeddedSurface();

  void ReadMeshFile(XXX);
  void GetMesh











};

#endif
