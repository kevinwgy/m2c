#include <TriangulatedSurface.h>
#include <Utils.h>

//----------------------------------------------------

void TriangulatedSurface::BuildConnectivity()
{
  int nNodes = X.size();
  int nElems = elems.size();

  node2node.clear();
  node2node.resize(nNodes);
  node2elem.clear();
  node2elem.resize(nNodes);
  
  for(int i=0; i<nElems; i++)
    for(int j=0; j<3; j++) {
      node2node[elems[i][j]].insert(elems[i][(j+1)%3]);
      node2node[elems[i][j]].insert(elems[i][(j+2)%3]);
      node2elem[elems[i][j]].insert(i);
    }
}

//----------------------------------------------------

void TriangulatedSurface::BuildElementNormals()
{
  int nElems = elems.size();

  elemNorm.clear();
  elemNorm.resize(nElems);

  int n1, n2, n3;
  Vec3D dx2, dx3;

  // Also look to determine a point inside the solid but away from the structure.
  for (int i=0; i<nElems; i++) {
    n1 = elems[i][0];
    n2 = elems[i][1];
    n3 = elems[i][2];

    dx2 = X[n2] - X[n1];
    if(n2==n3) {//assuming this is a set of line segments in 2D (x-y)
      degenerate = true;
      elemNorm[i][0] = -dx2[1];
      elemNorm[i][1] = dx2[0];
      elemNorm[i][2] = 0.0;
    } else {
      dx3 = X[n3] - X[n1];
      elemNorm[i] = dx2^dx3; // cross product ==> normal dir.
    }

    double nrm = elemNorm[i].norm();
    // normalize the normal.
    if(nrm > 0.0)
      elemNorm[i] /= nrm;
    else {
      print_error("Error: area (length) of triangle (edge) %d is %e.\n", i+1, nrm);
      exit_mpi();
    }
  }

  //verify the degenerate case
  if(degenerate) {
    for(int i=0; i<nElems; i++)
      if(n2!=n3) {
        print_error("Error: Cannot handle a mix of triangles and line segments.\n");
        exit_mpi();
      }
    for(int i=0; i<X.size(); i++)
      if(X[i][2] != 0) {
        print_error("Error: A degenerated triangulated surface must be on the x-y plane. Found z = %e.\n", 
                     X[i][2]);
        exit_mpi();
      }
  }
}

//----------------------------------------------------

