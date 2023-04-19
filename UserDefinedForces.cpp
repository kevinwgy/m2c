/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include "UserDefinedForces.h"

//------------------------------------------------------------
// This is a template for users to fill. Not compiled w/ M2C.
// If multiple UserDefinedForces need to be defined, they should
// have different names (e.g., MyForceCalculator1, MyForceCalculator2, etc.)
// Compilation script (an example):
// g++ -O3 -fPIC -I/path/to/folder/that/contains/UserDefinedForces.h -c UserDefinedForces.cpp; g++ -shared UserDefinedForces.o -o UserDefinedForces.so; rm UserDefinedForces.o
//------------------------------------------------------------

class MyForcesCalculator : public UserDefinedForces{

  // Calculates a force (fx,fy,fz) for each node; elements are triangles.
  // The dimension of X0, X, and nodal_forces is 3*nNodes. For example, X0[3*i+j] is the 
  // j-th coordinate (j = 0,1,2) of node i.
  void GetUserDefinedForces(double time, int nNodes, double *X0, double *X,
                            int nElems, int *elems, double* nodal_forces/*output*/);
};

//------------------------------------------------------------

void
MyForcesCalculator::GetUserDefinedForces(double time, int nNodes, double *X0, double *X,
                                         int nElems, int *elems, double* nodal_forces/*output*/)
{
  //TODO: The user should complete this function
  for(int i=0; i<nNodes; i++) {
    for(int j=0; j<3; j++) {
      nodal_forces[3*i+j] = 0.0;
    }
  }

}

//------------------------------------------------------------
// The class factory (Note: Do NOT change these functions except the word "MyForcesCalculator".)
extern "C" UserDefinedForces* Create() {
  return new MyForcesCalculator; //TODO: If you've changed the name of the derived class, this needs to be updated.
}

extern "C" void Destroy(UserDefinedForces* udf) {
  delete udf;
}

//------------------------------------------------------------

