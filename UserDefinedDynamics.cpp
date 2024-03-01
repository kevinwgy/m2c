/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include "UserDefinedDynamics.h"

//------------------------------------------------------------
// This is a template for users to fill. Not compiled w/ M2C.
// If multiple UserDefinedDynamics need to be defined, they should
// have different names (e.g., MyDynamicsCalculator1, MyDynamicsCalculator2, etc.)
// Compilation script (an example):
// g++ -O3 -fPIC -I/path/to/folder/that/contains/UserDefinedDynamics.h -c UserDefinedDynamics.cpp; g++ -shared UserDefinedDynamics.o -o UserDefinedDynamics.so; rm UserDefinedDynamics.o
//------------------------------------------------------------

class MyDynamicsCalculator : public UserDefinedDynamics{

public:
  // The dimension of X0, X, disp, and velo is 3*nNodes. For example, X0[3*i+j] is the 
  // j-th coordinate (j = 0,1,2) of node i.
  void GetUserDefinedDynamics(double time, double dt, int nNodes, double *X0, double *X,//inputs
                              double *disp, double *velo/*output*/);
};

//------------------------------------------------------------

void
MyDynamicsCalculator::GetUserDefinedDynamics(double time, double dt, int nNodes, double *X0, double *X,//inputs
                                             double *disp, double *velo)
{
  //TODO: The user should complete this function
  for(int i=0; i<nNodes; i++) {
    for(int j=0; j<3; j++) {
      disp[3*i+j] = 0.0;
      velo[3*i+j] = 0.0;
    }
  }

}

//------------------------------------------------------------
// The class factory (Note: Do NOT change these functions except the word "MyDynamicsCalculator".)
extern "C" UserDefinedDynamics* Create() {
  return new MyDynamicsCalculator; //TODO: If you've changed the name of the derived class, this needs to be updated.
}

extern "C" void Destroy(UserDefinedDynamics* udd) {
  delete udd;
}

//------------------------------------------------------------

