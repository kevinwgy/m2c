/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include "UserDefinedState.h"
#include "Vector3D.h"
#include "Vector5D.h"

//------------------------------------------------------------
// This is a template for users to fill. Not compiled w/ M2C.
// If multiple UserDefinedState need to be defined, they should
// have different names (e.g., MyStateCalculator1, MyStateCalculator2, etc.)
// Compilation script (an example):
// g++ -O3 -fPIC -I/path/to/folder/that/contains/UserDefinedState.h -c UserDefinedState.cpp; g++ -shared UserDefinedState.o -o UserDefinedState.so; rm UserDefinedState.o
//------------------------------------------------------------

class MyStateCalculator : public UserDefinedState{

public:
  void GetUserDefinedState(int i0, int j0, int k0, int imax, int jmax, int kmax,
                           Vec3D*** coords, Vec5D*** v, double*** id);
};

//------------------------------------------------------------

void
MyStateCalculator::GetUserDefinedState(int i0, int j0, int k0, int imax, int jmax, int kmax,
                                       Vec3D*** coords, Vec5D*** v, double*** id);
{
  // Note:
  // 1. Each processor core calls this functio with different inputs.
  // 2. Override the original values of v and id. 
  // 4. The state variables are each node are: rho,u,v,w,p.
  // 3. Do NOT change coords.

  //TODO: The user should complete this function
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        if(false) {
          v[k][j][i]  = 0.0; //Vec5D
          id[k][j][i] = 0.0;
        }
      }

}

//------------------------------------------------------------
// The class factory (Note: Do NOT change these functions except the word "MyStateCalculator".)
extern "C" UserDefinedState* Create() {
  return new MyStateCalculator; //TODO: If you've changed the name of the derived class, this needs to be updated.
}

extern "C" void Destroy(UserDefinedState* uds) {
  delete uds;
}

//------------------------------------------------------------

