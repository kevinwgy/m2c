/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <vector>
#include "UserDefinedSolution.h"
#include "Vector3D.h"
#include "Vector5D.h"

//------------------------------------------------------------
// This is a template for users to fill. Not compiled w/ M2C.
// If multiple UserDefinedSolutions are to be defined, they should
// have different names (e.g., MySolutionCalculator1, MySolutionCalculator2, etc.)
// Compilation script (an example):
// g++ -O3 -fPIC -I/path/to/folder/that/contains/UserDefinedSolution.h -c UserDefinedSolution.cpp; g++ -shared UserDefinedSolution.o -o UserDefinedSolution.so; rm UserDefinedSolution.o
//------------------------------------------------------------

class MySolutionCalculator : public UserDefinedSolution{

public:
  void GetUserDefinedSolution(int i0, int j0, int k0, int imax, int jmax, int kmax,
                              Vec3D*** coords, Vec5D*** v, double*** id, std::vector<double***> phi,
                              double time);
};

//------------------------------------------------------------

void
MySolutionCalculator::GetUserDefinedSolution(int i0, int j0, int k0, int imax, int jmax, int kmax,
                                             Vec3D*** coords, Vec5D*** v, double*** id,
                                             std::vector<double***> phi, double time)
{
  // Note:
  // 1. Each processor core calls this function with different inputs.
  // 2. Override the original values of v, id, and possibly phi.
  // 3. If phi is updated, make sure it is consistent with id.
  // 4. The state variables (v) are each node are: rho,u,v,w,p.
  // 5. Do NOT change coords.

  //TODO: The user should complete this function
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        if(time>0.0 && false) {
          v[k][j][i]  = 0.0; //Vec5D
          id[k][j][i] = 0.0;
        }
      }

}

//------------------------------------------------------------
// The class factory (Note: Do NOT change these functions except the word "MySolutionCalculator".)
extern "C" UserDefinedSolution* Create() {
  return new MySolutionCalculator; //TODO: If you've changed the name of the derived class, this needs to be updated.
}

extern "C" void Destroy(UserDefinedSolution* udsl) {
  delete udsl;
}

//------------------------------------------------------------

