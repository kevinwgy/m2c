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
bool already_called = false;
//------------------------------------------------------------
void
MySolutionCalculator::GetUserDefinedSolution(int i0, int j0, int k0, int imax, int jmax, int kmax,
                                             Vec3D*** coords, Vec5D*** v, double*** id,
                                             std::vector<double***> phi, double time)
{
  // This function is called by the M2C solver at the end of each time step.

  // Note:
  // 1. Each processor core calls this function with different inputs.
  // 2. Override the original values of v, id, and possibly phi.
  // 3. If phi is updated, make sure it is consistent with id.
  // 4. The state variables (v) are each node are: rho,u,v,w,p.
  // 5. Do NOT change coords.

  if(already_called)
    return; //make sure that this solution is prescribed ONLY at the very beginning of the simulation

  // [IMPORTANT] Parameters values must be the same as those in input.st
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0-1; i<imax; i++) {
        if(coords[k][j][i][0]>0.2 && coords[k][j][i][0]<0.3) {// for x in (0.2, 0.3)
          v[k][j][i][0] = 0.8*v[k][j][i][0]; //drop density by 20%
          //v[k][j][i][1] = u;
          //v[k][j][i][2] = 0.0;
          //v[k][j][i][3] = 0.0;
          //v[k][j][i][4] = p;
        }
      }

  already_called = true;
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

