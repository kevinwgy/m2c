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

  // [IMPORTANT] Parameters values must be the same as those in input.st
  //             Assumed units: mm, s, g
  double gamma = 6.59;
  double Pstiff = 4.1e8; //Pa
  double rho0 = 1.0e-3; //g/mm3
  double p0 = 1.01e5;; //Pa
  double soundspeed = sqrt(gamma*(p0+Pstiff)/rho0);

  // pressure signal
  double ps = 3.5e7; //Pa;
  double alpha = 1.48e6; //1/s
  double omega = 1.21e6; //1/s
  double PI = acos(-1.0);

  double p = p0 + 2.0*ps*exp(-1.0*alpha*time)*cos(omega*time+PI/3.0);
  double u = (p - p0)/(rho0*soundspeed);

  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0-1; i<imax; i++) {
        if(i<3) { //populate i=-1, 0, 1, 2 --- 4 layers, more than enough
          v[k][j][i][0] = rho0;
          v[k][j][i][1] = u;
          v[k][j][i][2] = 0.0;
          v[k][j][i][3] = 0.0;
          v[k][j][i][4] = p;
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

