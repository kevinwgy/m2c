/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <vector>
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
                           Vec3D*** coords, Vec5D*** v, double*** id, std::vector<double***> phi);
};

//------------------------------------------------------------

void
MyStateCalculator::GetUserDefinedState(int i0, int j0, int k0, int imax, int jmax, int kmax,
                                       Vec3D*** coords, Vec5D*** v, double*** id,
                                       std::vector<double***> phi)
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
  double soundspeed = gamma*(p0+Pstiff)/rho0;

  // pressure signal
  double ps = 3.5e7; //Pa;
  double alpha = 1.48e6; //1/s
  double omega = 1.21e6; //1/s
  double PI = acos(-1.0);

  double x, t, p; //calculate p(x)
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        x = -(coords[k][j][i][0] + 0.18);
        if(x>0 && x<3.7) { //within the wave form
          t = x/soundspeed; 
          p = 2.0*ps*exp(-1.0*alpha*t)*cos(omega*t+PI/3.0);
          v[k][j][i][0] = rho0;
          v[k][j][i][1] = (p - p0)/(rho0*soundspeed);
          v[k][j][i][2] = 0.0;
          v[k][j][i][3] = 0.0;
          v[k][j][i][4] = p;
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

