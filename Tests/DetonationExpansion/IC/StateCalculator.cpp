/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <climits>
#include <cassert>
#include "UserDefinedState.h"
#include "Vector3D.h"
#include "Vector5D.h"

using namespace std;

//------------------------------------------------------------
// This is a template for users to fill. Not compiled w/ M2C.
// If multiple UserDefinedState need to be defined, they should
// have different names (e.g., StateCalculator1, StateCalculator2, etc.)
// Compilation script (an example):
// g++ -O3 -fPIC -I/path/to/folder/that/contains/UserDefinedState.h -c UserDefinedState.cpp; g++ -shared UserDefinedState.o -o UserDefinedState.so; rm UserDefinedState.o
//------------------------------------------------------------

class StateCalculator : public UserDefinedState{

public:
  void GetUserDefinedState(int i0, int j0, int k0, int imax, int jmax, int kmax,
                           Vec3D*** coords, Vec5D*** v, double*** id, vector<double***> phi);
};

//------------------------------------------------------------

void
StateCalculator::GetUserDefinedState(int i0, int j0, int k0, int imax, int jmax, int kmax,
                                     Vec3D*** coords, Vec5D*** v, double*** id, vector<double***> phi)
{
  // Note:
  // 1. Each processor core calls this functio with different inputs.
  // 2. Override the original values of v and id. 
  // 4. The state variables are each node are: rho,u,v,w,p.
  // 3. Do NOT change coords.

  ifstream input("./IC/SphericalShock.txt", ios::in); //Be careful with relative path
  if(!input) {
    fprintf(stderr, "***Error: Could not open SphericalShock.txt.\n");
    exit(-1);
  }
  
  // This is the radius specified in GeometricEntities
  double bubble_radius = 32.14; 

  // common string for reading input file
  string line;

  // ignore headers
  for(int i=0; i<4; ++i) getline(input, line);

  // initialize vectors
  int initial_size = 1000;
  vector<double> xi(initial_size, 0.0);
  vector<double> material(initial_size, 0.0);
  vector<double> density(initial_size, 0.0);
  vector<double> velocity(initial_size, 0.0);
  vector<double> pressure(initial_size, 0.0);

  int m;
  for(m=0; m<INT_MAX; ++m) {
    double xi_m, material_m, density_m, velocity_m, pressure_m;

    getline(input, line);
    istringstream is(line);

    is >> xi_m >> material_m >> density_m >> velocity_m >> pressure_m;    
    //fprintf(stdout,"%e %e.\n", x_m, material_m);

    // eof check
    if(is.fail()) break;

    if(m>=xi.size()) {
      xi.resize(initial_size + xi.size(), 0.0);
      material.resize(initial_size + material.size(), 0.0);
      density.resize(initial_size + density.size(), 0.0);
      pressure.resize(initial_size + pressure.size(), 0.0);
      velocity.resize(initial_size + velocity.size(), 0.0);
    }

    xi[m] = xi_m; // Note: These are self-similar coordinates.
    material[m] = material_m;
    density[m] = density_m;
    pressure[m] = pressure_m;
    velocity[m] = velocity_m;
  }
  xi.resize(m);
  material.resize(m);
  density.resize(m);
  pressure.resize(m);
  velocity.resize(m);

  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        // calculate radial distance from the origin.
        Vec3D dir(coords[k][j][i][0], coords[k][j][i][1], coords[k][j][i][2]);
        double radius = dir.norm();

        dir /= radius;
        double xi_r = radius / bubble_radius;

        // For ids check input file under Equations.
        if(xi_r < xi.back()) { // inside bubble
          // find the interval
          int p1=-1, p2=-1;
          for(int p=1; p<xi.size(); ++p) {
            if(xi_r < xi[p]) {
              p1 = p-1;
              p2 = p;
              break;
            }
          }

          assert(p1>=0 && p2>=0);

          id[k][j][i] = 1.0; // JWL
          //fprintf(stdout,"[%d][%d][%d] setting id = 1.\n", k,j,i);

          double c1 = (xi_r - xi[p1])/(xi[p2] - xi[p1]);
          double c2 = 1.0 - c1;

          // assign state variables
          // The state variables are each node are: rho,u,v,w,p.
          v[k][j][i][0] = c1*density[p2] + c2*density[p1];
          v[k][j][i][1] = (c1*velocity[p2] + c2*velocity[p1]) * dir[0];
          v[k][j][i][2] = (c1*velocity[p2] + c2*velocity[p1]) * dir[1];
          v[k][j][i][3] = (c1*velocity[p2] + c2*velocity[p1]) * dir[2];
          v[k][j][i][4] = c1*pressure[p2] + c2*pressure[p1];

        }
        else {

          id[k][j][i] = 0.0; // PG

          // assign ambient state -- check input file.
          v[k][j][i][0] = 1.177e-6;
          v[k][j][i][1] = 0.0;
          v[k][j][i][2] = 0.0;
          v[k][j][i][3] = 0.0;
          v[k][j][i][4] = 1e5;

        }

      }

}

//------------------------------------------------------------
// The class factory (Note: Do NOT change these functions except the word "StateCalculator".)
extern "C" UserDefinedState* Create() {
  return new StateCalculator;
}

extern "C" void Destroy(UserDefinedState* uds) {
  delete uds;
}

//------------------------------------------------------------

