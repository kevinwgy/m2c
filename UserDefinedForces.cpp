/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include "UserDefinedForces.h"
#include "Vector3D.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <cassert>
using namespace std;

//------------------------------------------------------------
// This is a template for users to fill. Not compiled w/ M2C.
// If multiple UserDefinedForces need to be defined, they should
// have different names (e.g., MyForceCalculator1, MyForceCalculator2, etc.)
// Compilation script (an example):
// g++ -O3 -fPIC -I/path/to/folder/that/contains/UserDefinedForces.h -c UserDefinedForces.cpp; g++ -shared UserDefinedForces.o -o UserDefinedForces.so; rm UserDefinedForces.o
//------------------------------------------------------------

class MyForcesCalculator : public UserDefinedForces{

public:

  // Calculates a force (fx,fy,fz) for each node; elements are triangles.
  // The dimension of X0, X, and nodal_forces is 3*nNodes. For example, X0[3*i+j] is the 
  // j-th coordinate (j = 0,1,2) of node i.
  void GetUserDefinedForces(double time, int nNodes, double *X0, double *X,
                            int nElems, int *elems, double* nodal_forces/*output*/);


// The user may define private functions & data if needed. The following is just an example
private:

  double GetPressure(Vec3D &point, double t);

  void LoadData();
  vector<double> ptime;
  vector<vector<double> > psensor;
  
};

//------------------------------------------------------------

void
MyForcesCalculator::GetUserDefinedForces(double time, int nNodes, double *X0, double *X,
                                         int nElems, int *elems, double* nodal_forces/*output*/)
{
  //TODO: The user should complete this function
  //      The following is just an example

  if(psensor.empty())
    LoadData();


  Vec3D *xyz0 = (Vec3D*)X0;
  Vec3D *xyz  = (Vec3D*)X;
  Int3  *elem = (Int3*)elems;
  Vec3D *F    = (Vec3D*)nodal_forces;

  for(int i=0; i<nNodes; i++)
    F[i] = 0.0; //clear


  // Loop through elements
  Vec3D g0; //location of gauss point IN THE ORIGINAL CONFIGURATION
  int n1, n2, n3;
  double area, p;
  Vec3D normal; //unit normal, POINTING INWARD
  Vec3D nodal_load;
  double experiment_time;
  for(int i=0; i<nElems; i++) {
    n1 = elem[i][0]; 
    n2 = elem[i][1]; 
    n3 = elem[i][2]; 
    normal = (xyz[n2]-xyz[n1])^(xyz[n3]-xyz[n1]);
    area = 0.5*normal.norm();
    normal /= (2.0*area); //normalization
    if(normal[2]<0)
      normal *= -1.0; //POINTING IN +Z DIRECTION (for this specific structure)

    g0 = (xyz0[n1]+xyz0[n2]+xyz0[n3])/3.0; //center of triangle in the original configuration

    experiment_time = time; //TODO: CHECK WITH JASON

    p = GetPressure(g0, experiment_time);
 
    nodal_load = 1.0/3.0*p*area*normal;
    F[n1] += nodal_load;
    F[n2] += nodal_load;
    F[n3] += nodal_load;
    
  }

}

//------------------------------------------------------------

void
MyForcesCalculator::LoadData()
{

  ifstream infile;
  infile.open("ForceCalculator/FlexAl_ptime_DH1.txt");
  if(!infile.is_open()) {
    fprintf(stdout,"*** Error: Cannot open ptime file.\n");
    exit(-1);
  }
  double time;
  while(infile >> time)
    ptime.push_back(time);
  infile.close();
  fprintf(stdout, "- Read pressure time stamps. Size = %d, Last = %e.\n", ptime.size(), ptime.back());

  infile.open("ForceCalculator/FlexAl_pressure_DH1.txt");
  if(!infile.is_open()) {
    fprintf(stdout,"*** Error: Cannot open pressure data file.\n");
    exit(-1);
  }
  int Nsensor = 7;
  double p;
  string line, word;
  while(getline(infile,line)) {
    stringstream linestream(line);
    psensor.push_back(vector<double>());
    while(getline(linestream,word,','))
      psensor.back().push_back(stod(word));
    assert(psensor.back().size()==Nsensor);
  }
  assert(psensor.size() == ptime.size());
  fprintf(stdout, "- Read pressure measurements. Size = %d, Last = %e.\n",
          psensor.size(), psensor.back().back());
  infile.close();
 
}

//------------------------------------------------------------

double
MyForcesCalculator::GetPressure(Vec3D& point, double t)
{

  if(point[0]<=1e-2 || point[0]>=286.0 || point[1]>=317.5-1e-2)
    return 0.0; //not on the bottom surface

  double PI = 2.0*acos(0.0);
  double c = cos(20.0/180.0*PI);
  double xs[7];
  //TODO: CHECK THIS WITH JASON
  xs[0] = 0.5*49.265*c;
  for(int i=1; i<7; i++) 
    xs[i] = xs[i-1] + 33.9*c;

  
  int id1, id2;
  double c1, c2;  //two sensors & associated interpolation weights

  double x = point[0]; //extract the x coordinate
  
  if(x<=xs[0]) {
    id1 = id2 = 0;
    c1 = 1.0; c2 = 0.0;
  } 
  else if(x>=xs[6]) {
    id1 = id2 = 6;
    c1 = 1.0; c2 = 0;
  }
  else {
    for(int i=1; i<7; i++) {
      if(x<xs[i]) {
        id1 = i-1;
        id2 = i;
        c1 = (xs[i]-x)/(xs[i]-xs[i-1]);
        c2 = (x-xs[i-1])/(xs[i]-xs[i-1]);
        break;
      }
    }
  }
  

  // get sensor data at the right time
  int ta, tb;
  double cta, ctb;
  if(t<=ptime.front()) {
    ta = tb = 0;
    cta = 1.0; ctb = 0.0;
  }
  else if(t>=ptime.back()) {
    ta = tb = ptime.size()-1;
    cta = 1.0; ctb = 0.0;
  }
  else { 
    tb = lower_bound(ptime.begin(), ptime.end(), t) - ptime.begin();
    ta = tb - 1;
    cta = (ptime[tb]-t)/(ptime[tb]-ptime[ta]);
    ctb = (t-ptime[ta])/(ptime[tb]-ptime[ta]);
  }


  double p1 = cta*psensor[ta][id1] + ctb*psensor[tb][id1];
  double p2 = cta*psensor[ta][id2] + ctb*psensor[tb][id2];

  return c1*p1 + c2*p2;

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

