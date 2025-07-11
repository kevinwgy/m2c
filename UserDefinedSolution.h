/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _USER_DEFINED_SOLUTION_H_
#define _USER_DEFINED_SOLUTION_H_

#include<vector>

struct Vec3D;
struct Vec5D;

/*****************************************************************
 * class UserDefinedSolution is an interface for user to modify 
 * state variable field (V), material ID, and level sets (Phi)
 * at the end of each time step.
 ****************************************************************/
class UserDefinedSolution{

public:

  UserDefinedSolution() {}
  virtual ~UserDefinedSolution() {}

  virtual void GetUserDefinedSolution(int i0, int j0, int k0, int imax, int jmax, int kmax,
                                   Vec3D*** coords, Vec5D*** v, //!<  converted from Vec5D (dof = 5)
                                   double*** id, std::vector<double***> phi, double time) = 0;
};

typedef UserDefinedSolution* CreateUDSL();
typedef void                 DestroyUDSL(UserDefinedSolution*);

#endif
