/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _USER_DEFINED_FORCES_H_
#define _USER_DEFINED_FORCES_H_

/*****************************************************************
 * class UserDefinedForces is an interface for user to provide a
 * force calculator through dynamic loading (see, e.g., dlopen)
 ****************************************************************/
class UserDefinedForces {

public:

  UserDefinedForces() {}
  virtual ~UserDefinedForces() {}

  //! Calculates a force (fx,fy,fz) for each node; elements are triangles
  virtual void GetUserDefinedForces(double time, int nNodes, double *X0, double *X,
                                    int nElems, int *elems, double* nodal_forces) = 0;

};

typedef UserDefinedForces* CreateUDF();
typedef void               DestroyUDF(UserDefinedForces*);

#endif
