/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _USER_DEFINED_DYNAMICS_H_
#define _USER_DEFINED_DYNAMICS_H_

/*****************************************************************
 * class UserDefinedDynamics is an interface for user to provide a
 * dynamics calculator through dynamic loading (see, e.g., dlopen)
 ****************************************************************/
class UserDefinedDynamics {

public:

  UserDefinedDynamics() {}
  virtual ~UserDefinedDynamics() {}

  virtual void GetUserDefinedDynamics(double time, double dt, int nNodes, double *X0, double *X,//these are inputs
                                       double *disp, double *velo/*output*/) = 0;

};

typedef UserDefinedDynamics* CreateUDD();
typedef void                 DestroyUDD(UserDefinedDynamics*);

#endif
