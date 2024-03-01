/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _USER_DEFINED_STATE_H_
#define _USER_DEFINED_STATE_H_

struct Vec3D;
struct Vec5D;
/*****************************************************************
 * class UserDefinedState is an interface for user to provide a
 * state variable field (V) and material ID through dynamic loading (see, e.g., dlopen)
 ****************************************************************/
class UserDefinedState {

public:

  UserDefinedState() {}
  virtual ~UserDefinedState() {}

  virtual void GetUserDefinedState(int i0, int j0, int k0, int imax, int jmax, int kmax,
                                   Vec3D*** coords, Vec5D*** v, //!<  converted from Vec5D (dof = 5)
                                   double*** id) = 0;
};

typedef UserDefinedState* CreateUDS();
typedef void              DestroyUDS(UserDefinedState*);

#endif
