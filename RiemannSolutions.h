/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _RIEMANN_SOLUTIONS_H_
#define _RIEMANN_SOLUTIONS_H_

#include <Vector3D.h>
#include <Vector5D.h>

class RiemannSolutions {

public: //for the moment, keep everything public

  //Solutions of exact Riemann problems
  std::map<Int3/*k,j,i*/, std::pair<Vec5D/*state*/, int/*ID*/> > left;   //x-
  std::map<Int3/*k,j,i*/, std::pair<Vec5D/*state*/, int/*ID*/> > right;  //x+
  std::map<Int3/*k,j,i*/, std::pair<Vec5D/*state*/, int/*ID*/> > bottom; //y-
  std::map<Int3/*k,j,i*/, std::pair<Vec5D/*state*/, int/*ID*/> > top;    //y+
  std::map<Int3/*k,j,i*/, std::pair<Vec5D/*state*/, int/*ID*/> > back;   //z-
  std::map<Int3/*k,j,i*/, std::pair<Vec5D/*state*/, int/*ID*/> > front;  //z+

  RiemannSolutions() {}
  ~RiemannSolutions() {}

  void Clear() {
    left.clear(); right.clear(); bottom.clear(); top.clear(); back.clear(); front.clear();
  }
};


#endif
