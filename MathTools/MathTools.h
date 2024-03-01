/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#pragma once

namespace MathTools {

//! Get the k-th binary digit (0 or 1) of integer n
inline bool GetBit(int n, int k) {
  return (n >> k) & 1;
}


}
