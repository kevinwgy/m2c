#pragma once

namespace MathTools {

//! Get the k-th binary digit (0 or 1) of integer n
inline bool GetBit(int n, int k) {
  return (n >> k) & 1;
}


}
