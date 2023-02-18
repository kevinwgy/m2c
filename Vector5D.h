/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _VECTOR5D_H_
#define _VECTOR5D_H_

#include <cstdio>
#include <cmath>
#include <algorithm>

//------------------------------------------------------------------------------
// real number 
struct Vec5D {

  double v[5];

  Vec5D() {v[0] = v[1] = v[2] = v[3] = v[4] = 0.0; }
  Vec5D(double x[5]) { v[0] = x[0]; v[1] = x[1]; v[2] = x[2]; v[3] = x[3]; v[4] = x[4];}
  Vec5D(double x, double y, double z, double zz, double zzz) { v[0] = x; v[1] = y; v[2] = z; v[3] = zz; v[4] = zzz;}
  Vec5D(const Vec5D &v2) { v[0] = v2.v[0]; v[1] = v2.v[1]; v[2] = v2.v[2]; v[3] = v2.v[3]; v[4] = v2.v[4];}
  Vec5D(double x) { v[0] = v[1] = v[2] = v[3] = v[4] = x; }
  ~Vec5D() {}

  Vec5D &operator=(const double);
  Vec5D &operator=(const Vec5D &);
  Vec5D &operator+=(const Vec5D &);
  Vec5D &operator+=(const double &);
  Vec5D &operator-=(const Vec5D &);
  Vec5D &operator-=(const double &);
  Vec5D &operator*=(double);
  Vec5D &operator/=(double);
  Vec5D operator/(double) const;

  Vec5D operator+(const Vec5D &) const;
  Vec5D operator-(const Vec5D &) const;
  Vec5D operator-() const;

  double operator*(const Vec5D &) const;

  operator double*() { return v; }  //convert Vec5D to double*

  double &operator[](int i) { return v[i]; }
  double operator[](int i) const { return v[i]; }

  void print(const char *msg = "") { fprintf(stdout, "%s(%e %e %e %e %e)\n", msg, v[0], v[1], v[2], v[3], v[4]); }

  double norm() { return(sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]+v[3]*v[3]+v[4]*v[4])); }

  double norm1() { return fabs(v[0])+fabs(v[1])+fabs(v[2])+fabs(v[3])+fabs(v[4]); }
  double norm2() { return norm(); }
  double norminf() { double vmax = fabs(v[0]); 
                     for(int i=1; i<5; i++) vmax = std::max(vmax, fabs(v[i]));
                     return vmax; }

};

//------------------------------------------------------------------------------

inline 
Vec5D &Vec5D::operator=(const double v2)
{

  v[0] = v2;
  v[1] = v2;
  v[2] = v2;
  v[3] = v2;
  v[4] = v2;

  return *this;

}

//------------------------------------------------------------------------------

inline 
Vec5D &Vec5D::operator=(const Vec5D &v2)
{

  v[0] = v2.v[0];
  v[1] = v2.v[1];
  v[2] = v2.v[2];
  v[3] = v2.v[3];
  v[4] = v2.v[4];

  return *this;

}

//------------------------------------------------------------------------------

inline 
Vec5D &Vec5D::operator+=(const Vec5D &v2)
{

  v[0] += v2.v[0];
  v[1] += v2.v[1];
  v[2] += v2.v[2];
  v[3] += v2.v[3];
  v[4] += v2.v[4];

  return *this;

}

//------------------------------------------------------------------------------

inline
Vec5D &Vec5D::operator+=(const double &c)
{

  v[0] += c;
  v[1] += c;
  v[2] += c;
  v[3] += c;
  v[4] += c;

  return *this;

}

//------------------------------------------------------------------------------

inline 
Vec5D &Vec5D::operator-=(const Vec5D &v2)
{

  v[0] -= v2.v[0];
  v[1] -= v2.v[1];
  v[2] -= v2.v[2];
  v[3] -= v2.v[3];
  v[4] -= v2.v[4];

  return *this;

}

//------------------------------------------------------------------------------

inline
Vec5D &Vec5D::operator-=(const double &c)
{

  v[0] -= c;
  v[1] -= c;
  v[2] -= c;
  v[3] -= c;
  v[4] -= c;

  return *this;

}

//------------------------------------------------------------------------------

inline 
Vec5D &Vec5D::operator*=(double cst)
{

  v[0] *= cst;
  v[1] *= cst;
  v[2] *= cst;
  v[3] *= cst;
  v[4] *= cst;

  return *this;

}

//------------------------------------------------------------------------------

inline
Vec5D &Vec5D::operator/=(double cst)
{
  cst = 1.0/cst;

  v[0] *= cst;
  v[1] *= cst;
  v[2] *= cst;
  v[3] *= cst;
  v[4] *= cst;

  return *this;
}

//------------------------------------------------------------------------------

inline
Vec5D Vec5D::operator/(double cst)  const
{
  Vec5D vnew;
  vnew[0] = v[0] / cst;
  vnew[1] = v[1] / cst;
  vnew[2] = v[2] / cst;
  vnew[3] = v[3] / cst;
  vnew[4] = v[4] / cst;

  return vnew;
}

//------------------------------------------------------------------------------

inline 
Vec5D Vec5D::operator+(const Vec5D &v2) const
{

  Vec5D res;

  res.v[0] = v[0]+v2.v[0];
  res.v[1] = v[1]+v2.v[1];
  res.v[2] = v[2]+v2.v[2];
  res.v[3] = v[3]+v2.v[3];
  res.v[4] = v[4]+v2.v[4];

  return res;

}

//------------------------------------------------------------------------------

inline 
Vec5D Vec5D::operator-(const Vec5D &v2) const
{

  Vec5D res;

  res.v[0] = v[0]-v2.v[0];
  res.v[1] = v[1]-v2.v[1];
  res.v[2] = v[2]-v2.v[2];
  res.v[3] = v[3]-v2.v[3];
  res.v[4] = v[4]-v2.v[4];

  return res;

}

//------------------------------------------------------------------------------

inline 
Vec5D Vec5D::operator-() const
{
  Vec5D res;

  res.v[0] = -v[0];
  res.v[1] = -v[1];
  res.v[2] = -v[2];
  res.v[3] = -v[3];
  res.v[4] = -v[4];

  return res;
}

//------------------------------------------------------------------------------

inline 
double Vec5D::operator*(const Vec5D &v2) const
{

  return v[0]*v2.v[0] + v[1]*v2.v[1] + v[2]*v2.v[2] + v[3]*v2.v[3] + v[4]*v2.v[4];

}

//------------------------------------------------------------------------------

inline 
Vec5D operator*(double c, const Vec5D &v)
{

  Vec5D res;

  res.v[0] = c*v.v[0];
  res.v[1] = c*v.v[1];
  res.v[2] = c*v.v[2];
  res.v[3] = c*v.v[3];
  res.v[4] = c*v.v[4];

  return res;

}

inline 
Vec5D operator*(const Vec5D &v,double c )
{

  Vec5D res;

  res.v[0] = c*v.v[0];
  res.v[1] = c*v.v[1];
  res.v[2] = c*v.v[2];
  res.v[3] = c*v.v[3];
  res.v[4] = c*v.v[4];

  return res;

}

//------------------------------------------------------------------------------

#endif
