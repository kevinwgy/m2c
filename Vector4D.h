#ifndef _VECTOR4D_H_
#define _VECTOR4D_H_

#include <cstdio>
#include <cmath>
#include <algorithm>

//------------------------------------------------------------------------------
// real number 
struct Vec4D {

  double v[4];

  Vec4D() {v[0] = v[1] = v[2] = v[3] = 0.0; }
  Vec4D(double x[4]) { v[0] = x[0]; v[1] = x[1]; v[2] = x[2]; v[3] = x[3];}
  Vec4D(double x, double y, double z, double zz) { v[0] = x; v[1] = y; v[2] = z; v[3] = zz;}
  Vec4D(const Vec4D &v2) { v[0] = v2.v[0]; v[1] = v2.v[1]; v[2] = v2.v[2]; v[3] = v2.v[3]}
  Vec4D(double x) { v[0] = v[1] = v[2] = v[3] = x; }
  ~Vec4D() {}

  Vec4D &operator=(const double);
  Vec4D &operator=(const Vec4D &);
  Vec4D &operator+=(const Vec4D &);
  Vec4D &operator+=(const double &);
  Vec4D &operator-=(const Vec4D &);
  Vec4D &operator-=(const double &);
  Vec4D &operator*=(double);
  Vec4D &operator/=(double);
  Vec4D operator/(double) const;

  Vec4D operator+(const Vec4D &) const;
  Vec4D operator-(const Vec4D &) const;
  Vec4D operator-() const;
  Vec4D operator^(const Vec4D &) const;

  double operator*(const Vec4D &) const;

  operator double*() { return v; }

  double &operator[](int i) { return v[i]; }
  double operator[](int i) const { return v[i]; }

  void print(const char *msg = "") { fprintf(stdout, "%s(%e %e %e %e)\n", msg, v[0], v[1], v[2], v[3]); }

  double norm() { return(sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]+v[3]*v[3])); }
};

//------------------------------------------------------------------------------

inline 
Vec4D &Vec4D::operator=(const double v2)
{

  v[0] = v2;
  v[1] = v2;
  v[2] = v2;
  v[3] = v2;

  return *this;

}

//------------------------------------------------------------------------------

inline 
Vec4D &Vec4D::operator=(const Vec4D &v2)
{

  v[0] = v2.v[0];
  v[1] = v2.v[1];
  v[2] = v2.v[2];
  v[3] = v2.v[3];

  return *this;

}

//------------------------------------------------------------------------------

inline 
Vec4D &Vec4D::operator+=(const Vec4D &v2)
{

  v[0] += v2.v[0];
  v[1] += v2.v[1];
  v[2] += v2.v[2];
  v[3] += v2.v[3];

  return *this;

}

//------------------------------------------------------------------------------

inline
Vec4D &Vec4D::operator+=(const double &c)
{

  v[0] += c;
  v[1] += c;
  v[2] += c;
  v[3] += c;

  return *this;

}

//------------------------------------------------------------------------------

inline 
Vec4D &Vec4D::operator-=(const Vec4D &v2)
{

  v[0] -= v2.v[0];
  v[1] -= v2.v[1];
  v[2] -= v2.v[2];
  v[3] -= v2.v[3];

  return *this;

}

//------------------------------------------------------------------------------

inline
Vec4D &Vec4D::operator-=(const double &c)
{

  v[0] -= c;
  v[1] -= c;
  v[2] -= c;
  v[3] -= c;

  return *this;

}

//------------------------------------------------------------------------------

inline 
Vec4D &Vec4D::operator*=(double cst)
{

  v[0] *= cst;
  v[1] *= cst;
  v[2] *= cst;
  v[3] *= cst;

  return *this;

}

//------------------------------------------------------------------------------

inline
Vec4D &Vec4D::operator/=(double cst)
{
  cst = 1.0/cst;

  v[0] *= cst;
  v[1] *= cst;
  v[2] *= cst;
  v[3] *= cst;

  return *this;
}

//------------------------------------------------------------------------------

inline
Vec4D Vec4D::operator/(double cst)  const
{
  Vec4D vnew;
  vnew[0] = v[0] / cst;
  vnew[1] = v[1] / cst;
  vnew[2] = v[2] / cst;
  vnew[3] = v[3] / cst;

  return vnew;
}

//------------------------------------------------------------------------------

inline 
Vec4D Vec4D::operator+(const Vec4D &v2) const
{

  Vec4D res;

  res.v[0] = v[0]+v2.v[0];
  res.v[1] = v[1]+v2.v[1];
  res.v[2] = v[2]+v2.v[2];
  res.v[3] = v[3]+v2.v[3];

  return res;

}

//------------------------------------------------------------------------------

inline 
Vec4D Vec4D::operator-(const Vec4D &v2) const
{

  Vec4D res;

  res.v[0] = v[0]-v2.v[0];
  res.v[1] = v[1]-v2.v[1];
  res.v[2] = v[2]-v2.v[2];
  res.v[3] = v[3]-v2.v[3];

  return res;

}

//------------------------------------------------------------------------------

inline 
Vec4D Vec4D::operator-() const
{
  Vec4D res;

  res.v[0] = -v[0];
  res.v[1] = -v[1];
  res.v[2] = -v[2];
  res.v[3] = -v[3];

  return res;
}

//------------------------------------------------------------------------------

inline 
double Vec4D::operator*(const Vec4D &v2) const
{

  return v[0]*v2.v[0] + v[1]*v2.v[1] + v[2]*v2.v[2] + v[3]*v2.v[3];

}

//------------------------------------------------------------------------------

inline 
Vec4D operator*(double c, const Vec4D &v)
{

  Vec4D res;

  res.v[0] = c*v.v[0];
  res.v[1] = c*v.v[1];
  res.v[2] = c*v.v[2];
  res.v[3] = c*v.v[3];

  return res;

}

inline 
Vec4D operator*(const Vec4D &v,double c )
{

  Vec4D res;

  res.v[0] = c*v.v[0];
  res.v[1] = c*v.v[1];
  res.v[2] = c*v.v[2];
  res.v[3] = c*v.v[3];

  return res;

}

//------------------------------------------------------------------------------

#endif
