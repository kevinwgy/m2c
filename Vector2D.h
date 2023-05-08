/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _VECTOR2D_H_
#define _VECTOR2D_H_

#include <cstdio>
#include <cmath>
#include <algorithm>

//------------------------------------------------------------------------------
// integer
struct Int2 {
  int v[2];
  Int2() {v[0] = v[1] = 0;}
  Int2(int i1, int i2) {v[0] = i1; v[1] = i2;}
  Int2(const Int2& v2) {for(int i=0; i<2; i++) v[i] = v2.v[i];}
  Int2 &operator=(const Int2& v2) {for(int i=0; i<2; i++) v[i] = v2.v[i]; return *this;}
  Int2 operator+(const Int2& v2) {Int2 res; for(int i=0; i<2; i++) res.v[i] = v[i] + v2.v[i]; return res;}
  Int2 operator-(const Int2& v2) {Int2 res; for(int i=0; i<2; i++) res.v[i] = v[i] - v2.v[i]; return res;}
  int &operator[](int i) {return v[i];}
  int operator[](int i) const {return v[i];}
};

//------------------------------------------------------------------------------
// real number 
struct Vec2D {

  double v[2];
  Vec2D() {v[0] = v[1] = 0.0; }
  Vec2D(double x[2]) { v[0] = x[0]; v[1] = x[1];}
  Vec2D(double x, double y) { v[0] = x; v[1] = y;}
  Vec2D(const Vec2D &v2) { v[0] = v2.v[0]; v[1] = v2.v[1]; }
  Vec2D(const Int2 &v2) { v[0] = v2.v[0]; v[1] = v2.v[1]; }
  Vec2D(double x) { v[0] = v[1] = x; }
  ~Vec2D() {}

  Vec2D &operator=(const double);
  Vec2D &operator=(const Vec2D &);
  Vec2D &operator+=(const Vec2D &);
  Vec2D &operator+=(const double &);
  Vec2D &operator-=(const Vec2D &);
  Vec2D &operator-=(const double &);
  Vec2D &operator*=(double);
  Vec2D &operator/=(double);
  Vec2D operator/(double) const;

  Vec2D operator+(const Vec2D &) const;
  Vec2D operator-(const Vec2D &) const;
  Vec2D operator-() const;

  double operator*(const Vec2D &) const;

  operator double*() { return v; }

  double &operator[](int i) { return v[i]; }
  double operator[](int i) const { return v[i]; }

  void print(const char *msg = "") { fprintf(stdout, "%s(%e %e)\n", msg, v[0], v[1]); }

  double norm() { return(sqrt(v[0]*v[0]+v[1]*v[1])); }

  double norm1() { return fabs(v[0])+fabs(v[1]); }
  double norm2() { return norm(); }
  double norminf() { return std::max(fabs(v[0]), fabs(v[1])); }
};

//------------------------------------------------------------------------------

inline 
Vec2D &Vec2D::operator=(const double v2)
{

  v[0] = v2;
  v[1] = v2;

  return *this;

}

//------------------------------------------------------------------------------

inline 
Vec2D &Vec2D::operator=(const Vec2D &v2)
{

  v[0] = v2.v[0];
  v[1] = v2.v[1];

  return *this;

}

//------------------------------------------------------------------------------

inline 
Vec2D &Vec2D::operator+=(const Vec2D &v2)
{

  v[0] += v2.v[0];
  v[1] += v2.v[1];

  return *this;

}

//------------------------------------------------------------------------------

inline
Vec2D &Vec2D::operator+=(const double &c)
{

  v[0] += c;
  v[1] += c;

  return *this;

}

//------------------------------------------------------------------------------

inline 
Vec2D &Vec2D::operator-=(const Vec2D &v2)
{

  v[0] -= v2.v[0];
  v[1] -= v2.v[1];

  return *this;

}

//------------------------------------------------------------------------------

inline
Vec2D &Vec2D::operator-=(const double &c)
{

  v[0] -= c;
  v[1] -= c;

  return *this;

}

//------------------------------------------------------------------------------

inline 
Vec2D &Vec2D::operator*=(double cst)
{

  v[0] *= cst;
  v[1] *= cst;

  return *this;

}

//------------------------------------------------------------------------------

inline
Vec2D &Vec2D::operator/=(double cst)
{
  cst = 1.0/cst;

  v[0] *= cst;
  v[1] *= cst;

  return *this;
}

//------------------------------------------------------------------------------

inline
Vec2D Vec2D::operator/(double cst)  const
{
  Vec2D vnew;
  vnew[0] = v[0] / cst;
  vnew[1] = v[1] / cst;

  return vnew;
}

//------------------------------------------------------------------------------

inline 
Vec2D Vec2D::operator+(const Vec2D &v2) const
{

  Vec2D res;

  res.v[0] = v[0]+v2.v[0];
  res.v[1] = v[1]+v2.v[1];

  return res;

}

//------------------------------------------------------------------------------

inline 
Vec2D Vec2D::operator-(const Vec2D &v2) const
{

  Vec2D res;

  res.v[0] = v[0]-v2.v[0];
  res.v[1] = v[1]-v2.v[1];

  return res;

}

//------------------------------------------------------------------------------

inline 
Vec2D Vec2D::operator-() const
{
  Vec2D res;

  res.v[0] = -v[0];
  res.v[1] = -v[1];

  return res;
}

//------------------------------------------------------------------------------

inline 
double Vec2D::operator*(const Vec2D &v2) const
{

  return v[0]*v2.v[0] + v[1]*v2.v[1];

}

//------------------------------------------------------------------------------

inline 
Vec2D operator*(double c, const Vec2D &v)
{

  Vec2D res;

  res.v[0] = c*v.v[0];
  res.v[1] = c*v.v[1];

  return res;

}

inline 
Vec2D operator*(const Vec2D &v,double c )
{

  Vec2D res;

  res.v[0] = c*v.v[0];
  res.v[1] = c*v.v[1];

  return res;

}

//inline double min(double x,double y) { return (x<y?x:y); }
//inline double max(double x,double y) { return (x>y?x:y); }

inline
Vec2D min( const Vec2D& a, const Vec2D& b) {

  return Vec2D( std::min(a[0],b[0]),
                std::min(a[1],b[1]));
}

inline
Vec2D max( const Vec2D& a, const Vec2D& b) {

  return Vec2D( std::max(a[0],b[0]),
                std::max(a[1],b[1]));
}

//------------------------------------------------------------------------------

// An instantiation of "Obj" in KDTree.h
class PointIn2D {
public:
  int id;
  Vec2D x;
public:
  PointIn2D() {}
  PointIn2D(int i, Vec2D xin) {id = i; x = xin;}
  double val(int i) const {return x[i];}
  double width([[maybe_unused]] int i) const {return 0.0;}
  int pid() const {return id;}
};



#endif
