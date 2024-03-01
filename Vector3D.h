/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _VECTOR3D_H_
#define _VECTOR3D_H_

#include <cstdio>
#include <cmath>
#include <algorithm>

//------------------------------------------------------------------------------
// integer
struct Int3 {
  int v[3];
  Int3() {v[0] = v[1] = v[2] = 0;}
  Int3(int x[3]) { v[0] = x[0]; v[1] = x[1]; v[2] = x[2]; }
  Int3(int i1, int i2, int i3) {v[0] = i1; v[1] = i2; v[2] = i3;}
  Int3(int i) {v[0] = v[1] = v[2] = i;}
  Int3(const Int3& v2) {for(int i=0; i<3; i++) v[i] = v2.v[i];}
  Int3 &operator=(const int v2) {v[0] = v[1] = v[2] = v2; return *this;}
  Int3 &operator=(const Int3& v2) {for(int i=0; i<3; i++) v[i] = v2.v[i]; return *this;}
  Int3 &operator+=(const Int3 &v2) {for(int i=0; i<3; i++) v[i] += v2.v[i]; return *this;}
  Int3 operator+(const Int3& v2) {Int3 res; for(int i=0; i<3; i++) res.v[i] = v[i] + v2.v[i]; return res;}
  Int3 operator-(const Int3& v2) {Int3 res; for(int i=0; i<3; i++) res.v[i] = v[i] - v2.v[i]; return res;}
  int &operator[](int i) {return v[i];}
  int operator[](int i) const {return v[i];}
  operator int*() { return v; } //convert Int3D to int*

  bool operator<(const Int3& v2) const {
    if     (v[0]<v2.v[0]) {return true;} else if(v[0]>v2.v[0]) {return false;} 
    else if(v[1]<v2.v[1]) {return true;} else if(v[1]>v2.v[1]) {return false;} 
    else if(v[2]<v2.v[2]) {return true;} else return false;
  }

  bool operator==(const Int3& v2) const {
    if     (v[0] != v2.v[0])   return false;
    else if(v[1] != v2.v[1])   return false;
    else if(v[2] != v2.v[2])   return false;
    else                       return true;
  }
};

//------------------------------------------------------------------------------
// real number 
struct Vec3D {

  double v[3];

  Vec3D() {v[0] = v[1] = v[2] = 0.0; }
  Vec3D(double x[3]) { v[0] = x[0]; v[1] = x[1]; v[2] = x[2]; }
  Vec3D(double x, double y, double z) { v[0] = x; v[1] = y; v[2] = z; }
  Vec3D(const Vec3D &v2) { v[0] = v2.v[0]; v[1] = v2.v[1]; v[2] = v2.v[2]; }
  Vec3D(const Int3 &v2) { v[0] = v2.v[0]; v[1] = v2.v[1]; v[2] = v2.v[2]; }
  Vec3D(double x) { v[0] = v[1] = v[2] = x; }
  ~Vec3D() {}

  Vec3D &operator=(const double);
  Vec3D &operator=(const Vec3D &);
  Vec3D &operator+=(const Vec3D &);
  Vec3D &operator+=(const double &);
  Vec3D &operator-=(const Vec3D &);
  Vec3D &operator-=(const double &);
  Vec3D &operator*=(double);
  Vec3D &operator/=(double);
  Vec3D operator/(double) const;

  Vec3D operator+(const Vec3D &) const;
  Vec3D operator-(const Vec3D &) const;
  Vec3D operator-() const;
  Vec3D operator^(const Vec3D &) const;

  double operator*(const Vec3D &) const;

  operator double*() { return v; } //convert Vec3D to double*

  double &operator[](int i) { return v[i]; }
  double operator[](int i) const { return v[i]; }

  void print(const char *msg = "") { fprintf(stdout, "%s(%e %e %e)\n", msg, v[0], v[1], v[2]); }

  double norm() { return(sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])); }

  double norm1() { return fabs(v[0])+fabs(v[1])+fabs(v[2]); }
  double norm2() { return norm(); }
  double norminf() { return std::max(fabs(v[0]), std::max(fabs(v[1]), fabs(v[2]))); }
};

//------------------------------------------------------------------------------

inline 
Vec3D &Vec3D::operator=(const double v2)
{

  v[0] = v2;
  v[1] = v2;
  v[2] = v2;

  return *this;

}

//------------------------------------------------------------------------------

inline 
Vec3D &Vec3D::operator=(const Vec3D &v2)
{

  v[0] = v2.v[0];
  v[1] = v2.v[1];
  v[2] = v2.v[2];

  return *this;

}

//------------------------------------------------------------------------------

inline 
Vec3D &Vec3D::operator+=(const Vec3D &v2)
{

  v[0] += v2.v[0];
  v[1] += v2.v[1];
  v[2] += v2.v[2];

  return *this;

}

//------------------------------------------------------------------------------

inline
Vec3D &Vec3D::operator+=(const double &c)
{

  v[0] += c;
  v[1] += c;
  v[2] += c;

  return *this;

}

//------------------------------------------------------------------------------

inline 
Vec3D &Vec3D::operator-=(const Vec3D &v2)
{

  v[0] -= v2.v[0];
  v[1] -= v2.v[1];
  v[2] -= v2.v[2];

  return *this;

}

//------------------------------------------------------------------------------

inline
Vec3D &Vec3D::operator-=(const double &c)
{

  v[0] -= c;
  v[1] -= c;
  v[2] -= c;

  return *this;

}

//------------------------------------------------------------------------------

inline 
Vec3D &Vec3D::operator*=(double cst)
{

  v[0] *= cst;
  v[1] *= cst;
  v[2] *= cst;

  return *this;

}

//------------------------------------------------------------------------------

inline
Vec3D &Vec3D::operator/=(double cst)
{
  cst = 1.0/cst;

  v[0] *= cst;
  v[1] *= cst;
  v[2] *= cst;

  return *this;
}

//------------------------------------------------------------------------------

inline
Vec3D Vec3D::operator/(double cst)  const
{
  Vec3D vnew;
  vnew[0] = v[0] / cst;
  vnew[1] = v[1] / cst;
  vnew[2] = v[2] / cst;

  return vnew;
}

//------------------------------------------------------------------------------

inline 
Vec3D Vec3D::operator+(const Vec3D &v2) const
{

  Vec3D res;

  res.v[0] = v[0]+v2.v[0];
  res.v[1] = v[1]+v2.v[1];
  res.v[2] = v[2]+v2.v[2];

  return res;

}

//------------------------------------------------------------------------------

inline 
Vec3D Vec3D::operator-(const Vec3D &v2) const
{

  Vec3D res;

  res.v[0] = v[0]-v2.v[0];
  res.v[1] = v[1]-v2.v[1];
  res.v[2] = v[2]-v2.v[2];

  return res;

}

//------------------------------------------------------------------------------

inline 
Vec3D Vec3D::operator-() const
{
  Vec3D res;

  res.v[0] = -v[0];
  res.v[1] = -v[1];
  res.v[2] = -v[2];

  return res;
}

//------------------------------------------------------------------------------
// define vector cross product

inline 
Vec3D Vec3D::operator^(const Vec3D &v2) const
{

  Vec3D res;

  res.v[0] = v[1]*v2.v[2]-v[2]*v2.v[1];
  res.v[1] = v[2]*v2.v[0]-v[0]*v2.v[2];
  res.v[2] = v[0]*v2.v[1]-v[1]*v2.v[0];

  return res;

}

//------------------------------------------------------------------------------

inline 
double Vec3D::operator*(const Vec3D &v2) const
{

  return v[0]*v2.v[0] + v[1]*v2.v[1] + v[2]*v2.v[2];

}

//------------------------------------------------------------------------------

inline 
Vec3D operator*(double c, const Vec3D &v)
{

  Vec3D res;

  res.v[0] = c*v.v[0];
  res.v[1] = c*v.v[1];
  res.v[2] = c*v.v[2];

  return res;

}

inline 
Vec3D operator*(const Vec3D &v,double c )
{

  Vec3D res;

  res.v[0] = c*v.v[0];
  res.v[1] = c*v.v[1];
  res.v[2] = c*v.v[2];

  return res;

}

//inline double min(double x,double y) { return (x<y?x:y); }
//inline double max(double x,double y) { return (x>y?x:y); }

inline
Vec3D min( const Vec3D& a, const Vec3D& b) {

  return Vec3D( std::min(a[0],b[0]),
                std::min(a[1],b[1]),
                std::min(a[2],b[2]));
}

inline
Vec3D max( const Vec3D& a, const Vec3D& b) {

  return Vec3D( std::max(a[0],b[0]),
                std::max(a[1],b[1]),
                std::max(a[2],b[2]));
}

//------------------------------------------------------------------------------

// An instantiation of "Obj" in KDTree.h
class PointIn3D {
public:
  int id;
  Vec3D x;
public:
  PointIn3D() {}
  PointIn3D(int i, Vec3D &xin) {id = i; x = xin;}
  double val(int i) const {return x[i];}
  double width([[maybe_unused]] int i) const {return 0.0;}
  int pid() const {return id;}
};



#endif
