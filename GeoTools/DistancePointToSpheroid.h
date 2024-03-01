/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#pragma once

#include <assert.h>
#include <Vector2D.h>
#include <Vector3D.h>

namespace GeoTools {

/**********************************************************
 * Based on the algorithms described in David Eberly, 
 * Distance from a Point to an Ellipse, an Ellipsoid, or a
 * Hyperellipsoid, a documentation of Geometric Tools
 * https://www.geometrictools.com
 * *******************************************************/
 
struct Ellipse {

  //! following the notations in Eberly's document
  //! convention: e0 >= e1
  Vec2D C; //!< center
  Vec2D U0, U1; //!< axes (must be normalized, and orthogonal)
  double e0, e1; //!< semi axis length

  Vec2D GetPointCoordinates(Vec2D Q) { //!< coords of Q using C, U0, and U1
    Vec2D y;
    y[0] = U0*(Q-C);
    y[1] = U1*(Q-C);
    return y;
  }

  int GetPointLocation(Vec2D Q) { //!< 0 = on the ellipse, -1 = inside, 1 = outside
    Vec2D y = GetPointCoordinates(Q); 
    double p = y[0]*y[0]/(e0*e0) + y[1]*y[1]/(e1*e1) - 1.0;
    if(p<0) return -1;
    if(p>0) return 1;
    return 0;
  }
};

//***************************************************************************

struct Ellipsoid {

  //! following the notations in Eberly's document
  //! convention: e0 >= e1 >= e2
  Vec3D C; //!< center
  Vec3D U0, U1, U2; //!< axes (must be normalized, and orthogonal)
  double e0, e1, e2; //!< semi axis length

  void SetupSpheroid(Vec3D C_, Vec3D Umain, double semi_length, double radius) {
    C = C_;
    U0 = Umain/Umain.norm();
    e0 = semi_length;
    e1 = e2 = radius;

    U1 = 0.0;
    U2 = 0.0;
    // calculate U1 and U2
    bool done = false;
    for(int i=0; i<3; i++) {
      if(U0[i]==0) {
        U1[i] = 1.0; //got U1
        bool gotU2 = false;
        for(int j=i+1; j<3; j++) {
          if(U0[j]==0) {
            U2[j] = 1.0; //got U2;
            gotU2 = true;
            break;
          }
        }
        if(!gotU2) {
          int i1 = (i+1) % 3;
          int i2 = (i+2) % 3;
          U2[i1] = -U0[i2];
          U2[i2] = U0[i1];
          U2 /= U2.norm();
        }
        done = true;
        break;
      }
    }
    if(!done) { //!< all the three components of U0 are nonzero
      U1[0] = 1.0;  
      U1[1] = 0.0;
      U1[2] = -U0[0]/U0[2];
      U1 /= U1.norm();
      U2[0] = 1.0;
      U2[1] = -(U0[2]*U0[2]/U0[0] + U0[0])/U0[1];
      U2[2] = U0[2]/U0[0];
      U2 /= U2.norm();
    }

    if(e0<e2) { //!< make sure e0 >= e1 >= e2
      e2 = e0;
      e0 = e1;
      Vec3D tmp(U0);
      U0 = U2;
      U2 = tmp;
    }
  }

  Vec3D GetPointCoordinates(Vec3D Q) { //!< coords of Q using C, U0, U1, U2
    Vec3D y;
    y[0] = U0*(Q-C);
    y[1] = U1*(Q-C);
    y[2] = U2*(Q-C);
    return y;
  }

  int GetPointLocation(Vec3D Q) { //!< 0 = on the ellipsoid, -1 = inside, 1 = outside
    Vec3D y = GetPointCoordinates(Q); 
    double p = y[0]*y[0]/(e0*e0) + y[1]*y[1]/(e1*e1) + y[2]*y[2]/(e2*e2) - 1.0;
    if(p<0) return -1;
    if(p>0) return 1;
    return 0;
  }

};

//***************************************************************************

class DistanceFromPointToEllipse {

  Ellipse ellipse;

public:

  /****************************************************************************
   * Constructor
   ***************************************************************************/
  DistanceFromPointToEllipse() {}
  ~DistanceFromPointToEllipse() {}

  // U0_ and U1_ should have been normalized!
  void Setup(double *C_, double *U0_, double *U1_, double e0_, double e1_) {

    for(int i=0; i<2; i++)
      ellipse.C[i] = C_[i]; 

    // make sure e0 and e1 are non-negative
    if(e0_<0) e0_ = -e0_;
    if(e1_<0) e1_ = -e1_;

    if(e0_<e1_) {//swap the two to make sure e0 >= e1 (convention)
      ellipse.e0 = e1_; 
      ellipse.e1 = e0_;
      for(int i=0; i<2; i++) {
        ellipse.U0[i] = U1_[i];
        ellipse.U1[i] = U0_[i];
      }
    } else {
      ellipse.e0 = e0_; 
      ellipse.e1 = e1_;
      for(int i=0; i<2; i++) {
        ellipse.U0[i] = U0_[i];
        ellipse.U1[i] = U1_[i];
      }
    }
  }

  /****************************************************************************
   * calculate signed distance (>0 outside) from an arbitrary point Q to the ellipse, 
   * including a closest point (P) (may not be unique)
   ***************************************************************************/
  double Calculate(double *Q_, double *P = NULL) {

    Vec2D Q(Q_);

    Vec2D y = ellipse.GetPointCoordinates(Q); //coords of Q using C, U0, and U1

    // consider the first quadrant
    Vec2D factor(1.0,1.0);
    for(int i=0; i<2; i++) {
      if(y[i]<0) {
        y[i] = -y[i];
        factor[i] = -1.0;
      }
    }

    // Now, implement the algorithm in "list 1" of Eberly
    double distance;
    Vec2D x;
    double e0 = ellipse.e0, e1 = ellipse.e1;
    if(y[1]>0) {
      if(y[0]>0) {
        //Compute the unique root tbar of F(t) on (-e1*e1, +inf)
        double t0 = -e1*e1 + e1*y[1];
        double t1 = -e1*e1 + sqrt(e0*e0*y[0]*y[0] + e1*e1*y[1]*y[1]);
        double f0 = F(t0, y);
        double f1 = F(t1, y);
        assert(f0*f1<=0); //initial bracketing interval
        double tbar = GetRoot(y, t0, t1, f0, f1, 1.0e-8);
        x[0] = e0*e0*y[0]/(tbar + e0*e0);
        x[1] = e1*e1*y[1]/(tbar + e1*e1);
        distance = (x-y).norm();
      } else {//y[0] == 0
        x[0] = 0;
        x[1] = e1;
        distance = fabs(y[1] - e1);
      }
    } else {//y[1] == 0
      if(y[0] < (e0*e0 - e1*e1)/e0) {
        x[0] = e0*e0*y[0]/(e0*e0-e1*e1);
        x[1] = e1*sqrt(1-(x[0]/e0)*(x[0]/e0));
        distance = (x-y).norm(); 
      } else {
        x[0] = e0;
        x[1] = 0;
        distance = fabs(y[0] - e0);
      }
    } 

    // Find the closest point in original (x,y,z) coordinates
    if(P != NULL) {
      for(int i=0; i<2; i++)
        x[i] *= factor[i];
      for(int i=0; i<2; i++)
        P[i] = ellipse.C[i] + x[0]*ellipse.U0[i] + x[1]*ellipse.U1[i];
    }

    if(ellipse.GetPointLocation(Q) == -1) //inside
      distance *= -1.0;

    return distance;
  }

/****************************************************************************
 * Internal functions
 ****************************************************************************/
private:
  /****************************************************************************
   * returns F(t), Eq. (11) of Eberly
   ***************************************************************************/
  double F(double t, Vec2D y) {
    double p1 = ellipse.e0*y[0]/(t + ellipse.e0*ellipse.e0);
    double p2 = ellipse.e1*y[1]/(t + ellipse.e1*ellipse.e1);
    return p1*p1 + p2*p2 - 1.0;
  }

  /****************************************************************************
   *Brent method for root-finding (bisection and Secant)
   ***************************************************************************/
  double GetRoot(Vec2D y, double t0, double t1, double f0, double f1, double tol) {

    int maxIts = 10000;
    double root;

    if(f0==0.0)
      root = t0;
    else if(f1==0.0)
      root = t1;
    else {
      double t2 = t1; //t2 is always the latest one
      double f2 = f1;
      int it;
      for(it = 0; it<maxIts; it++) {
        t2 = t2 - f2*(t1 - t0)/(f1 - f0); //secant method
        if(t2 >= t1 || t2 <= t0) //discard and switch to bisection
          t2 = 0.5*(t0+t1);
        f2 = F(t2, y);
        if(f2==0.0) {
          root = t2;
          break;
        }
        if(f2*f0<0) {
          t1 = t2;
          f1 = f2;
        } else {
          t0 = t2;
          f0 = f2;
        }
        if(t1-t0<tol) {
          root = 0.5*(t0+t1);
          break;
        }
      }

      if(it==maxIts) {
        fprintf(stdout,"*** Warning: Root-finding method failed to converge after %d iterations.\n", it);
        root = 0.5*(t0+t1);
      }
    }

    return root;
  }

};

//***************************************************************************

class DistanceFromPointToSpheroid {

  Ellipsoid ellipsoid;
  enum Type {NONE = 0, SPHERE = 1, PROLATE = 2, OBLATE = 3, SIZE = 4} type;

  DistanceFromPointToEllipse d2ellipse; 

public:

  /****************************************************************************
   * Constructor 
   ***************************************************************************/
  DistanceFromPointToSpheroid(double *C_, double *U0_, double *U1_, double *U2_,
                              double e0_, double e1_, double e2_) : type(NONE)
  {

    for(int i=0; i<3; i++)
      ellipsoid.C[i] = C_[i]; 

    // make sure e0, e1, e2 are non-negative
    if(e0_<0) e0_ = -e0_;
    if(e1_<0) e1_ = -e1_;
    if(e2_<0) e2_ = -e2_;

    // make sure e0 >= e1 >= e2;
    Vec3D e(e0_, e1_, e2_);
    Vec3D U[3];
    for(int i=0; i<3; i++) {
      U[0][i] = U0_[i];
      U[1][i] = U1_[i];
      U[2][i] = U2_[i];
    }

    int imax = 0, imin = 2;
    double emax = e[imax], emin = e[imin]; 
    for(int i=0; i<3; i++) {
      if(e[i]>emax) {emax = e[i]; imax = i;}
      if(e[i]<emin) {emin = e[i]; imin = i;}
    }
    int imid = 3 - imax - imin;
    ellipsoid.e0 = e[imax];
    ellipsoid.e1 = e[imid];
    ellipsoid.e2 = e[imin];
    ellipsoid.U0 = U[imax];
    ellipsoid.U1 = U[imid];
    ellipsoid.U2 = U[imin];

    // figure out the type
    if(ellipsoid.e0 == ellipsoid.e1 && ellipsoid.e1 == ellipsoid.e2)
      type = SPHERE;
    else if(ellipsoid.e0 == ellipsoid.e1)
      type = OBLATE;
    else if(ellipsoid.e1 == ellipsoid.e2)
      type = PROLATE;
    else {
      fprintf(stdout,"*** Error: Found a general ellipsoid instead of a spheroid (%e, %e, %e).\n", 
              ellipsoid.e0, ellipsoid.e1, ellipsoid.e2);
      type = NONE;
    }

    //setup d2ellipse
    double C_2D[2] = {0.0, 0.0};
    double U0_2D[2] = {1.0, 0.0};
    double U1_2D[2] = {0.0, 1.0};
    d2ellipse.Setup(C_2D, U0_2D, U1_2D, ellipsoid.e0, ellipsoid.e2); 
  }



  /****************************************************************************
   * Constructor
   ***************************************************************************/
  DistanceFromPointToSpheroid(double *C_, double *U0_, double semi_length, double radius) : type(NONE)
  {
    ellipsoid.SetupSpheroid(Vec3D(C_), Vec3D(U0_), semi_length, radius);

    // figure out the type
    if(ellipsoid.e0 == ellipsoid.e1 && ellipsoid.e1 == ellipsoid.e2)
      type = SPHERE;
    else if(ellipsoid.e0 == ellipsoid.e1)
      type = OBLATE;
    else if(ellipsoid.e1 == ellipsoid.e2)
      type = PROLATE;
    else {
      fprintf(stdout,"*** Error: Found a general ellipsoid instead of a spheroid (%e, %e, %e).\n", 
              ellipsoid.e0, ellipsoid.e1, ellipsoid.e2);
      type = NONE;
    }

    //setup d2ellipse
    double C_2D[2] = {0.0, 0.0};
    double U0_2D[2] = {1.0, 0.0};
    double U1_2D[2] = {0.0, 1.0};
    d2ellipse.Setup(C_2D, U0_2D, U1_2D, ellipsoid.e0, ellipsoid.e2); 

  }

  /***************************************************************************
   * calculate signed distance (>0 outside) from an arbitrary point Q to the ellipsoid, 
   * including a closest point (P) (may not be unique)
   ***************************************************************************/
  double Calculate(double *Q_, double *P = NULL) {

    Vec3D Q(Q_);
    Vec3D y = ellipsoid.GetPointCoordinates(Q);
 
    if(type == SPHERE) {
      double distance = y.norm() - ellipsoid.e0;
      if(P != NULL) {
        double r = ellipsoid.e0;
        Vec3D x = (r/y.norm())*y;
        for(int i=0; i<3; i++)
          P[i] = ellipsoid.C[i] + x[0]*ellipsoid.U0[i] + x[1]*ellipsoid.U1[i] 
               + x[2]*ellipsoid.U2[i];
      }
      return distance;
    }

    Vec2D y2;
    if(type == OBLATE) {
      y2[0] = sqrt(y[0]*y[0] + y[1]*y[1]);
      y2[1] = y[2];
    } else if(type == PROLATE) {
      y2[0] = y[0];
      y2[1] = sqrt(y[1]*y[1] + y[2]*y[2]);
    }

    double distance;
    if(P != NULL) {
      double P2[2];
      distance = d2ellipse.Calculate(y2, P2);
      Vec3D x; //coords of the closest point
      if(type == OBLATE) {
        x[0] = y2[0]==0 ? P2[0] : P2[0]*(y[0]/y2[0]);
        x[1] = y2[0]==0 ? 0.0 : P2[0]*(y[1]/y2[0]);
        x[2] = P2[1]; 
      } else if(type == PROLATE) {
        x[0] = P2[0]; 
        x[1] = y2[1]==0 ? P2[1] : P2[1]*(y[1]/y2[1]);
        x[2] = y2[1]==0 ? 0.0 : P2[1]*(y[2]/y2[1]);
      }
      for(int i=0; i<3; i++)
        P[i] = ellipsoid.C[i] + x[0]*ellipsoid.U0[i] + x[1]*ellipsoid.U1[i] 
             + x[2]*ellipsoid.U2[i];
    } 
    else
      distance = d2ellipse.Calculate(y2, NULL);
  
    return distance;

  }

};

}//end of namespace
