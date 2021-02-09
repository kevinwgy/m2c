#ifndef _RECONSTRUCTOR_H_
#define _RECONSTRUCTOR_H_
#include <IoData.h>
#include <SpaceVariable.h>
using std::min;
using std::max;

/*******************************************
 * class Reconstructor handles spatial
 * reconstruction and slope limiters 
 * (for variables of any dimension)
 ******************************************/
class Reconstructor
{
  ReconstructionData& iod_rec;

  //! Mesh info
  SpaceVariable3D &delta_xyz;
  
  /** Internal variables
   *  CoeffA, CoeffB, and CoeffK are only defined within the real domain
   */
  SpaceVariable3D CoeffA;  //!< constant coefficient for non-uniform meshes, dim=3 for 3D 
  SpaceVariable3D CoeffB;  //!< constant coefficient for non-uniform meshes, dim=3 for 3D 
  SpaceVariable3D CoeffK;  //!< the integer k in Zeng's 2016 paper for Van Albada (stored as real numbers)


public:
  Reconstructor(MPI_Comm &comm_, DataManagers3D &dm_all_, ReconstructionData &iod_rec_, 
                SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_);

  ~Reconstructor();

  void Setup(); //!< compute AB and K

  /** Linear reconstruction in 3D 
    * Although the input argument is denoted by U, this function can be applied to
    * any other form (e.g., primitive state variables, characteristic varaibles, etc.*/
  void Reconstruct(SpaceVariable3D &U, SpaceVariable3D &Ul, SpaceVariable3D &Ur, 
           SpaceVariable3D &Ub, SpaceVariable3D &Ut, SpaceVariable3D &Uk, SpaceVariable3D &Uf);

  void ReconstructIn1D(int dir/*0~x,1~y,2~z*/, SpaceVariable3D &U, SpaceVariable3D &Um, SpaceVariable3D &Up
                       SpaceVariable3D *Slope = NULL);

  void Destroy(); //!< destroy the SpaceVariable3D variables

private:

  int CalculateSlopeLimiterCoefficientK(double A, double B); //calculate k in Zeng's 2016 paper

  
  /** slope limiter functions
   * (frequently called; need to be as fast as possible)
   */
  inline double GeneralizedMinMod(double A, double B, double alpha, double dq0, double dq1) {
    if(dq0<0.0 && dq1<0.0) 
      return max(alpha*dq0, max(B/(A+1.0)*(dq0+dq1), alpha*dq1));
    else if(dq0>0.0 && dq1>0.0) 
      return min(alpha*dq0, min(B/(A+1.0)*(dq0+dq1), alpha*dq1));
    else
      return 0.0;
  }

  inline double VanAlbada(double A, double B, int k, double dq0, double dq1) {
    if(dq0*dq1<=0) return 0.0;
    if(fabs(dq0*dq1)<1e-16) return 0.0;
    double dq0_power_k = pow(dq0,k);
    double dq1_power_k = pow(dq1,k);
    double denom = dq0_power_k + A*dq1_power_k;
    return B*(dq1*dq0_power_k + dq0*dq1_power_k)/(dq0_power_k + A*dq1_power_k);
  }

};

#endif
