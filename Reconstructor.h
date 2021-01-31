#ifndef _RECONSTRUCTOR_H_
#define _RECONSTRUCTOR_H_
#include <IoData.h>
#include <SpaceVariable.h>
using std::min;
using std::max;

/*******************************************
 * class Reconstructor handles spatial
 * reconstruction and slope limiters 
 ******************************************/
class Reconstructor
{
  IoData&   iod;

  //! Mesh info
  SpaceVariable3D &delta_xyz;
  
  /** Internal variables
   *  CoeffA, CoeffB, and CoeffK are only defined within the real domain
   */
  SpaceVariable3D CoeffA;  //!< constant coefficient for non-uniform meshes, dim=3 for 3D 
  SpaceVariable3D CoeffB;  //!< constant coefficient for non-uniform meshes, dim=3 for 3D 
  SpaceVariable3D CoeffK;  //!< the integer k in Zeng's 2016 paper for Van Albada (stored as real numbers)


public:
  Reconstructor(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_, 
                SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_);

  ~Reconstructor();

  void Setup(); //!< compute AB and K

  /** Linear reconstruction in 3D 
    * Although the input argument is denoted by U, this function can be applied to
    * any other form (e.g., primitive state variables, characteristic varaibles, etc.*/
  void Reconstruct(SpaceVariable3D &U, SpaceVariable3D &Ul, SpaceVariable3D &Ur, 
           SpaceVariable3D &Ub, SpaceVariable3D &Ut, SpaceVariable3D &Uk, SpaceVariable3D &Uf);

  void Destroy(); //!< destroy the SpaceVariable3D variables

private:

  int CalculateSlopeLimiterCoefficientK(double A, double B); //calculate k in Zeng's 2016 paper

  
  /** slope limiter functions
   * (frequently called; need to be as fast as possible)
   */
  inline double GeneralizedMinMod(double A, double B, double alpha, double theta) {
    return max(0.0, min(alpha*theta, min(B/(A+1.0)*(theta+1.0), alpha) ) );
  }

  inline double VanAlbada(double A, double B, int k, double theta) {
    double theta_power_k = pow(theta,k);
    return B*(theta_power_k + theta)/(theta_power_k + A);
  }

  inline double ModifiedVanAlbada(double A, double B, int k, double theta) {
    if(theta <= 0) return 0.0;
    else {
      double theta_power_k = pow(theta,k);
      return B*(theta_power_k + theta)/(theta_power_k + A);
    }
  }
  
};

#endif
