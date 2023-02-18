/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _ATOMIC_IONIZATION_DATA_H_
#define _ATOMIC_IONIZATION_DATA_H_

#include<VarFcnBase.h>
#include<tuple>
#include<boost/math/interpolators/cubic_b_spline.hpp>  //spline interpolation
#include<mpi.h>

/***********************************************************************
 * Class AtomicIonizationData stores the atomic ionization parameters
 * (excitation energy, ionization energy, etc.) for a chemical
 * element. It is also responsible for computing the energy partition 
 * function, by interpolation (faster, but requires storage and setup) 
 * or on-the-fly (slower, w/o overhead)
***********************************************************************/

class AtomicIonizationData { 

public: //!< data is public

  double molar_fraction;
  double molar_mass;
  int atomic_number;
  std::vector<double> I; //!< ionization energy
  std::vector<std::vector<int> > g; //!< degeneracy (=2l+1, l: angular momentum)
  std::vector<std::vector<double> > E; //!< excitation energy
  int rmax; //max charge state 

private:

  //! Whether the ideal or non-ideal ionization model is applied
  bool ideal;

  //! Constants
  double h; //planck_constant;
  double e; //electron_charge;
  double me; //electron_mass;
  double kb; //boltzmann_constant;

  //! Method for partition function evaluation
  int interpolation; //!< 0~not used, 1~cubic spline, 2~piecewise linear
  std::vector<std::vector<std::vector<double> > > Us; /**< partition function at sample Temperatures \n
                                                        *  Us[r][k][i] --> r: charge, k: reduction of summation upper \n
                                                        *  bound due to non-ideal depression (for the ideal case, only
                                                        *  k = 0 is stored, i: i-th sample temperature) */
  std::vector<std::tuple<double,double,double> > UsCoeffs; //!< for each r: ("factor", "expmin", "delta_exp")
  std::vector<int> max_terms; //!< number of terms in summation, up to E[r][n] <= I[r]

  std::vector<std::vector<boost::math::cubic_b_spline<double>* > > spline; 
                                         //!< spline[r][k] --> r: charge, k: reduction of summation upperbound
  int sample_size; //!< the size of Us[r], r = 0,1,...,rmax
  double sample_Tmin, sample_Tmax;

public:

  AtomicIonizationData();
  ~AtomicIonizationData();

  void Setup(AtomicIonizationModel* iod_aim, double h_, double e_, double me_, double kb_,
             bool ideal_, int interp, double sample_Tmin_, double sample_Tmax_, int sample_size_,
             MPI_Comm* comm);

  double CalculatePartitionFunction(int r, double T, double deltaI = 0.0);

private: 

  template<class T>
  void GetDataInFile(std::fstream& file, vector<T> &X, int MaxCount, bool non_negative);

  //! Interpolation for both ideal and non-ideal ionization models 
  void InitializeInterpolation(MPI_Comm& comm);
  void InitializeInterpolationForCharge(int r, MPI_Comm &comm);


  double CalculatePartitionFunctionOnTheFly(int r, double T, double deltaI);

  double CalculatePartitionFunctionOnTheFly2(int r, double T, int nmax);

  double CalculatePartitionFunctionByInterpolation(int r, double T, double deltaI);

};

#endif
