/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<NonIdealSahaEquationSolver.h>
#include<fstream>
#include<boost/math/tools/roots.hpp>
using namespace boost::math::tools;
using std::vector;
using std::map;
using std::pair;

extern double avogadro_number;

//--------------------------------------------------------------------------

NonIdealSahaEquationSolver::NonIdealSahaEquationSolver(MaterialIonizationModel& iod_ion_mat_, 
                                IoData& iod_, VarFcnBase* vf_, MPI_Comm* comm)
                          : SahaEquationSolver(iod_ion_mat_, iod_, vf_, comm),
                            eps0(iod_.ion.vacuum_permittivity)
{
  pi = 2.0*acos(0);
  factor_deltaI = e*e/(4.0*pi*eps0);
  factor_LambB  = h/sqrt(2.0*pi*me*kb);

  f.resize(elem.size(), vector<double>());
  alpha.resize(elem.size(), vector<double>());
  for(int j=0; j<(int)f.size(); j++) {
    f[j].resize(elem[j].rmax+1, 0.0);
    alpha[j].resize(elem[j].rmax+1, 0.0);
  }
}

//--------------------------------------------------------------------------

NonIdealSahaEquationSolver::~NonIdealSahaEquationSolver()
{ }

//--------------------------------------------------------------------------

double 
NonIdealSahaEquationSolver::ComputeDeltaI(int r, [[maybe_unused]] int j, double T, double nh, double zav,
                                          double one_over_lambD)
{

  if(iod_ion_mat->depression == MaterialIonizationModel::NO_DEPRESSION)
    return 0.0;

  double dI(0.0);
  if(iod_ion_mat->depression == MaterialIonizationModel::GRIEM) {
    dI = (r+1.0)*factor_deltaI*one_over_lambD;
  } 
  else if(iod_ion_mat->depression == MaterialIonizationModel::EBELING) {//Ebeling
    if(one_over_lambD==0.0)
      dI = 0.0;
    else 
      dI = (r+1.0)*factor_deltaI/(1.0/one_over_lambD + 0.125*factor_LambB/sqrt(T));
  }
  else if(iod_ion_mat->depression == MaterialIonizationModel::GRIEM_FLETCHER) {
    if(one_over_lambD==0.0)
      dI = 0.0;
    else {
      double Rr = 2.0/3.0*pow((3.0*(r+1)/(4.0*pi*nh*(1+zav))), 1.0/3.0);
      dI = (r+1.0)*factor_deltaI/sqrt(1.0/(one_over_lambD*one_over_lambD) + Rr*Rr);
    }
  }

  assert(std::isfinite(dI));
  return dI;
}

//--------------------------------------------------------------------------

double 
NonIdealSahaEquationSolver::ComputeStateForElement(int j, double T, double nh, double zav, 
                                                   double one_over_lambD, bool compute_alpha)
{
  assert(zav>=0.0 && nh>0.0);

  int rmax = elem[j].rmax;

  if(zav==0.0) {
    double zej = elem[j].molar_fraction*rmax;
    if(compute_alpha) {
      alpha[j][0] = elem[j].molar_fraction;
      for(int r=1; r<=rmax; r++)
        alpha[j][r] = 0.0;
    }
    return zej;
  }


  // Now, zav must be positive
  
  double kbT = kb*T;
  double fcore = pow( (2.0*pi*(me/h)*(kbT/h)), 1.5)/nh;

  // compute f_{r,j}, r = 1, ..., rmax
  double Ur0(0.0), Ur1(0.0), deltaI0(0.0), deltaI1(0.0);
  for(int r=0; r<=rmax; r++) {

    deltaI0 = deltaI1;
    Ur0     = Ur1;

    deltaI1 = ComputeDeltaI(r,j,T,nh,zav,one_over_lambD);
    if(r<rmax && deltaI1>iod_ion_mat->depression_max*elem[j].I[r]) {
      fprintf(stdout,"\033[0;35mWarning: Depression energy (deltaI) is truncated at %e: [j, r] = [%d, %d], "
                     "deltaI before truncation = %e, I = %e.\033[0m\n", iod_ion_mat->depression_max*elem[j].I[r],
                     j, r, deltaI1, elem[j].I[r]);
      deltaI1 = iod_ion_mat->depression_max*elem[j].I[r];
    }

    Ur1     = elem[j].CalculatePartitionFunction(r, T, deltaI1);

    f[j][r] = r==0 ? 0.0 //f[j][0] is not used anyway
                   : 2.0*Ur1/Ur0*fcore*exp(-(elem[j].I[r-1]-deltaI0)/kbT);
  }

  vector<double> zav_power(rmax+1, 0.0);
  zav_power[0] = 1.0;
  for(int r=1; r<=rmax; r++)
    zav_power[r] = zav_power[r-1]*zav;

  double denominator = 0.0; 
  double numerator = 0.0;
  double fprod = 1.0;
  double newterm(0.0);
  for(int r=1; r<=rmax; r++) {
    fprod *= f[j][r];
    newterm = fprod*zav_power[rmax-r];
    numerator += (double)r*newterm;
    denominator += newterm;
  }

  denominator += zav_power[rmax];

  //check denominator and compute zej
  double zej(0.0);
  if(denominator==0.0) //zav = 0 or extremely close to 0)
    zej = 0.0;
  else { 
    zej = numerator/denominator;
    if(!std::isfinite(zej)) {
      fprintf(stdout,"\033[0;31m*** Error: Non-Ideal Saha equation solver failed. "
                     "Z for element %d is not finite.\n\033[0m", j);
      exit(-1);
    }
    if(zej<0.0) {
      fprintf(stdout, "\033[0;35mWarning: Found negative Z (%e) for element %d. Setting it to 0.\n\033[0m",
              zej, j);
      zej = 0.0;
    } else if (zej>rmax){
      fprintf(stdout, "\033[0;35mWarning: Found Z greater than rmax (%e vs. %d) for element %d. Setting it to %d.\n\033[0m",
              zej, rmax, j, rmax);
      zej = (double)rmax;
    }
  }
  zej *= elem[j].molar_fraction;
  //fprintf(stdout,"Zej = %e for j = %d.\n", zej, j);

  //compute molar fractions if needed
  if(compute_alpha) {
    assert(std::isfinite(numerator)); 
    if(numerator==0.0) {//zej must be 0
      assert(zej==0.0);
      alpha[j][0] = elem[j].molar_fraction;
      for(int r=1; r<=rmax; r++)
        alpha[j][r] = 0.0;
    }

    alpha[j][0] = std::min(zej*zav_power[rmax]/numerator, elem[j].molar_fraction);
    double summation = alpha[j][0];
    for(int r=1; r<=rmax; r++) {
      alpha[j][r] = alpha[j][r-1]/zav*f[j][r];
      summation += alpha[j][r];
      if(summation>elem[j].molar_fraction) {
        alpha[j][r] -= (summation - elem[j].molar_fraction);
        summation = elem[j].molar_fraction;
        for(int rr=r+1; rr<=rmax; rr++)
          alpha[j][rr] = 0.0;
        break;
      }
    } 

    //make sure the sum of alphas is equal to molar_fraction
    double ratio = elem[j].molar_fraction/summation;
    if(ratio != 1.0) {
      for(auto&& alph : alpha[j]) 
        alph *= ratio;
    }
  }

  return zej;
} 

//--------------------------------------------------------------------------


void
NonIdealSahaEquationSolver::Solve(double* v, double& zav, double& nh, double& ne, 
                                  map<int, vector<double> >& alpha_rj, double* lambD)
{

  if(!iod_ion_mat) { //dummy solver 
    zav = 0.0;
    ne = 0.0;
    nh = 0.0;
    for(auto it = alpha_rj.begin(); it != alpha_rj.end(); it++) {
      vector<double> &alpha = it->second;
      for(int r=0; r<(int)alpha.size(); r++)
        alpha[r] = (it==alpha_rj.begin() && r==0) ? 1.0 : 0.0;
    }
    return;
  }


  double T = vf->GetTemperature(v[0], vf->GetInternalEnergyPerUnitMass(v[0], v[4]));
  nh = v[0]/molar_mass*avogadro_number;

  if(T<=Tmin) { //no ionization
    zav = 0.0;
    ne = 0.0;
    for(auto it = alpha_rj.begin(); it != alpha_rj.end(); it++) {
      vector<double> &alpha = it->second;
      for(int r=0; r<(int)alpha.size(); r++)
        alpha[r] = (r==0) ? elem[it->first].molar_fraction : 0.0;
    }
    return;
  }


  // ------------------------------------------
  // Solve the master equation for one_over_lambD
  // ------------------------------------------
  LambDEquation fun(*this, T, nh, &zav);

  //Find initial bracketing interval (one_over_lambD_0, one_over_lambD_1)
  double one_over_lambD_0, one_over_lambD_1, f0, f1;
  one_over_lambD_0 = 0.0;
  f0 = fun(one_over_lambD_0);

  //find the upper bound 
  one_over_lambD_1 = 4.0*fun.ComputeRHS(0.0, max_mean_atomic_number);
  double tmp = one_over_lambD_1; //store for possible use later
  bool found_initial_interval = false;
  for(int i=0; i<(int)(0.5*iod_ion_mat->maxIts); i++) {
    f1 = fun(one_over_lambD_1);
    if(f0*f1<=0.0) {
      found_initial_interval = true;
      break;
    }
    one_over_lambD_1 /= 2.0; 
  }
  if(!found_initial_interval) { //go the other way
    one_over_lambD_0 = 0.5*tmp;
    f0 = fun(one_over_lambD_0);
    one_over_lambD_1 = tmp; 
    for(int i=0; i<(int)(0.5*iod_ion_mat->maxIts); i++) {
      f1 = fun(one_over_lambD_1);
      if(f0*f1<=0.0) {
        found_initial_interval = true;
        break;
      }
      one_over_lambD_1 *= 2.0; 
    }
  }
  if(!found_initial_interval) {
    fprintf(stdout,"\033[0;31m*** Error: Non-ideal Saha equation solver failed. "
            "Cannot find an initial bracketing interval. (T = %e, nh = %e)\n\033[0m",
            T, nh);
    exit(-1);
  }


  // Calling boost function for root-finding
  // Warning: "maxit" is BOTH AN INPUT AND AN OUTPUT  
  double one_over_lambD(0.0); //solution
  boost::uintmax_t maxit = iod_ion_mat->maxIts;
  double tol = iod_ion_mat->convergence_tol;
  if(f0==0.0) {
    one_over_lambD = one_over_lambD_0; maxit = 0;
  } else if(f1==0.0) {
    one_over_lambD = one_over_lambD_1; maxit = 0;
  } else {
    for(int trial = 0; trial < 2; trial++) {
      pair<double,double> sol; 
      sol = toms748_solve(fun, one_over_lambD_0, one_over_lambD_1, f0, f1,
                          [=](double r0, double r1)
                          {return r1-r0<std::min(tol,0.001*(one_over_lambD_1-one_over_lambD_0));},
                          maxit);
      one_over_lambD = 0.5*(sol.first + sol.second);
      if(one_over_lambD>=0) break;

      // fail-safe
      if(trial>0) break;

      if(!isfinite(sol.first) || sol.first<0)
        sol.first = one_over_lambD_0;
      if(!isfinite(sol.second) || sol.second<sol.first)
        sol.second = one_over_lambD_1;
      one_over_lambD_0 = sol.first;
      one_over_lambD_1 = sol.second;
      f0 = fun(one_over_lambD_0);
      f1 = fun(one_over_lambD_1); 
      if(f0*f1>0)
        break;
    }
  }

  //*******************************************************************

#if DEBUG_SAHA_SOLVER == 1
  fprintf(stdout,"-- Non-Ideal Saha equation solver terminated in %d iterations, Zav = %.12e,"
          " 1/lambD = %.12e.\n",
          (int)maxit, zav, one_over_lambD);
#endif


  if(one_over_lambD<=0.0) {
    one_over_lambD = 0.0;
    zav = 0.0;
    ne = 0.0;
    for(auto it = alpha_rj.begin(); it != alpha_rj.end(); it++) {
      vector<double> &alpha = it->second;
      for(int r=0; r<(int)alpha.size(); r++)
        alpha[r] = (r==0) ? elem[it->first].molar_fraction : 0.0;
    }
    return;
  }


  // post-processing.
  fun(one_over_lambD); //updates zav and alpha's
  ne = zav*nh;
  for(auto it = alpha_rj.begin(); it != alpha_rj.end(); it++) {
    int j = it->first; //element id
    vector<double> &my_alpha = it->second; //alpha_r

    if(j>=(int)elem.size()) {//this material does not have element j
      for(int r=0; r<(int)my_alpha.size(); r++)
        my_alpha[r] = 0.0;
      continue;
    }

    int max_size = std::min((int)my_alpha.size()-1, elem[j].rmax);
    for(int r=0; r<max_size; r++) {
      my_alpha[r] = alpha[j][r];
    }

    my_alpha[max_size] = elem[j].molar_fraction;
    for(int r=0; r<max_size; r++)
      my_alpha[max_size] -= my_alpha[r];

    //allow some roundoff error
    //assert(alpha[last_one]>=-1.0e-4);
    if(my_alpha[max_size]<0)
      my_alpha[max_size] = 0;

    //too many slots? put 0
    for(int r=max_size+1; r<(int)my_alpha.size(); r++)
      my_alpha[r] = 0.0;
  }

  if(lambD) {
    if(one_over_lambD==0 || !std::isfinite(one_over_lambD))
      *lambD = 0.0;
    else
      *lambD = 1.0/one_over_lambD;
  }
 
}

//--------------------------------------------------------------------------

NonIdealSahaEquationSolver::
LambDEquation::LambDEquation(NonIdealSahaEquationSolver &saha_, double T_, double nh_,
                             double *zav_ptr_)
             : saha(saha_), T(T_), nh(nh_), zav_ptr(zav_ptr_)
{
  assert(nh>0.0 && T>0.0); 

  factor_lambD = saha.e*saha.e/(saha.eps0*saha.kb*T);
}

//--------------------------------------------------------------------------

double NonIdealSahaEquationSolver::
LambDEquation::operator()(double one_over_lambD)
{

  // ------------------------------
  // Step 1: Solve for Zav
  // ------------------------------
  ZavEquation fun(saha, T, nh, one_over_lambD);

  double zav(0.0);

  //Find initial bracketing interval (zav0, zav1)
  double zav0, zav1, f0, f1;
  zav0 = 0.0;
  zav1 = saha.max_mean_atomic_number; //zav1>zav0
  f0 = fun(zav0);
  bool found_initial_interval = false;
  for(int i=0; i<saha.iod_ion_mat->maxIts; i++) {
    f1 = fun(zav1);
    if(f0*f1<=0.0) {
      found_initial_interval = true;
      break;
    }
    zav1 /= 2.0;
  }
  if(!found_initial_interval) {
    fprintf(stdout,"\033[0;31m*** Error: Non-ideal Saha equation solver failed (1). "
            "Cannot find an initial bracketing interval. (T = %e, nh = %e, one_over_lambD = %e)\n\033[0m",
            T, nh, one_over_lambD);
    exit(-1);
  }

  // Calling boost function for root-finding
  // Warning: "maxit" is BOTH AN INPUT AND AN OUTPUT
  boost::uintmax_t maxit = saha.iod_ion_mat->maxIts;
  double tol = saha.iod_ion_mat->convergence_tol;
  if(f0==0.0) {
    zav = zav0; maxit = 0;
  } else if(f1==0.0) {
    zav = zav1; maxit = 0;
  } else {
    for(int trial = 0; trial < 2; trial++) {
      pair<double,double> sol;
      sol = toms748_solve(fun, zav0, zav1, f0, f1,
                          [=](double r0, double r1){return r1-r0<std::min(tol,0.001*(zav1-zav0));},
                          maxit);
      zav = 0.5*(sol.first + sol.second);
      if(zav>=0) break;

      // fail-safe
      if(trial>0) break;

      if(!std::isfinite(sol.first) || sol.first<0)
        sol.first = 0.0;
      if(!std::isfinite(sol.second) || sol.second<sol.first)
        sol.second = std::max(zav1, sol.first);
      zav0 = sol.first;
      zav1 = sol.second;
      f0 = fun(zav0);
      f1 = fun(zav1);
      if(f0*f1>0) //not a bracketing interval...
        break;
    }
  }

  if(zav<0.0)
    zav = 0.0;

  // Store zav
  if(zav_ptr)
    *zav_ptr = zav;

  // ------------------------------
  // Step 2: Compute RHS (and updates alphas)
  // ------------------------------

  double rhs = ComputeRHS(one_over_lambD, zav);

  return one_over_lambD - rhs;

}

//--------------------------------------------------------------------------

double NonIdealSahaEquationSolver::
LambDEquation::ComputeRHS(double one_over_lambD, double zav)
{
  double summation = zav; 
  for(int j=0; j<(int)saha.elem.size(); j++) {
    saha.ComputeStateForElement(j, T, nh, zav, one_over_lambD, true);  //we need alphas 
    for(int r=1; r<=saha.elem[j].rmax; r++)
      summation += r*r*saha.alpha[j][r];
  }

  return sqrt(factor_lambD*nh*summation);
}

//--------------------------------------------------------------------------

NonIdealSahaEquationSolver::LambDEquation::
ZavEquation::ZavEquation(NonIdealSahaEquationSolver &saha_, double T_, double nh_,
                         double one_over_lambD_)
           : saha(saha_), T(T_), nh(nh_), one_over_lambD(one_over_lambD_)
{
  assert(T>0.0 && nh>0.0);
}

//--------------------------------------------------------------------------

double NonIdealSahaEquationSolver::LambDEquation::
ZavEquation::ComputeRHS(double zav)
{
  double rhs = 0.0;
  for(int j=0; j<(int)saha.elem.size(); j++)
    rhs += saha.ComputeStateForElement(j, T, nh, zav, one_over_lambD, false); //we don't need alpha here

  return rhs;
}

//--------------------------------------------------------------------------






//--------------------------------------------------------------------------

