#include<SahaEquationSolver.h>
#include<fstream>
#include<boost/math/tools/roots.hpp>
using namespace boost::math::tools;
using std::vector;
using std::map;
using std::pair;

extern double avogadro_number;

//--------------------------------------------------------------------------

SahaEquationSolver::SahaEquationSolver(IoData& iod_, VarFcnBase* vf_)
                  : iod_ion_mat(NULL), vf(vf_),
                    h(iod_.ion.planck_constant),
                    e(iod_.ion.electron_charge),
                    me(iod_.ion.electron_mass),
                    kb(iod_.ion.boltzmann_constant)
{ }

//--------------------------------------------------------------------------

SahaEquationSolver::SahaEquationSolver(MaterialIonizationModel& iod_ion_mat_, IoData& iod_, 
                                       VarFcnBase* vf_, MPI_Comm* comm)
                  : iod_ion_mat(&iod_ion_mat_), vf(vf_),
                    h(iod_.ion.planck_constant),
                    e(iod_.ion.electron_charge),
                    me(iod_.ion.electron_mass),
                    kb(iod_.ion.boltzmann_constant)
{
  if(iod_ion_mat->type != MaterialIonizationModel::SAHA_IDEAL) {
    print_error("*** Error: Currently, only the ideal Saha equation is implemented in the solver.\n");
    exit_mpi();
  }

  Tmin = iod_ion_mat->ionization_Tmin;
  if(Tmin<=0) {
    print_error("*** Error: Detected Tmin = %e in material ionization model. Must be positive.\n", Tmin);
    exit_mpi();
  }

  // Read element data
  int numElems = iod_ion_mat->elementMap.dataMap.size();
  elem.resize(numElems, AtomicIonizationData());
  double total_molar = 0.0;
  for(auto it = iod_ion_mat->elementMap.dataMap.begin(); 
      it != iod_ion_mat->elementMap.dataMap.end(); it++) {
    int element_id = it->first;
    if(element_id<0 || element_id>=numElems) {
      print_error("*** Error: Chemical element id should be between 0 and %d. Found %d.\n",
                  numElems-1, element_id);
      exit_mpi();
    }

    total_molar += it->second->molar_fraction;

    AtomicIonizationData& data(elem[it->first]);
    if(iod_ion_mat->partition_evaluation == MaterialIonizationModel::CUBIC_SPLINE_INTERPOLATION)
      data.Setup(it->second, h, e, me, kb, 1, iod_ion_mat->sample_Tmin, iod_ion_mat->sample_Tmax, 
                 iod_ion_mat->sample_size, comm);
    else if(iod_ion_mat->partition_evaluation == MaterialIonizationModel::LINEAR_INTERPOLATION)
      data.Setup(it->second, h, e, me, kb, 2, iod_ion_mat->sample_Tmin, iod_ion_mat->sample_Tmax, 
                 iod_ion_mat->sample_size, comm);
    else if(iod_ion_mat->partition_evaluation == MaterialIonizationModel::ON_THE_FLY)
      data.Setup(it->second, h, e, me, kb, 0, 0, 0, 0, NULL);
    else {
      print_error("*** Error: Detected an unknown method for partition function evaluation.\n");
      exit_mpi();
    }
  }

  if(total_molar >= 1.0 + 1e-12) {
    print_error("*** Error: Sum of molar fractions (%e) exceeds 1.\n", total_molar);
    exit_mpi();
  }
  if(total_molar <= 1.0 - 1e-12) //throw out a warning, but continue to run
    print("Warning: Sum of molar fractions (%e) is less than 1.\n", total_molar);


  // find molar mass and max atomic number among all the species/elements
  molar_mass = 0.0;
  max_atomic_number = 0;
  for(auto it = elem.begin(); it != elem.end(); it++) {
    molar_mass += (it->molar_fraction)*(it->molar_mass);
    max_atomic_number = std::max(max_atomic_number, it->atomic_number);
  }

}

//--------------------------------------------------------------------------

SahaEquationSolver::~SahaEquationSolver()
{ }

//--------------------------------------------------------------------------

void
SahaEquationSolver::Solve(double* v, double& zav, double& nh, double& ne, 
                          map<int, vector<double> >& alpha_rj)
{
  //nh = T>0 ? v[4]/(kb*T) : 0; //for dummy solver, there may not be a temperature law --> T = 0

  if(!iod_ion_mat) { //dummy solver 
    zav = 0.0;
    ne = 0.0;
    nh = 0.0;
    for(auto it = alpha_rj.begin(); it != alpha_rj.end(); it++) {
      vector<double> &alpha = it->second;
      for(int r=0; r<alpha.size(); r++)
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
      for(int r=0; r<alpha.size(); r++)
        alpha[r] = (r==0) ? elem[it->first].molar_fraction : 0.0;
    }
    return;
  }

  // ------------------------------
  // Step 1: Solve for Zav 
  // ------------------------------
  ZavEquation fun(kb, T, nh, me, h, elem);

  //Find initial bracketing interval (zav0, zav1)
  double zav0, zav1, f0, f1; 
  zav0 = 0.0;
  zav1 = max_atomic_number; //zav1>zav0
  f0 = fun(zav0);
  bool found_initial_interval = false;
  for(int i=0; i<iod_ion_mat->maxIts; i++) {
    f1 = fun(zav1);
    if(f0*f1<=0.0) {
      found_initial_interval = true;
      break;
    }
    zav1 /= 2.0;
  }
  if(!found_initial_interval) {
    fprintf(stderr,"\033[0;31m*** Error: Saha equation solver failed. "
            "Cannot find an initial bracketing interval. (p = %e, T = %e)\n\033[0m", v[4], T);
    exit(-1);
  }

  // Calling boost function for root-finding
  // Warning: "maxit" is BOTH AN INPUT AND AN OUTPUT  
  boost::uintmax_t maxit = iod_ion_mat->maxIts;
  double tol = iod_ion_mat->convergence_tol;
  if(f0==0.0) {
    zav = zav0; maxit = 0;
  } else if(f1==0.0) {
    zav = zav1; maxit = 0;
  } else {
    pair<double,double> sol; 
    sol = toms748_solve(fun, zav0, zav1, f0, f1,
                        [=](double r0, double r1){return r1-r0<std::min(tol,0.001*(zav1-zav0));},
                        maxit);
    zav = 0.5*(sol.first + sol.second);
  }
  assert(zav>0.0);
  //*******************************************************************

#if DEBUG_SAHA_SOLVER == 1
  fprintf(stderr,"-- Saha equation solver converged in %d iterations, Zav = %e.\n", (int)maxit, zav);
#endif

  //post-processing.
  ne = zav*nh;

  for(auto it = alpha_rj.begin(); it != alpha_rj.end(); it++) {
    int j = it->first; //element id
    vector<double> &alpha = it->second; //alpha_r

    double zej = fun.GetZej(zav, j);
    double denom = 0.0;
    double zav_power = 1.0;
    for(int i=1; i<=elem[j].rmax; i++) {
      zav_power *= zav;
      denom += (double)i/zav_power*fun.GetFProd(i,j);
    }
    assert(denom>0.0);
    alpha[0] = zej/denom;
 
    for(int r=1; r<alpha.size()-1; r++)
      alpha[r] = (r<=elem[j].rmax) ? alpha[r-1]/zav*fun.GetFProd(r,j)/fun.GetFProd(r-1,j) : 0.0;

    int last_one = alpha.size()-1; 
    alpha[last_one] = elem[j].molar_fraction;
    for(int r=0; r<last_one; r++)
      alpha[last_one] -= alpha[r];

    assert(alpha[last_one]>=-1.0e-10); //allow some roundoff error
    if(alpha[last_one]<0)
      alpha[last_one] = 0;
  }

}

//--------------------------------------------------------------------------

SahaEquationSolver::ZavEquation::ZavEquation(double kb, double T, double nh, double me, double h, 
                                             vector<AtomicIonizationData>& elem_)
                               : elem(elem_)
{

  double kbT = kb*T;

  double pi = 2.0*acos(0);
  double fcore = pow( (2.0*pi*(me/h)*(kbT/h)), 1.5)/nh;

  // compute fprod
  fprod.resize(elem.size(), vector<double>());
  double Ur0, Ur1, f1;
  for(int j=0; j<fprod.size(); j++) {

    fprod[j].resize(elem[j].rmax+1, 0);

    fprod[j][0] = 1.0; //must be set to 1.0, s.t. fprod[j][1]/fprod[j][0] = f[j][1]

    Ur0 = elem[j].CalculatePartitionFunction(0, T);

    for(int r=0; r<elem[j].rmax; r++) {
      Ur1 = elem[j].CalculatePartitionFunction(r+1, T);
      f1 = 2.0*Ur1/Ur0*fcore*exp(-elem[j].I[r]/kbT);
      //if(j==1) fprintf(stderr,"f(%d,%d) = %e, Ur1 = %e, Ur0 = %e, exp = %e.\n", r+1, j, f1, Ur1, Ur0, exp(-elem[j].I[r]/kbT));
      fprod[j][r+1] = fprod[j][r]*f1;

      if(fprod[j][r+1]<0) {
        print_error("*** Error: Negative partition function (%e) in Saha equation solver. "
                    "Input error or small sample size.\n", fprod[j][r+1]); 
        exit_mpi();
      }

      Ur0 = Ur1;
    }
  } 
}

//--------------------------------------------------------------------------

double
SahaEquationSolver::ZavEquation::ComputeRHS(double zav)
{
  double rhs = 0.0;

  for(int j=0; j<elem.size(); j++)
    rhs += ComputeRHS_ElementJ(zav, j);

  return rhs;
}

//--------------------------------------------------------------------------

double
SahaEquationSolver::ZavEquation::ComputeRHS_ElementJ(double zav, int j)
{

  if(zav == 0.0) //special case
    return elem[j].molar_fraction*elem[j].rmax;

  // zav != 0

  double denominator = 0.0;
  double numerator = 0.0;
  double zav_power = 1.0;

  double ith_term;

  for(int i=elem[j].rmax; i>=1; i--) {
    ith_term     = fprod[j][i]*zav_power;
    numerator   += (double)i*ith_term;
    denominator += ith_term;

    zav_power *= zav;
  } 

  denominator += zav_power;


  return elem[j].molar_fraction*numerator/denominator;

}

//--------------------------------------------------------------------------







//--------------------------------------------------------------------------

