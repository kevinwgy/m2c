#include<SahaEquationSolver.h>
#include<fstream>
#include<boost/math/tools/roots.hpp>
using namespace boost::math::tools;
using std::vector;
using std::map;
using std::pair;


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
                                       VarFcnBase* vf_)
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

  // Read element data
  int numElems = iod_ion_mat->elementMap.dataMap.size();
  elem.resize(numElems, ElementData());
  for(auto it = iod_ion_mat->elementMap.dataMap.begin(); 
      it != iod_ion_mat->elementMap.dataMap.end(); it++) {
    int element_id = it->first;
    if(element_id<0 || element_id>=numElems) {
      print_error("*** Error: Chemical element id should be between 0 and %d. Found %d.\n",
                  numElems-1, element_id);
      exit_mpi();
    }


    ElementData& data(elem[it->first]);

    // Step 1: Read molar fraction and atomic number
    data.molar_fraction = it->second->molar_fraction;
    data.atomic_number = it->second->atomic_number;
    if(data.molar_fraction<=0 || data.atomic_number<=0) {
      print_error("*** Error: Detected molar fraction %e, atomic number %d. (Should be greater than 0.)\n",
                  data.molar_fraction, data.atomic_number);
      exit_mpi();
    }
    print("  o Found chemical element with atomic number %d and molar fraction %e.\n", data.atomic_number,
          data.molar_fraction);


    // Step 2: Read the ionization energies
    if(!strcmp(it->second->ionization_energy_filename, "")) {//filename undefined
      print_error("*** Error: Missing the ionization energy file for the element with atomic number %d.\n",
                  data.atomic_number);
      exit_mpi();
    }
    std::fstream file;
    file.open(it->second->ionization_energy_filename, std::fstream::in);
    if(!file.is_open()) {
      print_error("*** Error: Cannot open ionization energy file %s.\n", it->second->ionization_energy_filename);
      exit_mpi();
    }
    GetDataInFile(file, data.I, 1000, true);
    file.close();
    print("    * Read ionization energy file %s: %d energy values.\n",
          it->second->ionization_energy_filename, data.I.size());
    

    // Step 3: Read excitation energy files
    data.E.resize(data.atomic_number);
    int max_size = 0;
    for(int i=0; i<data.atomic_number; i++) {
      std::string filename = std::string(it->second->excitation_energy_files_prefix) + "_" 
                           + std::to_string(data.atomic_number) + "_" + std::to_string(i)
                           + std::string(it->second->excitation_energy_files_suffix);
      file.open(filename.c_str(), std::fstream::in);
      if(!file.is_open()) {
        print_error("*** Error: Cannot open excitation energy file %s.\n", filename.c_str());
        exit_mpi();
      }
      GetDataInFile(file, data.E[i], 10000, true);
      file.close();
      if(max_size<data.E[i].size())
        max_size = data.E[i].size();
    }
    print("    * Read %d excitation energy files %s_%d_X%s. Max excited state: %d.\n", data.atomic_number,
          it->second->excitation_energy_files_prefix, data.atomic_number,
          it->second->excitation_energy_files_suffix, max_size);


    // Step 4: Read angular momentum energy files
    data.g.resize(data.atomic_number);
    max_size = 0;
    for(int i=0; i<data.atomic_number; i++) {
      std::string filename = std::string(it->second->angular_momentum_files_prefix) + "_" 
                           + std::to_string(data.atomic_number) + "_" + std::to_string(i)
                           + std::string(it->second->angular_momentum_files_suffix);
      file.open(filename.c_str(), std::fstream::in);
      if(!file.is_open()) {
        print_error("*** Error: Cannot open angular momentum file %s.\n", filename.c_str());
        exit_mpi();
      }
      GetDataInFile(file, data.g[i], 10000, true);
      file.close();
      if(max_size<data.g[i].size())
        max_size = data.g[i].size();
    }
    print("    * Read %d angular momentum files %s_%d_X%s. Max excited state: %d.\n", data.atomic_number,
          it->second->angular_momentum_files_prefix, data.atomic_number,
          it->second->angular_momentum_files_suffix, max_size);

  }
      
  // find rmax for each element/species
  for(auto it = elem.begin(); it != elem.end(); it++) {
    it->rmax = std::min(it->I.size(), std::min(it->g.size(), it->E.size()));
    it->rmax--; //rmax should be "size" - 1
  }
  // find max atomic number among all the species/elements
  max_atomic_number = 0;
  for(auto it = elem.begin(); it != elem.end(); it++)
    max_atomic_number = std::max(max_atomic_number, it->atomic_number);

}

//--------------------------------------------------------------------------

SahaEquationSolver::~SahaEquationSolver()
{ }

//--------------------------------------------------------------------------

void
SahaEquationSolver::GetDataInFile(std::fstream& file, vector<double> &X, int MaxCount, bool non_negative)
{
  double tmp;
  for(int i=0; i<MaxCount; i++) {
    file >> tmp;
    if(non_negative && tmp<=0) {
      print_error("*** Error: Detected negative number (%e) in an ionization data file.\n", tmp);
      exit_mpi();
    }
    X.push_back(tmp);
    if(file.eof())
      break;
  } 
} 

//--------------------------------------------------------------------------

void
SahaEquationSolver::Solve(double* v, double& zav, double& nh, double& ne, 
                          map<int, vector<double> >& alpha_rj)
{

  double T = vf->GetTemperature(v[0], vf->GetInternalEnergyPerUnitMass(v[0], v[4]));
  nh = v[4]/(kb*T);

  if(!iod_ion_mat) { //dummy solver
    zav = 0.0;
    ne = 0.0;
    for(auto it = alpha_rj.begin(); it != alpha_rj.end(); it++) {
      vector<double> &alpha = it->second;
      for(int r=0; r<alpha.size(); r++)
        alpha[r] = 0.0;
    }
    return;
  }

 
  ZavEquation fun(kb, T, v[4], me, h, elem);

  //Find initial bracketing interval (zav0, zav1)
  double zav0, zav1, f0, f1; 
  zav0 = 0.0;
  zav1 = max_atomic_number; //zav1>zav0
  f0 = fun(zav0);
  f1 = fun(zav1);

  //*******************************************************************
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
  //*******************************************************************

  fprintf(stderr,"-- Saha equation solver converged in %d iterations, Zav = %e.\n", (int)maxit, zav);

  //post-processing.
  ne = zav*nh;

  for(auto it = alpha_rj.begin(); it != alpha_rj.end(); it++) {
    int j = it->first; //element id
    vector<double> &alpha = it->second; //alpha_r

    for(int r=0; r<alpha.size(); r++) {
      alpha[r] = 0.0; //TODO 
    }

  }

}

//--------------------------------------------------------------------------

SahaEquationSolver::ZavEquation::ZavEquation(double kb, double T, double p, double me, double h, 
                                             vector<ElementData>& elem_)
                               : elem(elem_)
{
  kbT = kb*T;
  nh = p/kbT;
  double pi = 2.0*acos(0);
  fcore = pow( (2.0*pi*me*kbT)/(h*h), 1.5);
}

//--------------------------------------------------------------------------

double
SahaEquationSolver::ZavEquation::ComputeRHS(double zav)
{
  double rhs = 0.0;

  for(int j=0; j<elem.size(); j++) {

    // calculate Ur, r = 0, 1, 2, ..., rmax
    vector<double> U(elem[j].rmax+1, 0.0);
    for(int r=0; r<=elem[j].rmax; r++) {
      

    }

  }

  return rhs;
}

//--------------------------------------------------------------------------






