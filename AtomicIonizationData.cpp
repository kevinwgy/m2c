/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<AtomicIonizationData.h>
#include<fstream>
using std::vector;

//--------------------------------------------------------------------------

AtomicIonizationData::AtomicIonizationData() : interpolation(0), molar_fraction(-1.0), atomic_number(-1)
{ }

//--------------------------------------------------------------------------

AtomicIonizationData::~AtomicIonizationData()
{ 
  for(auto it = spline.begin(); it != spline.end(); it++)
    for(auto it2 = it->begin(); it2 != it->end(); it2++)
      if(*it2)
        delete *it2;
}

//--------------------------------------------------------------------------

void
AtomicIonizationData::Setup(AtomicIonizationModel* iod_aim, double h_, double e_, double me_, double kb_,
                            bool ideal_, int interp, double sample_Tmin_, double sample_Tmax_, int sample_size_, 
                            MPI_Comm* comm)
{

  ideal = ideal_;

  // Step 0: Get constants
  h = h_;
  e = e_;
  me = me_;
  kb = kb_;

  // Step 1: Read molar fraction, atomic number, and max charge
  molar_fraction = iod_aim->molar_fraction;
  molar_mass     = iod_aim->molar_mass;
  atomic_number = iod_aim->atomic_number;
  rmax = iod_aim->max_charge;
  if(molar_fraction<=0 || molar_mass<=0 || atomic_number<=0) {
    print_error("*** Error: Detected molar fraction %e, molar mass %e, atomic number %d. (Should be greater than 0.)\n",
                molar_fraction, molar_mass, atomic_number);
    exit_mpi();
  }
  print("  o Found chemical element with atomic number %d, molar mass %e, and molar fraction %e.\n", atomic_number,
        molar_mass, molar_fraction);


  // Step 2: Read the ionization energies
  if(!strcmp(iod_aim->ionization_energy_filename, "")) {//filename undefined
    print_error("*** Error: Missing the ionization energy file for the element with atomic number %d.\n",
                atomic_number);
    exit_mpi();
  }
  std::fstream file;
  file.open(iod_aim->ionization_energy_filename, std::fstream::in);
  if(!file.is_open()) {
    print_error("*** Error: Cannot open ionization energy file %s.\n", iod_aim->ionization_energy_filename);
    exit_mpi();
  }
  GetDataInFile(file, I, 1000, true);
  file.close();
  print("    * Read ionization energy file %s: %d energy values.\n",
        iod_aim->ionization_energy_filename, I.size());
  

  // Step 3: Read excitation energy files
  E.resize(atomic_number);
  int max_size = 0, i;
  for(i=0; i<atomic_number; i++) {
    std::string filename = std::string(iod_aim->excitation_energy_files_prefix) 
                         + std::to_string(i)
                         + std::string(iod_aim->excitation_energy_files_suffix);
    file.open(filename.c_str(), std::fstream::in);
    if(!file.is_open()) {
      if(i<rmax) 
        print_warning("Warning: Cannot open excitation energy file %s. Skipping higher "
                      "states of ionization.\n", filename.c_str());
      //if i>=rmax, no need to throw a warning.
      E.resize(i);
      break;
    }
    GetDataInFile(file, E[i], 10000, true);
    file.close();
    if(max_size<(int)E[i].size())
      max_size = E[i].size();
  }
  print("    * Read %d excitation energy files: %sX%s. Max excited state: %d.\n", i,
        iod_aim->excitation_energy_files_prefix, 
        iod_aim->excitation_energy_files_suffix, max_size);


  // Step 4: Read degeneracy files ("g") 
  g.resize(atomic_number);
  max_size = 0;
  for(i=0; i<atomic_number; i++) {
    std::string filename = std::string(iod_aim->degeneracy_files_prefix)
                         + std::to_string(i)
                         + std::string(iod_aim->degeneracy_files_suffix);
    file.open(filename.c_str(), std::fstream::in);
    if(!file.is_open()) {
      if(i<rmax) 
        print_warning("Warning: Cannot open degeneracy file %s. Skipping higher "
                    "states of ionization.\n", filename.c_str());
      //if i>=rmax, no need to throw a warning.
      g.resize(i);
      break;
    }
    GetDataInFile(file, g[i], 10000, true);
    file.close();
    if(max_size<(int)g[i].size())
      max_size = g[i].size();
  }
  print("    * Read %d degeneracy files: %sX%s. Max excited state: %d.\n", i,
        iod_aim->degeneracy_files_prefix, 
        iod_aim->degeneracy_files_suffix, max_size);

      
  // check to see if rmax should be reduced based on data files
  int data_size = std::min(I.size(), std::min(g.size(), E.size()));
  if(rmax > data_size) //rmax should be smaller than or equal to "size"
    rmax = data_size;



  // find out the max number of terms involved in the summation (up to E[r][n] <= I[r])
  max_terms.resize(rmax, 0); 
  for(int r=0; r<rmax; r++)
    max_terms[r] = std::min(int(std::upper_bound(E[r].begin(), E[r].end(), I[r]) - E[r].begin()),
                            (int)g[r].size());
 

  // if interpolation is true, create the interpolation database
  interpolation = interp;
  if(interpolation && sample_Tmax_<=sample_Tmin_){
    print_error("*** Error: Incorrect Tmin and/or Tmax for samples in atomic ionization model (atomic no. %d). "
                "Tmin must be less than Tmax.\n", atomic_number);
    exit_mpi();
  } 
  sample_Tmin = sample_Tmin_;
  sample_Tmax = sample_Tmax_;
  sample_size = sample_size_;
  assert(isfinite(sample_size));

  if(interpolation) {
    assert(comm);
    InitializeInterpolation(*comm);
  } 

}

//--------------------------------------------------------------------------

template<class T>
void
AtomicIonizationData::GetDataInFile(std::fstream& file, vector<T> &X, int MaxCount, bool non_negative)
{
  double tmp;
  for(int i=0; i<MaxCount; i++) {
    file >> tmp;
    if(file.eof())
      break;
    if(non_negative && tmp<0) {
      print_error("*** Error: Detected negative number (%e) in an ionization data file.\n", tmp);
      exit_mpi();
    }
    X.push_back(tmp);
  }
}

//--------------------------------------------------------------------------

void
AtomicIonizationData::InitializeInterpolation(MPI_Comm &comm)
{

  for(int r=0; r<rmax; r++) {
    Us.push_back(vector<vector<double> >(ideal ? 1 : max_terms[r], vector<double>(sample_size)));
  }
  UsCoeffs.assign(rmax, std::tuple<double,double,double>(0,0,0));

  spline.clear();
  if(interpolation == 1) {//cubic spline interpolation
    for(int r=0; r<rmax; r++) 
      spline.push_back(vector<boost::math::cubic_b_spline<double>*>(ideal ? 1 : max_terms[r], NULL));
  }

  for(int r=0; r<rmax; r++)
    InitializeInterpolationForCharge(r, comm);

}

//--------------------------------------------------------------------------

void
AtomicIonizationData::InitializeInterpolationForCharge(int r, MPI_Comm &comm)
{

  vector<vector<double> >& U(Us[r]); 
  std::tuple<double,double,double>& coeffs(UsCoeffs[r]);

  double factor(0);
  for(int i=0; i<(int)E[r].size(); i++){
    factor = -E[r][i]/kb;
    if(factor!=0)
      break;
  }

  double expmin = exp(factor/sample_Tmin);
  double expmax = exp(factor/sample_Tmax);
  assert(expmin<expmax);
  
  double delta_exp = (expmax-expmin)/(sample_size-1.0);

  //store coefficients
  std::get<0>(coeffs) = factor;
  std::get<1>(coeffs) = expmin;
  std::get<2>(coeffs) = delta_exp;

  //each processor calculates some sample points
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size); 

  int block_size = floor((double)sample_size/(double)mpi_size);
  int remainder = sample_size - block_size*mpi_size;
  assert(remainder>=0 && remainder<mpi_size);

  //prepare for comm.
  int* counts = new int[mpi_size];
  int* displacements = new int[mpi_size];
  for(int i=0; i<mpi_size; i++) {
    counts[i] = (i<remainder) ? block_size + 1 : block_size;
    displacements[i] = (i<remainder) ? (block_size+1)*i : block_size*i + remainder;
  }
  assert(displacements[mpi_size-1]+counts[mpi_size-1] == sample_size);

  int my_start_id = displacements[mpi_rank];
  int my_block_size = counts[mpi_rank];

  vector<double> T(my_block_size,0.0);
  for(int i=0; i<my_block_size; i++)
    T[i] = (my_start_id+i==0) ? sample_Tmin : factor/log(expmin+(my_start_id+i)*delta_exp); 
            //expmin can be 0 if Tmin is small (but nonzero)
                
  for(int k=0; k<(int)U.size(); k++) {

    for(int i=0; i<my_block_size; i++)
      U[k][my_start_id+i] = CalculatePartitionFunctionOnTheFly2(r, T[i], max_terms[r] - k); 

    //communication
    MPI_Allgatherv(MPI_IN_PLACE, my_block_size, MPI_DOUBLE, U[k].data(), counts, displacements, MPI_DOUBLE, comm);
 
    //initialize spline interpolation
    if(interpolation == 1)
      spline[r][k] = new boost::math::cubic_b_spline<double>(U[k].begin(), U[k].end(), expmin, delta_exp);

  }

  delete [] counts;
  delete [] displacements;

}

//--------------------------------------------------------------------------

double
AtomicIonizationData::CalculatePartitionFunction(int r, double T, double deltaI)
{
  assert(r<=rmax);

  //if(r==atomic_number)
  if(r==rmax) //if rmax<atomic_number, we need to avoid accessing non-existent E & g...
    return 1.0; //per Shafquat

  if(interpolation && T>=sample_Tmin && T<=sample_Tmax)
    return CalculatePartitionFunctionByInterpolation(r, T, deltaI);

  return CalculatePartitionFunctionOnTheFly(r, T, deltaI);
}

//--------------------------------------------------------------------------

double
AtomicIonizationData::CalculatePartitionFunctionOnTheFly(int r, double T, double deltaI)
{

  double Ur = 0.0;

  //upper_bound returns an iterator to the first element that is greater than I[r]-deltaI
  int nsize = std::min(int(std::upper_bound(E[r].begin(), E[r].end(), I[r]-deltaI) - E[r].begin()),
                       (int)g[r].size());

/*
  //TEST: TRUNCATION IS IGNORED!
  int nsize = std::min(g[r].size(), E[r].size());
  //END OF TEST
  
*/


  //fprintf(stdout,"I[%d] = %e, nsize = %d, E[%d][%d] = %e, Ieff-E = %e | %e.\n", r, I[r], nsize, r, 
  //        nsize-1, E[r][nsize-1],
  //        I[r]-deltaI-E[r][nsize-1], I[r]-deltaI-E[r][std::min(nsize, (int)E[r].size()-1)]);

  for(int n=0; n<nsize; n++)
    Ur += g[r][n]*exp(-E[r][n]/(kb*T));

  return Ur;

}

//--------------------------------------------------------------------------

double
AtomicIonizationData::CalculatePartitionFunctionOnTheFly2(int r, double T, int nterms)
{

  assert(nterms>=1);

  double Ur = 0.0;
  int nsize = std::min(nterms, std::min((int)g[r].size(), (int)E[r].size()));
  for(int n=0; n<nsize; n++)
    Ur += g[r][n]*exp(-E[r][n]/(kb*T));

  return Ur;

}

//--------------------------------------------------------------------------

double
AtomicIonizationData::CalculatePartitionFunctionByInterpolation(int r, double T, double deltaI)
{

  double factor = std::get<0>(UsCoeffs[r]);
  double expT = exp(factor/T);

  int k;
  if(ideal) {
    assert(deltaI==0.0);
    k = 0;
  }
  else {
    int nsize = std::min(int(std::upper_bound(E[r].begin(), E[r].end(), I[r]-deltaI) - E[r].begin()),
                         (int)g[r].size());
    k = max_terms[r] - nsize;
    assert(k>=0 && k<(int)spline[r].size());
  }

  if(interpolation==1) 

    return (*spline[r][k])(expT); 

  else if(interpolation==2) { //linear

    double expmin = std::get<1>(UsCoeffs[r]);
    double delta_exp = std::get<2>(UsCoeffs[r]);

    int i = floor((expT - expmin)/delta_exp);
    assert(i>=0 && i<=sample_size-1);

    if(i==sample_size-1)
      return Us[r][k][sample_size-1];

    double d = ((expT-expmin) - i*delta_exp)/delta_exp;

    return (1.0-d)*Us[r][k][i] + d*Us[r][k][i+1];

  } else {
    fprintf(stdout,"\033[0;31m*** Error: Encountered unknown interpolation code (%d).\n\033[0m", 
            interpolation);
    exit(-1);
  }

  return 0.0;

}

//--------------------------------------------------------------------------




