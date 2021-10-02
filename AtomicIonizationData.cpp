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
    if(*it)
      delete *it;
}

//--------------------------------------------------------------------------

void
AtomicIonizationData::Setup(AtomicIonizationModel* iod_aim, double h_, double e_, double me_, double kb_,
                            int interp, double Tmin_, double Tmax_, int sample_size_, MPI_Comm* comm)
{

  // Step 0: Get constants
  h = h_;
  e = e_;
  me = me_;
  kb = kb_;

  // Step 1: Read molar fraction, atomic number, and max charge
  molar_fraction = iod_aim->molar_fraction;
  atomic_number = iod_aim->atomic_number;
  rmax = iod_aim->max_charge - 1;
  if(molar_fraction<=0 || atomic_number<=0) {
    print_error("*** Error: Detected molar fraction %e, atomic number %d. (Should be greater than 0.)\n",
                molar_fraction, atomic_number);
    exit_mpi();
  }
  print("  o Found chemical element with atomic number %d and molar fraction %e.\n", atomic_number,
        molar_fraction);


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
  int max_size = 0;
  for(int i=0; i<atomic_number; i++) {
    std::string filename = std::string(iod_aim->excitation_energy_files_prefix) 
                         + std::to_string(i)
                         + std::string(iod_aim->excitation_energy_files_suffix);
    file.open(filename.c_str(), std::fstream::in);
    if(!file.is_open()) {
      print_error("*** Error: Cannot open excitation energy file %s.\n", filename.c_str());
      exit_mpi();
    }
    GetDataInFile(file, E[i], 10000, true);
    file.close();
    if(max_size<E[i].size())
      max_size = E[i].size();
  }
  print("    * Read %d excitation energy files: %sX%s. Max excited state: %d.\n", atomic_number,
        iod_aim->excitation_energy_files_prefix, 
        iod_aim->excitation_energy_files_suffix, max_size);


  // Step 4: Read angular momentum energy files ("l") and convert to the degeneracy, g = 2l + 1
  g.resize(atomic_number);
  max_size = 0;
  for(int i=0; i<atomic_number; i++) {
    std::string filename = std::string(iod_aim->angular_momentum_files_prefix)
                         + std::to_string(i)
                         + std::string(iod_aim->angular_momentum_files_suffix);
    file.open(filename.c_str(), std::fstream::in);
    if(!file.is_open()) {
      print_error("*** Error: Cannot open angular momentum file %s.\n", filename.c_str());
      exit_mpi();
    }
    GetDataInFile(file, g[i], 10000, true);
    //convert g to l
    for(auto it = g[i].begin(); it != g[i].end(); it++)
      *it = 2.0*(*it) + 1.0;
    file.close();
    if(max_size<g[i].size())
      max_size = g[i].size();
  }
  print("    * Read %d angular momentum files: %sX%s. Max excited state: %d.\n", atomic_number,
        iod_aim->angular_momentum_files_prefix, 
        iod_aim->angular_momentum_files_suffix, max_size);

      
  // check to see if rmax should be reduced based on data files
  int data_size = std::min(I.size(), std::min(g.size(), E.size()));
  if(rmax > data_size-1) //rmax should be smaller than or equal to "size" - 1
    rmax = data_size-1;


  // if interpolation is true, create the interpolation database
  interpolation = interp;
  if(interpolation && Tmax_<=Tmin_){
    print_error("*** Error: Incorrect Tmin and/or Tmax in atomic ionization model (atomic no. %d). "
                "Tmin must be less than Tmax.\n", atomic_number);
    exit_mpi();
  } 
  Tmin = Tmin_;
  Tmax = Tmax_;
  sample_size = sample_size_;
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
  for(int r=0; r<=rmax; r++) {
    Us.push_back(vector<double>(sample_size));
  }
  UsCoeffs.resize(rmax+1, std::tuple<double,double,double>(0,0,0));

  if(interpolation == 1) //cubic spline interpolation
    spline.resize(rmax+1, NULL);

  for(int r=0; r<=rmax; r++)
    InitializeInterpolationForCharge(r, comm);

}

//--------------------------------------------------------------------------

void
AtomicIonizationData::InitializeInterpolationForCharge(int r, MPI_Comm &comm)
{

  vector<double>& U(Us[r]); 
  std::tuple<double,double,double>& coeffs(UsCoeffs[r]);

  double factor(0);
  for(int i=0; i<E[r].size(); i++){
    factor = -E[r][i]/kb;
    if(factor!=0)
      break;
  }

  double expmin = exp(factor/Tmin);
  double expmax = exp(factor/Tmax);
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

  //prepare for comm.
  int* counts = new int[mpi_size];
  int* displacements = new int[mpi_size];
  for(int i=0; i<mpi_size; i++) {
    counts[i] = (i<remainder) ? block_size + 1 : block_size;
    displacements[i] = (i<remainder) ? (block_size+1)*i : block_size*i + remainder;
  }

  int my_start_id = displacements[mpi_rank];
  int my_block_size = counts[mpi_rank];
  double T;
  for(int i=0; i<my_block_size; i++) {
    T = (my_start_id+i==0) ? Tmin : factor/log(expmin+(my_start_id+i)*delta_exp); //expmin can be 0 if Tmin is small (but nonzero)
    U[my_start_id+i] = CalculatePartitionFunctionOnTheFly(r, T);
  }


  //communication
  MPI_Allgatherv(MPI_IN_PLACE, my_block_size, MPI_DOUBLE, U.data(), counts, displacements, MPI_DOUBLE, comm);
 

  //initialize spline interpolation
  if(interpolation == 1)
    spline[r] = new boost::math::cubic_b_spline<double>(U.begin(), U.end(), expmin, delta_exp);


  delete [] counts;
  delete [] displacements;

  exit_mpi();
}

//--------------------------------------------------------------------------

double
AtomicIonizationData::CalculatePartitionFunction(int r, double T)
{
  assert(r<=rmax);
  if(interpolation && T>=Tmin && T<=Tmax)
    return CalculatePartitionFunctionByInterpolation(r, T);

  return CalculatePartitionFunctionOnTheFly(r, T);
}

//--------------------------------------------------------------------------

double
AtomicIonizationData::CalculatePartitionFunctionOnTheFly(int r, double T)
{

  if(r==atomic_number)
    return 1.0; //per Shafquat

  double Ur = 0.0;
  int nsize = std::max(g[r].size(), E[r].size());
  for(int n=0; n<nsize; n++) {
    if(I[r]<E[r][n]) {
      assert(n>0);
      break;
    }
    Ur += g[r][n]*exp(-E[r][n]/(kb*T));
  }

  return Ur;

}

//--------------------------------------------------------------------------

double
AtomicIonizationData::CalculatePartitionFunctionByInterpolation(int r, double T)
{

  double factor = std::get<0>(UsCoeffs[r]);
  double expT = exp(factor/T);

  if(interpolation==1) 

    return (*spline[r])(expT); 

  else if(interpolation==2) { //linear

    double expmin = std::get<1>(UsCoeffs[r]);
    double delta_exp = std::get<2>(UsCoeffs[r]);

    int i = floor((expT - expmin)/delta_exp);
    assert(i>=0 && i<=sample_size-1);

    if(i==sample_size-1)
      return Us[r][sample_size-1];

    double d = ((expT-expmin) - i*delta_exp)/delta_exp;

    return (1.0-d)*Us[r][i] + d*Us[r][i+1];

  } else {
    fprintf(stderr,"\033[0;31m*** Error: Encountered unknown interpolation code (%d).\n\033[0m", 
            interpolation);
    exit(-1);
  }

  return 0.0;

}

//--------------------------------------------------------------------------




