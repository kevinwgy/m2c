/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <SteadyStateOperator.h>
#include <IoData.h>
#include <Vector5D.h>
#include <Utils.h>
#include <cassert>

extern int INACTIVE_MATERIAL_ID;

//--------------------------------------------------------------------------

SteadyStateOperator::SteadyStateOperator(MPI_Comm &comm_, TsData &iod_ts_) :
                     comm(comm_), ref_calculated(false), converged(false)
{
  Rtol = iod_ts_.convergence_tolerance;

  R1_init = R2_init = Rinf_init = -1.0; //a negative number means the residual has not been computed
  R1 = R2 = Rinf = -1.0;

  for(int i=0; i<5; i++)
    Rref[i] = 0.0;
}

//--------------------------------------------------------------------------

SteadyStateOperator::~SteadyStateOperator()
{ }

//--------------------------------------------------------------------------

void
SteadyStateOperator::Destroy()
{ }

//--------------------------------------------------------------------------

void
SteadyStateOperator::MonitorConvergence(SpaceVariable3D &R, SpaceVariable3D &ID)
{

  if(!ref_calculated) {
    CalculateReferenceResidual(R,ID);
    ref_calculated = true;
  }

  int i0, j0, k0, imax, jmax, kmax;
  R.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);

  assert(R.NumDOF()==5);
  
  Vec5D*** res = (Vec5D***)R.GetDataPointer();
  double*** id = ID.GetDataPointer();  

  Vec5D r1(0.0), r2(0.0), rinf(0.0);
  double s;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {  

        if(id[k][j][i] == INACTIVE_MATERIAL_ID)
          continue; //skip inactive nodes/cells

        for(int p=0; p<5; p++) {
          s = fabs(res[k][j][i][p]); 
          r1[p] += s;
          r2[p] += s*s;
          rinf[p] = std::max(rinf[p], s);
        }
      }

  MPI_Allreduce(MPI_IN_PLACE, r1, 5, MPI_DOUBLE, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, r2, 5, MPI_DOUBLE, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, rinf, 5, MPI_DOUBLE, MPI_MAX, comm);

  for(int i=0; i<5; i++)
    r2[i] = sqrt(r2[i]);

  for(int i=0; i<5; i++) { // divide each component by ref.
    r1[i]   /= Rref[i];
    r2[i]   /= Rref[i];
    rinf[i] /= Rref[i];
  }

  R1 = r1.norm1();
  R2 = r2.norm2();
  Rinf = rinf.norminf();

  if(R1_init<0 || R2_init<0 || Rinf_init<0) { // first time-step
    R1_init = R1; 
    R2_init = R2;
    Rinf_init = Rinf_init;
    bool found_zero = false;
    if(R1_init == 0.0) {
      found_zero = true;
      R1_init = 1.0;
    }
    if(R2_init == 0.0) {
      found_zero = true;
      R2_init = 1.0;
    }
    if(Rinf_init == 0.0) {
      found_zero = true;
      Rinf_init = 1.0;
    }
    if(found_zero)
      print_warning(comm, "Warning: Found zero residual. Replaced by 1.\n"); 

    print(comm, "- Initial residual: %e (2-norm), %e (inf-norm).\n", R2_init, Rinf_init);
  }

  
  // check for convergence (only consider 2-norm and inf-norm at this point)
  if(R2/R2_init<Rtol || Rinf/Rinf_init<Rtol)
    converged = true;
  else
    converged = false;


  R.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();

}

//--------------------------------------------------------------------------

void
SteadyStateOperator::CalculateReferenceResidual(SpaceVariable3D &R, SpaceVariable3D &ID)
{
  int i0, j0, k0, imax, jmax, kmax;
  R.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  assert(R.NumDOF()==5);

  Vec5D*** res = (Vec5D***)R.GetDataPointer();
  double*** id = ID.GetDataPointer();  

  for(int i=0; i<5; i++)
    Rref[i] = 0.0;

  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {  

        if(id[k][j][i] == INACTIVE_MATERIAL_ID)
          continue; // skip inactive nodes/cells

        for(int p=0; p<5; p++)
          Rref[p] = std::max(Rref[p], fabs(res[k][j][i][p]));
      }

  // find global max
  MPI_Allreduce(MPI_IN_PLACE, Rref, 5, MPI_DOUBLE, MPI_MAX, comm);
  
  // use the same ref for x-, y-, and z-momentum
  double Rm = std::max(Rref[1], std::max(Rref[2], Rref[3]));
  Rref[1] = Rref[2] = Rref[3] = Rm;

  bool found_zero_ref(false);
  for(int i=0; i<5; i++)
    if(Rref[i]==0) {
      found_zero_ref = true;
      Rref[i] = 1.0;
    }
  if(found_zero_ref)
    print_warning(comm, "Warning: Found zero residual component. Replaced by 1.\n");
  
  print(comm, "- N-S residual normalization factors: %e %e %e.\n", Rref[0], Rref[1], Rref[4]);

  R.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();
}

//--------------------------------------------------------------------------









