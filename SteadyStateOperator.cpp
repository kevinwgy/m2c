/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <SteadyStateOperator.h>
#include <GlobalMeshInfo.h>
#include <IoData.h>
#include <Utils.h>
#include <cassert>

extern int INACTIVE_MATERIAL_ID;

//--------------------------------------------------------------------------

SteadyStateOperator::SteadyStateOperator(MPI_Comm &comm_, TsData &iod_ts_, GlobalMeshInfo &global_mesh_) :
                     comm(comm_), global_mesh(global_mesh_), ref_calculated(false), converged(false)
{
  Rtol = iod_ts_.convergence_tolerance;

  R1_init = R2_init = Rinf_init = -1.0; //a negative number means the residual has not been computed
  R1 = R2 = Rinf = -1.0;
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

  assert(R.NumDOF()==(int)Rref.size());

  vector<double> r1, r2, rinf;
  R.CalculateFunctionNormsConRec(ID, global_mesh, r1, r2, rinf);

//  fprintf(stdout, "r1 = %e %e %e, r2 = %e %e %e, rinf = %e %e %e.\n", r1[0], r1[1], r1[2], r2[0], r2[1], r2[2], rinf[0], rinf[1], rinf[2]);

  for(int i=0; i<R.NumDOF(); i++) { // divide each component by ref.
    r1[i]   /= Rref[i];
    r2[i]   /= Rref[i];
    rinf[i] /= Rref[i];
  }

  R1 = r1[0];
  R2 = r2[0]*r2[0];
  Rinf = rinf[0];
  for(int i=1; i<(int)r1.size(); i++) {
    R1 += r1[i];
    R2 += r2[i]*r2[i];
    Rinf = std::max(Rinf, rinf[i]);
  }
  R2 = sqrt(R2);

  if(R1_init<0 || R2_init<0 || Rinf_init<0) { // first time-step
    R1_init = R1; 
    R2_init = R2;
    Rinf_init = Rinf;
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

    print(comm, "- Initial residual: %e (Func. L1 norm), %e (L2), %e (inf).\n", R1_init, R2_init, Rinf_init);
  }

  
  // check for convergence 
  if(R1/R1_init<Rtol || R2/R2_init<Rtol || Rinf/Rinf_init<Rtol)
    converged = true;
  else
    converged = false;


}

//--------------------------------------------------------------------------

void
SteadyStateOperator::CalculateReferenceResidual(SpaceVariable3D &R, SpaceVariable3D &ID)
{
  int dof = R.NumDOF();
  Rref.assign(dof, 0.0);

  int i0, j0, k0, imax, jmax, kmax;
  R.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);

  double*** res = R.GetDataPointer();
  double*** id  = ID.GetDataPointer();  

  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {  

        if(id[k][j][i] == INACTIVE_MATERIAL_ID)
          continue; // skip inactive nodes/cells

        for(int p=0; p<dof; p++)
          Rref[p] = std::max(Rref[p], fabs(res[k][j][i*dof+p]));
      }

  // find global max
  MPI_Allreduce(MPI_IN_PLACE, Rref.data(), dof, MPI_DOUBLE, MPI_MAX, comm);
  
  for(int i=0; i<dof; i++)
    if(Rref[i]==0)
      Rref[i] = 1.0;
  
  print(comm, "- Calculated residual normalization factors:");
  for(auto&& rref : Rref)
    print(comm, " %e", rref);
  print(comm, ".\n");

  R.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();
}

//--------------------------------------------------------------------------









