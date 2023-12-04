/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<IncompressibleOperator.h>
using std::vector;

//--------------------------------------------------------------------------

IncompressibleOperator::IncompressibleOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_,
                                               vector<VarFcnBase*> &varFcn_, SpaceOperator &spo_, 
                                               InterpolatorBase &interpolator_)
                      : vf(varFcn_), spo(spo_), interpolator(interp_), gfo(NULL)
{
  // Get i0, j0, etc.
  SpaceVariable3D &coordinates(spo.GetMeshCoordinates());
  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);
  coordinates.GetGlobalSize(&NX,&NY,&NZ);

  CheckInputs(iod_);
}

//--------------------------------------------------------------------------

IncompressibleOperator::~IncompressibleOperator()
{
  if(gfo) delete gfo;
}

//--------------------------------------------------------------------------

void
IncompressibleOperator::Destroy()
{
  if(gfo)
    gfo->Destroy();
}

//--------------------------------------------------------------------------

void
IncompressibleOperator::CheckInputs(IoData &iod)
{
  // Note: This function checks some inputs, but CANNOT detect all possible error...
  //       Therefore, passing this check does not mean there is no error.

  // -----------------------------
  // Check material models 
  // -----------------------------
  for(auto&& mat : iod.eqs.materials.dataMap()) {
    if(mat.second->viscosity.type != ViscosityModelData::NONE &&
       mat.second->viscosity.type != ViscosityModelData::CONSTANT) {
      print_error("*** Error: The incompressible flow solver only supports a constant diffusivity for "
                  "each materia.\n");
      exit_mpi();
    }
    if(mat.second->viscosity.bulkViscosity != 0.0) {
      print_error("*** Error: For incompressible flows, bulk viscosity is irrelevant. "
                  "Detected non-zero bulk viscosity coefficient.\n");
      exit_mpi();
    }
  } 

  // -----------------------------
  // Check solver
  // -----------------------------
  // TODO  




  double default_density = 1.0e-6; //based on IoData (StateVaraibles::StateVariable())
  // -----------------------------
  // Check initial conditions
  // -----------------------------
  int ic_error = 0;
  if(iod.ic.default_ic.density  != default_density)      ic_error++;
  if(iod.ic.default_ic.pressure != 0.0)                  ic_error++;
  if(iod.ic.default_ic.internal_energy_per_mass != 0.0)  ic_error++;

  MultiInitialConditionsData &ic(iod.ic.multiInitialConditions);
  for(auto&& obj : ic.planeMap.dataMap) {
    if(obj.second->initialConditions.density != default_density)       ic_error++;
    if(obj.second->initialConditions.pressure != 0.0)                  ic_error++;
    if(obj.second->initialConditions.internal_energy_per_mass != 0.0)  ic_error++;
  }
  for(auto&& obj : ic.cylinderconeMap.dataMap) {
    if(obj.second->initialConditions.density != default_density)       ic_error++;
    if(obj.second->initialConditions.pressure != 0.0)                  ic_error++;
    if(obj.second->initialConditions.internal_energy_per_mass != 0.0)  ic_error++;
  }
  for(auto&& obj : ic.cylindersphereMap.dataMap) {
    if(obj.second->initialConditions.density != default_density)       ic_error++;
    if(obj.second->initialConditions.pressure != 0.0)                  ic_error++;
    if(obj.second->initialConditions.internal_energy_per_mass != 0.0)  ic_error++;
  }
  for(auto&& obj : ic.sphereMap.dataMap) {
    if(obj.second->initialConditions.density != default_density)       ic_error++;
    if(obj.second->initialConditions.pressure != 0.0)                  ic_error++;
    if(obj.second->initialConditions.internal_energy_per_mass != 0.0)  ic_error++;
  }
  for(auto&& obj : ic.parallelepipedMap.dataMap) {
    if(obj.second->initialConditions.density != default_density)       ic_error++;
    if(obj.second->initialConditions.pressure != 0.0)                  ic_error++;
    if(obj.second->initialConditions.internal_energy_per_mass != 0.0)  ic_error++;
  }
  for(auto&& obj : ic.spheroidMap.dataMap) {
    if(obj.second->initialConditions.density != default_density)       ic_error++;
    if(obj.second->initialConditions.pressure != 0.0)                  ic_error++;
    if(obj.second->initialConditions.internal_energy_per_mass != 0.0)  ic_error++;
  }
  for(auto&& obj : ic.enclosureMap.dataMap) {
    if(obj.second->initialConditions.density != default_density)       ic_error++;
    if(obj.second->initialConditions.pressure != 0.0)                  ic_error++;
    if(obj.second->initialConditions.internal_energy_per_mass != 0.0)  ic_error++;
  }
   
  if(ic_error>0) {
    print_error("*** Error: The incompressible flow solver does not accept initial values for "
                "density, pressure, or internal energy. Detected %d violations.\n", ic_error);
    exit_mpi();
  }


  // -----------------------------
  // Check boundary conditions
  // -----------------------------
  int bc_error = 0;
  if(iod.bc.inlet.density != default_density)       bc_error++;
  if(iod.bc.inlet.pressure != 0.0)                  bc_error++;
  if(iod.bc.inlet.internal_energy_per_mass != 0.0)  bc_error++;
  if(iod.bc.outlet.density != default_density)      bc_error++;
  if(iod.bc.outlet.pressure != 0.0)                 bc_error++;
  if(iod.bc.outlet.internal_energy_per_mass != 0.0) bc_error++;

  for(auto&& obj : iod.bc.multiBoundaryConditions.diskMap.dataMap)  {
    if(obj.second->state.density != default_density)       bc_error++;
    if(obj.second->state.pressure != 0.0)                  bc_error++;
    if(obj.second->state.internal_energy_per_mass != 0.0)  bc_error++; 
  }
  for(auto&& obj : iod.bc.multiBoundaryConditions.rectangleMap.dataMap)  {
    if(obj.second->state.density != default_density)       bc_error++;
    if(obj.second->state.pressure != 0.0)                  bc_error++;
    if(obj.second->state.internal_energy_per_mass != 0.0)  bc_error++; 
  }

  if(ic_error>0) {
    print_error("*** Error: The incompressible flow solver does not accept boundary values for "
                "density, pressure, or internal energy. Detected %d violations.\n", ic_error);
    exit_mpi();
  }

}

//--------------------------------------------------------------------------

void
IncompressibleOperator::FinalizeInitialCondition(SpaceVariable3D &V, SpaceVariable3D &ID)
{
  Vec5D***   v = (Vec5D***) V.GetDataPointer();
  double*** id = ID.GetDataPointer();

  // only update the domain interior
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        v[k][j][i][0] = vf[id[k][j][i]]->GetDensity(0.0,0.0); //get rho0
        v[k][j][i][4] = 0.0; //initialize p = 0
      }

  V.RestoreDataPointerAndInsert();
  ID.RestoreDataPointerToLocalVector();
}

//--------------------------------------------------------------------------

void
IncompressibleOperator::ApplyBoundaryConditions(SpaceVariable3D &V)
{
  spo.ApplyBoundaryConditions(V);
}

//--------------------------------------------------------------------------





