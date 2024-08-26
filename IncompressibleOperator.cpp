/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<IncompressibleOperator.h>
#include<LinearOperator.h>
#include<GeoTools.h>
using std::vector;

using namespace GeoTools;

extern int verbose;
extern int INACTIVE_MATERIAL_ID;

//--------------------------------------------------------------------------

IncompressibleOperator::IncompressibleOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_,
                                               vector<VarFcnBase*> &varFcn_, SpaceOperator &spo_, 
                                               InterpolatorBase &interp_)
                      : comm(comm_), iod(iod_), vf(varFcn_), spo(spo_), interpolator(interp_), gfo(NULL),
                        V3(comm_, &(dm_all_.ghosted1_3dof))
{
  // Get i0, j0, etc.
  SpaceVariable3D &coordinates(spo.GetMeshCoordinates());
  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);
  coordinates.GetGlobalSize(&NX,&NY,&NZ);

  CheckInputs(iod_);

  // Fill mu (dynamic viscosity)
  Mu.assign(vf.size(), 0.0); //by default, mu = 0
  ViscosityOperator *visco_ptr = spo.GetPointerToViscosityOperator();
  if(visco_ptr) {
    vector<ViscoFcnBase*>& visFcns(visco_ptr->GetViscoFcns());
    for(int i=0; i<(int)visFcns.size(); i++) {
      if(visFcns[i]->type == ViscoFcnBase::CONSTANT)
        Mu[i] = visFcns[i]->GetMu();
      else if(visFcns[i]->type != ViscoFcnBase::NONE) {
        print_error("*** Error: Currently, IncompressibleOperator only supports CONSTANT viscosity.\n");
        exit_mpi();
      }
    }
  }

  // Calculate deltas
  GlobalMeshInfo& global_mesh(spo.GetGlobalMeshInfo());
  Dx.resize(3); //dir = 0, 1, 2
  Dy.resize(3);
  Dz.resize(3);
  dx_l.resize(3);
  dx_r.resize(3);
  dy_b.resize(3);
  dy_t.resize(3);
  dz_k.resize(3);
  dz_f.resize(3);
  for(int d=0; d<3; d++) {
    for(int i=i0; i<imax; i++) {
      Dx[d].push_back(d==0 ? 0.5*(global_mesh.GetDx(i-1)+global_mesh.GetDx(i)) : global_mesh.GetDx(i));
      dx_l[d].push_back(d==0 ? global_mesh.GetDx(i-1) : 0.5*(global_mesh.GetDx(i-1)+global_mesh.GetDx(i)));
      dx_r[d].push_back(d==0 ? global_mesh.GetDx(i)   : 0.5*(global_mesh.GetDx(i)+global_mesh.GetDx(i+1)));
    }
    for(int j=j0; j<jmax; j++) {
      Dy[d].push_back(d==1 ? 0.5*(global_mesh.GetDy(j-1)+global_mesh.GetDy(j)) : global_mesh.GetDy(j));
      dy_b[d].push_back(d==1 ? global_mesh.GetDy(j-1) : 0.5*(global_mesh.GetDy(j-1)+global_mesh.GetDy(j)));
      dy_t[d].push_back(d==1 ? global_mesh.GetDy(j)   : 0.5*(global_mesh.GetDy(j)+global_mesh.GetDy(j+1)));
    }
    for(int k=k0; k<kmax; k++) {
      Dz[d].push_back(d==2 ? 0.5*(global_mesh.GetDz(k-1)+global_mesh.GetDz(k)) : global_mesh.GetDz(k));
      dz_k[d].push_back(d==2 ? global_mesh.GetDz(k-1) : 0.5*(global_mesh.GetDz(k-1)+global_mesh.GetDz(k)));
      dz_f[d].push_back(d==2 ? global_mesh.GetDz(k)   : 0.5*(global_mesh.GetDz(k)+global_mesh.GetDz(k+1)));
    }
  }

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
  V3.Destroy();
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
  for(auto&& mat : iod.eqs.materials.dataMap) {
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
  // Check mesh 
  // -----------------------------

  if(iod.mesh.type == MeshData::CYLINDRICAL || iod.mesh.type == MeshData::SPHERICAL) {
    print_error("*** Error: The incompressible flow solver does not support cylindrical or spherical "
                "symmetry at this moment.\n");
    exit_mpi();
  }
  if(iod.mesh.bc_x0 == MeshData::OVERSET || iod.mesh.bc_xmax == MeshData::OVERSET ||
     iod.mesh.bc_y0 == MeshData::OVERSET || iod.mesh.bc_ymax == MeshData::OVERSET ||
     iod.mesh.bc_z0 == MeshData::OVERSET || iod.mesh.bc_zmax == MeshData::OVERSET) {
    print_error("*** Error: The incompressible flow solver does not support the overset method "
                "at this moment.\n");
    exit_mpi();
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
  if(iod.bc.inlet2.density != default_density)      bc_error++;
  if(iod.bc.inlet2.pressure != 0.0)                 bc_error++;
  if(iod.bc.inlet2.internal_energy_per_mass != 0.0) bc_error++;

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

  // Calculate coefficients for interpolating velocity (i.e. shifting to cell boundaries)
  GlobalMeshInfo& global_mesh(spo.GetGlobalMeshInfo());
  vector<double> cx0(imax-i0, 0.0), cx1(imax-i0, 0.0);
  vector<double> cy0(jmax-j0, 0.0), cy1(jmax-j0, 0.0);
  vector<double> cz0(kmax-k0, 0.0), cz1(kmax-k0, 0.0);
  double sumD;
  for(int i=i0; i<imax; i++) {
    sumD = global_mesh.GetDx(i-1) + global_mesh.GetDx(i);
    cx0[i-i0] = global_mesh.GetDx(i)/sumD;
    cx1[i-i0] = global_mesh.GetDx(i-1)/sumD;
  }
  for(int j=j0; j<jmax; j++) {
    sumD = global_mesh.GetDy(j-1) + global_mesh.GetDy(j);
    cy0[j-j0] = global_mesh.GetDy(j)/sumD;
    cy1[j-j0] = global_mesh.GetDy(j-1)/sumD;
  }
  for(int k=k0; k<kmax; k++) {
    sumD = global_mesh.GetDz(k-1) + global_mesh.GetDz(k);
    cz0[k-k0] = global_mesh.GetDz(k)/sumD;
    cz1[k-k0] = global_mesh.GetDz(k-1)/sumD;
  }

  
  // Backup velocity at cell centers
  Vec5D***   v = (Vec5D***) V.GetDataPointer();
  Vec3D***  v3 = (Vec3D***) V3.GetDataPointer();
  double*** id = ID.GetDataPointer();

  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++) {
        v3[k][j][i][0] = v[k][j][i][1];  
        v3[k][j][i][1] = v[k][j][i][2];  
        v3[k][j][i][2] = v[k][j][i][3];  
      }


  // Only update the domain interior
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        v[k][j][i][0] = vf[id[k][j][i]]->GetDensity(0.0,0.0); //get rho0
        v[k][j][i][1] = cx0[i-i0]*v3[k][j][i-1][0] + cx1[i-i0]*v3[k][j][i][0];
        v[k][j][i][2] = cy0[j-j0]*v3[k][j-1][i][1] + cy1[j-j0]*v3[k][j][i][1];
        v[k][j][i][3] = cz0[k-k0]*v3[k-1][j][i][2] + cz1[k-k0]*v3[k][j][i][2];
        v[k][j][i][4] = 0.0; //initialize p = 0
      }

  //Note: At this point some velocity values ON THE DOMAIN BOUNDARIES ARE INCORRECT.
  //      They NEED TO BE CORRECTED by calling ApplyBoundaryConditions.

  V.RestoreDataPointerAndInsert();
  V3.RestoreDataPointerToLocalVector(); //no need to exchange data
  ID.RestoreDataPointerToLocalVector();

/*
  //debug
  v = (Vec5D***) V.GetDataPointer();
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jj0+3; j++)
      for(int i=ii0; i<iimax; i++)
        fprintf(stdout,"v[%d][%d][%d] = %e %e %e %e %e.\n", k,j,i, v[k][j][i][0], v[k][j][i][1], v[k][j][i][2], v[k][j][i][3], v[k][j][i][4]);
*/  
}

//--------------------------------------------------------------------------

void
IncompressibleOperator::ApplyBoundaryConditions(SpaceVariable3D &V)
{
  // Note: Modify both ghost and non-ghost entries due to the use of MAC / staggered grids  
  //       Only modify velocity components

  Vec5D*** v = (Vec5D***) V.GetDataPointer();

  GlobalMeshInfo &global_mesh(spo.GetGlobalMeshInfo());
  vector<GhostPoint>* ghost_nodes_outer(spo.GetPointerToOuterGhostNodes());

  for(auto it = ghost_nodes_outer->begin(); it != ghost_nodes_outer->end(); it++) {

    if(it->type_projection != GhostPoint::FACE)
      continue; //corner (edge or vertex) nodes are not populated

    // preparation
    int i(it->ijk[0]), j(it->ijk[1]), k(it->ijk[2]);
    int im_i(it->image_ijk[0]), im_j(it->image_ijk[1]), im_k(it->image_ijk[2]);

    Vec3D v0(0.0); //local velocity at the boundary
    if(it->bcType == MeshData::INLET) {
      v0[0] = iod.bc.inlet.velocity_x;
      v0[1] = iod.bc.inlet.velocity_y;
      v0[2] = iod.bc.inlet.velocity_z;
    }
    else if(it->bcType == MeshData::INLET2) {
      v0[0] = iod.bc.inlet2.velocity_x;
      v0[1] = iod.bc.inlet2.velocity_y;
      v0[2] = iod.bc.inlet2.velocity_z;
    }


    // Impose boundary conditions for velocity
    // Treating the six sides separately might seem clumsy. But in some cases imposing a boundary
    // condition requires populating two entries. 
    if(it->side == GhostPoint::LEFT) {
      if(it->bcType == MeshData::INLET || it->bcType == MeshData::INLET2) {
        v[k][j][im_i][1] = v0[0];
        v[k][j][i][2]    = v0[1];
        v[k][j][i][3]    = v0[2];
      }
      else if(it->bcType == MeshData::OUTLET) {
        v[k][j][im_i][1] = v[k][j][im_i+1][1];
        v[k][j][i][2]    = v[k][j][im_i][2];
        v[k][j][i][3]    = v[k][j][im_i][3];
      }
      else if(it->bcType == MeshData::SLIPWALL || it->bcType == MeshData::SYMMETRY) {
        v[k][j][im_i][1] = 0.0;
        v[k][j][i][2]    = v[k][j][im_i][2];
        v[k][j][i][3]    = v[k][j][im_i][3];
      }
      else if(it->bcType == MeshData::STICKWALL) {
        v[k][j][im_i][1] = 0.0;
        v[k][j][i][2]    = -v[k][j][im_i][2];
        v[k][j][i][3]    = -v[k][j][im_i][3];
      }
      else if(it->bcType == MeshData::OVERSET) {
        // nothing to be done here. ghost nodes will be populated below
      } else {
        fprintf(stdout,"*** Error: Detected unknown boundary condition type (%d).\n", (int)it->bcType);
        exit(-1);
      }
    }
    else if(it->side == GhostPoint::RIGHT) {
      if(it->bcType == MeshData::INLET || it->bcType == MeshData::INLET2) {
        v[k][j][i][1]    = v0[0];
        v[k][j][i][2]    = v0[1];
        v[k][j][i][3]    = v0[2];
      }
      else if(it->bcType == MeshData::OUTLET) {
        v[k][j][i][1]    = v[k][j][im_i][1];
        v[k][j][i][2]    = v[k][j][im_i][2];
        v[k][j][i][3]    = v[k][j][im_i][3];
      }
       else if(it->bcType == MeshData::SLIPWALL || it->bcType == MeshData::SYMMETRY) {
        v[k][j][i][1]    = 0.0;
        v[k][j][i][2]    = v[k][j][im_i][2];
        v[k][j][i][3]    = v[k][j][im_i][3];
      }
      else if(it->bcType == MeshData::STICKWALL) {
        v[k][j][i][1]    = 0.0;
        v[k][j][i][2]    = -v[k][j][im_i][2];
        v[k][j][i][3]    = -v[k][j][im_i][3];
      }
      else if(it->bcType == MeshData::OVERSET) {
        // nothing to be done here. ghost nodes will be populated below
      } else {
        fprintf(stdout,"*** Error: Detected unknown boundary condition type (%d).\n", (int)it->bcType);
        exit(-1);
      }
    }
    else if(it->side == GhostPoint::BOTTOM) {
      if(it->bcType == MeshData::INLET || it->bcType == MeshData::INLET2) {
        v[k][j][i][1]    = v0[0];
        v[k][im_j][i][2] = v0[1];
        v[k][j][i][3]    = v0[2];
      }
      else if(it->bcType == MeshData::OUTLET) {
        v[k][j][i][1]    = v[k][im_j][i][1];
        v[k][im_j][i][2] = v[k][im_j+1][i][2];
        v[k][j][i][3]    = v[k][im_j][i][3];
      }
      else if(it->bcType == MeshData::SLIPWALL || it->bcType == MeshData::SYMMETRY) {
        v[k][j][i][1]    = v[k][im_j][i][1];
        v[k][im_j][i][2] = 0.0;
        v[k][j][i][3]    = v[k][im_j][i][3];
      }
      else if(it->bcType == MeshData::STICKWALL) {
        v[k][j][i][1]    = -v[k][im_j][i][1];
        v[k][im_j][i][2] = 0.0;
        v[k][j][i][3]    = -v[k][im_j][i][3];
      }
      else if(it->bcType == MeshData::OVERSET) {
        // nothing to be done here. ghost nodes will be populated below
      } else {
        fprintf(stdout,"*** Error: Detected unknown boundary condition type (%d).\n", (int)it->bcType);
        exit(-1);
      }

      // special treatment for flat plate
      if(iod.rans.example == RANSTurbulenceModelData::FLAT_PLATE) {
        double x = global_mesh.GetX(i);
        if(x>=iod.rans.example_param_1) {
          if(x-0.5*global_mesh.GetDx(i) >= iod.rans.example_param_1) //check x-velo grid point
            v[k][j][i][1] = -v[k][im_j][i][1];
          v[k][im_j][i][2] = 0.0; //y-velocity is stored on the bottom facet of the cell
          v[k][j][i][3]    = -v[k][im_j][i][3]; 
        }
      }

    }
    else if(it->side == GhostPoint::TOP) {
      if(it->bcType == MeshData::INLET || it->bcType == MeshData::INLET2) {
        v[k][j][i][1]    = v0[0];
        v[k][j][i][2]    = v0[1];
        v[k][j][i][3]    = v0[2];
      }
      else if(it->bcType == MeshData::OUTLET) {
        v[k][j][i][1]    = v[k][im_j][i][1];
        v[k][j][i][2]    = v[k][im_j][i][2];
        v[k][j][i][3]    = v[k][im_j][i][3];
      }
      else if(it->bcType == MeshData::SLIPWALL || it->bcType == MeshData::SYMMETRY) {
        v[k][j][i][1]    = v[k][im_j][i][1];
        v[k][j][i][2]    = 0.0;
        v[k][j][i][3]    = v[k][im_j][i][3];
      }
      else if(it->bcType == MeshData::STICKWALL) {
        v[k][j][i][1]    = -v[k][im_j][i][1];
        v[k][j][i][2]    = 0.0;
        v[k][j][i][3]    = -v[k][im_j][i][3];
      }
      else if(it->bcType == MeshData::OVERSET) {
        // nothing to be done here. ghost nodes will be populated below
      } else {
        fprintf(stdout,"*** Error: Detected unknown boundary condition type (%d).\n", (int)it->bcType);
        exit(-1);
      }
    }
    else if(it->side == GhostPoint::BACK) {
      if(it->bcType == MeshData::INLET || it->bcType == MeshData::INLET2) {
        v[k][j][i][1]    = v0[0];
        v[k][i][i][2]    = v0[1];
        v[im_k][j][i][3] = v0[2];
      }
      else if(it->bcType == MeshData::OUTLET) {
        v[k][j][i][1]    = v[im_k][j][i][1];
        v[k][j][i][2]    = v[im_k][j][i][2];
        v[im_k][j][i][3] = v[im_k+1][j][i][3];
      }
      else if(it->bcType == MeshData::SLIPWALL || it->bcType == MeshData::SYMMETRY) {
        v[k][j][i][1]    = v[im_k][j][i][1];
        v[k][j][i][2]    = v[im_k][j][i][2];
        v[im_k][j][i][3] = 0.0;
      }
      else if(it->bcType == MeshData::STICKWALL) {
        v[k][j][i][1]    = -v[im_k][j][i][1];
        v[k][j][i][2]    = -v[im_k][j][i][2];
        v[im_k][j][i][3] = 0.0;
      }
      else if(it->bcType == MeshData::OVERSET) {
        // nothing to be done here. ghost nodes will be populated below
      } else {
        fprintf(stdout,"*** Error: Detected unknown boundary condition type (%d).\n", (int)it->bcType);
        exit(-1);
      }
    }
    else if(it->side == GhostPoint::FRONT) {
      if(it->bcType == MeshData::INLET || it->bcType == MeshData::INLET2) {
        v[k][j][i][1]    = v0[0];
        v[k][j][i][2]    = v0[1];
        v[k][j][i][3]    = v0[2];
      }
      else if(it->bcType == MeshData::OUTLET) {
        v[k][j][i][1]    = v[im_k][j][i][1];
        v[k][j][i][2]    = v[im_k][j][i][2];
        v[k][j][i][3]    = v[im_k][j][i][3];
      }
      else if(it->bcType == MeshData::SLIPWALL || it->bcType == MeshData::SYMMETRY) {
        v[k][j][i][1]    = v[im_k][j][i][1];
        v[k][j][i][2]    = v[im_k][j][i][2];
        v[k][j][i][3]    = 0.0;
      }
      else if(it->bcType == MeshData::STICKWALL) {
        v[k][j][i][1]    = -v[im_k][j][i][1];
        v[k][j][i][2]    = -v[im_k][j][i][2];
        v[k][j][i][3]    = 0.0;
      }
      else if(it->bcType == MeshData::OVERSET) {
        // nothing to be done here. ghost nodes will be populated below
      } else {
        fprintf(stdout,"*** Error: Detected unknown boundary condition type (%d).\n", (int)it->bcType);
        exit(-1);
      }
    }

  }



/*
  // update overset ghosts (if any)
  for(auto&& g : ghost_overset)
    v[g.first[2]][g.first[1]][g.first[0]] = g.second;
*/

  ApplyBoundaryConditionsGeometricEntities(v);




  V.RestoreDataPointerAndInsert();

}

//--------------------------------------------------------------------------

void
IncompressibleOperator::ApplyBoundaryConditionsGeometricEntities(Vec5D*** v)
{

  // Very similar to the same function in SpaceOperator, but handles staggered grids
  // Only applies velocity boundary conditions

  map<int, DiskData* >& disks(iod.bc.multiBoundaryConditions.diskMap.dataMap);
  map<int, RectangleData* >& rectangles(iod.bc.multiBoundaryConditions.rectangleMap.dataMap);

  if(!disks.size() && !rectangles.size())
    return;

  GlobalMeshInfo &global_mesh(spo.GetGlobalMeshInfo());

  if(ii0==-1) { 
    if (iod.mesh.bc_x0 == MeshData::INLET    || iod.mesh.bc_x0 == MeshData::INLET2 ||
        iod.mesh.bc_x0 == MeshData::SLIPWALL || iod.mesh.bc_x0 == MeshData::STICKWALL) {

      vector<DiskData* > mydisks;
      for(auto it=disks.begin(); it!=disks.end(); it++)
        if(it->second->cen_x == iod.mesh.x0) {
          Vec3D n(it->second->normal_x, it->second->normal_y, it->second->normal_z);
          if(fabs(n[0])/n.norm()>1-1e-8)
            mydisks.push_back(it->second); 
        }
      vector<RectangleData* > myrects;
      for(auto it=rectangles.begin(); it!=rectangles.end(); it++)       
        if(it->second->cen_x == iod.mesh.x0) {
          Vec3D n(it->second->normal_x, it->second->normal_y, it->second->normal_z);
          if(fabs(n[0])/n.norm()>1-1e-8)
            myrects.push_back(it->second); 
         }
      if(mydisks.size() || myrects.size()) {

        for(int k=kk0; k<kkmax; k++)
          for(int j=jj0; j<jjmax; j++) {

            if(global_mesh.OutsidePhysicalDomainAndUnpopulated(ii0,j,k))
              continue; //skip corner nodes

            for(int p=0; p<(int)mydisks.size(); p++) {
              if(IsPointInDisk(global_mesh.GetY(j), global_mesh.GetZ(k),
                               mydisks[p]->cen_y, mydisks[p]->cen_z, mydisks[p]->radius)){
                v[k][j][ii0+1][1] = mydisks[p]->state.velocity_x;
                v[k][j][ii0][2]   = mydisks[p]->state.velocity_y;
                v[k][j][ii0][3]   = mydisks[p]->state.velocity_z;
              }
            }

            for(int p=0; p<(int)myrects.size(); p++) {
              if(IsPointInRectangle(global_mesh.GetY(j), global_mesh.GetZ(k),
                                    myrects[p]->cen_y, myrects[p]->cen_z, myrects[p]->a, myrects[p]->b)){
                v[k][j][ii0+1][1] = myrects[p]->state.velocity_x;
                v[k][j][ii0][2]   = myrects[p]->state.velocity_y;
                v[k][j][ii0][3]   = myrects[p]->state.velocity_z;
              }
            }

          }
      }

    }
  }

  if(iimax==NX+1) { 
    if (iod.mesh.bc_xmax == MeshData::INLET    || iod.mesh.bc_xmax == MeshData::INLET2 ||
        iod.mesh.bc_xmax == MeshData::SLIPWALL || iod.mesh.bc_xmax == MeshData::STICKWALL) {

      vector<DiskData* > mydisks;
      for(auto it=disks.begin(); it!=disks.end(); it++)
        if(it->second->cen_x == iod.mesh.xmax) {
          Vec3D n(it->second->normal_x, it->second->normal_y, it->second->normal_z);
          if(fabs(n[0])/n.norm()>1-1e-8)
            mydisks.push_back(it->second); 
        } 
      vector<RectangleData* > myrects;
      for(auto it=rectangles.begin(); it!=rectangles.end(); it++)       
        if(it->second->cen_x == iod.mesh.xmax) {
          Vec3D n(it->second->normal_x, it->second->normal_y, it->second->normal_z);
          if(fabs(n[0])/n.norm()>1-1e-8)
            myrects.push_back(it->second); 
         }

      if(mydisks.size() || myrects.size()) {

        for(int k=kk0; k<kkmax; k++)
          for(int j=jj0; j<jjmax; j++) {

            if(global_mesh.OutsidePhysicalDomainAndUnpopulated(iimax-1,j,k))
              continue; //skip corner nodes

            for(int p=0; p<(int)mydisks.size(); p++) {
              if(IsPointInDisk(global_mesh.GetY(j), global_mesh.GetZ(k),
                               mydisks[p]->cen_y, mydisks[p]->cen_z, mydisks[p]->radius)){
                v[k][j][iimax-1][1] = mydisks[p]->state.velocity_x;
                v[k][j][iimax-1][2] = mydisks[p]->state.velocity_y;
                v[k][j][iimax-1][3] = mydisks[p]->state.velocity_z;
              }
            }

            for(int p=0; p<(int)myrects.size(); p++) {
              if(IsPointInRectangle(global_mesh.GetY(j), global_mesh.GetZ(k),
                                    myrects[p]->cen_y, myrects[p]->cen_z, myrects[p]->a, myrects[p]->b)){
                v[k][j][iimax-1][1] = myrects[p]->state.velocity_x;
                v[k][j][iimax-1][2] = myrects[p]->state.velocity_y;
                v[k][j][iimax-1][3] = myrects[p]->state.velocity_z;
              }
            }

          }
      }

    }
  }

  
  if(jj0==-1) { 
    if (iod.mesh.bc_y0 == MeshData::INLET    || iod.mesh.bc_y0 == MeshData::INLET2 ||
        iod.mesh.bc_y0 == MeshData::SLIPWALL || iod.mesh.bc_y0 == MeshData::STICKWALL) {

      vector<DiskData* > mydisks;
      for(auto it=disks.begin(); it!=disks.end(); it++)
        if(it->second->cen_y == iod.mesh.y0) {
          Vec3D n(it->second->normal_x, it->second->normal_y, it->second->normal_z);
          if(fabs(n[1])/n.norm()>1-1e-8)
            mydisks.push_back(it->second); 
        }
      vector<RectangleData* > myrects;
      for(auto it=rectangles.begin(); it!=rectangles.end(); it++)       
        if(it->second->cen_y == iod.mesh.y0) {
          Vec3D n(it->second->normal_x, it->second->normal_y, it->second->normal_z);
          if(fabs(n[1])/n.norm()>1-1e-8)
            myrects.push_back(it->second); 
        }
      if(mydisks.size() || myrects.size()) {

        for(int k=kk0; k<kkmax; k++)
          for(int i=ii0; i<iimax; i++) {

            if(global_mesh.OutsidePhysicalDomainAndUnpopulated(i,jj0,k))
              continue; //skip corner nodes

            for(int p=0; p<(int)mydisks.size(); p++) {
              if(IsPointInDisk(global_mesh.GetZ(k), global_mesh.GetX(i),
                               mydisks[p]->cen_z, mydisks[p]->cen_x, mydisks[p]->radius)){
                v[k][jj0][i][1]   = mydisks[p]->state.velocity_x;
                v[k][jj0+1][i][2] = mydisks[p]->state.velocity_y;
                v[k][jj0][i][3]   = mydisks[p]->state.velocity_z;
              }
            }

            for(int p=0; p<(int)myrects.size(); p++) {
              if(IsPointInRectangle(global_mesh.GetZ(k), global_mesh.GetX(i),
                                    myrects[p]->cen_z, myrects[p]->cen_x, myrects[p]->a, myrects[p]->b)){
                v[k][jj0][i][1]   = myrects[p]->state.velocity_x;
                v[k][jj0+1][i][2] = myrects[p]->state.velocity_y;
                v[k][jj0][i][3]   = myrects[p]->state.velocity_z;
              }
            }

          }
      }

    }
  }

  if(jjmax==NY+1) { 
    if (iod.mesh.bc_ymax == MeshData::INLET    || iod.mesh.bc_ymax == MeshData::INLET2 ||
        iod.mesh.bc_ymax == MeshData::SLIPWALL || iod.mesh.bc_ymax == MeshData::STICKWALL) {

      vector<DiskData* > mydisks;
      for(auto it=disks.begin(); it!=disks.end(); it++)
        if(it->second->cen_y == iod.mesh.ymax) {
          Vec3D n(it->second->normal_x, it->second->normal_y, it->second->normal_z);
          if(fabs(n[1])/n.norm()>1-1e-8)
            mydisks.push_back(it->second); 
        }
      vector<RectangleData* > myrects;
      for(auto it=rectangles.begin(); it!=rectangles.end(); it++)       
        if(it->second->cen_y == iod.mesh.ymax) {
          Vec3D n(it->second->normal_x, it->second->normal_y, it->second->normal_z);
          if(fabs(n[1])/n.norm()>1-1e-8)
            myrects.push_back(it->second); 
        }

      if(mydisks.size() || myrects.size()) {

        for(int k=kk0; k<kkmax; k++)
          for(int i=ii0; i<iimax; i++) {

            if(global_mesh.OutsidePhysicalDomainAndUnpopulated(i,jjmax-1,k))
              continue; //skip corner nodes

            for(int p=0; p<(int)mydisks.size(); p++) {
              if(IsPointInDisk(global_mesh.GetZ(k), global_mesh.GetX(i),
                               mydisks[p]->cen_z, mydisks[p]->cen_x, mydisks[p]->radius)){
                v[k][jjmax-1][i][1] = mydisks[p]->state.velocity_x;
                v[k][jjmax-1][i][2] = mydisks[p]->state.velocity_y;
                v[k][jjmax-1][i][3] = mydisks[p]->state.velocity_z;
              }
            }

            for(int p=0; p<(int)myrects.size(); p++) {
              if(IsPointInRectangle(global_mesh.GetZ(k), global_mesh.GetX(i),
                                    myrects[p]->cen_z, myrects[p]->cen_x, myrects[p]->a, myrects[p]->b)){
                v[k][jjmax-1][i][1] = myrects[p]->state.velocity_x;
                v[k][jjmax-1][i][2] = myrects[p]->state.velocity_y;
                v[k][jjmax-1][i][3] = myrects[p]->state.velocity_z;
              }
            }

          }
      }

    }
  }

  
  if(kk0==-1) { 
    if (iod.mesh.bc_z0 == MeshData::INLET    || iod.mesh.bc_z0 == MeshData::INLET2 ||
        iod.mesh.bc_z0 == MeshData::SLIPWALL || iod.mesh.bc_z0 == MeshData::STICKWALL) {

      vector<DiskData* > mydisks;
      for(auto it=disks.begin(); it!=disks.end(); it++)
        if(it->second->cen_z == iod.mesh.z0) {
          Vec3D n(it->second->normal_x, it->second->normal_y, it->second->normal_z);
          if(fabs(n[2])/n.norm()>1-1e-8)
            mydisks.push_back(it->second); 
        }
      vector<RectangleData* > myrects;
      for(auto it=rectangles.begin(); it!=rectangles.end(); it++)       
        if(it->second->cen_z == iod.mesh.z0) {
          Vec3D n(it->second->normal_x, it->second->normal_y, it->second->normal_z);
          if(fabs(n[2])/n.norm()>1-1e-8)
            myrects.push_back(it->second); 
        }

      if(mydisks.size() || myrects.size()) {

        for(int j=jj0; j<jjmax; j++)
          for(int i=ii0; i<iimax; i++) {

            if(global_mesh.OutsidePhysicalDomainAndUnpopulated(i,j,kk0))
              continue; //skip corner nodes

            for(int p=0; p<(int)mydisks.size(); p++) {
              if(IsPointInDisk(global_mesh.GetX(i), global_mesh.GetY(j),
                               mydisks[p]->cen_x, mydisks[p]->cen_y, mydisks[p]->radius)){
                v[kk0][j][i][1]   = mydisks[p]->state.velocity_x;
                v[kk0][j][i][2]   = mydisks[p]->state.velocity_y;
                v[kk0+1][j][i][3] = mydisks[p]->state.velocity_z;
              }
            }

            for(int p=0; p<(int)myrects.size(); p++) {
              if(IsPointInRectangle(global_mesh.GetX(i), global_mesh.GetY(j),
                                    myrects[p]->cen_x, myrects[p]->cen_y, myrects[p]->a, myrects[p]->b)){
                v[kk0][j][i][1]   = myrects[p]->state.velocity_x;
                v[kk0][j][i][2]   = myrects[p]->state.velocity_y;
                v[kk0+1][j][i][3] = myrects[p]->state.velocity_z;
              }
            }

          }
      }

    }
  }

  if(kkmax==NZ+1) { 
    if (iod.mesh.bc_zmax == MeshData::INLET    || iod.mesh.bc_zmax == MeshData::INLET2 ||
        iod.mesh.bc_zmax == MeshData::SLIPWALL || iod.mesh.bc_zmax == MeshData::STICKWALL) {

      vector<DiskData* > mydisks;
      for(auto it=disks.begin(); it!=disks.end(); it++)
        if(it->second->cen_z == iod.mesh.zmax) {
          Vec3D n(it->second->normal_x, it->second->normal_y, it->second->normal_z);
          if(fabs(n[2])/n.norm()>1-1e-8)
            mydisks.push_back(it->second); 
        }

      vector<RectangleData* > myrects;
      for(auto it=rectangles.begin(); it!=rectangles.end(); it++)       
        if(it->second->cen_z == iod.mesh.zmax) {
          Vec3D n(it->second->normal_x, it->second->normal_y, it->second->normal_z);
          if(fabs(n[2])/n.norm()>1-1e-8)
            myrects.push_back(it->second); 
        }

      if(mydisks.size() || myrects.size()) {

        for(int j=jj0; j<jjmax; j++)
          for(int i=ii0; i<iimax; i++) {

            if(global_mesh.OutsidePhysicalDomainAndUnpopulated(i,j,kkmax-1))
              continue; //skip corner nodes

            for(int p=0; p<(int)mydisks.size(); p++) {
              if(IsPointInDisk(global_mesh.GetX(i), global_mesh.GetY(j),
                               mydisks[p]->cen_x, mydisks[p]->cen_y, mydisks[p]->radius)){
                v[kkmax-1][j][i][1] = mydisks[p]->state.velocity_x;
                v[kkmax-1][j][i][2] = mydisks[p]->state.velocity_y;
                v[kkmax-1][j][i][3] = mydisks[p]->state.velocity_z;
              }
            }

            for(int p=0; p<(int)myrects.size(); p++) {
              if(IsPointInRectangle(global_mesh.GetX(i), global_mesh.GetY(j),
                                    myrects[p]->cen_x, myrects[p]->cen_y, myrects[p]->a, myrects[p]->b)){
                v[kkmax-1][j][i][1] = myrects[p]->state.velocity_x;
                v[kkmax-1][j][i][2] = myrects[p]->state.velocity_y;
                v[kkmax-1][j][i][3] = myrects[p]->state.velocity_z;
              }
            }

          }
      }

    }
  }

}

//--------------------------------------------------------------------------

void
IncompressibleOperator::ComputeTimeStepSize(SpaceVariable3D &V, SpaceVariable3D &ID,
                                            double &dt, double &cfl, SpaceVariable3D *LocalDt)
{

  if(LocalDt) { //local time-stepping, dealt with separately.
    assert(iod.ts.timestep<=0.0); //shouldn't have constant time-step size.
    ComputeLocalTimeStepSizes(V,ID,dt,cfl,*LocalDt);
    return;
  }

  GlobalMeshInfo &global_mesh(spo.GetGlobalMeshInfo());

  Vec5D***   v = (Vec5D***)V.GetDataPointer();
  double*** id = ID.GetDataPointer();

  double vel_over_dx_max = 0.0;
  double dx, dy, dz;
  for(int k=k0; k<kmax; k++) {
    dz = global_mesh.GetDz(k);
    for(int j=j0; j<jmax; j++) {
      dy = global_mesh.GetDy(j);
      for(int i=i0; i<imax; i++) {
        if(id[k][j][i] == INACTIVE_MATERIAL_ID)
          continue;
        dx = global_mesh.GetDx(i);
        if(global_mesh.x_glob.size()>1) vel_over_dx_max = max(vel_over_dx_max, fabs(v[k][j][i][1])/dx);
        if(global_mesh.y_glob.size()>1) vel_over_dx_max = max(vel_over_dx_max, fabs(v[k][j][i][2])/dy);
        if(global_mesh.z_glob.size()>1) vel_over_dx_max = max(vel_over_dx_max, fabs(v[k][j][i][3])/dz);
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &vel_over_dx_max, 1, MPI_DOUBLE, MPI_MAX, comm);

  if(vel_over_dx_max == 0.0) {//velocity is zero everywhere in the domain ==> check ghost layers
    for(int k=kk0; k<kkmax; k++) {
      dz = global_mesh.GetDz(k);
      for(int j=jj0; j<jjmax; j++) {
        dy = global_mesh.GetDy(j);
        for(int i=ii0; i<iimax; i++) {
          if(!global_mesh.OutsidePhysicalDomain(i,j,k) || global_mesh.OutsidePhysicalDomainAndUnpopulated(i,j,k)) 
            continue;
          // do not check id == INACTIVE_MATERIAL_ID (ghost cells may not have valid ID)
          dx = global_mesh.GetDx(i);
          if(global_mesh.x_glob.size()>1) vel_over_dx_max = max(vel_over_dx_max, fabs(v[k][j][i][1])/dx);
          if(global_mesh.y_glob.size()>1) vel_over_dx_max = max(vel_over_dx_max, fabs(v[k][j][i][2])/dy);
          if(global_mesh.z_glob.size()>1) vel_over_dx_max = max(vel_over_dx_max, fabs(v[k][j][i][3])/dz);
        }
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &vel_over_dx_max, 1, MPI_DOUBLE, MPI_MAX, comm);
  if(vel_over_dx_max == 0.0) //still zero... give it an arbitrary value
    vel_over_dx_max = 1.0e8;
 

  V.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();

  if(iod.ts.timestep > 0) {
    dt = iod.ts.timestep;
    cfl = dt*vel_over_dx_max;
  } else {//apply the CFL number
    cfl = iod.ts.cfl;
    dt = cfl/vel_over_dx_max;
  }

  // Check surface tension (written by Wentao)
  if(iod.exact_riemann.surface_tension == ExactRiemannSolverData::YES) {
    double dt_surfaceTension = 0.9*spo.ComputeTimeStepSizeSurfaceTension(V, ID); // 0.9 is a safety factor
    if(dt>dt_surfaceTension) {
      dt = dt_surfaceTension;
      cfl = dt*vel_over_dx_max;
    }
  }

}

//--------------------------------------------------------------------------

void
IncompressibleOperator::ComputeLocalTimeStepSizes(SpaceVariable3D &V, SpaceVariable3D &ID, double &dt, 
                                                  double &cfl, SpaceVariable3D &LocalDt)
{

  cfl = iod.ts.cfl;

  GlobalMeshInfo &global_mesh(spo.GetGlobalMeshInfo());

  Vec5D*** v    = (Vec5D***)V.GetDataPointer();
  double*** id  = ID.GetDataPointer();
  double*** dtl = LocalDt.GetDataPointer();


  // Step 1: Get a global dt (in case some cells have zero velocity --- see below)
  double vel_over_dx_max = 0.0;
  double dx, dy, dz;
  for(int k=k0; k<kmax; k++) {
    dz = global_mesh.GetDz(k);
    for(int j=j0; j<jmax; j++) {
      dy = global_mesh.GetDy(j);
      for(int i=i0; i<imax; i++) {
        if(id[k][j][i] == INACTIVE_MATERIAL_ID)
          continue;
        dx = global_mesh.GetDx(i);
        if(global_mesh.x_glob.size()>1) vel_over_dx_max = max(vel_over_dx_max, fabs(v[k][j][i][1])/dx);
        if(global_mesh.y_glob.size()>1) vel_over_dx_max = max(vel_over_dx_max, fabs(v[k][j][i][2])/dy);
        if(global_mesh.z_glob.size()>1) vel_over_dx_max = max(vel_over_dx_max, fabs(v[k][j][i][3])/dz);
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &vel_over_dx_max, 1, MPI_DOUBLE, MPI_MAX, comm);

  if(vel_over_dx_max == 0.0) {//velocity is zero everywhere in the domain ==> check ghost layers
    for(int k=kk0; k<kkmax; k++) {
      dz = global_mesh.GetDz(k);
      for(int j=jj0; j<jjmax; j++) {
        dy = global_mesh.GetDy(j);
        for(int i=ii0; i<iimax; i++) {
          if(!global_mesh.OutsidePhysicalDomain(i,j,k) || global_mesh.OutsidePhysicalDomainAndUnpopulated(i,j,k)) 
            continue;
          // do not check id == INACTIVE_MATERIAL_ID (ghost cells may not have valid ID)
          dx = global_mesh.GetDx(i);
          if(global_mesh.x_glob.size()>1) vel_over_dx_max = max(vel_over_dx_max, fabs(v[k][j][i][1])/dx);
          if(global_mesh.y_glob.size()>1) vel_over_dx_max = max(vel_over_dx_max, fabs(v[k][j][i][2])/dy);
          if(global_mesh.z_glob.size()>1) vel_over_dx_max = max(vel_over_dx_max, fabs(v[k][j][i][3])/dz);
        }
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &vel_over_dx_max, 1, MPI_DOUBLE, MPI_MAX, comm);
  if(vel_over_dx_max == 0.0) //still zero... give it an arbitrary value
    vel_over_dx_max = 1.0e8;
 
  // global min dt
  dt = cfl/vel_over_dx_max;


  //Step 2: Determine local dt
  double vel_over_dx_max_loc;
  // Loop through the real domain (excluding the ghost layer)
  for(int k=k0; k<kmax; k++) {
    dz = global_mesh.GetDz(k);
    for(int j=j0; j<jmax; j++) {
      dy = global_mesh.GetDy(j);
      for(int i=i0; i<imax; i++) {
        if(id[k][j][i] == INACTIVE_MATERIAL_ID)
          continue;
        dx = global_mesh.GetDx(i);
        
        vel_over_dx_max_loc = 0.0;
        if(global_mesh.x_glob.size()>1) vel_over_dx_max_loc = max(vel_over_dx_max_loc, fabs(v[k][j][i][1])/dx);
        if(global_mesh.y_glob.size()>1) vel_over_dx_max_loc = max(vel_over_dx_max_loc, fabs(v[k][j][i][2])/dy);
        if(global_mesh.z_glob.size()>1) vel_over_dx_max_loc = max(vel_over_dx_max_loc, fabs(v[k][j][i][3])/dz);

        // calculates local time-step size
        if(vel_over_dx_max_loc>0) {
          if(vel_over_dx_max_loc<0.0025*vel_over_dx_max)
            dtl[k][j][i] = cfl/vel_over_dx_max_loc;
          else
            dtl[k][j][i] = cfl/vel_over_dx_max*400.0; //at most, 400 times than dt_min
        } else
          dtl[k][j][i] = dt;
      }
    }
  }

  V.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();

  LocalDt.RestoreDataPointerAndInsert();

}

//--------------------------------------------------------------------------

void
IncompressibleOperator::BuildVelocityEquationSIMPLE(int dir, Vec5D*** v0, Vec5D*** v, double*** id,
                                                    double*** vturb, //turbulent working term
                                                    double*** homo, //node is in a homogeneous region?
                                                    vector<RowEntries> &vlin_rows, SpaceVariable3D &B,
                                                    SpaceVariable3D &Ddiag, bool SIMPLEC, double Efactor,
                                                    double dt, SpaceVariable3D *LocalDt)
{
  assert(dir==0 || dir==1 || dir==2);

  //local utility function for evaluating dynamic turbulence eddy viscosity (mu_T)
  auto GetMut = [&] (int i, int j, int k) {
    return GetDynamicEddyViscosity(v[k][j][i][0], Mu[id[k][j][i]], vturb[k][j][i]);
  };


  double*** dtloc = NULL;
  if(LocalDt) { //local time-stepping, dealt with separately.
    assert(iod.ts.timestep<=0.0); //shouldn't have constant time-step size.
    dtloc = LocalDt->GetDataPointer();
  }

  GlobalMeshInfo& global_mesh(spo.GetGlobalMeshInfo());

  double*** bb = B.GetDataPointer();
  double*** diag = Ddiag.GetDataPointer();

  vlin_rows.clear(); //clear existing data (to be safe)
  int row_counter = 0;

  double dx, dy, dz, dxl, dxr, dyb, dyt, dzk, dzf, dxdy, dydz, dxdz;
  double cm(1.0), cp(1.0), cm_plus_cp(0.0), rhou1, rhou2;
  double a, ap, ap0, F, D, mu, mu1, mu2, anb;
  double nut, nut1, nut2; //turbulent eddy viscosity (kinematic)
  double rho, rho1, rho2;

  for(int k=k0; k<kmax; k++) {
    dz  = Dz[dir][k-k0];
    dzk = dz_k[dir][k-k0];
    dzf = dz_f[dir][k-k0];
    for(int j=j0; j<jmax; j++) {
      dy   = Dy[dir][j-j0];
      dyb  = dy_b[dir][j-j0];
      dyt  = dy_t[dir][j-j0];
      dydz = dy*dz;
      for(int i=i0; i<imax; i++) {

        // insert the row
        vlin_rows.push_back(RowEntries(7)); // at most 7 non-zero entries on each row
        RowEntries &row(vlin_rows[row_counter]);
        row.SetRow(i,j,k);
        row_counter++;

        if((dir==0 && i==0) || (dir==1 && j==0) || (dir==2 && k==0)) {
          row.PushEntry(i,j,k, 1.0); 
          bb[k][j][i] = v[k][j][i][1+dir];
          diag[k][j][i] = 0.0; // not really used
          continue;
        }

        dx   = Dx[dir][i-i0];
        dxl  = dx_l[dir][i-i0];
        dxr  = dx_r[dir][i-i0];
        dxdy = dx*dy;
        dxdz = dx*dz;


        //---------------------------------------------------
        // Reference: Patankar's book, Eqs. (5.61) - (5.64)
        //---------------------------------------------------

        ap = 0.0; //diagonal
        bb[k][j][i] = 0.0; //rhs

        //-----------
        // LEFT
        //-----------
        // Calculate F at the correct location (different for dir=0,1,2)
        if(dir==0)
          F = v[k][j][i-1][0]*0.5*(v[k][j][i-1][1] + v[k][j][i][1]);
        else {
          cm = global_mesh.GetDx(i-1);
          cp = global_mesh.GetDx(i);
          cm_plus_cp = cm + cp;
          if(homo[k][j][i]) { //this "cell" (i,j,k) is in a neighborhood w/ constant rho and mu
            rhou1 = dir==1 ? v[k][j-1][i][1]*v[k][j][i][0] : v[k-1][j][i][1]*v[k][j][i][0];
            rhou2 = v[k][j][i][1]*v[k][j][i][0];
          } else {
            rhou1 = dir==1 ? v[k][j-1][i][1]*(cp*v[k][j-1][i-1][0] + cm*v[k][j-1][i][0])/cm_plus_cp
                           : v[k-1][j][i][1]*(cp*v[k-1][j][i-1][0] + cm*v[k-1][j][i][0])/cm_plus_cp;
            rhou2 = v[k][j][i][1]*(cp*v[k][j][i-1][0] + cm*v[k][j][i][0])/(cm+cp);
          }
          F = dir== 1 ? (dyt*rhou1 + dyb*rhou2)/(dyb+dyt) : (dzf*rhou1 + dzk*rhou2)/(dzk+dzf);
        }
        F *= dydz;
        // Calculate D and the "a" coefficient
        a = std::max(F, 0.0);
        if(dir==0) {
          mu = Mu[id[k][j][i-1]];
          if(vturb)
            mu += GetDynamicEddyViscosity(v[k][j][i-1][0], mu, vturb[k][j][i-1]);
        } else {
          // cm and cp have been calculated!   
          if(i==0) {
            mu1 = dir==1 ? Mu[id[k][j-1][i]] : Mu[id[k-1][j][i]];
            mu2 = Mu[id[k][j][i]];
            if(vturb) {
              rho1 = dir==1 ? v[k][j-1][i][0] : v[k-1][j][i][0];
              rho2 = v[k][j][i][0];
              nut1 = dir==1 ? vturb[k][j-1][i] : vturb[k-1][j][i];
              nut2 = vturb[k][j][i];
            }
          } else {
            mu1 = dir==1 ? (cp*Mu[id[k][j-1][i-1]] + cm*Mu[id[k][j-1][i]])/cm_plus_cp
                         : (cp*Mu[id[k-1][j][i-1]] + cm*Mu[id[k-1][j][i]])/cm_plus_cp;
            mu2 = (cp*Mu[id[k][j][i-1]] + cm*Mu[id[k][j][i]])/cm_plus_cp;
            if(vturb) {
              rho1 = dir==1 ? (cp*v[k][j-1][i-1][0] + cm*v[k][j-1][i][0])/cm_plus_cp :
                              (cp*v[k-1][j][i-1][0] + cm*v[k-1][j][i][0])/cm_plus_cp;
              rho2 = (cp*v[k][j][i-1][0] + cm*v[k][j][i][0])/cm_plus_cp;
              nut1 = dir==1 ? (cp*vturb[k][j-1][i-1] + cm*vturb[k][j-1][i])/cm_plus_cp :
                              (cp*vturb[k-1][j][i-1] + cm*vturb[k-1][j][i])/cm_plus_cp;
              nut2 = (cp*vturb[k][j][i-1] + cm*vturb[k][j][i])/cm_plus_cp;
            }
          }
          mu = dir==1 ? (dyt*mu1 + dyb*mu2)/(dyb+dyt) : (dzf*mu1 + dzk*mu2)/(dzk+dzf);
          if(vturb) {
            rho = dir==1 ? (dyt*rho1 + dyb*rho2)/(dyb+dyt) : (dzf*rho1 + dzk*rho2)/(dzk+dzf);
            nut = dir==1 ? (dyt*nut1 + dyb*nut2)/(dyb+dyt) : (dzf*nut1 + dzk*nut2)/(dzk+dzf);
            mu += GetDynamicEddyViscosity(rho, mu, nut);
          }
        }
        if(mu>0.0) {
          D  = mu*dydz/dxl;
          a += D*PowerLaw(F/D);
        }
        ap += a;

        // Add entry (or add to bb or ap, if i-1 is outside boundary)
        if(i-1>=0)
          row.PushEntry(i-1,j,k, -a);  //on the left hand side
        else { //i-1 is outside domain boundary (if it gets here, dir must be y or z (1 or 2))
          if(iod.mesh.bc_x0 == MeshData::INLET || iod.mesh.bc_x0 == MeshData::INLET2 ||
             iod.mesh.bc_x0 == MeshData::OVERSET)
            bb[k][j][i] += a*v[k][j][i-1][dir+1]; //+a*v or +a*w to the RHS 
          else if(iod.mesh.bc_x0 == MeshData::SLIPWALL || iod.mesh.bc_x0 == MeshData::SYMMETRY ||
                  iod.mesh.bc_x0 == MeshData::OUTLET)
            ap -= a;
          else if(iod.mesh.bc_x0 == MeshData::STICKWALL)
            ap += a; //bb[k][j][i] -= a*v[k][j][i][dir+1]; //v[k][j][i-1][dir+1] should be -v[k][j][i][dir+1]
                                                //to have zero velocity on the wall
          else {
            fprintf(stdout,"*** Error: Detected unknown boundary condition type (%d).\n",
                    (int)iod.mesh.bc_x0);
            exit(-1);
          }
        }


        //-----------
        // RIGHT 
        //-----------
        // Calculate F at the correct location (different for dir=0,1,2)
        if(dir==0)
          F = v[k][j][i][0]*0.5*(v[k][j][i][1] + v[k][j][i+1][1]);
        else {
          cm = global_mesh.GetDx(i);
          cp = global_mesh.GetDx(i+1);
          cm_plus_cp = cm + cp;
          if(homo[k][j][i]) {
            rhou1 = dir==1 ? v[k][j-1][i+1][1]*v[k][j][i][0] : v[k-1][j][i+1][1]*v[k][j][i][0];
            rhou2 = v[k][j][i+1][1]*v[k][j][i][0];
          } else {
            rhou1 = dir==1 ? v[k][j-1][i+1][1]*(cp*v[k][j-1][i][0] + cm*v[k][j-1][i+1][0])/cm_plus_cp
                           : v[k-1][j][i+1][1]*(cp*v[k-1][j][i][0] + cm*v[k-1][j][i+1][0])/cm_plus_cp;
            rhou2 = v[k][j][i+1][1]*(cp*v[k][j][i][0] + cm*v[k][j][i+1][0])/(cm+cp);
          }
          F = dir== 1 ? (dyt*rhou1 + dyb*rhou2)/(dyb+dyt) : (dzf*rhou1 + dzk*rhou2)/(dzk+dzf);
        }
        F *= dydz;
        // Calculate D and the "a" coefficient
        a = std::max(-F, 0.0);
        if(dir==0) {
          mu = Mu[id[k][j][i]];
          if(vturb)
            mu += GetDynamicEddyViscosity(v[k][j][i][0], mu, vturb[k][j][i]);
        } else {
          // cm and cp have been calculated!   
          if(i==NX-1) {
            mu1 = dir==1 ? Mu[id[k][j-1][i]] : Mu[id[k-1][j][i]];
            mu2 = Mu[id[k][j][i]];
            if(vturb) {
              rho1 = dir==1 ? v[k][j-1][i][0] : v[k-1][j][i][0];
              rho2 = v[k][j][i][0];
              nut1 = dir==1 ? vturb[k][j-1][i] : vturb[k-1][j][i];
              nut2 = vturb[k][j][i];
            }
          } else {
            mu1 = dir==1 ? (cp*Mu[id[k][j-1][i]] + cm*Mu[id[k][j-1][i+1]])/cm_plus_cp
                         : (cp*Mu[id[k-1][j][i]] + cm*Mu[id[k-1][j][i+1]])/cm_plus_cp;
            mu2 = (cp*Mu[id[k][j][i]] + cm*Mu[id[k][j][i+1]])/cm_plus_cp;
            if(vturb) {
              rho1 = dir==1 ? (cp*v[k][j-1][i][0] + cm*v[k][j-1][i+1][0])/cm_plus_cp :
                              (cp*v[k-1][j][i][0] + cm*v[k-1][j][i+1][0])/cm_plus_cp;
              rho2 = (cp*v[k][j][i][0] + cm*v[k][j][i+1][0])/cm_plus_cp;
              nut1 = dir==1 ? (cp*vturb[k][j-1][i] + cm*vturb[k][j-1][i+1])/cm_plus_cp :
                              (cp*vturb[k-1][j][i] + cm*vturb[k-1][j][i+1])/cm_plus_cp;
              nut2 = (cp*vturb[k][j][i] + cm*vturb[k][j][i+1])/cm_plus_cp;
            }
          }
          mu = dir==1 ? (dyt*mu1 + dyb*mu2)/(dyb+dyt) : (dzf*mu1 + dzk*mu2)/(dzk+dzf);
          if(vturb) {
            rho = dir==1 ? (dyt*rho1 + dyb*rho2)/(dyb+dyt) : (dzf*rho1 + dzk*rho2)/(dzk+dzf);
            nut = dir==1 ? (dyt*nut1 + dyb*nut2)/(dyb+dyt) : (dzf*nut1 + dzk*nut2)/(dzk+dzf);
            mu += GetDynamicEddyViscosity(rho, mu, nut);
          }
        }
        if(mu>0.0) {
          D  = mu*dydz/dxr;
          a += D*PowerLaw(F/D);
        }
        ap += a;
        // Add entry (or add to bb or ap, if i+1 is outside boundary)
        if(i+1<NX)
          row.PushEntry(i+1,j,k, -a);  //on the left hand side
        else { //i+1 is outside domain boundary
          if(iod.mesh.bc_xmax == MeshData::INLET || iod.mesh.bc_xmax == MeshData::INLET2 ||
             iod.mesh.bc_xmax == MeshData::OVERSET)
            bb[k][j][i] += a*v[k][j][i+1][dir+1];
          else if(iod.mesh.bc_xmax == MeshData::SLIPWALL || iod.mesh.bc_xmax == MeshData::SYMMETRY) {
            if(dir != 0) //otherwise, v[k][j][i+1][1] = 0
              ap -= a;
          } else if(iod.mesh.bc_xmax == MeshData::OUTLET)
            ap -= a;
          else if(iod.mesh.bc_xmax == MeshData::STICKWALL) {
            if(dir != 0) //otherwise, v[k][j][i+1][1] should be 0 as it is on the wall
              ap += a; //bb[k][j][i] -= a*v[k][j][i][dir+1]; //v[k][j][i+1][dir+1] should be -v[k][j][i][dir+1]
          } else {
            fprintf(stdout,"*** Error: Detected unknown boundary condition type (%d).\n",
                    (int)iod.mesh.bc_xmax);
            exit(-1);
          }
        }

    
        //-----------
        // BOTTOM 
        //-----------
        // Calculate F at the correct location (different for dir=0,1,2)
        if(dir==1)
          F = v[k][j-1][i][0]*0.5*(v[k][j-1][i][2] + v[k][j][i][2]);
        else {
          cm = global_mesh.GetDy(j-1);
          cp = global_mesh.GetDy(j);
          cm_plus_cp = cm + cp;
          if(homo[k][j][i]) {
            rhou1 = dir==0 ? v[k][j][i-1][2]*v[k][j][i][0] : v[k-1][j][i][2]*v[k][j][i][0];
            rhou2 = v[k][j][i][2]*v[k][j][i][0];
          } else {
            rhou1 = dir==0 ? v[k][j][i-1][2]*(cp*v[k][j-1][i-1][0] + cm*v[k][j][i-1][0])/cm_plus_cp
                           : v[k-1][j][i][2]*(cp*v[k-1][j-1][i][0] + cm*v[k-1][j][i][0])/cm_plus_cp;
            rhou2 = v[k][j][i][2]*(cp*v[k][j-1][i][0] + cm*v[k][j][i][0])/(cm+cp);
          }
          F = dir== 0 ? (dxr*rhou1 + dxl*rhou2)/(dxl+dxr) : (dzf*rhou1 + dzk*rhou2)/(dzk+dzf);
        }
        F *= dxdz;
        // Calculate D and the "a" coefficient
        a = std::max(F, 0.0);
        if(dir==1) {
          mu = Mu[id[k][j-1][i]];
          if(vturb)
            mu += GetDynamicEddyViscosity(v[k][j-1][i][0], mu, vturb[k][j-1][i]);
        } else {
          // cm and cp have been calculated!   
          if(j==0) {
            mu1 = dir==0 ? Mu[id[k][j][i-1]] : Mu[id[k-1][j][i]];
            mu2 = Mu[id[k][j][i]];
            if(vturb) {
              rho1 = dir==0 ? v[k][j][i-1][0] : v[k-1][j][i][0];
              rho2 = v[k][j][i][0];
              nut1 = dir==0 ? vturb[k][j][i-1] : vturb[k-1][j][i];
              nut2 = vturb[k][j][i];
            }
          } else {
            mu1 = dir==0 ? (cp*Mu[id[k][j-1][i-1]] + cm*Mu[id[k][j][i-1]])/cm_plus_cp
                         : (cp*Mu[id[k-1][j-1][i]] + cm*Mu[id[k-1][j][i]])/cm_plus_cp;
            mu2 = (cp*Mu[id[k][j-1][i]] + cm*Mu[id[k][j][i]])/cm_plus_cp;
            if(vturb) {
              rho1 = dir==0 ? (cp*v[k][j-1][i-1][0] + cm*v[k][j][i-1][0])/cm_plus_cp :
                              (cp*v[k-1][j-1][i][0] + cm*v[k-1][j][i][0])/cm_plus_cp;
              rho2 = (cp*v[k][j-1][i][0] + cm*v[k][j][i][0])/cm_plus_cp;
              nut1 = dir==0 ? (cp*vturb[k][j-1][i-1] + cm*vturb[k][j][i-1])/cm_plus_cp :
                              (cp*vturb[k-1][j-1][i] + cm*vturb[k-1][j][i])/cm_plus_cp;
              nut2 = (cp*vturb[k][j-1][i] + cm*vturb[k][j][i])/cm_plus_cp;
            }
          }
          mu = dir==0 ? (dxr*mu1 + dxl*mu2)/(dxl+dxr) : (dzf*mu1 + dzk*mu2)/(dzk+dzf);
          if(vturb) {
            rho = dir==0 ? (dxr*rho1 + dxl*rho2)/(dxl+dxr) : (dzf*rho1 + dzk*rho2)/(dzk+dzf);
            nut = dir==0 ? (dxr*nut1 + dxl*nut2)/(dxl+dxr) : (dzf*nut1 + dzk*nut2)/(dzk+dzf);
            mu += GetDynamicEddyViscosity(rho, mu, nut);
          }
        }
        if(mu>0.0) {
          D  = mu*dxdz/dyb;
          a += D*PowerLaw(F/D);
        }
        ap += a;
        // Add entry (or add to bb or ap, if j-1 is outside boundary)
        if(j-1>=0)
          row.PushEntry(i,j-1,k, -a);  //on the left hand side
        else { //j-1 is outside domain boundary (if it gets here, dir must be x or z (0 or 2))
          // special treatment for flat plate
          if(iod.rans.example == RANSTurbulenceModelData::FLAT_PLATE &&
             ((dir == 0 && global_mesh.GetX(i) - 0.5*global_mesh.GetDx(i) >= iod.rans.example_param_1) ||
              (dir == 2 && global_mesh.GetX(i) >= iod.rans.example_param_1))) { //stick wall
            ap += a; //bb[k][j][i] -= a*v[k][j][i][dir+1]; //v[k][j-1][i][dir+1] should be -v[k][j][i][dir+1]
                                                //to have zero velocity on the wall
          }
          else if(iod.mesh.bc_y0 == MeshData::INLET || iod.mesh.bc_y0 == MeshData::INLET2 ||
             iod.mesh.bc_y0 == MeshData::OVERSET)
            bb[k][j][i] += a*v[k][j-1][i][dir+1]; //+a*u or +a*w to the RHS 
          else if(iod.mesh.bc_y0 == MeshData::SLIPWALL || iod.mesh.bc_y0 == MeshData::SYMMETRY ||
                  iod.mesh.bc_y0 == MeshData::OUTLET)
            ap -= a;
          else if(iod.mesh.bc_y0 == MeshData::STICKWALL)
            ap += a; //bb[k][j][i] -= a*v[k][j][i][dir+1]; //v[k][j-1][i][dir+1] should be -v[k][j][i][dir+1]
                                                //to have zero velocity on the wall
          else {
            fprintf(stdout,"*** Error: Detected unknown boundary condition type (%d).\n",
                    (int)iod.mesh.bc_y0);
            exit(-1);
          }
        }


        //-----------
        // TOP 
        //-----------
        // Calculate F at the correct location (different for dir=0,1,2)
        if(dir==1)
          F = v[k][j][i][0]*0.5*(v[k][j][i][2] + v[k][j+1][i][2]);
        else {
          cm = global_mesh.GetDy(j);
          cp = global_mesh.GetDy(j+1);
          cm_plus_cp = cm + cp;
          if(homo[k][j][i]) {
            rhou1 = dir==0 ? v[k][j+1][i-1][2]*v[k][j][i][0] : v[k-1][j+1][i][2]*v[k][j][i][0];
            rhou2 = v[k][j+1][i][2]*v[k][j][i][0];
          } else {
            rhou1 = dir==0 ? v[k][j+1][i-1][2]*(cp*v[k][j][i-1][0] + cm*v[k][j+1][i-1][0])/cm_plus_cp
                           : v[k-1][j+1][i][2]*(cp*v[k-1][j][i][0] + cm*v[k-1][j+1][i][0])/cm_plus_cp;
            rhou2 = v[k][j+1][i][2]*(cp*v[k][j][i][0] + cm*v[k][j+1][i][0])/(cm+cp);
          }
          F = dir== 0 ? (dxr*rhou1 + dxl*rhou2)/(dxl+dxr) : (dzf*rhou1 + dzk*rhou2)/(dzk+dzf);
        }
        F *= dxdz;
        // Calculate D and the "a" coefficient
        a = std::max(-F, 0.0);
        if(dir==1) {
          mu = Mu[id[k][j][i]];
          if(vturb)
            mu += GetDynamicEddyViscosity(v[k][j][i][0], mu, vturb[k][j][i]);
        } else {
          // cm and cp have been calculated!   
          if(j==NY-1) {
            mu1 = dir==0 ? Mu[id[k][j][i-1]] : Mu[id[k-1][j][i]];
            mu2 = Mu[id[k][j][i]];
            if(vturb) {
              rho1 = dir==0 ? v[k][j][i-1][0] : v[k-1][j][i][0];
              rho2 = v[k][j][i][0];
              nut1 = dir==0 ? vturb[k][j][i-1] : vturb[k-1][j][i];
              nut2 = vturb[k][j][i];
            }
          } else {
            mu1 = dir==0 ? (cp*Mu[id[k][j][i-1]] + cm*Mu[id[k][j+1][i-1]])/cm_plus_cp
                         : (cp*Mu[id[k-1][j][i]] + cm*Mu[id[k-1][j+1][i]])/cm_plus_cp;
            mu2 = (cp*Mu[id[k][j][i]] + cm*Mu[id[k][j+1][i]])/cm_plus_cp;
            if(vturb) {
              rho1 = dir==0 ? (cp*v[k][j][i-1][0] + cm*v[k][j+1][i-1][0])/cm_plus_cp :
                              (cp*v[k-1][j][i][0] + cm*v[k-1][j+1][i][0])/cm_plus_cp;
              rho2 = (cp*v[k][j][i][0] + cm*v[k][j+1][i][0])/cm_plus_cp;
              nut1 = dir==0 ? (cp*vturb[k][j][i-1] + cm*vturb[k][j+1][i-1])/cm_plus_cp :
                              (cp*vturb[k-1][j][i] + cm*vturb[k-1][j+1][i])/cm_plus_cp;
              nut2 = (cp*vturb[k][j][i] + cm*vturb[k][j+1][i])/cm_plus_cp;
            }
          }
          mu = dir==0 ? (dxr*mu1 + dxl*mu2)/(dxl+dxr) : (dzf*mu1 + dzk*mu2)/(dzk+dzf);
          if(vturb) {
            rho = dir==0 ? (dxr*rho1 + dxl*rho2)/(dxl+dxr) : (dzf*rho1 + dzk*rho2)/(dzk+dzf);
            nut = dir==0 ? (dxr*nut1 + dxl*nut2)/(dxl+dxr) : (dzf*nut1 + dzk*nut2)/(dzk+dzf);
            mu += GetDynamicEddyViscosity(rho, mu, nut);
          }
        }
        if(mu>0.0) {
          D  = mu*dxdz/dyt;
          a += D*PowerLaw(F/D);
        }
        ap += a;
        // Add entry (or add to bb or ap, if j+1 is outside boundary)
        if(j+1<NY)
          row.PushEntry(i,j+1,k, -a);  //on the left hand side
        else { //j+1 is outside domain boundary
          if(iod.mesh.bc_ymax == MeshData::INLET || iod.mesh.bc_ymax == MeshData::INLET2 ||
             iod.mesh.bc_ymax == MeshData::OVERSET)
            bb[k][j][i] += a*v[k][j+1][i][dir+1];
          else if(iod.mesh.bc_ymax == MeshData::SLIPWALL || iod.mesh.bc_ymax == MeshData::SYMMETRY) {
            if(dir != 1) //otherwise, v[k][j+1][i][2] = 0
              ap -= a;
          } else if(iod.mesh.bc_ymax == MeshData::OUTLET)
              ap -= a;
          else if(iod.mesh.bc_ymax == MeshData::STICKWALL) {
            if(dir != 1) //otherwise, v[k][j+1][i][2] should be 0 as it is on the wall
              ap += a; //bb[k][j][i] -= a*v[k][j][i][dir+1]; //v[k][j+1][i][dir+1] should be -v[k][j][i][dir+1]
          } else {
            fprintf(stdout,"*** Error: Detected unknown boundary condition type (%d).\n",
                    (int)iod.mesh.bc_ymax);
            exit(-1);
          }
        }
   
    
        //-----------
        // BACK
        //-----------
        // Calculate F at the correct location (different for dir=0,1,2)
        if(dir==2)
          F = v[k-1][j][i][0]*0.5*(v[k-1][j][i][3] + v[k][j][i][3]);
        else {
          cm = global_mesh.GetDz(k-1);
          cp = global_mesh.GetDz(k);
          cm_plus_cp = cm + cp;
          if(homo[k][j][i]) {
            rhou1 = dir==0 ? v[k][j][i-1][3]*v[k][j][i][0] : v[k][j-1][i][3]*v[k][j][i][0];
            rhou2 = v[k][j][i][3]*v[k][j][i][0];
          } else {
            rhou1 = dir==0 ? v[k][j][i-1][3]*(cp*v[k-1][j][i-1][0] + cm*v[k][j][i-1][0])/cm_plus_cp
                           : v[k][j-1][i][3]*(cp*v[k-1][j-1][i][0] + cm*v[k][j-1][i][0])/cm_plus_cp;
            rhou2 = v[k][j][i][3]*(cp*v[k-1][j][i][0] + cm*v[k][j][i][0])/(cm+cp);
          }
          F = dir==0 ? (dxr*rhou1 + dxl*rhou2)/(dxl+dxr) : (dyt*rhou1 + dyb*rhou2)/(dyb+dyt);
        }
        F *= dxdy;
        // Calculate D and the "a" coefficient
        a = std::max(F, 0.0);
        if(dir==2) {
          mu = Mu[id[k-1][j][i]];
          if(vturb)
            mu += GetDynamicEddyViscosity(v[k-1][j][i][0], mu, vturb[k-1][j][i]);
        } else {
          // cm and cp have been calculated!   
          if(k==0) {
            mu1 = dir==0 ? Mu[id[k][j][i-1]] : Mu[id[k][j-1][i]];
            mu2 = Mu[id[k][j][i]];
            if(vturb) {
              rho1 = dir==0 ? v[k][j][i-1][0] : v[k][j-1][i][0];
              rho2 = v[k][j][i][0];
              nut1 = dir==0 ? vturb[k][j][i-1] : vturb[k][j-1][i];
              nut2 = vturb[k][j][i];
            }
          } else {
            mu1 = dir==0 ? (cp*Mu[id[k-1][j][i-1]] + cm*Mu[id[k][j][i-1]])/cm_plus_cp
                         : (cp*Mu[id[k-1][j-1][i]] + cm*Mu[id[k][j-1][i]])/cm_plus_cp;
            mu2 = (cp*Mu[id[k-1][j][i]] + cm*Mu[id[k][j][i]])/cm_plus_cp;
            if(vturb) {
              rho1 = dir==0 ? (cp*v[k-1][j][i-1][0] + cm*v[k][j][i-1][0])/cm_plus_cp :
                              (cp*v[k-1][j-1][i][0] + cm*v[k][j-1][i][0])/cm_plus_cp;
              rho2 = (cp*v[k-1][j][i][0] + cm*v[k][j][i][0])/cm_plus_cp;
              nut1 = dir==0 ? (cp*vturb[k-1][j][i-1] + cm*vturb[k][j][i-1])/cm_plus_cp :
                              (cp*vturb[k-1][j-1][i] + cm*vturb[k][j-1][i])/cm_plus_cp;
              nut2 = (cp*vturb[k-1][j][i] + cm*vturb[k][j][i])/cm_plus_cp;
            }
          }
          mu = dir==0 ? (dxr*mu1 + dxl*mu2)/(dxl+dxr) : (dyt*mu1 + dyb*mu2)/(dyb+dyt);
          if(vturb) {
            rho = dir==0 ? (dxr*rho1 + dxl*rho2)/(dxl+dxr) : (dyt*rho1 + dyb*rho2)/(dyb+dyt);
            nut = dir==0 ? (dxr*nut1 + dxl*nut2)/(dxl+dxr) : (dyt*nut1 + dyb*nut2)/(dyb+dyt);
            mu += GetDynamicEddyViscosity(rho, mu, nut);
          }
        }
        if(mu>0.0) {
          D  = mu*dxdy/dzk;
          a += D*PowerLaw(F/D);
        }
        ap += a;
        // Add entry (or add to bb or ap, if k-1 is outside boundary)
        if(k-1>=0)
          row.PushEntry(i,j,k-1, -a);  //on the left hand side
        else { //k-1 is outside domain boundary (if it gets here, dir must be x or y (0 or 1))
          if(iod.mesh.bc_z0 == MeshData::INLET || iod.mesh.bc_z0 == MeshData::INLET2 ||
             iod.mesh.bc_z0 == MeshData::OVERSET)
            bb[k][j][i] += a*v[k-1][j][i][dir+1]; //+a*u or +a*v to the RHS 
          else if(iod.mesh.bc_z0 == MeshData::SLIPWALL || iod.mesh.bc_z0 == MeshData::SYMMETRY ||
                  iod.mesh.bc_z0 == MeshData::OUTLET)
            ap -= a;
          else if(iod.mesh.bc_z0 == MeshData::STICKWALL)
            ap += a; //bb[k][j][i] -= a*v[k][j][i][dir+1]; //v[k-1][j][i][dir+1] should be -v[k][j][i][dir+1]
                                                //to have zero velocity on the wall
          else {
            fprintf(stdout,"*** Error: Detected unknown boundary condition type (%d).\n",
                    (int)iod.mesh.bc_z0);
            exit(-1);
          }
        }


              

        //-----------
        // FRONT 
        //-----------
        // Calculate F at the correct location (different for dir=0,1,2)
        if(dir==2)
          F = v[k][j][i][0]*0.5*(v[k][j][i][3] + v[k+1][j][i][3]);
        else {
          cm = global_mesh.GetDz(k);
          cp = global_mesh.GetDz(k+1);
          cm_plus_cp = cm + cp;
          if(homo[k][j][i]) {
            rhou1 = dir==0 ? v[k+1][j][i-1][3]*v[k][j][i][0] : v[k+1][j-1][i][3]*v[k][j][i][0];
            rhou2 = v[k+1][j][i][3]*v[k][j][i][0];
          } else {
            rhou1 = dir==0 ? v[k+1][j][i-1][3]*(cp*v[k][j][i-1][0] + cm*v[k+1][j][i-1][0])/cm_plus_cp
                           : v[k+1][j-1][i][3]*(cp*v[k][j-1][i][0] + cm*v[k+1][j-1][i][0])/cm_plus_cp;
            rhou2 = v[k+1][j][i][3]*(cp*v[k][j][i][0] + cm*v[k+1][j][i][0])/(cm+cp);
            F = dir== 0 ? (dxr*rhou1 + dxl*rhou2)/(dxl+dxr) : (dyt*rhou1 + dyb*rhou2)/(dyb+dyt);
          }
        }
        F *= dxdy;
        // Calculate D and the "a" coefficient
        a = std::max(-F, 0.0);
        if(dir==2) {
          mu = Mu[id[k][j][i]];
          if(vturb)
            mu += GetDynamicEddyViscosity(v[k][j][i][0], mu, vturb[k][j][i]);
        } else {
          // cm and cp have been calculated!   
          if(k==NZ-1) {
            mu1 = dir==0 ? Mu[id[k][j][i-1]] : Mu[id[k][j-1][i]];
            mu2 = Mu[id[k][j][i]];
            if(vturb) {
              rho1 = dir==0 ? v[k][j][i-1][0] : v[k][j-1][i][0];
              rho2 = v[k][j][i][0];
              nut1 = dir==0 ? vturb[k][j][i-1] : vturb[k][j-1][i];
              nut2 = vturb[k][j][i];
            }
          } else {
            mu1 = dir==0 ? (cp*Mu[id[k][j][i-1]] + cm*Mu[id[k+1][j][i-1]])/cm_plus_cp
                         : (cp*Mu[id[k][j-1][i]] + cm*Mu[id[k+1][j-1][i]])/cm_plus_cp;
            mu2 = (cp*Mu[id[k][j][i]] + cm*Mu[id[k+1][j][i]])/cm_plus_cp;
            if(vturb) {
              rho1 = dir==0 ? (cp*v[k][j][i-1][0] + cm*v[k+1][j][i-1][0])/cm_plus_cp :
                              (cp*v[k][j-1][i][0] + cm*v[k+1][j-1][i][0])/cm_plus_cp;
              rho2 = (cp*v[k][j][i][0] + cm*v[k+1][j][i][0])/cm_plus_cp;
              nut1 = dir==0 ? (cp*vturb[k][j][i-1] + cm*vturb[k+1][j][i-1])/cm_plus_cp :
                              (cp*vturb[k][j-1][i] + cm*vturb[k+1][j-1][i])/cm_plus_cp;
              nut2 = (cp*vturb[k][j][i] + cm*vturb[k+1][j][i])/cm_plus_cp;
            }
          }
          mu = dir==0 ? (dxr*mu1 + dxl*mu2)/(dxl+dxr) : (dyt*mu1 + dyb*mu2)/(dyb+dyt);
          if(vturb) {
            rho = dir==0 ? (dxr*rho1 + dxl*rho2)/(dxl+dxr) : (dyt*rho1 + dyb*rho2)/(dyt+dyb);
            nut = dir==0 ? (dxr*nut1 + dxl*nut2)/(dxl+dxr) : (dyt*nut1 + dyb*nut2)/(dyt+dyb);
            mu += GetDynamicEddyViscosity(rho, mu, nut);
          }
        }
        if(mu>0.0) {
          D  = mu*dxdy/dzf;
          a += D*PowerLaw(F/D);
        }
        ap += a;
        // Add entry (or add to bb or ap, if j+1 is outside boundary)
        if(k+1<NZ)
          row.PushEntry(i,j,k+1, -a);  //on the left hand side
        else { //k+1 is outside domain boundary
          if(iod.mesh.bc_zmax == MeshData::INLET || iod.mesh.bc_zmax == MeshData::INLET2 ||
             iod.mesh.bc_zmax == MeshData::OVERSET)
            bb[k][j][i] += a*v[k+1][j][i][dir+1];
          else if(iod.mesh.bc_zmax == MeshData::SLIPWALL || iod.mesh.bc_zmax == MeshData::SYMMETRY) {
            if(dir != 2) //otherwise, v[k+1][j][i][3] = 0
              ap -= a;
          } else if(iod.mesh.bc_zmax == MeshData::OUTLET)
            ap -= a;
          else if(iod.mesh.bc_zmax == MeshData::STICKWALL) {
            if(dir != 2) //otherwise, v[k+1][j][i][3] should be 0 as it is on the wall
              ap += a; //bb[k][j][i] -= a*v[k][j][i][dir+1]; //v[k+1][j][i][dir+1] should be -v[k][j][i][dir+1]
          } else {
            fprintf(stdout,"*** Error: Detected unknown boundary condition type (%d).\n",
                    (int)iod.mesh.bc_zmax);
            exit(-1);
          }
        }
   
   
        //------------------------------------------------------
        // Calculate and add the diagonal entry and the RHS (b)
        // Ref: Eqs. (5.62) and (6.8) in Patankar's book
        //------------------------------------------------------
        anb = SIMPLEC ? ap : 0.0; //needed later
        ap0 = dir==0 ? (dxr*v0[k][j][i-1][0] + dxl*v0[k][j][i][0])/(dxl+dxr) :
              dir==1 ? (dyt*v0[k][j-1][i][0] + dyb*v0[k][j][i][0])/(dyb+dyt) :
                       (dzf*v0[k-1][j][i][0] + dzk*v0[k][j][i][0])/(dzk+dzf);
        ap0 *= LocalDt ? dxdy*dz/dtloc[k][j][i] : dxdy*dz/dt;
        ap += ap0; //!< -Sp*dx*dy*dz, for source terms

        bb[k][j][i] += ap0*v0[k][j][i][dir+1]; 


        //------------------------------------------------------
        // Evaluating the source term due to turbulent eddy viscosity
        //------------------------------------------------------
        if(vturb) {
          double Sc_turb = 0.0;
          if (dir == 0) { //X-momentum equation case (i>=1)

            //Velocity gradients
            double dudx = -dxr/(dxl*(dxl+dxr))*v[k][j][i-1][1] + (dxr-dxl)/(dxl*dxr)*v[k][j][i][1]
                        +  dxl/(dxr*(dxl+dxr))*v[k][j][i+1][1]; //2nd-order accuracy, see Kevin's notes
            double v_left  = (v[k][j+1][i-1][2] + v[k][j][i-1][2])/2.0;
            double v_right = (v[k][j+1][i][2]   + v[k][j][i][2])/2.0;
            double dvdx    = (v_right - v_left)/dx; //TODO: only 1st order for non-unif grids (also dwdx, dir=1,2)
            double w_left  = (v[k+1][j][i-1][3] + v[k][j][i-1][3])/2.0;
            double w_right = (v[k+1][j][i][3]   + v[k][j][i][3])/2.0;
            double dwdx    = (w_right - w_left)/dx;

            //mu_t gradients
            double dmutdx = (GetMut(i,j,k) - GetMut(i-1,j,k))/dx;
            double mut[3];
            for(int p=-1; p<2; p++)
              mut[p+1] = (dxr*GetMut(i-1,j+p,k) + dxl*GetMut(i,j+p,k))/(dxl+dxr);
            double dmutdy = -dyt/(dyb*(dyt+dyb))*mut[0] + (dyt-dyb)/(dyt*dyb)*mut[1] + dyb/(dyt*(dyt+dyb))*mut[2];
            for(int p=-1; p<2; p++)
              mut[p+1] = (dxr*GetMut(i-1,j,k+p) + dxl*GetMut(i,j,k+p))/(dxl+dxr);
            double dmutdz = -dzf/(dzk*(dzf+dzk))*mut[0] + (dzf-dzk)/(dzf*dzk)*mut[1] + dzk/(dzf*(dzf+dzk))*mut[2];

            Sc_turb = dudx*dmutdx + dvdx*dmutdy + dwdx*dmutdz;

          }
          else if (dir == 1) { //Y-Momentum equation case (j>=1)

            //Velocity gradients
            double u_btm = (v[k][j-1][i+1][1] + v[k][j-1][i][1])/2.0;
            double u_top = (v[k][j][i+1][1]   + v[k][j][i][1])/2.0;
            double dudy  = (u_top - u_btm)/dy;
            double dvdy  = -dyt/(dyb*(dyt+dyb))*v[k][j-1][i][2] + (dyt-dyb)/(dyt*dyb)*v[k][j][i][2]
                         +  dyb/(dyt*(dyt+dyb))*v[k][j+1][i][2]; //2nd-order accuracy, see Kevin's notes
            double w_btm = (v[k+1][j-1][i][3] + v[k][j-1][i][3])/2.0;
            double w_top = (v[k+1][j][i][3]   + v[k][j][i][3])/2.0;
            double dwdy  = (w_top - w_btm)/dy;

            //mu_t gradients
            double mut[3];
            for(int p=-1; p<2; p++)
              mut[p+1] = (dyt*GetMut(i+p,j-1,k) + dyb*GetMut(i+p,j,k))/(dyt+dyb);
            double dmutdx = -dxr/(dxl*(dxl+dxr))*mut[0] + (dxr-dxl)/(dxl*dxr)*mut[1] + dxl/(dxr*(dxl+dxr))*mut[2];
            double dmutdy = (GetMut(i,j,k) - GetMut(i,j-1,k))/dy;
            for(int p=-1; p<2; p++)
              mut[p+1] = (dyt*GetMut(i,j-1,k+p) + dyb*GetMut(i,j,k+p))/(dyt+dyb);
            double dmutdz = -dzf/(dzk*(dzf+dzk))*mut[0] + (dzf-dzk)/(dzf*dzk)*mut[1] + dzk/(dzf*(dzf+dzk))*mut[2];

            Sc_turb = dudy*dmutdx + dvdy*dmutdy + dwdy*dmutdz;

          }
          else if (dir == 2) {//Z-Momentum equation case (k>=1)

            //Velocity gradients
            double u_bk = (v[k-1][j][i+1][1] + v[k-1][j][i][1])/2.0;
            double u_ft = (v[k][j][i+1][1]   + v[k][j][i][1])/2.0;
            double dudz = (u_ft - u_bk)/dz;
            double v_bk = (v[k-1][j+1][i][2] + v[k-1][j][i][2])/2.0;
            double v_ft = (v[k][j+1][i][2] + v[k][j][i][2])/2.0;
            double dvdz = (v_ft - v_bk)/dz;
            double dwdz = -dzf/(dzk*(dzf+dzk))*v[k-1][j][i][3] + (dzf-dzk)/(dzf*dzk)*v[k][j][i][3]
                        +  dzk/(dzf*(dzf+dzk))*v[k+1][j][i][3];

            //mu_t gradients
            double mut[3];
            for(int p=-1; p<2; p++)
              mut[p+1] = (dzf*GetMut(i+p,j,k-1) + dzk*GetMut(i+p,j,k))/(dzf+dzk);
            double dmutdx = -dxr/(dxl*(dxl+dxr))*mut[0] + (dxr-dxl)/(dxl*dxr)*mut[1] + dxl/(dxr*(dxl+dxr))*mut[2];
            for(int p=-1; p<2; p++)
              mut[p+1] = (dzf*GetMut(i,j+p,k-1) + dzk*GetMut(i,j+p,k))/(dzf+dzk);
            double dmutdy = -dyt/(dyb*(dyt+dyb))*mut[0] + (dyt-dyb)/(dyt*dyb)*mut[1] + dyb/(dyt*(dyt+dyb))*mut[2];
            double dmutdz = (GetMut(i,j,k) - GetMut(i,j,k-1))/dz;
   
            Sc_turb = dudz*dmutdx + dvdz*dmutdy + dwdz*dmutdz;

          }

          bb[k][j][i] += v[k][j][i][0]*Sc_turb*dx*dy*dz; //updates bb
        }
        //------------------------------------------------------
        // END of evaluating the source term due to turbulent eddy viscosity
        //------------------------------------------------------


        bb[k][j][i] += dir==0 ? (v[k][j][i-1][4] - v[k][j][i][4])*dydz :
                       dir==1 ? (v[k][j-1][i][4] - v[k][j][i][4])*dxdz :
                                (v[k-1][j][i][4] - v[k][j][i][4])*dxdy;

        // Apply relaxation (Ref: Eq.(6) of Van Doormaal and Rathby, 1984)
        assert(Efactor>0.0);
        bb[k][j][i] += ap*v0[k][j][i][dir+1]/Efactor;

        ap *= 1.0 + 1.0/Efactor; 
        row.PushEntry(i,j,k, ap);

        // Store diagonal for use in pressure correction equation
        assert(ap-anb!=0.0);
        diag[k][j][i] = 1.0/(ap-anb);
        diag[k][j][i] *= dir==0 ? dydz : dir==1 ? dxdz : dxdy;

      }
    }
  }


  B.RestoreDataPointerAndInsert();
  Ddiag.RestoreDataPointerAndInsert();

  if(LocalDt)
    LocalDt->RestoreDataPointerToLocalVector();
}

//--------------------------------------------------------------------------

void
IncompressibleOperator::BuildPressureEquationSIMPLE(Vec5D*** v, double*** homo, SpaceVariable3D &VXstar,
                                                    SpaceVariable3D &VYstar, SpaceVariable3D &VZstar,
                                                    SpaceVariable3D &DX, SpaceVariable3D &DY, SpaceVariable3D &DZ,
                                                    vector<RowEntries> &plin_rows, SpaceVariable3D &B,
                                                    Int3 *ijk_zero_p)
{

  GlobalMeshInfo& global_mesh(spo.GetGlobalMeshInfo());

  double*** bb = B.GetDataPointer();
  double*** diagx = DX.GetDataPointer();
  double*** diagy = DY.GetDataPointer();
  double*** diagz = DZ.GetDataPointer();
  double*** ustar = VXstar.GetDataPointer();
  double*** vstar = VYstar.GetDataPointer();
  double*** wstar = VZstar.GetDataPointer();

  int row_counter = 0;
  plin_rows.clear(); //clear existing data

  double ap, a, rho, dx, dxl, dxr, dy, dyb, dyt, dz, dzk, dzf;

  for(int k=k0; k<kmax; k++) {
    dz  = global_mesh.GetDz(k);
    dzk = global_mesh.GetDz(k-1);
    dzf = global_mesh.GetDz(k+1);
    for(int j=j0; j<jmax; j++) {
      dy  = global_mesh.GetDy(j);
      dyb = global_mesh.GetDy(j-1);
      dyt = global_mesh.GetDy(j+1);
      for(int i=i0; i<imax; i++) {
        dx  = global_mesh.GetDx(i);
        dxl = global_mesh.GetDx(i-1);
        dxr = global_mesh.GetDx(i+1);

        // locate the row
        plin_rows.push_back(RowEntries(7)); // at most 7 non-zero entries on each row
        RowEntries &row(plin_rows[row_counter]);
        row.SetRow(i,j,k);
        row_counter++;

        // initialization
        ap = 0.0;
        bb[k][j][i] = 0.0;

        // set p = 0? (Otherwise, the linear system is singular, but *may* still be solvable by iterative methods)
        if(ijk_zero_p && i==(*ijk_zero_p)[0] && j==(*ijk_zero_p)[1] && k==(*ijk_zero_p)[2]) {
          row.PushEntry(i,j,k, 1.0);
          bb[k][j][i] = 0.0;
          continue;
        }

        //-------
        // LEFT
        //-------
        if(i>0) {
          rho = homo[k][j][i] ? v[k][j][i][0] : (dx*v[k][j][i-1][0] + dxl*v[k][j][i][0])/(dxl+dx);
          a = rho*diagx[k][j][i]/dx;
          if(!ijk_zero_p || *ijk_zero_p != Int3(i-1,j,k))
            row.PushEntry(i-1,j,k, -a);  //on the left hand side
          ap += a;
          bb[k][j][i] += rho*ustar[k][j][i]/dx;
        } else {
          // TODO: Only consider velocity b.c. for now (velocity is known)
          rho = v[k][j][i][0];
          bb[k][j][i] += rho*v[k][j][i][1]/dx; //boundary velocity
        }

        //-------
        // RIGHT 
        //-------
        if(i<NX-1) {
          rho = homo[k][j][i] ? v[k][j][i][0] : (dxr*v[k][j][i][0] + dx*v[k][j][i+1][0])/(dx+dxr);
          a = rho*diagx[k][j][i+1]/dx;
          if(!ijk_zero_p || *ijk_zero_p != Int3(i+1,j,k))
            row.PushEntry(i+1,j,k, -a);  //on the left hand side
          ap += a;
          bb[k][j][i] -= rho*ustar[k][j][i+1]/dx;
        } else {
          // TODO: Only consider velocity b.c. for now (velocity is known)
          rho = v[k][j][i][0];
          bb[k][j][i] -= rho*v[k][j][i+1][1]/dx; //boundary velocity
        }
 
        //-------
        // BOTTOM 
        //-------
        if(j>0) {
          rho = homo[k][j][i] ? v[k][j][i][0] : (dy*v[k][j-1][i][0] + dyb*v[k][j][i][0])/(dyb+dy);
          a = rho*diagy[k][j][i]/dy;
          if(!ijk_zero_p || *ijk_zero_p != Int3(i,j-1,k))
            row.PushEntry(i,j-1,k, -a);  //on the left hand side
          ap += a;
          bb[k][j][i] += rho*vstar[k][j][i]/dy;
        } else {
          // TODO: Only consider velocity b.c. for now (velocity is known)
          rho = v[k][j][i][0];
          bb[k][j][i] += rho*v[k][j][i][2]/dy; //boundary velocity
        }

        //-------
        // TOP 
        //-------
        if(j<NY-1) {
          rho = homo[k][j][i] ? v[k][j][i][0] : (dyt*v[k][j][i][0] + dy*v[k][j+1][i][0])/(dy+dyt);
          a = rho*diagy[k][j+1][i]/dy;
          if(!ijk_zero_p || *ijk_zero_p != Int3(i,j+1,k))
            row.PushEntry(i,j+1,k, -a);  //on the left hand side
          ap += a;
          bb[k][j][i] -= rho*vstar[k][j+1][i]/dy;
        } else {
          // TODO: Only consider velocity b.c. for now (velocity is known)
          rho = v[k][j][i][0];
          bb[k][j][i] -= rho*v[k][j+1][i][2]/dy; //boundary velocity
        }
  
        //-------
        // BACK 
        //-------
        if(k>0) {
          rho = homo[k][j][i] ? v[k][j][i][0] : (dz*v[k-1][j][i][0] + dzk*v[k][j][i][0])/(dzk+dz);
          a = rho*diagz[k][j][i]/dz;
          if(!ijk_zero_p || *ijk_zero_p != Int3(i,j,k-1))
            row.PushEntry(i,j,k-1, -a);  //on the left hand side
          ap += a;
          bb[k][j][i] += rho*wstar[k][j][i]/dz;
        } else {
          // TODO: Only consider velocity b.c. for now (velocity is known)
          rho = v[k][j][i][0];
          bb[k][j][i] += rho*v[k][j][i][3]/dz; //boundary velocity
        }

        //-------
        // FRONT 
        //-------
        if(k<NZ-1) {
          rho = homo[k][j][i] ? v[k][j][i][0] : (dzf*v[k][j][i][0] + dz*v[k+1][j][i][0])/(dz+dzf);
          a = rho*diagz[k+1][j][i]/dz;
          if(!ijk_zero_p || *ijk_zero_p != Int3(i,j,k+1))
            row.PushEntry(i,j,k+1, -a);  //on the left hand side
          ap += a;
          bb[k][j][i] -= rho*wstar[k+1][j][i]/dz;
        } else {
          // TODO: Only consider velocity b.c. for now (velocity is known)
          rho = v[k][j][i][0];
          bb[k][j][i] -= rho*v[k+1][j][i][3]/dz; //boundary velocity
        }
 
        //Push the diagonal
        row.PushEntry(i,j,k, ap);
      }
    }
  }


  B.RestoreDataPointerAndInsert();
  DX.RestoreDataPointerToLocalVector();
  DY.RestoreDataPointerToLocalVector();
  DZ.RestoreDataPointerToLocalVector();
  VXstar.RestoreDataPointerToLocalVector();
  VYstar.RestoreDataPointerToLocalVector();
  VZstar.RestoreDataPointerToLocalVector();

}

//--------------------------------------------------------------------------

void
IncompressibleOperator::CalculateCoefficientsSIMPLER(int dir, Vec5D*** v0, Vec5D*** v, double*** id,
                                                     double*** homo, vector<RowEntries> &vlin_rows,
                                                     SpaceVariable3D &Bv,
                                                     SpaceVariable3D &Vhat, SpaceVariable3D &Ddiag,
                                                     double Efactor, double dt, SpaceVariable3D *LocalDt)
{

  // very similar to BuildVelocityEquationSIMPLE

  assert(dir==0 || dir==1 || dir==2);

  double*** dtloc = NULL;
  if(LocalDt) { //local time-stepping, dealt with separately.
    assert(iod.ts.timestep<=0.0); //shouldn't have constant time-step size.
    dtloc = LocalDt->GetDataPointer();
  }

  GlobalMeshInfo& global_mesh(spo.GetGlobalMeshInfo());

  double*** vhat = Vhat.GetDataPointer();
  double*** diag = Ddiag.GetDataPointer();
  double*** bb   = Bv.GetDataPointer();

  int row_counter = 0;

  double dx, dy, dz, dxl, dxr, dyb, dyt, dzk, dzf, dxdy, dydz, dxdz;
  double cm(1.0), cp(1.0), cm_plus_cp(0.0), rhou1, rhou2;
  double a, ap, ap0, F, D, mu, mu1, mu2, vhat_denom;

  for(int k=k0; k<kmax; k++) {
    dz  = Dz[dir][k-k0];
    dzk = dz_k[dir][k-k0];
    dzf = dz_f[dir][k-k0];
    for(int j=j0; j<jmax; j++) {
      dy   = Dy[dir][j-j0];
      dyb  = dy_b[dir][j-j0];
      dyt  = dy_t[dir][j-j0];
      dydz = dy*dz;
      for(int i=i0; i<imax; i++) {

        // create the row
        vlin_rows.push_back(RowEntries(7)); // at most 7 non-zero entries on each row
        RowEntries &row(vlin_rows[row_counter]);
        row.SetRow(i,j,k);
        row_counter++;

        if((dir==0 && i==0) || (dir==1 && j==0) || (dir==2 && k==0)) {
          row.PushEntry(i,j,k, 1.0);
          bb[k][j][i] = v[k][j][i][1+dir];
          diag[k][j][i] = 0.0; // not really used
          vhat[k][j][i] = v[k][j][i][dir+1];
          continue;
        }

        dx   = Dx[dir][i-i0];
        dxl  = dx_l[dir][i-i0];
        dxr  = dx_r[dir][i-i0];
        dxdy = dx*dy;
        dxdz = dx*dz;


        //---------------------------------------------------
        // Reference: Patankar's book, Eqs. (5.61) - (5.64)
        //---------------------------------------------------

        ap = 0.0; //diagonal
        vhat[k][j][i] = 0.0;
        vhat_denom = 0.0; //denominator of vhat, usually same as ap, except near some boundaries
        bb[k][j][i] = 0.0;

        //-----------
        // LEFT
        //-----------
        // Calculate F at the correct location (different for dir=0,1,2)
        if(dir==0)
          F = v[k][j][i-1][0]*0.5*(v[k][j][i-1][1] + v[k][j][i][1]);
        else {
          if(homo[k][j][i]) { //this "cell" (i,j,k) is in a neighborhood w/ constant rho and mu
            rhou1 = dir==1 ? v[k][j-1][i][1]*v[k][j][i][0] : v[k-1][j][i][1]*v[k][j][i][0];
            rhou2 = v[k][j][i][1]*v[k][j][i][0];
          } else {
            cm = global_mesh.GetDx(i-1);
            cp = global_mesh.GetDx(i);
            cm_plus_cp = cm + cp;
            rhou1 = dir==1 ? v[k][j-1][i][1]*(cp*v[k][j-1][i-1][0] + cm*v[k][j-1][i][0])/cm_plus_cp
                           : v[k-1][j][i][1]*(cp*v[k-1][j][i-1][0] + cm*v[k-1][j][i][0])/cm_plus_cp;
            rhou2 = v[k][j][i][1]*(cp*v[k][j][i-1][0] + cm*v[k][j][i][0])/(cm+cp);
          }
          F = dir== 1 ? (dyt*rhou1 + dyb*rhou2)/(dyb+dyt) : (dzf*rhou1 + dzk*rhou2)/(dzk+dzf);
        }
        F *= dydz;
        // Calculate D and the "a" coefficient
        a = std::max(F, 0.0);
        if(dir==0)
          mu = Mu[id[k][j][i-1]];
        else {
          if(homo[k][j][i])
            mu = Mu[id[k][j][i]];
          else {
            // cm and cp have been calculated!   
            if(i==0) {
              mu1 = dir==1 ? Mu[id[k][j-1][i]] : Mu[id[k-1][j][i]];
              mu2 = Mu[id[k][j][i]];
            } else {
              mu1 = dir==1 ? (cp*Mu[id[k][j-1][i-1]] + cm*Mu[id[k][j-1][i]])/cm_plus_cp
                           : (cp*Mu[id[k-1][j][i-1]] + cm*Mu[id[k-1][j][i]])/cm_plus_cp;
              mu2 = (cp*Mu[id[k][j][i-1]] + cm*Mu[id[k][j][i]])/cm_plus_cp;
            }
            mu = dir==1 ? (dyt*mu1 + dyb*mu2)/(dyb+dyt) : (dzf*mu1 + dzk*mu2)/(dzk+dzf);
          }
        }
        if(mu>0.0) {
          D  = mu*dydz/dxl;
          a += D*PowerLaw(F/D);
        }
        ap += a;
        vhat[k][j][i] += a*v[k][j][i-1][dir+1];
        vhat_denom += a;

        // Add entry (or add to bb or ap, if i-1 is outside boundary)
        if(i-1>=0)
          row.PushEntry(i-1,j,k, -a);  //on the left hand side
        else { //i-1 is outside domain boundary (if it gets here, dir must be y or z (1 or 2))
          if(iod.mesh.bc_x0 == MeshData::INLET || iod.mesh.bc_x0 == MeshData::INLET2 ||
             iod.mesh.bc_x0 == MeshData::OVERSET)
            bb[k][j][i] += a*v[k][j][i-1][dir+1]; //+a*v or +a*w to the RHS
          else if(iod.mesh.bc_x0 == MeshData::SLIPWALL || iod.mesh.bc_x0 == MeshData::SYMMETRY ||
                  iod.mesh.bc_x0 == MeshData::OUTLET)
            ap -= a; // vhat_denom remains the same
          else if(iod.mesh.bc_x0 == MeshData::STICKWALL)
            ap += a; //bb[k][j][i] -= a*v[k][j][i][dir+1]; //v[k][j][i-1][dir+1] should be -v[k][j][i][dir+1]
                                                //to have zero velocity on the wall
          else {
            fprintf(stdout,"*** Error: Detected unknown boundary condition type (%d).\n",
                    (int)iod.mesh.bc_x0);
            exit(-1);
          }
        } 

        //-----------
        // RIGHT 
        //-----------
        // Calculate F at the correct location (different for dir=0,1,2)
        if(dir==0)
          F = v[k][j][i][0]*0.5*(v[k][j][i][1] + v[k][j][i+1][1]);
        else {
          if(homo[k][j][i]) {
            rhou1 = dir==1 ? v[k][j-1][i+1][1]*v[k][j][i][0] : v[k-1][j][i+1][1]*v[k][j][i][0];
            rhou2 = v[k][j][i+1][1]*v[k][j][i][0];
          } else {
            cm = global_mesh.GetDx(i);
            cp = global_mesh.GetDx(i+1);
            cm_plus_cp = cm + cp;
            rhou1 = dir==1 ? v[k][j-1][i+1][1]*(cp*v[k][j-1][i][0] + cm*v[k][j-1][i+1][0])/cm_plus_cp
                           : v[k-1][j][i+1][1]*(cp*v[k-1][j][i][0] + cm*v[k-1][j][i+1][0])/cm_plus_cp;
            rhou2 = v[k][j][i+1][1]*(cp*v[k][j][i][0] + cm*v[k][j][i+1][0])/(cm+cp);
          }
          F = dir== 1 ? (dyt*rhou1 + dyb*rhou2)/(dyb+dyt) : (dzf*rhou1 + dzk*rhou2)/(dzk+dzf);
        }
        F *= dydz;
        // Calculate D and the "a" coefficient
        a = std::max(-F, 0.0);
        if(dir==0)
          mu = Mu[id[k][j][i]];
        else {
          if(homo[k][j][i])
            mu = Mu[id[k][j][i]];
          else {
            // cm and cp have been calculated!   
            if(i==NX-1) {
              mu1 = dir==1 ? Mu[id[k][j-1][i]] : Mu[id[k-1][j][i]];
              mu2 = Mu[id[k][j][i]];
            } else {
              mu1 = dir==1 ? (cp*Mu[id[k][j-1][i]] + cm*Mu[id[k][j-1][i+1]])/cm_plus_cp
                           : (cp*Mu[id[k-1][j][i]] + cm*Mu[id[k-1][j][i+1]])/cm_plus_cp;
              mu2 = (cp*Mu[id[k][j][i]] + cm*Mu[id[k][j][i+1]])/cm_plus_cp;
            }
            mu = dir==1 ? (dyt*mu1 + dyb*mu2)/(dyb+dyt) : (dzf*mu1 + dzk*mu2)/(dzk+dzf);
          }
        }
        if(mu>0.0) {
          D  = mu*dydz/dxr;
          a += D*PowerLaw(F/D);
        }
        ap += a;
        vhat[k][j][i] += a*v[k][j][i+1][dir+1];
        vhat_denom += a;

        // Add entry (or add to bb or ap, if i+1 is outside boundary)
        if(i+1<NX)
          row.PushEntry(i+1,j,k, -a);  //on the left hand side
        else { //i+1 is outside domain boundary
          if(iod.mesh.bc_xmax == MeshData::INLET || iod.mesh.bc_xmax == MeshData::INLET2 ||
             iod.mesh.bc_xmax == MeshData::OVERSET)
            bb[k][j][i] += a*v[k][j][i+1][dir+1];
          else if(iod.mesh.bc_xmax == MeshData::SLIPWALL || iod.mesh.bc_xmax == MeshData::SYMMETRY) {
            if(dir != 0) //otherwise, v[k][j][i+1][1] = 0
              ap -= a;
          } else if(iod.mesh.bc_xmax == MeshData::OUTLET) {
            ap -= a;
          } else if(iod.mesh.bc_xmax == MeshData::STICKWALL) {
            if(dir != 0) //otherwise, v[k][j][i+1][1] should be 0 as it is on the wall
              ap += a; //bb[k][j][i] -= a*v[k][j][i][dir+1]; //v[k][j][i+1][dir+1] should be -v[k][j][i][dir+1]
          } else {
            fprintf(stdout,"*** Error: Detected unknown boundary condition type (%d).\n",
                    (int)iod.mesh.bc_xmax);
            exit(-1);
          }
        }
    

        //-----------
        // BOTTOM 
        //-----------
        // Calculate F at the correct location (different for dir=0,1,2)
        if(dir==1)
          F = v[k][j-1][i][0]*0.5*(v[k][j-1][i][2] + v[k][j][i][2]);
        else {
          if(homo[k][j][i]) {
            rhou1 = dir==0 ? v[k][j][i-1][2]*v[k][j][i][0] : v[k-1][j][i][2]*v[k][j][i][0];
            rhou2 = v[k][j][i][2]*v[k][j][i][0];
          } else {
            cm = global_mesh.GetDy(j-1);
            cp = global_mesh.GetDy(j);
            cm_plus_cp = cm + cp;
            rhou1 = dir==0 ? v[k][j][i-1][2]*(cp*v[k][j-1][i-1][0] + cm*v[k][j][i-1][0])/cm_plus_cp
                           : v[k-1][j][i][2]*(cp*v[k-1][j-1][i][0] + cm*v[k-1][j][i][0])/cm_plus_cp;
            rhou2 = v[k][j][i][2]*(cp*v[k][j-1][i][0] + cm*v[k][j][i][0])/(cm+cp);
          }
          F = dir== 0 ? (dxr*rhou1 + dxl*rhou2)/(dxl+dxr) : (dzf*rhou1 + dzk*rhou2)/(dzk+dzf);
        }
        F *= dxdz;
        // Calculate D and the "a" coefficient
        a = std::max(F, 0.0);
        if(dir==1)
          mu = Mu[id[k][j-1][i]];
        else {
          if(homo[k][j][i])
            mu = Mu[id[k][j][i]]; 
          else {
            // cm and cp have been calculated!   
            if(j==0) {
              mu1 = dir==0 ? Mu[id[k][j][i-1]] : Mu[id[k-1][j][i]];
              mu2 = Mu[id[k][j][i]];
            } else {
              mu1 = dir==0 ? (cp*Mu[id[k][j-1][i-1]] + cm*Mu[id[k][j][i-1]])/cm_plus_cp
                           : (cp*Mu[id[k-1][j-1][i]] + cm*Mu[id[k-1][j][i]])/cm_plus_cp;
              mu2 = (cp*Mu[id[k][j-1][i]] + cm*Mu[id[k][j][i]])/cm_plus_cp;
            }
            mu = dir==0 ? (dxr*mu1 + dxl*mu2)/(dxl+dxr) : (dzf*mu1 + dzk*mu2)/(dzk+dzf);
          }
        }
        if(mu>0.0) {
          D  = mu*dxdz/dyb;
          a += D*PowerLaw(F/D);
        }
        ap += a;
        vhat[k][j][i] += a*v[k][j-1][i][dir+1];
        vhat_denom += a;

        // Add entry (or add to bb or ap, if j-1 is outside boundary)
        if(j-1>=0)
          row.PushEntry(i,j-1,k, -a);  //on the left hand side
        else { //j-1 is outside domain boundary (if it gets here, dir must be x or z (0 or 2))
          if(iod.mesh.bc_y0 == MeshData::INLET || iod.mesh.bc_y0 == MeshData::INLET2 ||
             iod.mesh.bc_y0 == MeshData::OVERSET)
            bb[k][j][i] += a*v[k][j-1][i][dir+1]; //+a*u or +a*w to the RHS
          else if(iod.mesh.bc_y0 == MeshData::SLIPWALL || iod.mesh.bc_y0 == MeshData::SYMMETRY ||
                  iod.mesh.bc_y0 == MeshData::OUTLET)
            ap -= a;
          else if(iod.mesh.bc_y0 == MeshData::STICKWALL)
            ap += a; //bb[k][j][i] -= a*v[k][j][i][dir+1]; //v[k][j-1][i][dir+1] should be -v[k][j][i][dir+1]
                                                //to have zero velocity on the wall
          else {
            fprintf(stdout,"*** Error: Detected unknown boundary condition type (%d).\n",
                    (int)iod.mesh.bc_y0);
            exit(-1);
          }
        }             


        //-----------
        // TOP 
        //-----------
        // Calculate F at the correct location (different for 
        if(dir==1)
          F = v[k][j][i][0]*0.5*(v[k][j][i][2] + v[k][j+1][i][2]);
        else {
          if(homo[k][j][i]) {
            rhou1 = dir==0 ? v[k][j+1][i-1][2]*v[k][j][i][0] : v[k-1][j+1][i][2]*v[k][j][i][0];
            rhou2 = v[k][j+1][i][2]*v[k][j][i][0];
          } else {
            cm = global_mesh.GetDy(j);
            cp = global_mesh.GetDy(j+1);
            cm_plus_cp = cm + cp;
            rhou1 = dir==0 ? v[k][j+1][i-1][2]*(cp*v[k][j][i-1][0] + cm*v[k][j+1][i-1][0])/cm_plus_cp
                           : v[k-1][j+1][i][2]*(cp*v[k-1][j][i][0] + cm*v[k-1][j+1][i][0])/cm_plus_cp;
            rhou2 = v[k][j+1][i][2]*(cp*v[k][j][i][0] + cm*v[k][j+1][i][0])/(cm+cp);
          }
          F = dir== 0 ? (dxr*rhou1 + dxl*rhou2)/(dxl+dxr) : (dzf*rhou1 + dzk*rhou2)/(dzk+dzf);
        }
        F *= dxdz;
        // Calculate D and the "a" coefficient
        a = std::max(-F, 0.0);
        if(dir==1)
          mu = Mu[id[k][j][i]];
        else {
          if(homo[k][j][i])
            mu = Mu[id[k][j][i]];
          else {
            // cm and cp have been calculated!   
            if(j==NY-1) {
              mu1 = dir==0 ? Mu[id[k][j][i-1]] : Mu[id[k-1][j][i]];
              mu2 = Mu[id[k][j][i]];
            } else {
              mu1 = dir==0 ? (cp*Mu[id[k][j][i-1]] + cm*Mu[id[k][j+1][i-1]])/cm_plus_cp
                           : (cp*Mu[id[k-1][j][i]] + cm*Mu[id[k-1][j+1][i]])/cm_plus_cp;
              mu2 = (cp*Mu[id[k][j][i]] + cm*Mu[id[k][j+1][i]])/cm_plus_cp;
            }
            mu = dir==0 ? (dxr*mu1 + dxl*mu2)/(dxl+dxr) : (dzf*mu1 + dzk*mu2)/(dzk+dzf);
          }
        }
        if(mu>0.0) {
          D  = mu*dxdz/dyt;
          a += D*PowerLaw(F/D);
        }
        ap += a;
        vhat[k][j][i] += a*v[k][j+1][i][dir+1];
        vhat_denom += a;

        // Add entry (or add to bb or ap, if j+1 is outside boundary)
        if(j+1<NY)
          row.PushEntry(i,j+1,k, -a);  //on the left hand side
        else { //j+1 is outside domain boundary
          if(iod.mesh.bc_ymax == MeshData::INLET || iod.mesh.bc_ymax == MeshData::INLET2 ||
             iod.mesh.bc_ymax == MeshData::OVERSET)
            bb[k][j][i] += a*v[k][j+1][i][dir+1];
          else if(iod.mesh.bc_ymax == MeshData::SLIPWALL || iod.mesh.bc_ymax == MeshData::SYMMETRY) {
            if(dir != 1) //otherwise, v[k][j+1][i][2] = 0
              ap -= a;
          } else if(iod.mesh.bc_ymax == MeshData::OUTLET) {
              ap -= a;
          } else if(iod.mesh.bc_ymax == MeshData::STICKWALL) {
            if(dir != 1) //otherwise, v[k][j+1][i][2] should be 0 as it is on the wall
              ap += a; //bb[k][j][i] -= a*v[k][j][i][dir+1]; //v[k][j+1][i][dir+1] should be -v[k][j][i][dir+1]
          } else {
            fprintf(stdout,"*** Error: Detected unknown boundary condition type (%d).\n",
                    (int)iod.mesh.bc_ymax);
            exit(-1);
          }
        }     
    

        //-----------
        // BACK
        //-----------
        // Calculate F at the correct location (different for dir=0,1,2)
        if(dir==2)
          F = v[k-1][j][i][0]*0.5*(v[k-1][j][i][3] + v[k][j][i][3]);
        else {
          if(homo[k][j][i]) {
            rhou1 = dir==0 ? v[k][j][i-1][3]*v[k][j][i][0] : v[k][j-1][i][3]*v[k][j][i][0];
            rhou2 = v[k][j][i][3]*v[k][j][i][0];
          } else {
            cm = global_mesh.GetDz(k-1);
            cp = global_mesh.GetDz(k);
            cm_plus_cp = cm + cp;
            rhou1 = dir==0 ? v[k][j][i-1][3]*(cp*v[k-1][j][i-1][0] + cm*v[k][j][i-1][0])/cm_plus_cp
                           : v[k][j-1][i][3]*(cp*v[k-1][j-1][i][0] + cm*v[k][j-1][i][0])/cm_plus_cp;
            rhou2 = v[k][j][i][3]*(cp*v[k-1][j][i][0] + cm*v[k][j][i][0])/(cm+cp);
          }
          F = dir==0 ? (dxr*rhou1 + dxl*rhou2)/(dxl+dxr) : (dyt*rhou1 + dyb*rhou2)/(dyb+dyt);
        }
        F *= dxdy;
        // Calculate D and the "a" coefficient
        a = std::max(F, 0.0);
        if(dir==2)
          mu = Mu[id[k-1][j][i]];
        else {
          if(homo[k][j][i])
            mu = Mu[id[k][j][i]];
          else {
            // cm and cp have been calculated!   
            if(k==0) {
              mu1 = dir==0 ? Mu[id[k][j][i-1]] : Mu[id[k][j-1][i]];
              mu2 = Mu[id[k][j][i]];
            } else {
              mu1 = dir==0 ? (cp*Mu[id[k-1][j][i-1]] + cm*Mu[id[k][j][i-1]])/cm_plus_cp
                           : (cp*Mu[id[k-1][j-1][i]] + cm*Mu[id[k][j-1][i]])/cm_plus_cp;
              mu2 = (cp*Mu[id[k-1][j][i]] + cm*Mu[id[k][j][i]])/cm_plus_cp;
            }
            mu = dir==0 ? (dxr*mu1 + dxl*mu2)/(dxl+dxr) : (dyt*mu1 + dyb*mu2)/(dyb+dyt);
          }
        }
        if(mu>0.0) {
          D  = mu*dxdy/dzk;
          a += D*PowerLaw(F/D);
        }
        ap += a;
        vhat[k][j][i] += a*v[k-1][j][i][dir+1];
        vhat_denom += a;

        // Add entry (or add to bb or ap, if k-1 is outside boundary)
        if(k-1>=0)
          row.PushEntry(i,j,k-1, -a);  //on the left hand side
        else { //k-1 is outside domain boundary (if it gets here, dir must be x or y (0 or 1))
          if(iod.mesh.bc_z0 == MeshData::INLET || iod.mesh.bc_z0 == MeshData::INLET2 ||
             iod.mesh.bc_z0 == MeshData::OVERSET)
            bb[k][j][i] += a*v[k-1][j][i][dir+1]; //+a*u or +a*v to the RHS
          else if(iod.mesh.bc_z0 == MeshData::SLIPWALL || iod.mesh.bc_z0 == MeshData::SYMMETRY ||
                  iod.mesh.bc_z0 == MeshData::OUTLET)
            ap -= a;
          else if(iod.mesh.bc_z0 == MeshData::STICKWALL)
            ap += a; //bb[k][j][i] -= a*v[k][j][i][dir+1]; //v[k-1][j][i][dir+1] should be -v[k][j][i][dir+1]
                                                //to have zero velocity on the wall
          else {
            fprintf(stdout,"*** Error: Detected unknown boundary condition type (%d).\n",
                    (int)iod.mesh.bc_z0);
            exit(-1);
          }
        }
             

        //-----------
        // FRONT 
        //-----------
        // Calculate F at the correct location (different for dir=0,1,2)
        if(dir==2)
          F = v[k][j][i][0]*0.5*(v[k][j][i][3] + v[k+1][j][i][3]);
        else {
          if(homo[k][j][i]) {
            rhou1 = dir==0 ? v[k+1][j][i-1][3]*v[k][j][i][0] : v[k+1][j-1][i][3]*v[k][j][i][0];
            rhou2 = v[k+1][j][i][3]*v[k][j][i][0];
          } else {
            cm = global_mesh.GetDz(k);
            cp = global_mesh.GetDz(k+1);
            cm_plus_cp = cm + cp;
            rhou1 = dir==0 ? v[k+1][j][i-1][3]*(cp*v[k][j][i-1][0] + cm*v[k+1][j][i-1][0])/cm_plus_cp
                           : v[k+1][j-1][i][3]*(cp*v[k][j-1][i][0] + cm*v[k+1][j-1][i][0])/cm_plus_cp;
            rhou2 = v[k+1][j][i][3]*(cp*v[k][j][i][0] + cm*v[k+1][j][i][0])/(cm+cp);
            F = dir== 0 ? (dxr*rhou1 + dxl*rhou2)/(dxl+dxr) : (dyt*rhou1 + dyb*rhou2)/(dyb+dyt);
          }
        }
        F *= dxdy;
        // Calculate D and the "a" coefficient
        a = std::max(-F, 0.0);
        if(dir==2)
          mu = Mu[id[k][j][i]];
        else {
          if(homo[k][j][i]) 
            mu = Mu[id[k][j][i]];
          else {
            // cm and cp have been calculated!   
            if(k==NZ-1) {
              mu1 = dir==0 ? Mu[id[k][j][i-1]] : Mu[id[k][j-1][i]];
              mu2 = Mu[id[k][j][i]];
            } else {
              mu1 = dir==0 ? (cp*Mu[id[k][j][i-1]] + cm*Mu[id[k+1][j][i-1]])/cm_plus_cp
                           : (cp*Mu[id[k][j-1][i]] + cm*Mu[id[k+1][j-1][i]])/cm_plus_cp;
              mu2 = (cp*Mu[id[k][j][i]] + cm*Mu[id[k+1][j][i]])/cm_plus_cp;
            }
            mu = dir==0 ? (dxr*mu1 + dxl*mu2)/(dxl+dxr) : (dyt*mu1 + dyb*mu2)/(dyb+dyt);
          }
        }
        if(mu>0.0) {
          D  = mu*dxdy/dzf;
          a += D*PowerLaw(F/D);
        }
        ap += a;
        vhat[k][j][i] += a*v[k+1][j][i][dir+1];
        vhat_denom += a;

        // Add entry (or add to bb or ap, if j+1 is outside boundary)
        if(k+1<NZ)
          row.PushEntry(i,j,k+1, -a);  //on the left hand side
        else { //k+1 is outside domain boundary
          if(iod.mesh.bc_zmax == MeshData::INLET || iod.mesh.bc_zmax == MeshData::INLET2 ||
             iod.mesh.bc_zmax == MeshData::OVERSET)
            bb[k][j][i] += a*v[k+1][j][i][dir+1];
          else if(iod.mesh.bc_zmax == MeshData::SLIPWALL || iod.mesh.bc_zmax == MeshData::SYMMETRY) {
            if(dir != 2) //otherwise, v[k+1][j][i][3] = 0
              ap -= a;
          } else if(iod.mesh.bc_zmax == MeshData::OUTLET) {
            ap -= a;
          } else if(iod.mesh.bc_zmax == MeshData::STICKWALL) {
            if(dir != 2) //otherwise, v[k+1][j][i][3] should be 0 as it is on the wall
              ap += a; //bb[k][j][i] -= a*v[k][j][i][dir+1]; //v[k+1][j][i][dir+1] should be -v[k][j][i][dir+1]
          } else {
            fprintf(stdout,"*** Error: Detected unknown boundary condition type (%d).\n",
                    (int)iod.mesh.bc_zmax);
            exit(-1);
          }
        }
             

        //------------------------------------------------------
        // Calculate and add the diagonal entry and the RHS (b)
        // Ref: Eqs. (5.62) and (6.8) in Patankar's book
        //------------------------------------------------------
        ap0 = dir==0 ? (dxr*v0[k][j][i-1][0] + dxl*v0[k][j][i][0])/(dxl+dxr) :
              dir==1 ? (dyt*v0[k][j-1][i][0] + dyb*v0[k][j][i][0])/(dyb+dyt) :
                       (dzf*v0[k-1][j][i][0] + dzk*v0[k][j][i][0])/(dzk+dzf);
        ap0 *= LocalDt ? dxdy*dz/dtloc[k][j][i] : dxdy*dz/dt;
        ap += ap0; //!< -Sp*dx*dy*dz, for source terms

        vhat[k][j][i] += ap0*v0[k][j][i][dir+1]; //!< +Sc*dx*dy*dz for source terms
        vhat_denom += ap0;

        // no pressure here (SIMPLER)

        bb[k][j][i] += ap0*v0[k][j][i][dir+1]; //!< +Sc*dx*dy*dz for source terms

        // Apply relaxation (Ref: Eq.(6) of Van Doormaal and Rathby, 1984)
        assert(Efactor>0.0);
        bb[k][j][i]   += ap*v0[k][j][i][dir+1]/Efactor; 

        ap *= 1.0 + 1.0/Efactor; 
        row.PushEntry(i,j,k, ap);

        // Also apply relaxation to vhat
        vhat[k][j][i] += vhat_denom*v0[k][j][i][dir+1]/Efactor;
        vhat[k][j][i] *= (1.0 + 1.0/Efactor)/vhat_denom;

        // Store diagonal for use in pressure correction equation
        assert(ap!=0.0);
        diag[k][j][i] = 1.0/ap;
        diag[k][j][i] *= dir==0 ? dydz : dir==1 ? dxdz : dxdy;

      }
    }
  }

  Vhat.RestoreDataPointerAndInsert();
  Ddiag.RestoreDataPointerAndInsert();
  Bv.RestoreDataPointerAndInsert();

  if(LocalDt)
    LocalDt->RestoreDataPointerToLocalVector();
}

//--------------------------------------------------------------------------

void
IncompressibleOperator::UpdateVelocityEquationRHS_SIMPLER(int dir, SpaceVariable3D &P, SpaceVariable3D &B)
{

  assert(dir==0 || dir==1 || dir==2);

  double*** p  = P.GetDataPointer();
  double*** bb = B.GetDataPointer();

  double dx, dy, dz, dxdy, dydz, dxdz;

  for(int k=k0; k<kmax; k++) {
    dz  = Dz[dir][k-k0];
    for(int j=j0; j<jmax; j++) {
      dy   = Dy[dir][j-j0];
      dydz = dy*dz;
      for(int i=i0; i<imax; i++) {

        if((dir==0 && i==0) || (dir==1 && j==0) || (dir==2 && k==0))
          continue; 

        dx   = Dx[dir][i-i0];
        dxdy = dx*dy;
        dxdz = dx*dz;

        bb[k][j][i] += dir==0 ? (p[k][j][i-1] - p[k][j][i])*dydz :
                       dir==1 ? (p[k][j-1][i] - p[k][j][i])*dxdz :
                                (p[k-1][j][i] - p[k][j][i])*dxdy;
      }
    }
  }

  P.RestoreDataPointerToLocalVector();
  B.RestoreDataPointerAndInsert();

}

//--------------------------------------------------------------------------

void
IncompressibleOperator::BuildPressureEquationRHS_SIMPLER(Vec5D*** v, double*** homo,
                                                         SpaceVariable3D &VXstar, SpaceVariable3D &VYstar,
                                                         SpaceVariable3D &VZstar, SpaceVariable3D &B,
                                                         Int3 *ijk_zero_p)
{

  GlobalMeshInfo& global_mesh(spo.GetGlobalMeshInfo());

  double*** bb = B.GetDataPointer();
  double*** ustar = VXstar.GetDataPointer();
  double*** vstar = VYstar.GetDataPointer();
  double*** wstar = VZstar.GetDataPointer();

  double rho, dx, dxl, dxr, dy, dyb, dyt, dz, dzk, dzf;

  for(int k=k0; k<kmax; k++) {
    dz  = global_mesh.GetDz(k);
    dzk = global_mesh.GetDz(k-1);
    dzf = global_mesh.GetDz(k+1);
    for(int j=j0; j<jmax; j++) {
      dy  = global_mesh.GetDy(j);
      dyb = global_mesh.GetDy(j-1);
      dyt = global_mesh.GetDy(j+1);
      for(int i=i0; i<imax; i++) {
        dx  = global_mesh.GetDx(i);
        dxl = global_mesh.GetDx(i-1);
        dxr = global_mesh.GetDx(i+1);

        // initialization
        bb[k][j][i] = 0.0;

        // set p = 0? (Otherwise, the linear system is singular, but *may* still be solvable by iterative methods)
        if(ijk_zero_p && i==(*ijk_zero_p)[0] && j==(*ijk_zero_p)[1] && k==(*ijk_zero_p)[2]) {
          bb[k][j][i] = 0.0;
          continue;
        }

        //-------
        // LEFT
        //-------
        if(i>0) {
          rho = homo[k][j][i] ? v[k][j][i][0] : (dx*v[k][j][i-1][0] + dxl*v[k][j][i][0])/(dxl+dx);
          bb[k][j][i] += rho*ustar[k][j][i]/dx;
        } else {
          // TODO: Only consider velocity b.c. for now (velocity is known)
          rho = v[k][j][i][0];
          bb[k][j][i] += rho*v[k][j][i][1]/dx; //boundary velocity
        }

        //-------
        // RIGHT 
        //-------
        // if i==NX-1, do nothing: at a boundary where normal velocity is known, a = 0. (Chap 6.7-3 of Patankar)
        if(i<NX-1) {
          rho = homo[k][j][i] ? v[k][j][i][0] : (dxr*v[k][j][i][0] + dx*v[k][j][i+1][0])/(dx+dxr);
          bb[k][j][i] -= rho*ustar[k][j][i+1]/dx;
        } else {
          // TODO: Only consider velocity b.c. for now (velocity is known)
          rho = v[k][j][i][0];
          bb[k][j][i] -= rho*v[k][j][i+1][1]/dx; //boundary velocity
        }
 
        //-------
        // BOTTOM 
        //-------
        if(j>0) {
          rho = homo[k][j][i] ? v[k][j][i][0] : (dy*v[k][j-1][i][0] + dyb*v[k][j][i][0])/(dyb+dy);
          bb[k][j][i] += rho*vstar[k][j][i]/dy;
        } else {
          // TODO: Only consider velocity b.c. for now (velocity is known)
          rho = v[k][j][i][0];
          bb[k][j][i] += rho*v[k][j][i][2]/dy; //boundary velocity
        }

        //-------
        // TOP 
        //-------
        if(j<NY-1) {
          rho = homo[k][j][i] ? v[k][j][i][0] : (dyt*v[k][j][i][0] + dy*v[k][j+1][i][0])/(dy+dyt);
          bb[k][j][i] -= rho*vstar[k][j+1][i]/dy;
        } else {
          // TODO: Only consider velocity b.c. for now (velocity is known)
          rho = v[k][j][i][0];
          bb[k][j][i] -= rho*v[k][j+1][i][2]/dx; //boundary velocity
        }
  
        //-------
        // BACK 
        //-------
        if(k>0) {
          rho = homo[k][j][i] ? v[k][j][i][0] : (dz*v[k-1][j][i][0] + dzk*v[k][j][i][0])/(dzk+dz);
          bb[k][j][i] += rho*wstar[k][j][i]/dz;
        } else {
          // TODO: Only consider velocity b.c. for now (velocity is known)
          rho = v[k][j][i][0];
          bb[k][j][i] += rho*v[k][j][i][3]/dz; //boundary velocity
        }

        //-------
        // FRONT 
        //-------
        if(k<NZ-1) {
          rho = homo[k][j][i] ? v[k][j][i][0] : (dzf*v[k][j][i][0] + dz*v[k+1][j][i][0])/(dz+dzf);
          bb[k][j][i] -= rho*wstar[k+1][j][i]/dz;
        } else {
          // TODO: Only consider velocity b.c. for now (velocity is known)
          rho = v[k][j][i][0];
          bb[k][j][i] -= rho*v[k+1][j][i][3]/dz; //boundary velocity
        }
 
      }
    }
  }


  B.RestoreDataPointerAndInsert();
  VXstar.RestoreDataPointerToLocalVector();
  VYstar.RestoreDataPointerToLocalVector();
  VZstar.RestoreDataPointerToLocalVector();

}

//--------------------------------------------------------------------------

void
IncompressibleOperator::CalculateVelocityTildePISO(int dir, Vec5D*** v0, Vec5D*** v, double*** id, double*** homo,
                                                   SpaceVariable3D &Vprime, SpaceVariable3D &Vtilde,
                                                   double Efactor, double dt, SpaceVariable3D *LocalDt)
{

  // very similar to BuildVelocityEquationSIMPLE

  assert(dir==0 || dir==1 || dir==2);

  double*** dtloc = NULL;
  if(LocalDt) { //local time-stepping, dealt with separately.
    assert(iod.ts.timestep<=0.0); //shouldn't have constant time-step size.
    dtloc = LocalDt->GetDataPointer();
  }

  GlobalMeshInfo& global_mesh(spo.GetGlobalMeshInfo());

  double*** vprime = Vprime.GetDataPointer();
  double*** vtilde = Vtilde.GetDataPointer();

  double dx, dy, dz, dxl, dxr, dyb, dyt, dzk, dzf, dxdy, dydz, dxdz;
  double cm(1.0), cp(1.0), cm_plus_cp(0.0), rhou1, rhou2;
  double a, ap, ap0, F, D, mu, mu1, mu2;

  for(int k=k0; k<kmax; k++) {
    dz  = Dz[dir][k-k0];
    dzk = dz_k[dir][k-k0];
    dzf = dz_f[dir][k-k0];
    for(int j=j0; j<jmax; j++) {
      dy   = Dy[dir][j-j0];
      dyb  = dy_b[dir][j-j0];
      dyt  = dy_t[dir][j-j0];
      dydz = dy*dz;
      for(int i=i0; i<imax; i++) {

        if((dir==0 && i==0) || (dir==1 && j==0) || (dir==2 && k==0)) {
          vtilde[k][j][i] = 0.0;
          continue;
        }

        dx   = Dx[dir][i-i0];
        dxl  = dx_l[dir][i-i0];
        dxr  = dx_r[dir][i-i0];
        dxdy = dx*dy;
        dxdz = dx*dz;


        //---------------------------------------------------
        // Reference: Patankar's book, Eqs. (5.61) - (5.64)
        //---------------------------------------------------

        ap = 0.0; //diagonal
        vtilde[k][j][i] = 0.0;

        //-----------
        // LEFT
        //-----------
        // Calculate F at the correct location (different for dir=0,1,2)
        if(dir==0)
          F = v[k][j][i-1][0]*0.5*(v[k][j][i-1][1] + v[k][j][i][1]);
        else {
          if(homo[k][j][i]) { //this "cell" (i,j,k) is in a neighborhood w/ constant rho and mu
            rhou1 = dir==1 ? v[k][j-1][i][1]*v[k][j][i][0] : v[k-1][j][i][1]*v[k][j][i][0];
            rhou2 = v[k][j][i][1]*v[k][j][i][0];
          } else {
            cm = global_mesh.GetDx(i-1);
            cp = global_mesh.GetDx(i);
            cm_plus_cp = cm + cp;
            rhou1 = dir==1 ? v[k][j-1][i][1]*(cp*v[k][j-1][i-1][0] + cm*v[k][j-1][i][0])/cm_plus_cp
                           : v[k-1][j][i][1]*(cp*v[k-1][j][i-1][0] + cm*v[k-1][j][i][0])/cm_plus_cp;
            rhou2 = v[k][j][i][1]*(cp*v[k][j][i-1][0] + cm*v[k][j][i][0])/(cm+cp);
          }
          F = dir== 1 ? (dyt*rhou1 + dyb*rhou2)/(dyb+dyt) : (dzf*rhou1 + dzk*rhou2)/(dzk+dzf);
        }
        F *= dydz;
        // Calculate D and the "a" coefficient
        a = std::max(F, 0.0);
        if(dir==0)
          mu = Mu[id[k][j][i-1]];
        else {
          if(homo[k][j][i])
            mu = Mu[id[k][j][i]];
          else {
            // cm and cp have been calculated!   
            if(i==0) {
              mu1 = dir==1 ? Mu[id[k][j-1][i]] : Mu[id[k-1][j][i]];
              mu2 = Mu[id[k][j][i]];
            } else {
              mu1 = dir==1 ? (cp*Mu[id[k][j-1][i-1]] + cm*Mu[id[k][j-1][i]])/cm_plus_cp
                           : (cp*Mu[id[k-1][j][i-1]] + cm*Mu[id[k-1][j][i]])/cm_plus_cp;
              mu2 = (cp*Mu[id[k][j][i-1]] + cm*Mu[id[k][j][i]])/cm_plus_cp;
            }
            mu = dir==1 ? (dyt*mu1 + dyb*mu2)/(dyb+dyt) : (dzf*mu1 + dzk*mu2)/(dzk+dzf);
          }
        }
        if(mu>0.0) {
          D  = mu*dydz/dxl;
          a += D*PowerLaw(F/D);
        }
        ap += a;
        vtilde[k][j][i] += a*vprime[k][j][i-1];
        //TODO: BOUNDARY CONDITION!!!
 
        //-----------
        // RIGHT 
        //-----------
        // Calculate F at the correct location (different for dir=0,1,2)
        if(dir==0)
          F = v[k][j][i][0]*0.5*(v[k][j][i][1] + v[k][j][i+1][1]);
        else {
          if(homo[k][j][i]) {
            rhou1 = dir==1 ? v[k][j-1][i+1][1]*v[k][j][i][0] : v[k-1][j][i+1][1]*v[k][j][i][0];
            rhou2 = v[k][j][i+1][1]*v[k][j][i][0];
          } else {
            cm = global_mesh.GetDx(i);
            cp = global_mesh.GetDx(i+1);
            cm_plus_cp = cm + cp;
            rhou1 = dir==1 ? v[k][j-1][i+1][1]*(cp*v[k][j-1][i][0] + cm*v[k][j-1][i+1][0])/cm_plus_cp
                           : v[k-1][j][i+1][1]*(cp*v[k-1][j][i][0] + cm*v[k-1][j][i+1][0])/cm_plus_cp;
            rhou2 = v[k][j][i+1][1]*(cp*v[k][j][i][0] + cm*v[k][j][i+1][0])/(cm+cp);
          }
          F = dir== 1 ? (dyt*rhou1 + dyb*rhou2)/(dyb+dyt) : (dzf*rhou1 + dzk*rhou2)/(dzk+dzf);
        }
        F *= dydz;
        // Calculate D and the "a" coefficient
        a = std::max(-F, 0.0);
        if(dir==0)
          mu = Mu[id[k][j][i]];
        else {
          if(homo[k][j][i])
            mu = Mu[id[k][j][i]];
          else {
            // cm and cp have been calculated!   
            if(i==NX-1) {
              mu1 = dir==1 ? Mu[id[k][j-1][i]] : Mu[id[k-1][j][i]];
              mu2 = Mu[id[k][j][i]];
            } else {
              mu1 = dir==1 ? (cp*Mu[id[k][j-1][i]] + cm*Mu[id[k][j-1][i+1]])/cm_plus_cp
                           : (cp*Mu[id[k-1][j][i]] + cm*Mu[id[k-1][j][i+1]])/cm_plus_cp;
              mu2 = (cp*Mu[id[k][j][i]] + cm*Mu[id[k][j][i+1]])/cm_plus_cp;
            }
            mu = dir==1 ? (dyt*mu1 + dyb*mu2)/(dyb+dyt) : (dzf*mu1 + dzk*mu2)/(dzk+dzf);
          }
        }
        if(mu>0.0) {
          D  = mu*dydz/dxr;
          a += D*PowerLaw(F/D);
        }
        ap += a;
        vtilde[k][j][i] += a*vprime[k][j][i+1];
             
    
        //-----------
        // BOTTOM 
        //-----------
        // Calculate F at the correct location (different for dir=0,1,2)
        if(dir==1)
          F = v[k][j-1][i][0]*0.5*(v[k][j-1][i][2] + v[k][j][i][2]);
        else {
          if(homo[k][j][i]) {
            rhou1 = dir==0 ? v[k][j][i-1][2]*v[k][j][i][0] : v[k-1][j][i][2]*v[k][j][i][0];
            rhou2 = v[k][j][i][2]*v[k][j][i][0];
          } else {
            cm = global_mesh.GetDy(j-1);
            cp = global_mesh.GetDy(j);
            cm_plus_cp = cm + cp;
            rhou1 = dir==0 ? v[k][j][i-1][2]*(cp*v[k][j-1][i-1][0] + cm*v[k][j][i-1][0])/cm_plus_cp
                           : v[k-1][j][i][2]*(cp*v[k-1][j-1][i][0] + cm*v[k-1][j][i][0])/cm_plus_cp;
            rhou2 = v[k][j][i][2]*(cp*v[k][j-1][i][0] + cm*v[k][j][i][0])/(cm+cp);
          }
          F = dir== 0 ? (dxr*rhou1 + dxl*rhou2)/(dxl+dxr) : (dzf*rhou1 + dzk*rhou2)/(dzk+dzf);
        }
        F *= dxdz;
        // Calculate D and the "a" coefficient
        a = std::max(F, 0.0);
        if(dir==1)
          mu = Mu[id[k][j-1][i]];
        else {
          if(homo[k][j][i])
            mu = Mu[id[k][j][i]]; 
          else {
            // cm and cp have been calculated!   
            if(j==0) {
              mu1 = dir==0 ? Mu[id[k][j][i-1]] : Mu[id[k-1][j][i]];
              mu2 = Mu[id[k][j][i]];
            } else {
              mu1 = dir==0 ? (cp*Mu[id[k][j-1][i-1]] + cm*Mu[id[k][j][i-1]])/cm_plus_cp
                           : (cp*Mu[id[k-1][j-1][i]] + cm*Mu[id[k-1][j][i]])/cm_plus_cp;
              mu2 = (cp*Mu[id[k][j-1][i]] + cm*Mu[id[k][j][i]])/cm_plus_cp;
            }
            mu = dir==0 ? (dxr*mu1 + dxl*mu2)/(dxl+dxr) : (dzf*mu1 + dzk*mu2)/(dzk+dzf);
          }
        }
        if(mu>0.0) {
          D  = mu*dxdz/dyb;
          a += D*PowerLaw(F/D);
        }
        ap += a;
        vtilde[k][j][i] += a*vprime[k][j-1][i];
             

        //-----------
        // TOP 
        //-----------
        // Calculate F at the correct location (different for dir=0,1,2)
        if(dir==1)
          F = v[k][j][i][0]*0.5*(v[k][j][i][2] + v[k][j+1][i][2]);
        else {
          if(homo[k][j][i]) {
            rhou1 = dir==0 ? v[k][j+1][i-1][2]*v[k][j][i][0] : v[k-1][j+1][i][2]*v[k][j][i][0];
            rhou2 = v[k][j+1][i][2]*v[k][j][i][0];
          } else {
            cm = global_mesh.GetDy(j);
            cp = global_mesh.GetDy(j+1);
            cm_plus_cp = cm + cp;
            rhou1 = dir==0 ? v[k][j+1][i-1][2]*(cp*v[k][j][i-1][0] + cm*v[k][j+1][i-1][0])/cm_plus_cp
                           : v[k-1][j+1][i][2]*(cp*v[k-1][j][i][0] + cm*v[k-1][j+1][i][0])/cm_plus_cp;
            rhou2 = v[k][j+1][i][2]*(cp*v[k][j][i][0] + cm*v[k][j+1][i][0])/(cm+cp);
          }
          F = dir== 0 ? (dxr*rhou1 + dxl*rhou2)/(dxl+dxr) : (dzf*rhou1 + dzk*rhou2)/(dzk+dzf);
        }
        F *= dxdz;
        // Calculate D and the "a" coefficient
        a = std::max(-F, 0.0);
        if(dir==1)
          mu = Mu[id[k][j][i]];
        else {
          if(homo[k][j][i])
            mu = Mu[id[k][j][i]];
          else {
            // cm and cp have been calculated!   
            if(j==NY-1) {
              mu1 = dir==0 ? Mu[id[k][j][i-1]] : Mu[id[k-1][j][i]];
              mu2 = Mu[id[k][j][i]];
            } else {
              mu1 = dir==0 ? (cp*Mu[id[k][j][i-1]] + cm*Mu[id[k][j+1][i-1]])/cm_plus_cp
                           : (cp*Mu[id[k-1][j][i]] + cm*Mu[id[k-1][j+1][i]])/cm_plus_cp;
              mu2 = (cp*Mu[id[k][j][i]] + cm*Mu[id[k][j+1][i]])/cm_plus_cp;
            }
            mu = dir==0 ? (dxr*mu1 + dxl*mu2)/(dxl+dxr) : (dzf*mu1 + dzk*mu2)/(dzk+dzf);
          }
        }
        if(mu>0.0) {
          D  = mu*dxdz/dyt;
          a += D*PowerLaw(F/D);
        }
        ap += a;
        vtilde[k][j][i] += a*vprime[k][j+1][i];
             
    
        //-----------
        // BACK
        //-----------
        // Calculate F at the correct location (different for dir=0,1,2)
        if(dir==2)
          F = v[k-1][j][i][0]*0.5*(v[k-1][j][i][3] + v[k][j][i][3]);
        else {
          if(homo[k][j][i]) {
            rhou1 = dir==0 ? v[k][j][i-1][3]*v[k][j][i][0] : v[k][j-1][i][3]*v[k][j][i][0];
            rhou2 = v[k][j][i][3]*v[k][j][i][0];
          } else {
            cm = global_mesh.GetDz(k-1);
            cp = global_mesh.GetDz(k);
            cm_plus_cp = cm + cp;
            rhou1 = dir==0 ? v[k][j][i-1][3]*(cp*v[k-1][j][i-1][0] + cm*v[k][j][i-1][0])/cm_plus_cp
                           : v[k][j-1][i][3]*(cp*v[k-1][j-1][i][0] + cm*v[k][j-1][i][0])/cm_plus_cp;
            rhou2 = v[k][j][i][3]*(cp*v[k-1][j][i][0] + cm*v[k][j][i][0])/(cm+cp);
          }
          F = dir==0 ? (dxr*rhou1 + dxl*rhou2)/(dxl+dxr) : (dyt*rhou1 + dyb*rhou2)/(dyb+dyt);
        }
        F *= dxdy;
        // Calculate D and the "a" coefficient
        a = std::max(F, 0.0);
        if(dir==2)
          mu = Mu[id[k-1][j][i]];
        else {
          if(homo[k][j][i])
            mu = Mu[id[k][j][i]];
          else {
            // cm and cp have been calculated!   
            if(k==0) {
              mu1 = dir==0 ? Mu[id[k][j][i-1]] : Mu[id[k][j-1][i]];
              mu2 = Mu[id[k][j][i]];
            } else {
              mu1 = dir==0 ? (cp*Mu[id[k-1][j][i-1]] + cm*Mu[id[k][j][i-1]])/cm_plus_cp
                           : (cp*Mu[id[k-1][j-1][i]] + cm*Mu[id[k][j-1][i]])/cm_plus_cp;
              mu2 = (cp*Mu[id[k-1][j][i]] + cm*Mu[id[k][j][i]])/cm_plus_cp;
            }
            mu = dir==0 ? (dxr*mu1 + dxl*mu2)/(dxl+dxr) : (dyt*mu1 + dyb*mu2)/(dyb+dyt);
          }
        }
        if(mu>0.0) {
          D  = mu*dxdy/dzk;
          a += D*PowerLaw(F/D);
        }
        ap += a;
        vtilde[k][j][i] += a*vprime[k-1][j][i];
             

        //-----------
        // FRONT 
        //-----------
        // Calculate F at the correct location (different for dir=0,1,2)
        if(dir==2)
          F = v[k][j][i][0]*0.5*(v[k][j][i][3] + v[k+1][j][i][3]);
        else {
          if(homo[k][j][i]) {
            rhou1 = dir==0 ? v[k+1][j][i-1][3]*v[k][j][i][0] : v[k+1][j-1][i][3]*v[k][j][i][0];
            rhou2 = v[k+1][j][i][3]*v[k][j][i][0];
          } else {
            cm = global_mesh.GetDz(k);
            cp = global_mesh.GetDz(k+1);
            cm_plus_cp = cm + cp;
            rhou1 = dir==0 ? v[k+1][j][i-1][3]*(cp*v[k][j][i-1][0] + cm*v[k+1][j][i-1][0])/cm_plus_cp
                           : v[k+1][j-1][i][3]*(cp*v[k][j-1][i][0] + cm*v[k+1][j-1][i][0])/cm_plus_cp;
            rhou2 = v[k+1][j][i][3]*(cp*v[k][j][i][0] + cm*v[k+1][j][i][0])/(cm+cp);
            F = dir== 0 ? (dxr*rhou1 + dxl*rhou2)/(dxl+dxr) : (dyt*rhou1 + dyb*rhou2)/(dyb+dyt);
          }
        }
        F *= dxdy;
        // Calculate D and the "a" coefficient
        a = std::max(-F, 0.0);
        if(dir==2)
          mu = Mu[id[k][j][i]];
        else {
          if(homo[k][j][i]) 
            mu = Mu[id[k][j][i]];
          else {
            // cm and cp have been calculated!   
            if(k==NZ-1) {
              mu1 = dir==0 ? Mu[id[k][j][i-1]] : Mu[id[k][j-1][i]];
              mu2 = Mu[id[k][j][i]];
            } else {
              mu1 = dir==0 ? (cp*Mu[id[k][j][i-1]] + cm*Mu[id[k+1][j][i-1]])/cm_plus_cp
                           : (cp*Mu[id[k][j-1][i]] + cm*Mu[id[k+1][j-1][i]])/cm_plus_cp;
              mu2 = (cp*Mu[id[k][j][i]] + cm*Mu[id[k+1][j][i]])/cm_plus_cp;
            }
            mu = dir==0 ? (dxr*mu1 + dxl*mu2)/(dxl+dxr) : (dyt*mu1 + dyb*mu2)/(dyb+dyt);
          }
        }
        if(mu>0.0) {
          D  = mu*dxdy/dzf;
          a += D*PowerLaw(F/D);
        }
        ap += a;
        vtilde[k][j][i] += a*vprime[k+1][j][i];
             

        //------------------------------------------------------
        // Calculate and add the diagonal entry and the RHS (b)
        // Ref: Eqs. (5.62) and (6.8) in Patankar's book
        //------------------------------------------------------
        ap0 = dir==0 ? (dxr*v0[k][j][i-1][0] + dxl*v0[k][j][i][0])/(dxl+dxr) :
              dir==1 ? (dyt*v0[k][j-1][i][0] + dyb*v0[k][j][i][0])/(dyb+dyt) :
                       (dzf*v0[k-1][j][i][0] + dzk*v0[k][j][i][0])/(dzk+dzf);
        ap0 *= LocalDt ? dxdy*dz/dtloc[k][j][i] : dxdy*dz/dt;
        ap += ap0; //!< -Sp*dx*dy*dz, for source terms

        // Apply relaxation (Ref: Eq.(6) of Van Doormaal and Rathby, 1984)
        assert(Efactor>0.0);
        ap *= 1.0 + 1.0/Efactor; 

        // Store diagonal for use in pressure correction equation
        assert(ap!=0.0);
        vtilde[k][j][i] /= ap;

      }
    }
  }

  Vtilde.RestoreDataPointerAndInsert();
  Vprime.RestoreDataPointerToLocalVector();

  if(LocalDt)
    LocalDt->RestoreDataPointerToLocalVector();
}


//--------------------------------------------------------------------------

void
IncompressibleOperator::CalculateMomentumChanges(Vec5D*** v0, SpaceVariable3D &V, double*** id,
                                                 SpaceVariable3D &R3)
{

  Vec5D*** v = (Vec5D***)V.GetDataPointer();
  Vec3D*** r = (Vec3D***)R3.GetDataPointer();

  double rho0, rho;

  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        if(id[k][j][i] == INACTIVE_MATERIAL_ID) {
          r[k][j][i] = 0.0;
          continue;
        }

        rho0 = v0[k][j][i][0];
        rho  = v[k][j][i][0];

        for(int p=0; p<3; p++)
          r[k][j][i][p] = rho*v[k][j][i][p+1] - rho0*v0[k][j][i][p+1];
      }

  V.RestoreDataPointerToLocalVector();
  R3.RestoreDataPointerToLocalVector(); //!< no need of communication

}

//--------------------------------------------------------------------------

void
IncompressibleOperator::BuildSATurbulenceEquationSIMPLE(Vec5D*** v, double*** id, double*** vturb0,
                                                        double*** vturb,//SA eddy viscosity working var prev & current
                                                        vector<RowEntries> &vlin_rows, SpaceVariable3D &B,
                                                        double Efactor, double cw1_reduction,
                                                        double dt, SpaceVariable3D *LocalDt)
{
  if(iod.rans.model   != RANSTurbulenceModelData::SPALART_ALLMARAS ||
     iod.rans.example != RANSTurbulenceModelData::FLAT_PLATE) {
    print_error("*** Error: Detected unsupported turbulence model or example problem.\n");
    exit_mpi();
  }

  double*** dtloc = NULL;
  if(LocalDt) { //local time-stepping, dealt with separately.
    assert(iod.ts.timestep<=0.0); //shouldn't have constant time-step size.
    dtloc = LocalDt->GetDataPointer();
  }

  GlobalMeshInfo& global_mesh(spo.GetGlobalMeshInfo());

  double*** bb = B.GetDataPointer();

  vlin_rows.clear(); //clear existing data (to be safe)
  int row_counter = 0;

  double dx, dy, dz, dxl, dxr, dyb, dyt, dzk, dzf, dxdy, dydz, dxdz;
  double a, ap, ap0, F, D, nu, nustar;
  double dudy, dudz, dvdx, dvdz, dwdx, dwdy;
  Vec3D dnut;


  ///Spalart-Allmaras constants
  double cb1 = 0.1355; //production
  double sigma = 2.0/3.0, cb2 = 0.622; // diffusion
  double cw2 = 0.3, cw3 = 2.0, cv1 = 7.1, ct3 = 1.2, ct4 = 0.5, kappa = 0.41; //source
  double cw3_pow6 = pow(cw3, 6.0);
  double cw1 = cb1/(kappa*kappa) + (1.0+cb2)/sigma; //source
  cw1 /= cw1_reduction; //TODO: THIS IS A HACK
  double fv1,fv2,fw,ft2,fn,g,r,S,Sbar;
  double chi,Omega,d2w,k2d2;
  Vec3D vort;

  for(int k=k0; k<kmax; k++) {
    dz  = global_mesh.GetDz(k);
    dzk = 0.5*(global_mesh.GetDz(k-1) + dz);
    dzf = 0.5*(dz + global_mesh.GetDz(k+1));
    for(int j=j0; j<jmax; j++) {
      dy   = global_mesh.GetDy(j);
      dyb  = 0.5*(global_mesh.GetDy(j-1) + dy);
      dyt  = 0.5*(dy + global_mesh.GetDy(j+1));
      dydz = dy*dz;
      for(int i=i0; i<imax; i++) {
         
        nu = Mu[id[k][j][i]]/v[k][j][i][0]; //kinematic viscosity

        // insert the row
        vlin_rows.push_back(RowEntries(7)); // at most 7 non-zero entries on each row
        RowEntries &row(vlin_rows[row_counter]);
        row.SetRow(i,j,k);
        row_counter++;

/*
        // TEMPORARY FIX: Setting pts left of the wall to freestream value
        Vec3D coords = global_mesh.GetXYZ(i,j,k);
        if (coords[0] < 0.0){
          //fprintf(stdout,"\nClosest pt to leading edge is getting fixed");
          row.PushEntry(i,j,k, 1.0); 
          bb[k][j][i] = 8.0e-7; 
          continue;
        }
*/

        //---------------------------------------------------
        // Calculating d2w (distance to nearest wall)
        //---------------------------------------------------
        d2w = GetDistanceToWall(global_mesh.GetXYZ(i,j,k));
        if(d2w==0.0) { //at the wall, nu_t = 0.0
          row.PushEntry(i,j,k, 1.0); 
          bb[k][j][i] = 0.0; 
          continue;
        }
 
        dx   = global_mesh.GetDx(i);
        dxl  = 0.5*(global_mesh.GetDx(i-1) + dx);
        dxr  = 0.5*(dx + global_mesh.GetDx(i+1));
        dxdy = dx*dy;
        dxdz = dx*dz;
        
        bool flat_plate_wall = (iod.rans.example == RANSTurbulenceModelData::FLAT_PLATE) &&
                               (global_mesh.GetX(i) >= iod.rans.example_param_1); //stick wall

        //---------------------------------------------------
        // Calculate vorticity 
        // Note: Must avoid accessing edge/corner ghosts outside physical domain
        //---------------------------------------------------
        double vtmp[3];
        for(int p=-1; p<2; p++)
          vtmp[p+1] = (v[k][j+p][i][1] + v[k][j+p][i+1][1]) / 2.0;
        if(i+1>=global_mesh.NX) { //may need correction
          if(j-1<0) //edge
            vtmp[0] = ApplyVelocityBoundaryConditionLocal(0, GhostPoint::BOTTOM, vtmp[1], flat_plate_wall); //0~x-velo
          if(j+1>=global_mesh.NY) //edge
            vtmp[2] = ApplyVelocityBoundaryConditionLocal(0, GhostPoint::TOP, vtmp[1]);
        }
        dudy = -dyt/(dyb*(dyt+dyb))*vtmp[0] + (dyt-dyb)/(dyt*dyb)*vtmp[1] + dyb/(dyt*(dyt+dyb))*vtmp[2];

        for(int p=-1; p<2; p++)
          vtmp[p+1] = (v[k+p][j][i][1] + v[k+p][j][i+1][1]) / 2.0;

        if(i+1>=global_mesh.NX) { //may need correction
          if(k-1<0) //edge
            vtmp[0] = ApplyVelocityBoundaryConditionLocal(0, GhostPoint::BACK, vtmp[1]);
          if(k+1>=global_mesh.NZ) //edge
            vtmp[2] = ApplyVelocityBoundaryConditionLocal(0, GhostPoint::FRONT, vtmp[1]);
        }
        dudz = -dzf/(dzk*(dzf+dzk))*vtmp[0] + (dzf-dzk)/(dzf*dzk)*vtmp[1] + dzk/(dzf*(dzf+dzk))*vtmp[2];

        for(int p=-1; p<2; p++)
          vtmp[p+1] = (v[k][j][i+p][2] + v[k][j+1][i+p][2]) / 2.0;
        if(j+1>=global_mesh.NY) { //may need correction
          if(i-1<0) //edge
            vtmp[0] = ApplyVelocityBoundaryConditionLocal(1, GhostPoint::LEFT, vtmp[1]); //1~y-velo
          if(i+1>=global_mesh.NX) //edge
            vtmp[2] = ApplyVelocityBoundaryConditionLocal(1, GhostPoint::RIGHT, vtmp[1]);
        }
        dvdx = -dxr/(dxl*(dxl+dxr))*vtmp[0] + (dxr-dxl)/(dxl*dxr)*vtmp[1] + dxl/(dxr*(dxl+dxr))*vtmp[2];

        for(int p=-1; p<2; p++)
          vtmp[p+1] = (v[k+p][j][i][2] + v[k+p][j+1][i][2]) / 2.0;
        if(j+1>=global_mesh.NY) { //may need correction
          if(k-1<0) //edge
            vtmp[0] = ApplyVelocityBoundaryConditionLocal(1, GhostPoint::BACK, vtmp[1]); //1~y-velo
          if(k+1>=global_mesh.NZ) //edge
            vtmp[2] = ApplyVelocityBoundaryConditionLocal(1, GhostPoint::FRONT, vtmp[1]);
        }
        dvdz = -dzf/(dzk*(dzf+dzk))*vtmp[0] + (dzf-dzk)/(dzf*dzk)*vtmp[1] + dzk/(dzf*(dzf+dzk))*vtmp[2];

        for(int p=-1; p<2; p++)
          vtmp[p+1] = (v[k][j][i+p][3] + v[k+1][j][i+p][3]) / 2.0;
        if(k+1>=global_mesh.NZ) { //may need correction
          if(i-1<0) //edge
            vtmp[0] = ApplyVelocityBoundaryConditionLocal(2, GhostPoint::LEFT, vtmp[1]); //2~z-velo
          if(i+1>=global_mesh.NX) //edge
            vtmp[2] = ApplyVelocityBoundaryConditionLocal(2, GhostPoint::RIGHT, vtmp[1]);
        }
        dwdx = -dxr/(dxl*(dxl+dxr))*vtmp[0] + (dxr-dxl)/(dxl*dxr)*vtmp[1] + dxl/(dxr*(dxl+dxr))*vtmp[2];

        for(int p=-1; p<2; p++)
          vtmp[p+1] = (v[k][j+p][i][3] + v[k+1][j+p][i][3]) / 2.0;
        if(k+1>=global_mesh.NZ) { //may need correction
          if(j-1<0) //edge
            vtmp[0] = ApplyVelocityBoundaryConditionLocal(2, GhostPoint::BOTTOM, vtmp[1], flat_plate_wall); //2~z-velo
          if(j+1>=global_mesh.NY) //edge
            vtmp[2] = ApplyVelocityBoundaryConditionLocal(2, GhostPoint::TOP, vtmp[1]);
        }
        dwdy = -dyt/(dyb*(dyt+dyb))*vtmp[0] + (dyt-dyb)/(dyt*dyb)*vtmp[1] + dyb/(dyt*(dyt+dyb))*vtmp[2];

        vort[0] = dwdy-dvdz;
        vort[1] = dudz-dwdx;
        vort[2] = dvdx-dudy;
        Omega = vort.norm();

        //---------------------------------------------------
        // Spalart-Allmaras Additional Eqs. 
        //---------------------------------------------------
        bool neg_nut = vturb[k][j][i]<0.0;
        chi = vturb[k][j][i]/nu;
        fv1 = pow(chi,3)/(pow(chi,3) + pow(cv1,3));
        fv2 = 1.0 - chi/(1.0 + chi*fv1);
        ft2 = ct3*exp(-ct4*chi*chi);
        fn  = neg_nut ? (16.0 + pow(chi,3))/(16.0 - pow(chi,3)) : 1.0;
        k2d2 = kappa*kappa*d2w*d2w;
        Sbar = fv2*vturb[k][j][i]/k2d2; 

        //---------------------------------------------
        // Limiting S using Note 1(c) of https://turbmodels.larc.nasa.gov/spalart.html
        if(Sbar>=-0.7*Omega) //if(true) ---> original w/o limiting
          S = Omega + Sbar;
        else
          S = Omega - Omega*(0.49*Omega + 0.9*Sbar)/(0.5*Omega + Sbar);

        r = S==0 ? 10.0 : std::min(vturb[k][j][i]/(S*k2d2), 10.0); 
        //---------------------------------------------

        g = r + cw2*(pow(r,6.0)-r);
        fw = g * pow((1+cw3_pow6)/(pow(g,6.0)+cw3_pow6), 1.0/6.0);

        //---------------------------------------------------
        // Reference: Patankar's book, Eqs. (5.61) - (5.64)
        //---------------------------------------------------

        ap = 0.0; //diagonal
        bb[k][j][i] = 0.0; //rhs


        //-----------
        // LEFT
        //-----------
        // Calculate F at the correct location 
        F = v[k][j][i][1]*dydz;
        // Calculate D and the "a" coefficient
        a = std::max(F, 0.0);
        nustar = (nu + fn*vturb[k][j][i]) / sigma;
        if(nustar>0.0) {
          D  = nustar*dydz/dxl;
          a += D*PowerLaw(F/D);
        }
        ap += a;

        // Add entry (or add to bb or ap, if i-1 is outside boundary)
        if(i-1>=0)
          row.PushEntry(i-1,j,k, -a);  //on the left hand side
        else { //i-1 is outside domain boundary 
          if(iod.mesh.bc_x0 == MeshData::INLET || iod.mesh.bc_x0 == MeshData::INLET2 ||
             iod.mesh.bc_x0 == MeshData::OVERSET)
            bb[k][j][i] += a*vturb[k][j][i-1]; //+a*vturb to the RHS 
          else if(iod.mesh.bc_x0 == MeshData::SYMMETRY || iod.mesh.bc_x0 == MeshData::OUTLET)
            ap -= a;
          else if(iod.mesh.bc_x0 == MeshData::STICKWALL)
            ap += a; //bb[k][j][i] -= a*vturb[k][j][i]; //vturb[k][j][i-1] should be -vturb[k][j][i]
                                             //to have zero on the wall
          else {
            fprintf(stdout,"*** Error: Detected unsupported boundary condition type (%d).\n",
                    (int)iod.mesh.bc_x0);
            exit(-1);
          }
        }


        //-----------
        // RIGHT 
        //-----------
        // Calculate F at the correct location 
        F = v[k][j][i+1][1]*dydz;
        // Calculate D and the "a" coefficient
        a = std::max(-F, 0.0);
        nustar = (nu + fn*vturb[k][j][i]) / sigma;
        if(nustar>0.0) {
          D  = nustar*dydz/dxr;
          a += D*PowerLaw(F/D);
        }
        ap += a;
        // Add entry (or add to bb or ap, if i+1 is outside boundary)
        if(i+1<NX)
          row.PushEntry(i+1,j,k, -a);  //on the left hand side
        else { //i+1 is outside domain boundary
          if(iod.mesh.bc_xmax == MeshData::INLET || iod.mesh.bc_xmax == MeshData::INLET2 ||
             iod.mesh.bc_xmax == MeshData::OVERSET)
            bb[k][j][i] += a*vturb[k][j][i+1];
          else if(iod.mesh.bc_xmax == MeshData::SYMMETRY || iod.mesh.bc_xmax == MeshData::OUTLET) {
            ap -= a;
          } else if(iod.mesh.bc_xmax == MeshData::STICKWALL) {
            ap += a; //bb[k][j][i] -= a*vturb[k][j][i]; //vturb[k][j][i+1] should be -vturb[k][j][i]
          } else {
            fprintf(stdout,"*** Error: Detected unsupported boundary condition type (%d).\n",
                    (int)iod.mesh.bc_xmax);
            exit(-1);
          }
        }


        //-----------
        // BOTTOM 
        //-----------
        // Calculate F at the correct location
        F = v[k][j][i][2]*dxdz;
        // Calculate D and the "a" coefficient
        a = std::max(F, 0.0);
        nustar = (nu + fn*vturb[k][j][i]) / sigma;
        if(nustar>0.0) {
          D  = nustar*dxdz/dyb;
          a += D*PowerLaw(F/D);
        }
        ap += a;
        // Add entry (or add to bb or ap, if j-1 is outside boundary)
        if(j-1>=0)
          row.PushEntry(i,j-1,k, -a);  //on the left hand side
        else { //j-1 is outside domain boundary 
          if(flat_plate_wall) //stick wall
            ap += a; 
          else if(iod.mesh.bc_y0 == MeshData::INLET || iod.mesh.bc_y0 == MeshData::INLET2 ||
                  iod.mesh.bc_y0 == MeshData::OVERSET)
            bb[k][j][i] += a*vturb[k][j-1][i];
          else if(iod.mesh.bc_y0 == MeshData::SYMMETRY || iod.mesh.bc_y0 == MeshData::OUTLET)
            ap -= a;
          else if(iod.mesh.bc_y0 == MeshData::STICKWALL)
            ap += a; //bb[k][j][i] -= a*vturb[k][j][i]; 
          else {
            fprintf(stdout,"*** Error: Detected unsupported boundary condition type (%d).\n",
                    (int)iod.mesh.bc_y0);
            exit(-1);
          }
        }


        //-----------
        // TOP 
        //-----------
        // Calculate F at the correct location 
        F = v[k][j+1][i][2]*dxdz;
        // Calculate D and the "a" coefficient
        a = std::max(-F, 0.0);
        nustar = (nu + fn*vturb[k][j][i]) / sigma;
        if(nustar>0.0) {
          D  = nustar*dxdz/dyt;
          a += D*PowerLaw(F/D);
        }
        ap += a;
        // Add entry (or add to bb or ap, if j+1 is outside boundary)
        if(j+1<NY)
          row.PushEntry(i,j+1,k, -a);  //on the left hand side
        else { //j+1 is outside domain boundary
          if(iod.mesh.bc_ymax == MeshData::INLET || iod.mesh.bc_ymax == MeshData::INLET2 ||
             iod.mesh.bc_ymax == MeshData::OVERSET)
            bb[k][j][i] += a*vturb[k][j+1][i];
          else if(iod.mesh.bc_ymax == MeshData::SYMMETRY || iod.mesh.bc_ymax == MeshData::OUTLET) {
            ap -= a;
          } else if(iod.mesh.bc_ymax == MeshData::STICKWALL) {
            bb[k][j][i] -= a*vturb[k][j][i];
          } else {
            fprintf(stdout,"*** Error: Detected unsupported boundary condition type (%d).\n",
                    (int)iod.mesh.bc_ymax);
            exit(-1);
          }
        }
   
    
        //-----------
        // BACK
        //-----------
        // Calculate F at the correct location
        F = v[k][j][i][3]*dxdy;
        // Calculate D and the "a" coefficient
        a = std::max(F, 0.0);
        nustar = (nu + fn*vturb[k][j][i]) / sigma;
        if(nustar>0.0) {
          D  = nustar*dxdy/dzk;
          a += D*PowerLaw(F/D);
        }
        ap += a;
        // Add entry (or add to bb or ap, if k-1 is outside boundary)
        if(k-1>=0)
          row.PushEntry(i,j,k-1, -a);  //on the left hand side
        else { //k-1 is outside domain boundary 
          if(iod.mesh.bc_z0 == MeshData::INLET || iod.mesh.bc_z0 == MeshData::INLET2 ||
             iod.mesh.bc_z0 == MeshData::OVERSET)
            bb[k][j][i] += a*vturb[k-1][j][i]; 
          else if(iod.mesh.bc_z0 == MeshData::SYMMETRY || iod.mesh.bc_z0 == MeshData::OUTLET)
            ap -= a;
          else if(iod.mesh.bc_z0 == MeshData::STICKWALL)
            bb[k][j][i] -= a*vturb[k][j][i];
          else {
            fprintf(stdout,"*** Error: Detected unsupported boundary condition type (%d).\n",
                    (int)iod.mesh.bc_z0);
            exit(-1);
          }
        }
              

        //-----------
        // FRONT 
        //-----------
        // Calculate F at the correct location
        F = v[k+1][j][i][3]*dxdy;
        // Calculate D and the "a" coefficient
        a = std::max(-F, 0.0);
        nustar = (nu + fn*vturb[k][j][i]) / sigma;
        if(nustar>0.0) {
          D  = nustar*dxdy/dzf;
          a += D*PowerLaw(F/D);
        }
        ap += a;
        // Add entry (or add to bb or ap, if j+1 is outside boundary)
        if(k+1<NZ)
          row.PushEntry(i,j,k+1, -a);  //on the left hand side
        else { //k+1 is outside domain boundary
          if(iod.mesh.bc_zmax == MeshData::INLET || iod.mesh.bc_zmax == MeshData::INLET2 ||
             iod.mesh.bc_zmax == MeshData::OVERSET)
            bb[k][j][i] += a*vturb[k+1][j][i];
          else if(iod.mesh.bc_zmax == MeshData::SYMMETRY || iod.mesh.bc_zmax == MeshData::OUTLET) {
            ap -= a;
          } else if(iod.mesh.bc_zmax == MeshData::STICKWALL) {
            bb[k][j][i] -= a*vturb[k][j][i];
          } else {
            fprintf(stdout,"*** Error: Detected unsupported boundary condition type (%d).\n",
                    (int)iod.mesh.bc_zmax);
            exit(-1);
          }
        }
   
  
        //------------------------------------------------------
        // Calculate and add the diagonal entry and the RHS (b)
        // Ref: Eqs. (5.62) and (6.8) in Patankar's book
        //------------------------------------------------------
        ap0 = LocalDt ? dxdy*dz/dtloc[k][j][i] : dxdy*dz/dt;
        ap += ap0; //!< -Sp*dx*dy*dz, for source terms

        bb[k][j][i] += ap0*vturb0[k][j][i]; 


        //------------------------------------------------------
        // Addition of Source Terms
        // RHS of Spalart-Allmaras Model
        //------------------------------------------------------
        double Sc = 0.0; 
        //production term
        Sc += neg_nut ? cb1*(1-ct3)*Omega*vturb[k][j][i] : cb1*(1.0-ft2)*S*vturb[k][j][i];
        //destruction term
        Sc += neg_nut ? cw1*pow(vturb[k][j][i]/d2w,2) : -(cw1*fw - cb1/(kappa*kappa)*ft2)*pow(vturb[k][j][i]/d2w,2); //TODO: THIS IS A HACK

        //nonlinear diffusion term
        dnut[0] = -dxr/(dxl*(dxl+dxr))*vturb[k][j][i-1] + (dxr-dxl)/(dxl*dxr)*vturb[k][j][i]
                +  dxl/(dxr*(dxl+dxr))*vturb[k][j][i+1]; //2nd-order accuracy, see Kevin's notes
        dnut[1] = -dyt/(dyb*(dyt+dyb))*vturb[k][j-1][i] + (dyt-dyb)/(dyt*dyb)*vturb[k][j][i]
                +  dyb/(dyt*(dyt+dyb))*vturb[k][j+1][i];
        dnut[2] = -dzf/(dzk*(dzf+dzk))*vturb[k-1][j][i] + (dzf-dzk)/(dzf*dzk)*vturb[k][j][i]
                +  dzk/(dzf*(dzf+dzk))*vturb[k+1][j][i];
        Sc += cb2/sigma*(dnut*dnut);
        bb[k][j][i] += Sc*dxdy*dz; //multiply source term by cell volume;

        // Apply relaxation (Ref: Eq.(6) of Van Doormaal and Rathby, 1984)
        assert(Efactor>0.0);
        bb[k][j][i] += ap*vturb0[k][j][i]/Efactor;

        ap *= 1.0 + 1.0/Efactor; 
        row.PushEntry(i,j,k, ap);

      }
    }
  }
  

  B.RestoreDataPointerAndInsert();

  if(LocalDt)
    LocalDt->RestoreDataPointerToLocalVector();
}

//----------------------------------------------------------------------------

void
IncompressibleOperator::InitializeTurbulenceVariables(SpaceVariable3D &Vturb)
{
  assert(iod.rans.model == RANSTurbulenceModelData::SPALART_ALLMARAS &&
         iod.rans.example == RANSTurbulenceModelData::FLAT_PLATE); //this is all we can do right now

  //Vturb.SetConstantValue(0.0, true); //set Vturb = 0 (including ghost layer)
  Vturb.SetConstantValue(iod.rans.nu_tilde_farfield, true); //set Vturb = 0 (including ghost layer)

  ApplyBoundaryConditionsTurbulenceVariables(Vturb);
}

//--------------------------------------------------------------------------

void
IncompressibleOperator::ApplyBoundaryConditionsTurbulenceVariables(SpaceVariable3D &Vturb)
{
  double*** vturb = Vturb.GetDataPointer();

  vector<GhostPoint>* ghost_nodes_outer(spo.GetPointerToOuterGhostNodes());

  double nu_tilde_far = iod.rans.nu_tilde_farfield;
  assert(nu_tilde_far>=0.0);

  for(auto it = ghost_nodes_outer->begin(); it != ghost_nodes_outer->end(); it++) {

    if(it->type_projection != GhostPoint::FACE)
      continue; //corner (edge or vertex) nodes are not populated

    // preparation
    int i(it->ijk[0]), j(it->ijk[1]), k(it->ijk[2]);
    int im_i(it->image_ijk[0]), im_j(it->image_ijk[1]), im_k(it->image_ijk[2]);

    switch (it->bcType) {
      case MeshData::INLET :
      case MeshData::INLET2 :
        vturb[k][j][i] = nu_tilde_far;
        break;
      case MeshData::OUTLET :
      case MeshData::SYMMETRY :
        vturb[k][j][i] = vturb[im_k][im_j][im_i];
        break;
      case MeshData::STICKWALL :
        vturb[k][j][i] = -vturb[im_k][im_j][im_i];
        break;
      case MeshData::OVERSET :
        break; // nothing to be done here.
      default:
        fprintf(stdout,"*** Error: Detected unknown boundary condition type (%d).\n", (int)it->bcType);
        exit(-1);
    }
  }

  Vturb.RestoreDataPointerAndInsert();
}

//--------------------------------------------------------------------------

double
IncompressibleOperator::GetKinematicEddyViscosity(double rho, double mu, double nu_tilde)
{
  // Spalart-Allmaras for the moment
  if(nu_tilde<0)
    return 0.0; //if set nut to 0.0

  double cv1_3 = 357.911; //=7.1*7.1*7.1

  assert(rho>0 && mu>0);
  double chi  = rho*nu_tilde/mu;
  double chi3 = chi*chi*chi;
  double fv1  = chi3/(chi3 + cv1_3);

  return nu_tilde*fv1;
}

//--------------------------------------------------------------------------

double
IncompressibleOperator::GetDynamicEddyViscosity(double rho, double mu, double nu_tilde)
{
  // Spalart-Allmaras for the moment
  if(nu_tilde<0)
    return 0.0; //if set nut to 0.0

  double cv1_3 = 357.911; //=7.1*7.1*7.1

  assert(rho>0 && mu>0);
  double chi  = rho*nu_tilde/mu;
  double chi3 = chi*chi*chi;
  double fv1  = chi3/(chi3 + cv1_3);

  return rho*nu_tilde*fv1;
}

//--------------------------------------------------------------------------

double
IncompressibleOperator::GetDistanceToWall(Vec3D x)
{
  GlobalMeshInfo& global_mesh(spo.GetGlobalMeshInfo());

  if(iod.rans.example == RANSTurbulenceModelData::FLAT_PLATE) {
    double h = x[1] - global_mesh.GetYmin();
    assert(h>=0.0);
    double d = x[0] - iod.rans.example_param_1;
    if(d>=0.0)
      return h;
    else
      return sqrt(h*h + d*d);
  }

  double dist = DBL_MAX;
  if(iod.mesh.bc_x0 == MeshData::STICKWALL)
    dist = std::min(dist, x[0] - global_mesh.GetXmin());
  if(iod.mesh.bc_xmax == MeshData::STICKWALL)
    dist = std::min(dist, global_mesh.GetXmax() - x[0]);
  if(iod.mesh.bc_y0 == MeshData::STICKWALL)
    dist = std::min(dist, x[1] - global_mesh.GetYmin());
  if(iod.mesh.bc_ymax == MeshData::STICKWALL)
    dist = std::min(dist, global_mesh.GetYmax() - x[1]);
  if(iod.mesh.bc_z0 == MeshData::STICKWALL)
    dist = std::min(dist, x[2] - global_mesh.GetZmin());
  if(iod.mesh.bc_zmax == MeshData::STICKWALL)
    dist = std::min(dist, global_mesh.GetZmax() - x[2]);

  if(dist == DBL_MAX) {
    fprintf(stdout,"\033[0;31m*** Error: Unable to find a stick wall --- needed in turbulence model.\033[0m\n");
    exit(-1);
  }

  assert(dist>=0.0);
  return dist;
}    

//--------------------------------------------------------------------------

void
IncompressibleOperator::ComputeKinematicEddyViscosity(SpaceVariable3D &Vturb, SpaceVariable3D &V,
                                                      SpaceVariable3D &ID, SpaceVariable3D &NuT)
{
  double*** vturb = Vturb.GetDataPointer();
  Vec5D***  v     = (Vec5D***) V.GetDataPointer();
  double*** id     = ID.GetDataPointer();
  double*** nut    = NuT.GetDataPointer();

  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++)
        nut[k][j][i] = GetKinematicEddyViscosity(v[k][j][i][0], Mu[id[k][j][i]], vturb[k][j][i]);

  Vturb.RestoreDataPointerToLocalVector();
  V.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();

  NuT.RestoreDataPointerAndInsert();
}

//--------------------------------------------------------------------------

double
IncompressibleOperator::ApplyVelocityBoundaryConditionLocal(int dir, GhostPoint::Side side, double vim,
                                                            bool flat_plate_wall)
{
  // This function calculates velocity at ghost cell CENTER; and "vim" is the velocity at the image cell CENTER,
  // regardless of dir = 0,1,2.

  if(flat_plate_wall) //stick wall
    return -vim;

  double vghost(DBL_MAX);
  int wall_dir(0);
  MeshData::BcType type(MeshData::NONE);
  switch (side) {
    case GhostPoint::LEFT :
      type = iod.mesh.bc_x0;
      wall_dir = 0;
      break;
    case GhostPoint::RIGHT :
      type = iod.mesh.bc_xmax;
      wall_dir = 0;
      break;
    case GhostPoint::BOTTOM :
      type = iod.mesh.bc_y0;
      wall_dir = 1;
      break;
    case GhostPoint::TOP :
      type = iod.mesh.bc_ymax;
      wall_dir = 1;
      break;
    case GhostPoint::BACK :
      type = iod.mesh.bc_z0;
      wall_dir = 2;
      break;
    case GhostPoint::FRONT :
      type = iod.mesh.bc_zmax;
      wall_dir = 2;
      break;
    default :
      fprintf(stdout,"*** Error: Detected unsupported boundary code.\n");
      exit(-1);
  }

  switch (type) {
    case MeshData::INLET :
        vghost = dir==0 ? iod.bc.inlet.velocity_x : (dir==1 ? iod.bc.inlet.velocity_y : iod.bc.inlet.velocity_z);
        break;
      case MeshData::INLET2 :
        vghost = dir==0 ? iod.bc.inlet2.velocity_x : (dir==1 ? iod.bc.inlet2.velocity_y : iod.bc.inlet2.velocity_z);
        break;
      case MeshData::OUTLET :
        vghost = vim;
        break;
      case MeshData::SLIPWALL :
      case MeshData::SYMMETRY :
        if(dir==wall_dir) vghost = -vim; 
        else              vghost =  vim;
        break;
      case MeshData::STICKWALL :
        vghost = -vim;
        break;
      case MeshData::OVERSET :
        //nothing to be done here. TODO
        break;
      default :
        fprintf(stdout,"*** Error: Detected unsupported boundary condition type (%d).\n", (int)iod.mesh.bc_x0); 
        exit(-1);
  }

  return vghost;
}

//--------------------------------------------------------------------------

//--------------------------------------------------------------------------



