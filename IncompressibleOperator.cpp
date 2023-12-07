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
                        V3(comm_, &(dm_all_.ghosted1_3dof)), viscous(false)
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
      else if(visFcn[i]->type != ViscoFcnBase::NONE) {
        print_error("*** Error: Currently, IncompressibleOperator only supports CONSTANT viscosity.\n");
        exit_mpi();
      }
    }
  }

  // Calculate deltas
  Dx.resize(3); //dir = 0, 1, 2
  Dy.resize(3);
  Dz.resize(3);
  dx_l.resize(3);
  dx_r.resize(3);
  dy_b.resize(3);
  dy_t.resize(3);
  dz_k.resize(3);
  dz_f.resize(3);
  for(d=0; d<3; d++) {
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
        v[k][j][i][1] = cx0[i]*v3[k][j][i-1][0] + cx1[i]*v3[k][j][i][0];
        v[k][j][i][2] = cy0[j]*v3[k][j-1][i][1] + cy1[j]*v3[k][j][i][1];
        v[k][j][i][3] = cz0[k]*v3[k-1][j][i][2] + cz1[k]*v3[k][j][i][2];
        v[k][j][i][4] = 0.0; //initialize p = 0
      }

  //Note: At this point some velocity values ON THE DOMAIN BOUNDARIES ARE INCORRECT.
  //      They NEED TO BE CORRECTED by calling ApplyBoundaryConditions.

  V.RestoreDataPointerAndInsert();
  V3.RestoreDataPointerToLocalVector(); //no need to exchange data
  ID.RestoreDataPointerToLocalVector();
}

//--------------------------------------------------------------------------

void
IncompressibleOperator::ApplyBoundaryConditions(SpaceVariable3D &V)
{
  // Note: Modify both ghost and non-ghost entries due to the use of MAC / staggered grids  
  //       Only modify velocity components

  Vec5D*** v = (Vec5D***) V.GetDataPointer();

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
    else if(it->bcType == MeshData::OUTLET) {
      v0[0] = iod.bc.outlet.velocity_x;
      v0[1] = iod.bc.outlet.velocity_y;
      v0[2] = iod.bc.outlet.velocity_z;
    }


    // Impose boundary conditions for velocity
    // Treating the six sides separately might seem clumsy. But in some cases imposing a boundary
    // condition requires populating two entries. 
    if(it->side == GhostPoint::LEFT) {
      if(it->bcType == MeshData::INLET || it->bcType == MeshData::OUTLET) {
        v[k][j][im_i][1] = v0[0];
        v[k][j][i][2]    = v0[1];
        v[k][j][i][3]    = v0[2];
      }
      else if(it->bcType == MeshData::SLIPWALL || it->bcType == MeshData::SYMMETRY) {
        v[k][j][im_i][1] = 0.0;
        v[k][j][i][2]    = v[k][j][im_i][2];
        v[k][j][i][3]    = v[k][j][im_i][3];
      }
      else if(it->bcType == MeshData::STICKWALL) {
        v[k][j][im_i][0] = 0.0;
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
      if(it->bcType == MeshData::INLET || it->bcType == MeshData::OUTLET) {
        v[k][j][i][1]    = v0[0];
        v[k][j][i][2]    = v0[1];
        v[k][j][i][3]    = v0[2];
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
      if(it->bcType == MeshData::INLET || it->bcType == MeshData::OUTLET) {
        v[k][j][i][1]    = v0[0];
        v[k][im_j][i][2] = v0[1];
        v[k][j][i][3]    = v0[2];
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
    }
    else if(it->side == GhostPoint::TOP) {
      if(it->bcType == MeshData::INLET || it->bcType == MeshData::OUTLET) {
        v[k][j][i][1]    = v0[0];
        v[k][j][i][2]    = v0[1];
        v[k][j][i][3]    = v0[2];
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
      if(it->bcType == MeshData::INLET || it->bcType == MeshData::OUTLET) {
        v[k][j][i][1]    = v0[0];
        v[k][i][i][2]    = v0[1];
        v[im_k][j][i][3] = v0[2];
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
      if(it->bcType == MeshData::INLET || it->bcType == MeshData::OUTLET) {
        v[k][j][i][1]    = v0[0];
        v[k][j][i][2]    = v0[1];
        v[k][j][i][3]    = v0[2];
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
    if (iod.mesh.bc_x0 == MeshData::INLET    || iod.mesh.bc_x0 == MeshData::OUTLET ||
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
    if (iod.mesh.bc_xmax == MeshData::INLET    || iod.mesh.bc_xmax == MeshData::OUTLET ||
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
    if (iod.mesh.bc_y0 == MeshData::INLET    || iod.mesh.bc_y0 == MeshData::OUTLET ||
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
    if (iod.mesh.bc_ymax == MeshData::INLET    || iod.mesh.bc_ymax == MeshData::OUTLET ||
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
    if (iod.mesh.bc_z0 == MeshData::INLET    || iod.mesh.bc_z0 == MeshData::OUTLET ||
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
    if (iod.mesh.bc_zmax == MeshData::INLET    || iod.mesh.bc_zmax == MeshData::OUTLET ||
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
        vel_over_dx_max = max(vel_over_dx_max,
                              max(fabs(v[k][j][i][1])/dx, 
                                  max(fabs(v[k][j][i][2])/dy, fabs(v[k][j][i][3])/dz)));
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &vel_over_dx_max, 1, MPI_DOUBLE, MPI_MAX, comm);
  assert(vel_over_dx_max>0.0);

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

  double vel_over_dx_max = 0.0, vel_over_dx_max_loc;
  double dx, dy, dz;

  // Loop through the real domain (excluding the ghost layer)
  for(int k=k0; k<kmax; k++) {
    dz = global_mesh.GetDz(k);
    for(int j=j0; j<jmax; j++) {
      dy = global_mesh.GetDy(j);
      for(int i=i0; i<imax; i++) {
        if(id[k][j][i] == INACTIVE_MATERIAL_ID)
          continue;
        dx = global_mesh.GetDx(i);
        
        vel_over_dx_max_loc = max(fabs(v[k][j][i][1])/dx, 
                                  max(fabs(v[k][j][i][2])/dy, fabs(v[k][j][i][3])/dz));
        vel_over_dx_max = max(vel_over_dx_max, vel_over_dx_max_loc);

        // calculates local time-step size
        dtl[k][j][i] = cfl/vel_over_dx_max_loc;
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &vel_over_dx_max, 1, MPI_DOUBLE, MPI_MAX, comm);
  assert(vel_over_dx_max>0.0);

  // global min dt
  dt = cfl/vel_over_dx_max;

  V.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();

  LocalDt.RestoreDataPointerAndInsert();

}

//--------------------------------------------------------------------------

void
IncompressibleOperator::BuildVelocityEquationSIMPLE(int dir, Vec5D*** v, double*** id,
                                                    vector<RowEntries> &vlin_rows, SpaceVariable3D &B,
                                                    SemiImplicitTsData &iod_semi)
{

  assert(dir==0 || dir==1 || dir==2);

  GlobalMeshInfo& global_mesh(spo.GetGlobalMeshInfo());

  double*** bb = B.GetDataPointer();

  int row_counter = 0;
  int original_size = vlin_rows.size();
  int myid;

  double a_p, a_l, a_r, a_b, a_t, a_k, a_f, b, dp;
  double Dx, Dy, Dz, dx_l, dx_r, dy_b, dy_t, dz_k, dz_f;
  double mu;

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

        // locate the row
        if(row_counter>=original_size)
          vlin_rows.push_back(RowEntries(7)); // at most 7 non-zero entries on each row
        RowEntries &row(vlin_rows[row_counter]);
        row.SetRow(i,j,k);
        row.ClearEntries();
        row_counter++;

        if((dir==0 && i==0) || (dir==1 && j==0) || (dir==2 && k==0)) {
          row.PushEntry(i,j,k, 1.0); 
          bb[k][j][i] = v[k][j][i][1+dir];
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

        //-----------
        // LEFT
        //-----------
        // Calculate F at the correct location (different for dir=1,2,3)
        if(dir==0)
          F = v[k][j][i-1][0]*0.5*(v[k][j][i-1][1] + v[k][j][i][1]);
        else {
          cm = global_mesh.GetDx(i-1);
          cp = global_mesh.GetDx(i);
          cm_plus_cp = cm + cp;
          rhou1 = dir==1 ? v[k][j-1][i][1]*(cp*v[k][j-1][i-1][0] + cm*v[k][j-1][i][0])/cm_plus_cp
                         : v[k-1][j][i][1]*(cp*v[k-1][j][i-1][0] + cm*v[k-1][j][i][0])/cm_plus_cp;
          rhou2 = v[k][j][i][1]*(cp*v[k][j][i-1][0] + cm*v[k][j][i][0])/(cm+cp);
          F = dir== 1 ? (dyt*rhou1 + dyb*rhou2)/(dyb+dyt) : (dzf*rhou1 + dzk*rhou2)/(dzk+dzf);
        }
        F *= dydz;
        // Calculate D and the "a" coefficient
        a = std::max(F, 0.0);
        if(dir==0)
          mu = Mu[id[k][j][i-1]];
        else {
          // cm and cp has been calculated!   
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
        if(mu>0.0) {
          D  = mu*dydz/dxl;
          a += D*PowerLaw(F/D);
        }
        ap += a;
        // Add entry
        row.PushEntry(i-1,j,k, -a);  //on the left hand side
             

        //-----------
        // RIGHT 
        //-----------
        // Calculate F at the correct location (different for dir=1,2,3)
        if(dir==0)
          F = v[k][j][i][0]*0.5*(v[k][j][i][1] + v[k][j][i+1][1]);
        else {
          cm = global_mesh.GetDx(i);
          cp = global_mesh.GetDx(i+1);
          cm_plus_cp = cm + cp;
          rhou1 = dir==1 ? v[k][j-1][i+1][1]*(cp*v[k][j-1][i][0] + cm*v[k][j-1][i+1][0])/cm_plus_cp
                         : v[k-1][j][i+1][1]*(cp*v[k-1][j][i][0] + cm*v[k-1][j][i+1][0])/cm_plus_cp;
          rhou2 = v[k][j][i+1][1]*(cp*v[k][j][i][0] + cm*v[k][j][i+1][0])/(cm+cp);
          F = dir== 1 ? (dyt*rhou1 + dyb*rhou2)/(dyb+dyt) : (dzf*rhou1 + dzk*rhou2)/(dzk+dzf);
        }
        F *= dydz;
        // Calculate D and the "a" coefficient
        a = std::max(-F, 0.0);
        if(dir==0)
          mu = Mu[id[k][j][i]];
        else {
          // cm and cp has been calculated!   
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
        if(mu>0.0) {
          D  = mu*dydz/dxr;
          a += D*PowerLaw(F/D);
        }
        ap += a;
        // Add entry
        row.PushEntry(i+1,j,k, -a);  //on the left hand side
             
    
        //-----------
        // BOTTOM 
        //-----------
        // Calculate F at the correct location (different for dir=1,2,3)
        if(dir==1)
          F = v[k][j-1][i][0]*0.5*(v[k][j-1][i][2] + v[k][j][i][2]);
        else {
          cm = global_mesh.GetDy(j-1);
          cp = global_mesh.GetDy(j);
          cm_plus_cp = cm + cp;
          rhou1 = dir==0 ? v[k][j][i-1][2]*(cp*v[k][j-1][i-1][0] + cm*v[k][j][i-1][0])/cm_plus_cp
                         : v[k-1][j][i][2]*(cp*v[k-1][j-1][i][0] + cm*v[k-1][j][i][0])/cm_plus_cp;
          rhou2 = v[k][j][i][2]*(cp*v[k][j-1][i][0] + cm*v[k][j][i][0])/(cm+cp);
          F = dir== 0 ? (dxr*rhou1 + dxl*rhou2)/(dxl+dxr) : (dzf*rhou1 + dzk*rhou2)/(dzk+dzf);
        }
        F *= dxdz;
        // Calculate D and the "a" coefficient
        a = std::max(F, 0.0);
        if(dir==1)
          mu = Mu[id[k][j-1][i]];
        else {
          // cm and cp has been calculated!   
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
        if(mu>0.0) {
          D  = mu*dxdz/dyb;
          a += D*PowerLaw(F/D);
        }
        ap += a;
        // Add entry
        row.PushEntry(i,j-1,k, -a);  //on the left hand side
             

        //-----------
        // TOP 
        //-----------
        // Calculate F at the correct location (different for dir=1,2,3)
        if(dir==1)
          F = v[k][j][i][0]*0.5*(v[k][j][i][2] + v[k][j+1][i][2]);
        else {
          cm = global_mesh.GetDy(j);
          cp = global_mesh.GetDy(j+1);
          cm_plus_cp = cm + cp;
          rhou1 = dir==0 ? v[k][j+1][i-1][2]*(cp*v[k][j][i-1][0] + cm*v[k][j+1][i-1][0])/cm_plus_cp
                         : v[k-1][j+1][i][2]*(cp*v[k-1][j][i][0] + cm*v[k-1][j+1][i][0])/cm_plus_cp;
          rhou2 = v[k][j+1][i][2]*(cp*v[k][j][i][0] + cm*v[k][j+1][i][0])/(cm+cp);
          F = dir== 0 ? (dxr*rhou1 + dxl*rhou2)/(dxl+dxr) : (dzf*rhou1 + dzk*rhou2)/(dzk+dzf);
        }
        F *= dxdz;
        // Calculate D and the "a" coefficient
        a = std::max(-F, 0.0);
        if(dir==1)
          mu = Mu[id[k][j][i]];
        else {
          // cm and cp has been calculated!   
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
        if(mu>0.0) {
          D  = mu*dxdz/dyt;
          a += D*PowerLaw(F/D);
        }
        ap += a;
        // Add entry
        row.PushEntry(i,j+1,k, -a);  //on the left hand side
             
    
        //-----------
        // BACK 
        //-----------
        // Calculate F at the correct location (different for dir=1,2,3)
        if(dir==2)
          F = v[k-1][j][i][0]*0.5*(v[k-1][j][i][3] + v[k][j][i][3]);
        else {
          cm = global_mesh.GetDz(k-1);
          cp = global_mesh.GetDy(k);
          cm_plus_cp = cm + cp;
          rhou1 = dir==0 ? v[k][j][i-1][3]*(cp*v[k-1][j][i-1][0] + cm*v[k][j][i-1][0])/cm_plus_cp
                         : v[k][j-1][i][3]*(cp*v[k-1][j-1][i][0] + cm*v[k][j-1][i][0])/cm_plus_cp;
          rhou2 = v[k][j][i][3]*(cp*v[k-1][j][i][0] + cm*v[k][j][i][0])/(cm+cp);
          F = dir== 0 ? (dxr*rhou1 + dxl*rhou2)/(dxl+dxr) : (dyt*rhou1 + dyb*rhou2)/(dyb+dyt);
        }
        F *= dxdy;
        // Calculate D and the "a" coefficient
        a = std::max(F, 0.0);
        if(dir==2)
          mu = Mu[id[k-1][j][i]];
        else {
          // cm and cp has been calculated!   
          if(k==0) {
            mu1 = dir==0 ? Mu[id[k-1][j][i]] : Mu[id[k][j][i-1]];
            mu2 = Mu[id[k][j][i]];
I AM HERE
          } else {
            mu1 = dir==0 ? (cp*Mu[id[k][j-1][i-1]] + cm*Mu[id[k][j][i-1]])/cm_plus_cp
                         : (cp*Mu[id[k-1][j-1][i]] + cm*Mu[id[k-1][j][i]])/cm_plus_cp;
            mu2 = (cp*Mu[id[k][j-1][i]] + cm*Mu[id[k][j][i]])/cm_plus_cp;
          }
          mu = dir==0 ? (dxr*mu1 + dxl*mu2)/(dxl+dxr) : (dzf*mu1 + dzk*mu2)/(dzk+dzf);
        }
        if(mu>0.0) {
          D  = mu*dxdz/dyb;
          a += D*PowerLaw(F/D);
        }
        ap += a;
        // Add entry
        row.PushEntry(i,j-1,k, -a);  //on the left hand side
             

        //-----------
        // TOP 
        //-----------
        // Calculate F at the correct location (different for dir=1,2,3)
        if(dir==1)
          F = v[k][j][i][0]*0.5*(v[k][j][i][2] + v[k][j+1][i][2]);
        else {
          cm = global_mesh.GetDy(j);
          cp = global_mesh.GetDy(j+1);
          cm_plus_cp = cm + cp;
          rhou1 = dir==0 ? v[k][j+1][i-1][2]*(cp*v[k][j][i-1][0] + cm*v[k][j+1][i-1][0])/cm_plus_cp
                         : v[k-1][j+1][i][2]*(cp*v[k-1][j][i][0] + cm*v[k-1][j+1][i][0])/cm_plus_cp;
          rhou2 = v[k][j+1][i][2]*(cp*v[k][j][i][0] + cm*v[k][j+1][i][0])/(cm+cp);
          F = dir== 0 ? (dxr*rhou1 + dxl*rhou2)/(dxl+dxr) : (dzf*rhou1 + dzk*rhou2)/(dzk+dzf);
        }
        F *= dxdz;
        // Calculate D and the "a" coefficient
        a = std::max(-F, 0.0);
        if(dir==1)
          mu = Mu[id[k][j][i]];
        else {
          // cm and cp has been calculated!   
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
        if(mu>0.0) {
          D  = mu*dxdz/dyt;
          a += D*PowerLaw(F/D);
        }
        ap += a;
        // Add entry
        row.PushEntry(i,j+1,k, -a);  //on the left hand side
              
 
         
      }
    }
  }


  B.RestoreDataPointerAndInsert();
}

//--------------------------------------------------------------------------




//--------------------------------------------------------------------------



