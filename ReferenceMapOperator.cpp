/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<ReferenceMapOperator.h>
#include<GradientCalculatorFD3.h>
#include<GeoTools.h>
#include<Vector5D.h>

//-------------------------------------------------------------------

ReferenceMapOperator::ReferenceMapOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_,
                                           SpaceVariable3D &coordinates_,
                                           SpaceVariable3D &delta_xyz_, GlobalMeshInfo &global_mesh_,
                                           std::vector<GhostPoint> &ghost_nodes_inner_,
                                           std::vector<GhostPoint> &ghost_nodes_outer_)
                    : comm(comm_), iod(iod_), coordinates(coordinates_), global_mesh(global_mesh_),
                      ghost_nodes_inner(ghost_nodes_inner_), ghost_nodes_outer(ghost_nodes_outer_),
                      vectorG2(comm_, &(dm_all_.ghosted2_3dof)),    
                      Xil(comm_, &(dm_all_.ghosted1_3dof)),
                      Xir(comm_, &(dm_all_.ghosted1_3dof)),
                      Xib(comm_, &(dm_all_.ghosted1_3dof)),
                      Xit(comm_, &(dm_all_.ghosted1_3dof)),
                      Xik(comm_, &(dm_all_.ghosted1_3dof)),
                      Xif(comm_, &(dm_all_.ghosted1_3dof))
{

  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);

  //set up gradient calculator
  if(iod.refmap.fd == ReferenceMapData::UPWIND_CENTRAL_3) { //currently the only option
    grad_minus = new GradientCalculatorFD3(comm, dm_all_, coordinates, delta_xyz_, -1);
    grad_plus  = new GradientCalculatorFD3(comm, dm_all_, coordinates, delta_xyz_,  1);
  }

  TagExternalGhostNodes();
}

//-------------------------------------------------------------------

ReferenceMapOperator::~ReferenceMapOperator()
{
  if(grad_minus) delete grad_minus;
  if(grad_plus)  delete grad_plus;
}

//-------------------------------------------------------------------

void
ReferenceMapOperator::Destroy()
{
  vectorG2.Destroy();
  Xil.Destroy();
  Xir.Destroy();
  Xib.Destroy();
  Xit.Destroy();
  Xik.Destroy();
  Xif.Destroy();
}

//-------------------------------------------------------------------

void
ReferenceMapOperator::TagExternalGhostNodes()
{
  //tag ghost nodes outside open boundaries, which require special treatment when imposing b.c.
  ghost_nodes_outer_tag.assign(ghost_nodes_outer.size(), -1);

  map<int, DiskData* >& disks(iod.bc.multiBoundaryConditions.diskMap.dataMap);
  map<int, RectangleData* >& rectangles(iod.bc.multiBoundaryConditions.rectangleMap.dataMap);

  for(auto it = ghost_nodes_outer.begin(); it != ghost_nodes_outer.end();  it++) {

    if(it->type_projection != GhostPoint::FACE)
      continue; //corner (i.e. edge or vertex) nodes are not populated

    if(it->bcType == MeshData::SYMMETRY) {
      ghost_nodes_outer_tag[int(it-ghost_nodes_outer.begin())] = 0;
    }
    else if(it->bcType == MeshData::INLET || it->bcType == MeshData::INLET2 || it->bcType == MeshData::OUTLET ||
            it->bcType == MeshData::OVERSET) {
      ghost_nodes_outer_tag[int(it-ghost_nodes_outer.begin())] = 2;      
    }
    else if(it->bcType == MeshData::SLIPWALL || it->bcType == MeshData::STICKWALL) {
      
      ghost_nodes_outer_tag[int(it-ghost_nodes_outer.begin())] = (it->bcType == MeshData::SLIPWALL) ?
                                                                 0 : 1; 
      // Note that we allow user to make "windows" on a wall

      // Loop through user-specified disks (if any) and make correction if necessary
      for(auto&& disk : disks) {
        if((it->side == GhostPoint::LEFT  && disk.second->cen_x == iod.mesh.x0) ||
           (it->side == GhostPoint::RIGHT && disk.second->cen_x == iod.mesh.xmax)) {
          Vec3D n(disk.second->normal_x, disk.second->normal_y, disk.second->normal_z);
          if(fabs(n[0])/n.norm()>1-1e-8) {
            if(GeoTools::IsPointInDisk(global_mesh.GetY(it->ijk), global_mesh.GetZ(it->ijk),
                                       disk.second->cen_y, disk.second->cen_z, disk.second->radius)) {
              ghost_nodes_outer_tag[int(it-ghost_nodes_outer.begin())] = 2;
              break;
            }
          }
        }
        if((it->side == GhostPoint::BOTTOM && disk.second->cen_y == iod.mesh.y0) ||
           (it->side == GhostPoint::TOP    && disk.second->cen_y == iod.mesh.ymax)) {
          Vec3D n(disk.second->normal_x, disk.second->normal_y, disk.second->normal_z);
          if(fabs(n[1])/n.norm()>1-1e-8) {
            if(GeoTools::IsPointInDisk(global_mesh.GetZ(it->ijk), global_mesh.GetX(it->ijk),
                                       disk.second->cen_z, disk.second->cen_x, disk.second->radius)) {
              ghost_nodes_outer_tag[int(it-ghost_nodes_outer.begin())] = 2;
              break;
            }
          }
        }
        if((it->side == GhostPoint::BACK  && disk.second->cen_z == iod.mesh.z0) ||
           (it->side == GhostPoint::FRONT && disk.second->cen_z == iod.mesh.zmax)) {
          Vec3D n(disk.second->normal_x, disk.second->normal_y, disk.second->normal_z);
          if(fabs(n[2])/n.norm()>1-1e-8) {
            if(GeoTools::IsPointInDisk(global_mesh.GetX(it->ijk), global_mesh.GetY(it->ijk),
                                       disk.second->cen_x, disk.second->cen_y, disk.second->radius)) {
              ghost_nodes_outer_tag[int(it-ghost_nodes_outer.begin())] = 2;
              break;
            }
          }
        }
      } 

      if(ghost_nodes_outer_tag[int(it-ghost_nodes_outer.begin())] == 2) //already corrected
        continue;

      // Loop through user-specified rectangles (if any)
      for(auto&& rect : rectangles) {
        if((it->side == GhostPoint::LEFT  && rect.second->cen_x == iod.mesh.x0) ||
           (it->side == GhostPoint::RIGHT && rect.second->cen_x == iod.mesh.xmax)) {
          Vec3D n(rect.second->normal_x, rect.second->normal_y, rect.second->normal_z);
          if(fabs(n[0])/n.norm()>1-1e-8) {
            if(GeoTools::IsPointInRectangle(global_mesh.GetY(it->ijk), global_mesh.GetZ(it->ijk),
                                            rect.second->cen_y, rect.second->cen_z, 
                                            rect.second->a, rect.second->b)) {
              ghost_nodes_outer_tag[int(it-ghost_nodes_outer.begin())] = 2;
              break;
            }
          }
        }
        if((it->side == GhostPoint::BOTTOM && rect.second->cen_y == iod.mesh.y0) ||
           (it->side == GhostPoint::TOP    && rect.second->cen_y == iod.mesh.ymax)) {
          Vec3D n(rect.second->normal_x, rect.second->normal_y, rect.second->normal_z);
          if(fabs(n[1])/n.norm()>1-1e-8) {
            if(GeoTools::IsPointInRectangle(global_mesh.GetZ(it->ijk), global_mesh.GetX(it->ijk),
                                            rect.second->cen_z, rect.second->cen_x, 
                                            rect.second->a, rect.second->b)) {
              ghost_nodes_outer_tag[int(it-ghost_nodes_outer.begin())] = 2;
              break;
            }
          }
        }
        if((it->side == GhostPoint::BACK  && rect.second->cen_z == iod.mesh.z0) ||
           (it->side == GhostPoint::FRONT && rect.second->cen_z == iod.mesh.zmax)) {
          Vec3D n(rect.second->normal_x, rect.second->normal_y, rect.second->normal_z);
          if(fabs(n[2])/n.norm()>1-1e-8) {
            if(GeoTools::IsPointInRectangle(global_mesh.GetX(it->ijk), global_mesh.GetY(it->ijk),
                                            rect.second->cen_x, rect.second->cen_y, 
                                            rect.second->a, rect.second->b)) {
              ghost_nodes_outer_tag[int(it-ghost_nodes_outer.begin())] = 2;
              break;
            }
          }
        }
      }

    }
    else {
      fprintf(stdout,"\033[0;31m*** Error: Detected unknown boundary type (%d).\033[0m\n",(int)it->bcType);
      exit(-1);
    }

  }
}


//-------------------------------------------------------------------
// Initialize Xi to coords
void
ReferenceMapOperator::SetInitialCondition(SpaceVariable3D &Xi)
{

#if HYPERELASTICITY_TEST == 1
  // assume time t = 1.0
  double time = 1.0;
  double pi = acos(-1.0);
  double dmax = 0.5, Rmax = 1.0;
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
  Vec3D*** xi     = (Vec3D***)Xi.GetDataPointer();
  double z, r;

  if(iod.mesh.type == MeshData::CYLINDRICAL) {
    for(int k=kk0; k<kkmax; k++)
      for(int j=jj0; j<jjmax; j++)
        for(int i=ii0; i<iimax; i++) {
          z = coords[k][j][i][0];        
          r = coords[k][j][i][1];        
          xi[k][j][i][0] = r<=Rmax ? z - dmax*sin(pi*r/(2.0*Rmax))*time : z;
          xi[k][j][i][1] = r;
          xi[k][j][i][2] = coords[k][j][i][2];
        }
  }
  else { //true 3D
    for(int k=kk0; k<kkmax; k++)
      for(int j=jj0; j<jjmax; j++)
        for(int i=ii0; i<iimax; i++) {
          z = coords[k][j][i][2];
          r = sqrt(coords[k][j][i][0]*coords[k][j][i][0] +
                   coords[k][j][i][1]*coords[k][j][i][1]);        
          xi[k][j][i][0] = coords[k][j][i][0];
          xi[k][j][i][1] = coords[k][j][i][1];
          xi[k][j][i][2] = r<=Rmax ? z - dmax*sin(pi*r/(2.0*Rmax))*time : z;
        }
  }

  Xi.RestoreDataPointerAndInsert();
  coordinates.RestoreDataPointerToLocalVector();
  return;
#endif

  // ------------------------------------------
  Xi.AXPlusBY(0.0, 1.0, coordinates, true); //also populating ghost nodes
  // ------------------------------------------
}

//-------------------------------------------------------------------

// Apply boundary conditions by populating ghost cells of Xi
void
ReferenceMapOperator::ApplyBoundaryConditions(SpaceVariable3D &Xi)
{

  int NX, NY, NZ;
  coordinates.GetGlobalSize(&NX, &NY, &NZ);

  Vec3D*** xi = (Vec3D***)Xi.GetDataPointer();
  
  int i,j,k;
  double r, r1, r2;
  Vec3D xi1, xi2;

  for(auto it = ghost_nodes_outer.begin(); it != ghost_nodes_outer.end(); it++) {

    if(it->type_projection != GhostPoint::FACE)
      continue; //corner (i.e. edge or vertex) nodes are not populated

    i = it->ijk[0];
    j = it->ijk[1];
    k = it->ijk[2];

    Vec3D& xi_im(xi[it->image_ijk[2]][it->image_ijk[1]][it->image_ijk[0]]);
    Vec3D& normal(it->outward_normal);
    Vec3D& p(it->boundary_projection);

    switch(ghost_nodes_outer_tag[int(it-ghost_nodes_outer.begin())]) {

      case 0: //slip wall or symmetry
        xi[k][j][i] = 2.0*(p*normal - xi_im*normal)*normal + xi_im;
        break;

      case 1: //no-slip wall
        xi[k][j][i] = 2.0*p - xi_im;
        break;

      case 2: //"open"-->linear extrapolation
        if(it->side == GhostPoint::LEFT) {
          if(i+2<NX) {
            r  = global_mesh.GetX(i);
            r1 = global_mesh.GetX(i+1);  xi1 = xi[k][j][i+1];
            r2 = global_mesh.GetX(i+2);  xi2 = xi[k][j][i+2];
            xi[k][j][i] = xi1 + (xi2-xi1)/(r2-r1)*(r-r1);
          } else
            xi[k][j][i] = xi_im;
        }
        else if(it->side == GhostPoint::RIGHT) {
          if(i-2>=0) {
            r  = global_mesh.GetX(i);
            r1 = global_mesh.GetX(i-1);  xi1 = xi[k][j][i-1];
            r2 = global_mesh.GetX(i-2);  xi2 = xi[k][j][i-2];
            xi[k][j][i] = xi1 + (xi2-xi1)/(r2-r1)*(r-r1);
          } else
            xi[k][j][i] = xi_im;
        }
        else if(it->side == GhostPoint::BOTTOM) {
          if(j+2<NY) {
            r  = global_mesh.GetY(j);
            r1 = global_mesh.GetY(j+1);  xi1 = xi[k][j+1][i];
            r2 = global_mesh.GetY(j+2);  xi2 = xi[k][j+2][i];
            xi[k][j][i] = xi1 + (xi2-xi1)/(r2-r1)*(r-r1);
          } else
            xi[k][j][i] = xi_im;
        }
        else if(it->side == GhostPoint::TOP) {
          if(j-2>=0) {
            r  = global_mesh.GetY(j);
            r1 = global_mesh.GetY(j-1);  xi1 = xi[k][j-1][i];
            r2 = global_mesh.GetY(j-2);  xi2 = xi[k][j-2][i];
            xi[k][j][i] = xi1 + (xi2-xi1)/(r2-r1)*(r-r1);
          } else
            xi[k][j][i] = xi_im;
        }
        else if(it->side == GhostPoint::BACK) {
          if(k+2<NZ) {
            r  = global_mesh.GetZ(k);
            r1 = global_mesh.GetZ(k+1);  xi1 = xi[k+1][j][i];
            r2 = global_mesh.GetZ(k+2);  xi2 = xi[k+2][j][i];
            xi[k][j][i] = xi1 + (xi2-xi1)/(r2-r1)*(r-r1);
          } else
            xi[k][j][i] = xi_im;
        }
        else if(it->side == GhostPoint::FRONT) {
          if(k-2>=0) {
            r  = global_mesh.GetZ(k);
            r1 = global_mesh.GetZ(k-1);  xi1 = xi[k-1][j][i];
            r2 = global_mesh.GetZ(k-2);  xi2 = xi[k-2][j][i];
            xi[k][j][i] = xi1 + (xi2-xi1)/(r2-r1)*(r-r1);
          } else
            xi[k][j][i] = xi_im;
        }
        break;
      
      default:
        fprintf(stdout,"\033[0;31m*** Error: Unable to impose boundary conditions "
                       "for reference map.\033[0m\n");
        exit(-1);
    }

  }

  Xi.RestoreDataPointerAndInsert();
}

//-------------------------------------------------------------------

void
ReferenceMapOperator::ComputeResidual(SpaceVariable3D &V, SpaceVariable3D &Xi, SpaceVariable3D &R)
{
  Vec5D*** v   = (Vec5D***) V.GetDataPointer();
  Vec3D*** xi  = (Vec3D***) Xi.GetDataPointer();
  Vec3D*** res = (Vec3D***) R.GetDataPointer();

  //***************************************************************
  // Step 1: Calculate partial derivatives of xi 
  //***************************************************************
  vector<int> ind0{0,1,2};

  Vec3D*** s = (Vec3D***)vectorG2.GetDataPointer();
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++)
        s[k][j][i] = xi[k][j][i];
  vectorG2.RestoreDataPointerAndInsert(); //need to exchange

  grad_minus->CalculateFirstDerivativeAtNodes(0/*x*/, vectorG2, ind0, Xil, ind0);
  grad_plus->CalculateFirstDerivativeAtNodes(0/*x*/, vectorG2, ind0, Xir, ind0);
  grad_minus->CalculateFirstDerivativeAtNodes(1/*y*/, vectorG2, ind0, Xib, ind0);
  grad_plus->CalculateFirstDerivativeAtNodes(1/*y*/, vectorG2, ind0, Xit, ind0);
  grad_minus->CalculateFirstDerivativeAtNodes(2/*z*/, vectorG2, ind0, Xik, ind0);
  grad_plus->CalculateFirstDerivativeAtNodes(2/*z*/, vectorG2, ind0, Xif, ind0); 


  //***************************************************************
  // Step 2: Loop through active nodes and compute residual
  //***************************************************************
  int NX, NY, NZ;
  coordinates.GetGlobalSize(&NX, &NY, &NZ);

  Vec3D*** xil = (Vec3D***)Xil.GetDataPointer(); //d(Xi)/dx, left-biased
  Vec3D*** xir = (Vec3D***)Xir.GetDataPointer(); //d(Xi)/dx, right-biased
  Vec3D*** xib = (Vec3D***)Xib.GetDataPointer();
  Vec3D*** xit = (Vec3D***)Xit.GetDataPointer();
  Vec3D*** xik = (Vec3D***)Xik.GetDataPointer();
  Vec3D*** xif = (Vec3D***)Xif.GetDataPointer();
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();

  Vec3D a, b, c, d, e, f;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        a = (i-2>=-1) ? xil[k][j][i] : (xi[k][j][i]-xi[k][j][i-1])/(coords[k][j][i][0]-coords[k][j][i-1][0]);
        b = (i+2<=NX) ? xir[k][j][i] : (xi[k][j][i+1]-xi[k][j][i])/(coords[k][j][i+1][0]-coords[k][j][i][0]);
        c = (j-2>=-1) ? xib[k][j][i] : (xi[k][j][i]-xi[k][j-1][i])/(coords[k][j][i][1]-coords[k][j-1][i][1]);
        d = (j+2<=NY) ? xit[k][j][i] : (xi[k][j+1][i]-xi[k][j][i])/(coords[k][j+1][i][1]-coords[k][j][i][1]);
        e = (k-2>=-1) ? xik[k][j][i] : (xi[k][j][i]-xi[k-1][j][i])/(coords[k][j][i][2]-coords[k-1][j][i][2]);
        f = (k+2<=NZ) ? xif[k][j][i] : (xi[k+1][j][i]-xi[k][j][i])/(coords[k+1][j][i][2]-coords[k][j][i][2]);

        res[k][j][i] = (v[k][j][i][1]>=0 ? a*v[k][j][i][1] : b*v[k][j][i][1])
                     + (v[k][j][i][2]>=0 ? c*v[k][j][i][2] : d*v[k][j][i][2])
                     + (v[k][j][i][3]>=0 ? e*v[k][j][i][3] : f*v[k][j][i][3]);

        res[k][j][i] *= -1.0; //moves the residual to the RHS

      }


  Xil.RestoreDataPointerToLocalVector();
  Xir.RestoreDataPointerToLocalVector();
  Xib.RestoreDataPointerToLocalVector();
  Xit.RestoreDataPointerToLocalVector();
  Xik.RestoreDataPointerToLocalVector();
  Xif.RestoreDataPointerToLocalVector();
  coordinates.RestoreDataPointerToLocalVector();

  V.RestoreDataPointerToLocalVector();
  Xi.RestoreDataPointerToLocalVector();

  R.RestoreDataPointerAndInsert();

}

//-------------------------------------------------------------------


//-------------------------------------------------------------------




