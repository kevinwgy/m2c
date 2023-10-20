/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <LevelSetOperator.h>
#include <SpaceOperator.h>
#include <GeoTools.h>
#include <DistancePointToSpheroid.h>
#include <DistancePointToParallelepiped.h>
#include <GradientCalculatorFD3.h>
#include <EmbeddedBoundaryDataSet.h>
#include <EmbeddedBoundaryOperator.h>

#ifdef LEVELSET_TEST
  #include <Vector2D.h>
#endif

using std::min;
using std::max;
using std::pair;
using std::unique_ptr;

extern double domain_diagonal;
//-----------------------------------------------------

LevelSetOperator::LevelSetOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_,
                                   LevelSetSchemeData &iod_ls_, SpaceOperator &spo)
  : comm(comm_), iod(iod_), iod_ls(iod_ls_),
    dms(dm_all_),
    coordinates(spo.GetMeshCoordinates()),
    delta_xyz(spo.GetMeshDeltaXYZ()),
    volume(spo.GetMeshCellVolumes()),
    global_mesh(spo.GetGlobalMeshInfo()),
    spo_ghost_nodes_inner(*spo.GetPointerToInnerGhostNodes()),
    spo_ghost_nodes_outer(*spo.GetPointerToOuterGhostNodes()),
    reinit(NULL),
    scalar(comm_, &(dm_all_.ghosted1_1dof)),
    scalarG2(comm_, &(dm_all_.ghosted2_1dof)),
    ul(comm_, &(dm_all_.ghosted1_1dof)),
    ur(comm_, &(dm_all_.ghosted1_1dof)),
    vb(comm_, &(dm_all_.ghosted1_1dof)),
    vt(comm_, &(dm_all_.ghosted1_1dof)),
    wk(comm_, &(dm_all_.ghosted1_1dof)),
    wf(comm_, &(dm_all_.ghosted1_1dof)),
    dudx(comm_, &(dm_all_.ghosted1_1dof)),
    dvdy(comm_, &(dm_all_.ghosted1_1dof)),
    dwdz(comm_, &(dm_all_.ghosted1_1dof)),
    Phil(comm_, &(dm_all_.ghosted1_1dof)),
    Phir(comm_, &(dm_all_.ghosted1_1dof)),
    Phib(comm_, &(dm_all_.ghosted1_1dof)),
    Phit(comm_, &(dm_all_.ghosted1_1dof)),
    Phik(comm_, &(dm_all_.ghosted1_1dof)),
    Phif(comm_, &(dm_all_.ghosted1_1dof)),
    Level(comm_, &(dm_all_.ghosted1_1dof)),
    UsefulG2(comm_, &(dm_all_.ghosted2_1dof)),
    Active(comm_, &(dm_all_.ghosted1_1dof))
{

  materialid = iod_ls.materialid;

  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);

  CreateGhostNodeLists(); //create ghost_nodes_inner and ghost_nodes_outer

  // spatial discretization method
  if(iod_ls.solver == LevelSetSchemeData::FINITE_VOLUME) {
    rec = new Reconstructor(comm_, dm_all_, iod_ls_.rec, coordinates, delta_xyz);
    //TODO: currently, passing the ghost nodes in SPO. rec uses MesData::BoundaryType (needs to be updated!)
    rec->Setup(spo.GetPointerToInnerGhostNodes(),
               spo.GetPointerToOuterGhostNodes()); //this function requires mesh info (dxyz)
    grad_minus = NULL;
    grad_plus = NULL;
  }
  else if(iod_ls.solver == LevelSetSchemeData::FINITE_DIFFERENCE) {
    if(iod_ls.fd == LevelSetSchemeData::UPWIND_CENTRAL_3) { //currently this is the only option
      grad_minus = new GradientCalculatorFD3(comm, dm_all_, coordinates, delta_xyz, -1);
      grad_plus  = new GradientCalculatorFD3(comm, dm_all_, coordinates, delta_xyz,  1);
    }
    rec = NULL;
  }

  if(iod_ls.bandwidth < INT_MAX) {//user specified narrow-band level set method
    narrow_band = true;
    if(iod_ls.bandwidth <= 1) {
      print_error("*** Error: In the narrow-band level set method, bandwidth should be at least 3 (current: %d)\n",
                  iod_ls.bandwidth);
      exit_mpi();
    }
    if(iod_ls.reinit.frequency <= 0 || iod_ls.reinit.frequency >= iod_ls.bandwidth) {
      print_error("*** Error: The level set reinitialization "
                  "frequency (actually, time-step period) should be smaller than bandwidth (current: %d).\n",
                  iod_ls.reinit.frequency);
      exit_mpi();
    }
  }
  else
    narrow_band = false;

  if(iod_ls.reinit.frequency>0 || iod_ls.reinit.frequency_dt>0)
    reinit = new LevelSetReinitializer(comm, dm_all_, iod_ls, coordinates, delta_xyz,
                                       ghost_nodes_inner, ghost_nodes_outer);

}

//-----------------------------------------------------

LevelSetOperator::~LevelSetOperator()
{
  if(reinit) delete reinit;
  if(rec) delete rec;
  if(grad_minus) delete grad_minus;
  if(grad_plus) delete grad_plus;
}

//-----------------------------------------------------

void LevelSetOperator::Destroy()
{
  if(rec)         rec->Destroy();
  if(grad_minus)  grad_minus->Destroy();
  if(grad_plus)   grad_plus->Destroy();

  if(reinit)      reinit->Destroy();

  scalar.Destroy();
  scalarG2.Destroy();
  ul.Destroy();
  ur.Destroy();
  vb.Destroy();
  vt.Destroy();
  dudx.Destroy();
  dvdy.Destroy();
  dwdz.Destroy();
  wk.Destroy();
  wf.Destroy();
  Phil.Destroy();
  Phir.Destroy();
  Phib.Destroy();
  Phit.Destroy();
  Phik.Destroy();
  Phif.Destroy();
  Level.Destroy();
  UsefulG2.Destroy();
  Active.Destroy();
}

//-----------------------------------------------------

void 
LevelSetOperator::AXPlusBY(double a, SpaceVariable3D &X, double b, SpaceVariable3D &Y, bool workOnGhost)
{
  if(X.NumDOF() != Y.NumDOF()) {
    print_error("*** Error: Vector operation (AXPlusBY) failed due to inconsistent sizes (%d vs. %d)\n", 
                 X.NumDOF(), Y.NumDOF());
    exit_mpi();
  }

  int dof = X.NumDOF();

  double*** x = X.GetDataPointer();
  double*** y = Y.GetDataPointer();

  if(!narrow_band) {

    int myi0, myj0, myk0, myimax, myjmax, mykmax;

    if(workOnGhost)
      X.GetGhostedCornerIndices(&myi0, &myj0, &myk0, &myimax, &myjmax, &mykmax);
    else
      X.GetCornerIndices(&myi0, &myj0, &myk0, &myimax, &myjmax, &mykmax);

    for(int k=myk0; k<mykmax; k++)
      for(int j=myj0; j<myjmax; j++)
        for(int i=myi0; i<myimax; i++)
          for(int p=0; p<dof; p++)
            x[k][j][i*dof+p] = a*x[k][j][i*dof+p] + b*y[k][j][i*dof+p];

  }
  else { //narrow-band

    for(auto it = useful_nodes.begin(); it != useful_nodes.end(); it++) {
      int i((*it)[0]), j((*it)[1]), k((*it)[2]);
      if(!workOnGhost && !X.IsHere(i,j,k,false)) 
        continue;
      for(int p=0; p<dof; p++)
        x[k][j][i*dof+p] = a*x[k][j][i*dof+p] + b*y[k][j][i*dof+p];
    }    

  }

  X.RestoreDataPointerAndInsert();
  Y.RestoreDataPointerToLocalVector();
}

//-----------------------------------------------------

void LevelSetOperator::CreateGhostNodeLists()
{
  ghost_nodes_inner.clear();
  ghost_nodes_outer.clear();

  int NX, NY, NZ;
  coordinates.GetGlobalSize(&NX, &NY, &NZ);

  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();

  Int3 image;
  Vec3D proj(0.0), out_normal(0.0);
  LevelSetSchemeData::BcType bcType = LevelSetSchemeData::NONE;
  GhostPoint::Side side = GhostPoint::UNDEFINED;
  int counter;
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++) {

        if(k>=k0 && k<kmax && j>=j0 && j<jmax && i>=i0 && i<imax)
          continue; //interior of the subdomain

        if(coordinates.OutsidePhysicalDomain(i,j,k)) {//outside physical domain

          //determine the image, the projection point and boundary condition
          image      = 0;
          proj       = 0.0;
          out_normal = 0.0;
          counter    = 0;
          bcType     = LevelSetSchemeData::NONE;
          side       = GhostPoint::UNDEFINED;

          if(i<0)        {image[0] = -i-1;
                          proj[0]  = iod.mesh.x0;      out_normal[0] = -1.0;
                          bcType   = iod_ls.bc_x0;     side = GhostPoint::LEFT;     counter++;}
          else if(i>=NX) {image[0] = NX+(i-NX)-1;
                          proj[0]  = iod.mesh.xmax;    out_normal[0] =  1.0;
                          bcType   = iod_ls.bc_xmax;   side = GhostPoint::RIGHT;    counter++;}
          else           {image[0] = i;
                          proj[0]  = coords[k][j][i][0];}


          if(j<0)        {image[1] = -j-1;
                          proj[1]  = iod.mesh.y0;      out_normal[1] = -1.0;
                          bcType   = iod_ls.bc_y0;     side = GhostPoint::BOTTOM;   counter++;}
          else if(j>=NY) {image[1] = NY+(j-NY)-1;
                          proj[1]  = iod.mesh.ymax;    out_normal[1] =  1.0;
                          bcType   = iod_ls.bc_ymax;   side = GhostPoint::TOP;      counter++;}
          else           {image[1] = j;
                          proj[1]  = coords[k][j][i][1];}


          if(k<0)        {image[2] = -k-1;
                          proj[2]  = iod.mesh.z0;      out_normal[2] = -1.0;
                          bcType   = iod_ls.bc_z0;     side = GhostPoint::BACK;     counter++;}
          else if(k>=NZ) {image[2] = NZ+(k-NZ)-1;
                          proj[2]  = iod.mesh.zmax;    out_normal[2] =  1.0;
                          bcType   = iod_ls.bc_zmax;   side = GhostPoint::FRONT;    counter++;}
          else           {image[2] = k;
                          proj[2]  = coords[k][j][i][2];}

          out_normal /= out_normal.norm();

          assert(counter<=3 && counter>0);

          if(counter == 1)
            ghost_nodes_outer.push_back(GhostPoint(Int3(i,j,k), image, GhostPoint::FACE,
                                        proj, out_normal, (int)bcType, side));
          else if(counter == 2)
            ghost_nodes_outer.push_back(GhostPoint(Int3(i,j,k), image, GhostPoint::EDGE,
                                        proj, out_normal, 0));
          else
            ghost_nodes_outer.push_back(GhostPoint(Int3(i,j,k), image, GhostPoint::VERTEX,
                                        proj, out_normal, 0));

        }
        else //inside physical domain
          ghost_nodes_inner.push_back(GhostPoint(i,j,k));

      }

  int nInner = ghost_nodes_inner.size();
  int nOuter = ghost_nodes_outer.size();
  MPI_Allreduce(MPI_IN_PLACE, &nInner, 1, MPI_INT, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, &nOuter, 1, MPI_INT, MPI_SUM, comm);

  coordinates.RestoreDataPointerToLocalVector();

}

//-----------------------------------------------------

bool 
LevelSetOperator::ApplyInitialConditionWithinEnclosure(UserSpecifiedEnclosureData &enclosure,
                                                       SpaceVariable3D &Phi)
{

  bool applied_ic = false;
  double*** phi = Phi.GetDataPointer();

  // Create a reinitializer (if not available) for one-time use within this function
  // This is a full-domain reinitializer
  bool created_tmp_reinit = false;
  if(!reinit) {
    reinit = new LevelSetReinitializer(comm, dms, iod_ls, coordinates, delta_xyz,
                                       ghost_nodes_inner, ghost_nodes_outer);
    created_tmp_reinit = true;
  }

  // The construction and use of EmbeddedBoundaryOperator is a repetition of what has been done in
  // SpaceOperator::ApplyInitialConditionWithinEnclosure. However, this is a one-time effort. The
  // overhead is acceptable.

  //create an ad hoc EmbeddedSurfaceData
  EmbeddedSurfaceData esd;
  esd.provided_by_another_solver = EmbeddedSurfaceData::NO;
  esd.surface_thickness          = enclosure.surface_thickness;
  esd.filename                   = enclosure.surface_filename;

  //create the embedded operator to track the surface
  EmbeddedBoundaryOperator embed(comm, esd);
  embed.SetCommAndMeshInfo(dms, coordinates, spo_ghost_nodes_inner, spo_ghost_nodes_outer,
                           global_mesh);
  embed.SetupIntersectors();

  //track the surface, and get access to the results
  double max_dist = embed.TrackSurfaces(5); //compute "phi" for 5 layers
  auto EBDS = embed.GetPointerToEmbeddedBoundaryData(0);

  double*** color = EBDS->Color_ptr->GetDataPointer();
  double*** psi   = EBDS->Phi_ptr->GetDataPointer();
  double*** cpi   = EBDS->ClosestPointIndex_ptr->GetDataPointer();
  int nLayer = EBDS->Phi_nLayer;
  assert(nLayer>=2);

  //impose i.c. for phi inside and near enclosures 
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        if(cpi[k][j][i]<0 && color[k][j][i]>=0)
          continue;

        if(!applied_ic)
          applied_ic = true;

        // determine phi at [k][j][i]
        if(color[k][j][i]<0) { //"enclosed" (must be consistent with ID specified in space operator)
          if(phi[k][j][i]>=0)
            phi[k][j][i] = -psi[k][j][i];
          else
            phi[k][j][i] = std::max(phi[k][j][i], -psi[k][j][i]);
        } 
        else if(color[k][j][i] == 0) {//occluded
          phi[k][j][i] = 0.0; //note that psi may not be 0, because surface has finite thickness
        } 
        else { //"outside"
          if(phi[k][j][i]>=0)
            phi[k][j][i] = std::min(phi[k][j][i], psi[k][j][i]);
          else //this is really bad... the user should avoid this...
            phi[k][j][i] = std::max(phi[k][j][i], -psi[k][j][i]);
        }
      }

  int tmp_int = applied_ic ? 1 : 0;
  MPI_Allreduce(MPI_IN_PLACE, &tmp_int, 1, MPI_INT, MPI_MAX, comm);
  if(tmp_int>0)
    applied_ic = true;
 
  if(!applied_ic) {
    print_warning("Warning: The surface in %s may either have no enclosure, or"
                  " lie outside the computational domain.\n", enclosure.surface_filename);
  }

  // Give other enclosed nodes "phi_min" (to avoid large jumps)
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        if(color[k][j][i]<0 && phi[k][j][i]<-max_dist) //"enclosed" 
          phi[k][j][i] = -max_dist; 
        else if(color[k][j][i]>0 && phi[k][j][i]>max_dist)
          phi[k][j][i] = max_dist;
      }

  Phi.RestoreDataPointerAndInsert();

  // reinitialize phi here if tmp reinit is created. Otherwise, phi will be reinitialized near
  // the end of the function "SetInitialCondition"
  if(created_tmp_reinit) {
    ApplyBoundaryConditions(Phi); //need this before reinitialization!
    reinit->ReinitializeFullDomain(Phi, 600); //one-time use
    reinit->Destroy(); 
    delete reinit;
    reinit = NULL;
  }

  //clean-up
  EBDS->Color_ptr->RestoreDataPointerToLocalVector();
  EBDS->Phi_ptr->RestoreDataPointerToLocalVector();
  EBDS->ClosestPointIndex_ptr->RestoreDataPointerToLocalVector();

  embed.Destroy();

  return applied_ic;
}

//-----------------------------------------------------

// Apply boundary conditions by populating ghost cells of Phi
void LevelSetOperator::ApplyBoundaryConditions(SpaceVariable3D &Phi)
{

  double*** phi = (double***) Phi.GetDataPointer();
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();

  int NX, NY, NZ;
  Phi.GetGlobalSize(&NX, &NY, &NZ);

  double r, r1, r2, f1, f2;

  for(auto it = ghost_nodes_outer.begin(); it != ghost_nodes_outer.end();  it++) {

    if(it->type_projection != GhostPoint::FACE)
      continue; //corner (i.e. edge or vertex) nodes are not populated

    int i(it->ijk[0]), j(it->ijk[1]), k(it->ijk[2]);
    int im_i(it->image_ijk[0]), im_j(it->image_ijk[1]), im_k(it->image_ijk[2]);

    if(it->bcType == (int)LevelSetSchemeData::ZERO_NEUMANN) {

      phi[k][j][i] = phi[im_k][im_j][im_i];

    }
    else if ((it->bcType == (int)LevelSetSchemeData::LINEAR_EXTRAPOLATION) ||
             (it->bcType == (int)LevelSetSchemeData::NON_NEGATIVE)) {

      //make sure the width of the subdomain is big enough for linear extrapolation
      if(it->side == GhostPoint::LEFT) {
        if(i+2<NX) {
          r  = coords[k][j][i][0];
          r1 = coords[k][j][i+1][0];  f1 = phi[k][j][i+1];
          r2 = coords[k][j][i+2][0];  f2 = phi[k][j][i+2];
          phi[k][j][i] = f1 + (f2-f1)/(r2-r1)*(r-r1);
        } else
          phi[k][j][i] = phi[im_k][im_j][im_i];
      }
      else if(it->side == GhostPoint::RIGHT) {
        if(i-2>=0) {
          r  = coords[k][j][i][0];
          r1 = coords[k][j][i-1][0];  f1 = phi[k][j][i-1];
          r2 = coords[k][j][i-2][0];  f2 = phi[k][j][i-2];
          phi[k][j][i] = f1 + (f2-f1)/(r2-r1)*(r-r1);
        } else
          phi[k][j][i] = phi[im_k][im_j][im_i];
      }
      else if(it->side == GhostPoint::BOTTOM) {
        if(j+2<NY) {
          r  = coords[k][j][i][1];
          r1 = coords[k][j+1][i][1];  f1 = phi[k][j+1][i];
          r2 = coords[k][j+2][i][1];  f2 = phi[k][j+2][i];
          phi[k][j][i] = f1 + (f2-f1)/(r2-r1)*(r-r1);
        } else
          phi[k][j][i] = phi[im_k][im_j][im_i];
      }
      else if(it->side == GhostPoint::TOP) {
        if(j-2>=0) {
          r  = coords[k][j][i][1];
          r1 = coords[k][j-1][i][1];  f1 = phi[k][j-1][i];
          r2 = coords[k][j-2][i][1];  f2 = phi[k][j-2][i];
          phi[k][j][i] = f1 + (f2-f1)/(r2-r1)*(r-r1);
        } else
          phi[k][j][i] = phi[im_k][im_j][im_i];
      }
      else if(it->side == GhostPoint::BACK) {
        if(k+2<NZ) { 
          r  = coords[k][j][i][2];
          r1 = coords[k+1][j][i][2];  f1 = phi[k+1][j][i];
          r2 = coords[k+2][j][i][2];  f2 = phi[k+2][j][i];
          phi[k][j][i] = f1 + (f2-f1)/(r2-r1)*(r-r1);
        } else
          phi[k][j][i] = phi[im_k][im_j][im_i];
      }
      else if(it->side == GhostPoint::FRONT) {
        if(k-2>=0) {
          r  = coords[k][j][i][2];
          r1 = coords[k-1][j][i][2];  f1 = phi[k-1][j][i];
          r2 = coords[k-2][j][i][2];  f2 = phi[k-2][j][i];
          phi[k][j][i] = f1 + (f2-f1)/(r2-r1)*(r-r1);
        } else
          phi[k][j][i] = phi[im_k][im_j][im_i];
      }

      if(it->bcType == (int)LevelSetSchemeData::NON_NEGATIVE) {
        int ind(0);
        if(it->side == GhostPoint::LEFT || it->side == GhostPoint::RIGHT)  ind = 0;
        if(it->side == GhostPoint::BOTTOM || it->side == GhostPoint::TOP)  ind = 1;
        if(it->side == GhostPoint::BACK || it->side == GhostPoint::FRONT)  ind = 2;
        double d2wall = 0.5*fabs(coords[im_k][im_j][im_i][ind] - coords[k][j][i][ind]); //distance from this ghost node to the wall boundary
        phi[k][j][i] = std::max(-d2wall, phi[k][j][i]);
      }

    }

  }

  Phi.RestoreDataPointerAndInsert();

  coordinates.RestoreDataPointerToLocalVector();

}

//-----------------------------------------------------
// Apply boundary conditions by populating ghost cells of NPhi (unit normal vector of zero levelset contour)
void LevelSetOperator::ApplyBoundaryConditionsNPhi(SpaceVariable3D &NPhi)
{

  Vec3D*** normal = (Vec3D***) NPhi.GetDataPointer();
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();

  int NX, NY, NZ;
  NPhi.GetGlobalSize(&NX, &NY, &NZ);

  double r, r1, r2, f1, f2;

  for(auto it = ghost_nodes_outer.begin(); it != ghost_nodes_outer.end();  it++) {

    if(it->type_projection != GhostPoint::FACE)
      continue; //corner (i.e. edge or vertex) nodes are not populated

    int i(it->ijk[0]), j(it->ijk[1]), k(it->ijk[2]);
    int im_i(it->image_ijk[0]), im_j(it->image_ijk[1]), im_k(it->image_ijk[2]);

    if(it->bcType == (int)LevelSetSchemeData::ZERO_NEUMANN) {
      normal[k][j][i][0] = normal[im_k][im_j][im_i][0];
      normal[k][j][i][1] = normal[im_k][im_j][im_i][1];
      normal[k][j][i][2] = normal[im_k][im_j][im_i][2];

      if(it->side == GhostPoint::LEFT || it->side == GhostPoint::RIGHT) normal[k][j][i][0] = -1.*normal[im_k][im_j][im_i][0];
      if(it->side == GhostPoint::BOTTOM || it->side == GhostPoint::TOP) normal[k][j][i][1] = -1.*normal[im_k][im_j][im_i][1];
      if(it->side == GhostPoint::BACK || it->side == GhostPoint::FRONT) normal[k][j][i][2] = -1.*normal[im_k][im_j][im_i][2];
    }
    else if ((it->bcType == (int)LevelSetSchemeData::LINEAR_EXTRAPOLATION) ||
             (it->bcType == (int)LevelSetSchemeData::NON_NEGATIVE)) {

      //make sure the width of the subdomain is big enough for linear extrapolation
      if(it->side == GhostPoint::LEFT) {
        if(i+2<NX) {
          r  = coords[k][j][i][0];
          r1 = coords[k][j][i+1][0]; 
          r2 = coords[k][j][i+2][0];  
          for (int dir = 0; dir < 3; dir++) { //operating on three components of NPhi
            f1 = normal[k][j][i+1][dir];
            f2 = normal[k][j][i+2][dir];
            normal[k][j][i][dir] = f1 + (f2-f1)/(r2-r1)*(r-r1);
          }
        } else {
	  normal[k][j][i][0] = normal[im_k][im_j][im_i][0];
	  normal[k][j][i][1] = normal[im_k][im_j][im_i][1];
	  normal[k][j][i][2] = normal[im_k][im_j][im_i][2];
	} 
      }
      else if(it->side == GhostPoint::RIGHT) {
        if(i-2>=0) {
          r  = coords[k][j][i][0];
          r1 = coords[k][j][i-1][0];  
          r2 = coords[k][j][i-2][0]; 
          for (int dir = 0; dir < 3; dir++) { //operating on three components of NPhi
            f1 = normal[k][j][i-1][dir];
            f2 = normal[k][j][i-2][dir];
	    normal[k][j][i][dir] = f1 + (f2-f1)/(r2-r1)*(r-r1);
          } 
        } else {
	  normal[k][j][i][0] = normal[im_k][im_j][im_i][0];
	  normal[k][j][i][1] = normal[im_k][im_j][im_i][1];
	  normal[k][j][i][2] = normal[im_k][im_j][im_i][2];
	} 
      }
      else if(it->side == GhostPoint::BOTTOM) {
        if(j+2<NY) {
          r  = coords[k][j][i][1];
          r1 = coords[k][j+1][i][1]; 
          r2 = coords[k][j+2][i][1];
	  for (int dir = 0; dir < 3; dir++) { //operating on three components of NPhi
            f1 = normal[k][j+1][i][dir];
            f2 = normal[k][j+2][i][dir];
	    normal[k][j][i][dir] = f1 + (f2-f1)/(r2-r1)*(r-r1);
          } 
        } else {
	  normal[k][j][i][0] = normal[im_k][im_j][im_i][0];
	  normal[k][j][i][1] = normal[im_k][im_j][im_i][1];
	  normal[k][j][i][2] = normal[im_k][im_j][im_i][2];
	} 
      }
      else if(it->side == GhostPoint::TOP) {
        if(j-2>=0) {
          r  = coords[k][j][i][1];
          r1 = coords[k][j-1][i][1];
          r2 = coords[k][j-2][i][1];
	  for (int dir = 0; dir < 3; dir++) { //operating on three components of NPhi
            f1 = normal[k][j-1][i][dir];
            f2 = normal[k][j-2][i][dir];
	    normal[k][j][i][dir] = f1 + (f2-f1)/(r2-r1)*(r-r1);
          } 
        } else {
	  normal[k][j][i][0] = normal[im_k][im_j][im_i][0];
	  normal[k][j][i][1] = normal[im_k][im_j][im_i][1];
	  normal[k][j][i][2] = normal[im_k][im_j][im_i][2];
	} 
      }
      else if(it->side == GhostPoint::BACK) {
        if(k+2<NZ) { 
          r  = coords[k][j][i][2];
          r1 = coords[k+1][j][i][2]; 
          r2 = coords[k+2][j][i][2]; 
	  for (int dir = 0; dir < 3; dir++) { //operating on three components of NPhi
            f1 = normal[k+1][j][i][dir];
            f2 = normal[k+2][j][i][dir];
	    normal[k][j][i][dir] = f1 + (f2-f1)/(r2-r1)*(r-r1);
          } 
        } else {
	  normal[k][j][i][0] = normal[im_k][im_j][im_i][0];
	  normal[k][j][i][1] = normal[im_k][im_j][im_i][1];
	  normal[k][j][i][2] = normal[im_k][im_j][im_i][2];
	} 
      }
      else if(it->side == GhostPoint::FRONT) {
        if(k-2>=0) {
          r  = coords[k][j][i][2];
          r1 = coords[k-1][j][i][2]; 
          r2 = coords[k-2][j][i][2];
	  for (int dir = 0; dir < 3; dir++) { //operating on three components of NPhi
            f1 = normal[k-1][j][i][dir];
            f2 = normal[k-2][j][i][dir];
	    normal[k][j][i][dir] = f1 + (f2-f1)/(r2-r1)*(r-r1);
          } 
        } else {
	  normal[k][j][i][0] = normal[im_k][im_j][im_i][0];
	  normal[k][j][i][1] = normal[im_k][im_j][im_i][1];
	  normal[k][j][i][2] = normal[im_k][im_j][im_i][2];
	} 
      }

      if(it->bcType == (int)LevelSetSchemeData::NON_NEGATIVE) {
        std::cout << "Error: The NON_NEGATIVE boundary condition for NPhi (interface normal) has not been implemented at the moment. Execution ends" << std::endl;
        exit_mpi(); 
      }

    }

  }

  NPhi.RestoreDataPointerAndInsert();

  coordinates.RestoreDataPointerToLocalVector();

}

//-----------------------------------------------------

// Apply boundary conditions by populating ghost cells of KappaPhi (curvature of zero Levelset)
void LevelSetOperator::ApplyBoundaryConditionsKappaPhi(SpaceVariable3D &KappaPhi)
{

  double*** curvature = (double***) KappaPhi.GetDataPointer();
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();

  int NX, NY, NZ;
  KappaPhi.GetGlobalSize(&NX, &NY, &NZ);

  double r, r1, r2, f1, f2;

  for(auto it = ghost_nodes_outer.begin(); it != ghost_nodes_outer.end();  it++) {

    if(it->type_projection != GhostPoint::FACE)
      continue; //corner (i.e. edge or vertex) nodes are not populated

    int i(it->ijk[0]), j(it->ijk[1]), k(it->ijk[2]);
    int im_i(it->image_ijk[0]), im_j(it->image_ijk[1]), im_k(it->image_ijk[2]);

    if(it->bcType == (int)LevelSetSchemeData::ZERO_NEUMANN) {

      curvature[k][j][i] = curvature[im_k][im_j][im_i];

    }
    else if ((it->bcType == (int)LevelSetSchemeData::LINEAR_EXTRAPOLATION) ||
             (it->bcType == (int)LevelSetSchemeData::NON_NEGATIVE)) {

      //make sure the width of the subdomain is big enough for linear extrapolation
      if(it->side == GhostPoint::LEFT) {
        if(i+2<NX) {
          r  = coords[k][j][i][0];
          r1 = coords[k][j][i+1][0];  f1 = curvature[k][j][i+1];
          r2 = coords[k][j][i+2][0];  f2 = curvature[k][j][i+2];
          curvature[k][j][i] = f1 + (f2-f1)/(r2-r1)*(r-r1);
        } else
          curvature[k][j][i] = curvature[im_k][im_j][im_i];
      }
      else if(it->side == GhostPoint::RIGHT) {
        if(i-2>=0) {
          r  = coords[k][j][i][0];
          r1 = coords[k][j][i-1][0];  f1 = curvature[k][j][i-1];
          r2 = coords[k][j][i-2][0];  f2 = curvature[k][j][i-2];
          curvature[k][j][i] = f1 + (f2-f1)/(r2-r1)*(r-r1);
        } else
          curvature[k][j][i] = curvature[im_k][im_j][im_i];
      }
      else if(it->side == GhostPoint::BOTTOM) {
        if(j+2<NY) {
          r  = coords[k][j][i][1];
          r1 = coords[k][j+1][i][1];  f1 = curvature[k][j+1][i];
          r2 = coords[k][j+2][i][1];  f2 = curvature[k][j+2][i];
          curvature[k][j][i] = f1 + (f2-f1)/(r2-r1)*(r-r1);
        } else
          curvature[k][j][i] = curvature[im_k][im_j][im_i];
      }
      else if(it->side == GhostPoint::TOP) {
        if(j-2>=0) {
          r  = coords[k][j][i][1];
          r1 = coords[k][j-1][i][1];  f1 = curvature[k][j-1][i];
          r2 = coords[k][j-2][i][1];  f2 = curvature[k][j-2][i];
          curvature[k][j][i] = f1 + (f2-f1)/(r2-r1)*(r-r1);
        } else
          curvature[k][j][i] = curvature[im_k][im_j][im_i];
      }
      else if(it->side == GhostPoint::BACK) {
        if(k+2<NZ) { 
          r  = coords[k][j][i][2];
          r1 = coords[k+1][j][i][2];  f1 = curvature[k+1][j][i];
          r2 = coords[k+2][j][i][2];  f2 = curvature[k+2][j][i];
          curvature[k][j][i] = f1 + (f2-f1)/(r2-r1)*(r-r1);
        } else
          curvature[k][j][i] = curvature[im_k][im_j][im_i];
      }
      else if(it->side == GhostPoint::FRONT) {
        if(k-2>=0) {
          r  = coords[k][j][i][2];
          r1 = coords[k-1][j][i][2];  f1 = curvature[k-1][j][i];
          r2 = coords[k-2][j][i][2];  f2 = curvature[k-2][j][i];
          curvature[k][j][i] = f1 + (f2-f1)/(r2-r1)*(r-r1);
        } else
          curvature[k][j][i] = curvature[im_k][im_j][im_i];
      }

      if(it->bcType == (int)LevelSetSchemeData::NON_NEGATIVE) {
        std::cout << "Error: The NON_NEGATIVE boundary condition for KappaPhi (interface curvature) has not been implemented at the moment. Execution ends" << std::endl;
        exit_mpi(); 
 
        //int ind(0);
        //if(it->side == GhostPoint::LEFT || it->side == GhostPoint::RIGHT)  ind = 0;
        //if(it->side == GhostPoint::BOTTOM || it->side == GhostPoint::TOP)  ind = 1;
        //if(it->side == GhostPoint::BACK || it->side == GhostPoint::FRONT)  ind = 2;
        // double d2wall = 0.5*fabs(coords[im_k][im_j][im_i][ind] - coords[k][j][i][ind]); //distance from this ghost node to the wall boundary
        // phi[k][j][i] = std::max(-d2wall, phi[k][j][i]);
      }

    }

  }

  KappaPhi.RestoreDataPointerAndInsert();

  coordinates.RestoreDataPointerToLocalVector();

}



//-----------------------------------------------------

void LevelSetOperator::ComputeResidual(SpaceVariable3D &V, SpaceVariable3D &Phi, SpaceVariable3D &R,
                                       [[maybe_unused]] double time)
{

#ifdef LEVELSET_TEST
  PrescribeVelocityFieldForTesting(V, Phi, time);
#endif

  if(iod_ls.solver == LevelSetSchemeData::FINITE_VOLUME)
    ComputeResidualFVM(V,Phi,R);
  else if(iod_ls.solver == LevelSetSchemeData::FINITE_DIFFERENCE)
    ComputeResidualFDM(V,Phi,R);

}

//-----------------------------------------------------

void LevelSetOperator::ComputeResidualFDM(SpaceVariable3D &V, SpaceVariable3D &Phi, SpaceVariable3D &R)
{

  if(!narrow_band) 
    ComputeResidualFDM_FullDomain(V,Phi,R);
  else
    ComputeResidualFDM_NarrowBand(V,Phi,R);
}

//-----------------------------------------------------

void LevelSetOperator::ComputeResidualFDM_FullDomain(SpaceVariable3D &V, SpaceVariable3D &Phi, SpaceVariable3D &R)
{
  Vec5D***    v = (Vec5D***) V.GetDataPointer();
  double*** phi = Phi.GetDataPointer();
  double*** res = R.GetDataPointer(); //residual, on the right-hand-side of the ODE 

  //***************************************************************
  // Step 1: Calculate partial derivatives of phi
  //***************************************************************
  vector<int> ind0{0};
   
  double*** s = scalarG2.GetDataPointer();
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++)
        s[k][j][i] = phi[k][j][i]; 
  scalarG2.RestoreDataPointerAndInsert(); //need to exchange 

  grad_minus->CalculateFirstDerivativeAtNodes(0/*x*/, scalarG2, ind0, Phil, ind0);
  grad_plus->CalculateFirstDerivativeAtNodes(0/*x*/, scalarG2, ind0, Phir, ind0);
  grad_minus->CalculateFirstDerivativeAtNodes(1/*y*/, scalarG2, ind0, Phib, ind0);
  grad_plus->CalculateFirstDerivativeAtNodes(1/*y*/, scalarG2, ind0, Phit, ind0);
  grad_minus->CalculateFirstDerivativeAtNodes(2/*z*/, scalarG2, ind0, Phik, ind0);
  grad_plus->CalculateFirstDerivativeAtNodes(2/*z*/, scalarG2, ind0, Phif, ind0);


  //***************************************************************
  // Step 2: Loop through active nodes and compute residual
  //***************************************************************
  int NX, NY, NZ;
  coordinates.GetGlobalSize(&NX, &NY, &NZ);

  double*** phil = Phil.GetDataPointer(); //d(Phi)/dx, left-biased
  double*** phir = Phir.GetDataPointer(); //d(Phi)/dx, right-biased
  double*** phib = Phib.GetDataPointer();
  double*** phit = Phit.GetDataPointer();
  double*** phik = Phik.GetDataPointer();
  double*** phif = Phif.GetDataPointer();
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();

  double a, b, c, d, e, f;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        a = (i-2>=-1) ? phil[k][j][i] : (phi[k][j][i]-phi[k][j][i-1])/(coords[k][j][i][0]-coords[k][j][i-1][0]);
        b = (i+2<=NX) ? phir[k][j][i] : (phi[k][j][i+1]-phi[k][j][i])/(coords[k][j][i+1][0]-coords[k][j][i][0]);
        c = (j-2>=-1) ? phib[k][j][i] : (phi[k][j][i]-phi[k][j-1][i])/(coords[k][j][i][1]-coords[k][j-1][i][1]);
        d = (j+2<=NY) ? phit[k][j][i] : (phi[k][j+1][i]-phi[k][j][i])/(coords[k][j+1][i][1]-coords[k][j][i][1]);
        e = (k-2>=-1) ? phik[k][j][i] : (phi[k][j][i]-phi[k-1][j][i])/(coords[k][j][i][2]-coords[k-1][j][i][2]);
        f = (k+2<=NZ) ? phif[k][j][i] : (phi[k+1][j][i]-phi[k][j][i])/(coords[k+1][j][i][2]-coords[k][j][i][2]);

        res[k][j][i] = (v[k][j][i][1]>=0 ? a*v[k][j][i][1] : b*v[k][j][i][1])
                     + (v[k][j][i][2]>=0 ? c*v[k][j][i][2] : d*v[k][j][i][2])
                     + (v[k][j][i][3]>=0 ? e*v[k][j][i][3] : f*v[k][j][i][3]);

        res[k][j][i] *= -1.0; //moves the residual to the RHS

      }


  Phil.RestoreDataPointerToLocalVector();
  Phir.RestoreDataPointerToLocalVector();
  Phib.RestoreDataPointerToLocalVector();
  Phit.RestoreDataPointerToLocalVector();
  Phik.RestoreDataPointerToLocalVector();
  Phif.RestoreDataPointerToLocalVector();
  coordinates.RestoreDataPointerToLocalVector();

  V.RestoreDataPointerToLocalVector();
  Phi.RestoreDataPointerToLocalVector();

  R.RestoreDataPointerAndInsert();

}

//-----------------------------------------------------

void LevelSetOperator::ComputeResidualFDM_NarrowBand(SpaceVariable3D &V, SpaceVariable3D &Phi, SpaceVariable3D &R)
{
  Vec5D***    v = (Vec5D***) V.GetDataPointer();
  double*** phi = Phi.GetDataPointer();
  double*** res = R.GetDataPointer(); //residual, on the right-hand-side of the ODE 

  //***************************************************************
  // Step 1: Calculate partial derivatives of phi
  //***************************************************************
  vector<int> ind0{0};
   
  double*** s = scalarG2.GetDataPointer();
  for(auto it = useful_nodes.begin(); it != useful_nodes.end(); it++) {
    int i((*it)[0]), j((*it)[1]), k((*it)[2]);
    s[k][j][i] = phi[k][j][i]; 
  }
  scalarG2.RestoreDataPointerAndInsert(); //need to exchange (scalarG2 has 2 ghost layers, but useful_nodes only has 1)

  grad_minus->CalculateFirstDerivativeAtSelectedNodes(0/*x*/, active_nodes, scalarG2, ind0, Phil, ind0);
  grad_plus->CalculateFirstDerivativeAtSelectedNodes(0/*x*/, active_nodes, scalarG2, ind0, Phir, ind0);
  grad_minus->CalculateFirstDerivativeAtSelectedNodes(1/*y*/, active_nodes, scalarG2, ind0, Phib, ind0);
  grad_plus->CalculateFirstDerivativeAtSelectedNodes(1/*y*/, active_nodes, scalarG2, ind0, Phit, ind0);
  grad_minus->CalculateFirstDerivativeAtSelectedNodes(2/*z*/, active_nodes, scalarG2, ind0, Phik, ind0);
  grad_plus->CalculateFirstDerivativeAtSelectedNodes(2/*z*/, active_nodes, scalarG2, ind0, Phif, ind0);


  //***************************************************************
  // Step 2: Clear residual at useful_nodes (so that if a useful
  //         node becomes out-of-band next time-step, it does not
  //         carry a non-zero residual)
  //***************************************************************
  for(auto it = useful_nodes.begin(); it != useful_nodes.end(); it++)
    res[(*it)[2]][(*it)[1]][(*it)[0]] = 0.0;


  //***************************************************************
  // Step 3: Loop through active nodes inside subdomain and compute residual
  //***************************************************************
  double*** phil = Phil.GetDataPointer(); //d(Phi)/dx, left-biased
  double*** phir = Phir.GetDataPointer(); //d(Phi)/dx, right-biased
  double*** phib = Phib.GetDataPointer();
  double*** phit = Phit.GetDataPointer();
  double*** phik = Phik.GetDataPointer();
  double*** phif = Phif.GetDataPointer();
  double*** useful = UsefulG2.GetDataPointer();
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();

  double dm, dp; 
  for(auto it = active_nodes.begin(); it != active_nodes.end(); it++) {

    int i((*it)[0]), j((*it)[1]), k((*it)[2]);

    if(!coordinates.IsHere(i,j,k,false))
      continue; // only work on nodes in the subdomain interior

    dm = useful[k][j][i-2] ? phil[k][j][i] : (phi[k][j][i]-phi[k][j][i-1])/(coords[k][j][i][0]-coords[k][j][i-1][0]);
    dp = useful[k][j][i+2] ? phir[k][j][i] : (phi[k][j][i+1]-phi[k][j][i])/(coords[k][j][i+1][0]-coords[k][j][i][0]);

    res[k][j][i] = v[k][j][i][1]>=0 ? dm*v[k][j][i][1] : dp*v[k][j][i][1]; 

    dm = useful[k][j-2][i] ? phib[k][j][i] : (phi[k][j][i]-phi[k][j-1][i])/(coords[k][j][i][1]-coords[k][j-1][i][1]);
    dp = useful[k][j+2][i] ? phit[k][j][i] : (phi[k][j+1][i]-phi[k][j][i])/(coords[k][j+1][i][1]-coords[k][j][i][1]);

    res[k][j][i] += v[k][j][i][2]>=0 ? dm*v[k][j][i][2] : dp*v[k][j][i][2]; 

    dm = useful[k-2][j][i] ? phik[k][j][i] : (phi[k][j][i]-phi[k-1][j][i])/(coords[k][j][i][2]-coords[k-1][j][i][2]);
    dp = useful[k+2][j][i] ? phif[k][j][i] : (phi[k+1][j][i]-phi[k][j][i])/(coords[k+1][j][i][2]-coords[k][j][i][2]);

    res[k][j][i] += v[k][j][i][3]>=0 ? dm*v[k][j][i][3] : dp*v[k][j][i][3];

    res[k][j][i] *= -1.0; //moves the residual to the RHS

  }


  Phil.RestoreDataPointerToLocalVector();
  Phir.RestoreDataPointerToLocalVector();
  Phib.RestoreDataPointerToLocalVector();
  Phit.RestoreDataPointerToLocalVector();
  Phik.RestoreDataPointerToLocalVector();
  Phif.RestoreDataPointerToLocalVector();
  UsefulG2.RestoreDataPointerToLocalVector();
  coordinates.RestoreDataPointerToLocalVector();

  V.RestoreDataPointerToLocalVector();
  Phi.RestoreDataPointerToLocalVector();

  R.RestoreDataPointerAndInsert();

}

//-----------------------------------------------------

void LevelSetOperator::ComputeResidualFVM(SpaceVariable3D &V, SpaceVariable3D &Phi, SpaceVariable3D &R)
{

  if(!narrow_band) {

    Reconstruct(V, Phi); // => ul, ur, dudx, vb, vt, dvdy, wk, wf, dwdz, Phil, Phir, Phib, Phit, Phik, Phif

    ComputeAdvectionFlux(R);  //Advection flux, on the right-hand-side

    AddSourceTerm(Phi, R);
  }
  else {

    ReconstructInBand(V, Phi);

    ComputeAdvectionFluxInBand(R);

    AddSourceTermInBand(Phi, R);

  }

}

//-----------------------------------------------------

void LevelSetOperator::Reconstruct(SpaceVariable3D &V, SpaceVariable3D &Phi)
{
  Vec5D*** v = (Vec5D***) V.GetDataPointer();
   
  // Reconstruction: x-velocity
  double*** s = (double***) scalar.GetDataPointer();
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++)
        s[k][j][i] = v[k][j][i][1]; //x-velocity
  scalar.RestoreDataPointerToLocalVector(); //no need to communicate w/ neighbors
  rec->ReconstructIn1D(0/*x-dir*/, scalar, ul, ur, &dudx);

  // Reconstruction: y-velocity
  s = (double***) scalar.GetDataPointer();
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++)
        s[k][j][i] = v[k][j][i][2]; //y-velocity
  scalar.RestoreDataPointerToLocalVector(); //no need to communicate w/ neighbors
  rec->ReconstructIn1D(1/*y-dir*/, scalar, vb, vt, &dvdy);

  // Reconstruction: z-velocity
  s = (double***) scalar.GetDataPointer();
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++)
        s[k][j][i] = v[k][j][i][3]; //z-velocity
  scalar.RestoreDataPointerToLocalVector(); //no need to communicate w/ neighbors
  rec->ReconstructIn1D(2/*z-dir*/, scalar, wk, wf, &dwdz);

  // Reconstruction: Phi
  rec->Reconstruct(Phi, Phil, Phir, Phib, Phit, Phik, Phif);

  V.RestoreDataPointerToLocalVector(); //no changes made
}

//-----------------------------------------------------

void LevelSetOperator::ReconstructInBand(SpaceVariable3D &V, SpaceVariable3D &Phi)
{
  Vec5D*** v = (Vec5D***) V.GetDataPointer();
   
  // Reconstruction: x-velocity
  double*** s = (double***) scalar.GetDataPointer();
  for(auto it = useful_nodes.begin(); it != useful_nodes.end(); it++) {
    int i((*it)[0]), j((*it)[1]), k((*it)[2]);
    s[k][j][i] = v[k][j][i][1]; //x-velocity
  }
  scalar.RestoreDataPointerToLocalVector(); //no need to communicate w/ neighbors
  rec->ReconstructIn1D(0/*x-dir*/, scalar, ul, ur, &dudx, &Active);
        //KW noticed that turning on linear reconstruction at the boundary (i.e. useful but not active)
        //of the band would lead to oscillations. Consistent with what Hartmann et al. (2008) says

  // Reconstruction: y-velocity
  s = (double***) scalar.GetDataPointer();
  for(auto it = useful_nodes.begin(); it != useful_nodes.end(); it++) {
    int i((*it)[0]), j((*it)[1]), k((*it)[2]);
    s[k][j][i] = v[k][j][i][2]; //y-velocity
  }
  scalar.RestoreDataPointerToLocalVector(); //no need to communicate w/ neighbors
  rec->ReconstructIn1D(1/*y-dir*/, scalar, vb, vt, &dvdy, &Active);

  // Reconstruction: z-velocity
  s = (double***) scalar.GetDataPointer();
  for(auto it = useful_nodes.begin(); it != useful_nodes.end(); it++) {
    int i((*it)[0]), j((*it)[1]), k((*it)[2]);
    s[k][j][i] = v[k][j][i][3]; //z-velocity
  }
  scalar.RestoreDataPointerToLocalVector(); //no need to communicate w/ neighbors
  rec->ReconstructIn1D(2/*z-dir*/, scalar, wk, wf, &dwdz, &Active);

  // Reconstruction: Phi
  rec->Reconstruct(Phi, Phil, Phir, Phib, Phit, Phik, Phif, NULL/*ID*/, nullptr/*EBDS*/, &Active);

  V.RestoreDataPointerToLocalVector(); //no changes made
}

//-----------------------------------------------------

void LevelSetOperator::ComputeAdvectionFlux(SpaceVariable3D &R)
{
  Vec3D*** dxyz = (Vec3D***)delta_xyz.GetDataPointer();

  double*** phil   = Phil.GetDataPointer();
  double*** phir   = Phir.GetDataPointer();
  double*** phib   = Phib.GetDataPointer();
  double*** phit   = Phit.GetDataPointer();
  double*** phik   = Phik.GetDataPointer();
  double*** phif   = Phif.GetDataPointer();
  double*** uldata = ul.GetDataPointer();
  double*** urdata = ur.GetDataPointer();
  double*** vbdata = vb.GetDataPointer();
  double*** vtdata = vt.GetDataPointer();
  double*** wkdata = wk.GetDataPointer();
  double*** wfdata = wf.GetDataPointer();
  double*** res    = R.GetDataPointer(); //residual, on the right-hand-side of the ODE

  //initialize R to 0
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++)
          res[k][j][i] = 0.0;

  double localflux;

  // Loop through the domain interior, and the right and top ghost layers. For each cell, calculate the
  // numerical flux across the left and lower cell boundaries/interfaces
  for(int k=k0; k<kkmax; k++) {
    for(int j=j0; j<jjmax; j++) {
      for(int i=i0; i<iimax; i++) {

        //calculate F_{i-1/2,jk}
        if(k!=kkmax-1 && j!=jjmax-1) {
          localflux = ComputeLocalAdvectionFlux(phir[k][j][i-1], phil[k][j][i], 
                                                urdata[k][j][i-1], uldata[k][j][i]) / dxyz[k][j][i][0];
          res[k][j][i-1] -= localflux;
          res[k][j][i]   += localflux;
        }

        //calculate G_{i,j-1/2,k}
        if(k!=kkmax-1 && i!=iimax-1) {
          localflux = ComputeLocalAdvectionFlux(phit[k][j-1][i], phib[k][j][i], 
                                                vtdata[k][j-1][i], vbdata[k][j][i]) / dxyz[k][j][i][1];
          res[k][j-1][i] -= localflux;
          res[k][j][i]   += localflux;
        }

        //calculate H_{ij,k-1/2}
        if(j!=jjmax-1 && i!=iimax-1) {
          localflux = ComputeLocalAdvectionFlux(phif[k-1][j][i], phik[k][j][i], 
                                                wfdata[k-1][j][i], wkdata[k][j][i]) / dxyz[k][j][i][2];
          res[k-1][j][i] -= localflux;
          res[k][j][i]   += localflux;
        }

      }
    }
  }

  // Restore space variables
  delta_xyz.RestoreDataPointerToLocalVector();

  Phil.RestoreDataPointerToLocalVector();
  Phir.RestoreDataPointerToLocalVector();
  Phib.RestoreDataPointerToLocalVector();
  Phit.RestoreDataPointerToLocalVector();
  Phik.RestoreDataPointerToLocalVector();
  Phif.RestoreDataPointerToLocalVector();
  ul.RestoreDataPointerToLocalVector();
  ur.RestoreDataPointerToLocalVector();
  vb.RestoreDataPointerToLocalVector();
  vt.RestoreDataPointerToLocalVector();
  wk.RestoreDataPointerToLocalVector();
  wf.RestoreDataPointerToLocalVector();

  R.RestoreDataPointerAndInsert(); //insert
}

//-----------------------------------------------------

void LevelSetOperator::ComputeAdvectionFluxInBand(SpaceVariable3D &R)
{
  Vec3D*** dxyz = (Vec3D***)delta_xyz.GetDataPointer();

  double*** phil   = Phil.GetDataPointer();
  double*** phir   = Phir.GetDataPointer();
  double*** phib   = Phib.GetDataPointer();
  double*** phit   = Phit.GetDataPointer();
  double*** phik   = Phik.GetDataPointer();
  double*** phif   = Phif.GetDataPointer();
  double*** uldata = ul.GetDataPointer();
  double*** urdata = ur.GetDataPointer();
  double*** vbdata = vb.GetDataPointer();
  double*** vtdata = vt.GetDataPointer();
  double*** wkdata = wk.GetDataPointer();
  double*** wfdata = wf.GetDataPointer();
  double*** res    = R.GetDataPointer(); //residual, on the right-hand-side of the ODE
  double*** active = Active.GetDataPointer();
  //double*** useful = UsefulG2.GetDataPointer();

  //initialize R to 0
  for(auto it = useful_nodes.begin(); it != useful_nodes.end(); it++)
    res[(*it)[2]][(*it)[1]][(*it)[0]] = 0.0;

  double localflux;

  // Loop through the domain interior, and the right and top ghost layers. For each cell, calculate the
  // numerical flux across the left and lower cell boundaries/interfaces
  for(auto it = useful_nodes.begin(); it != useful_nodes.end(); it++) {

    int i((*it)[0]), j((*it)[1]), k((*it)[2]);

    if(i==ii0 || j==jj0 || k==kk0)
      continue;

    //calculate F_{i-1/2,jk}
    if(k!=kkmax-1 && j!=jjmax-1) {
      localflux = ComputeLocalAdvectionFlux(phir[k][j][i-1], phil[k][j][i], 
                                            urdata[k][j][i-1], uldata[k][j][i]) / dxyz[k][j][i][0];
      if(active[k][j][i-1]) 
        res[k][j][i-1] -= localflux;

      if(active[k][j][i])
        res[k][j][i]   += localflux;
    }

    //calculate G_{i,j-1/2,k}
    if(k!=kkmax-1 && i!=iimax-1) {
      localflux = ComputeLocalAdvectionFlux(phit[k][j-1][i], phib[k][j][i], 
                                            vtdata[k][j-1][i], vbdata[k][j][i]) / dxyz[k][j][i][1];
      if(active[k][j-1][i])
        res[k][j-1][i] -= localflux;

      if(active[k][j][i])
        res[k][j][i]   += localflux;
    }

    //calculate H_{ij,k-1/2}
    if(j!=jjmax-1 && i!=iimax-1) {
      localflux = ComputeLocalAdvectionFlux(phif[k-1][j][i], phik[k][j][i], 
                                            wfdata[k-1][j][i], wkdata[k][j][i]) / dxyz[k][j][i][2];
      if(active[k-1][j][i])
        res[k-1][j][i] -= localflux;

      if(active[k][j][i])
        res[k][j][i]   += localflux;
    }

  }

  // Restore space variables
  delta_xyz.RestoreDataPointerToLocalVector();

  Phil.RestoreDataPointerToLocalVector();
  Phir.RestoreDataPointerToLocalVector();
  Phib.RestoreDataPointerToLocalVector();
  Phit.RestoreDataPointerToLocalVector();
  Phik.RestoreDataPointerToLocalVector();
  Phif.RestoreDataPointerToLocalVector();
  ul.RestoreDataPointerToLocalVector();
  ur.RestoreDataPointerToLocalVector();
  vb.RestoreDataPointerToLocalVector();
  vt.RestoreDataPointerToLocalVector();
  wk.RestoreDataPointerToLocalVector();
  wf.RestoreDataPointerToLocalVector();
  Active.RestoreDataPointerToLocalVector();
  //UsefulG2.RestoreDataPointerToLocalVector();

  R.RestoreDataPointerAndInsert(); //insert
}

//-----------------------------------------------------

double LevelSetOperator::ComputeLocalAdvectionFlux(double phim, double phip, double um, double up)
{
  double flux = 0.0;
  double lam, tol, vmid;

  switch (iod_ls.flux) {

    case LevelSetSchemeData::UPWIND :
      vmid = 0.5*(um+up); 
      if(vmid>0)
        flux = phim*um;
      else if(vmid<0)
        flux = phip*up;
      else
        flux = 0.5*(phim*um + phip*up);
      break;

    case LevelSetSchemeData::LOCAL_LAX_FRIEDRICHS :
      flux = 0.5*(phim*um + phip*up) - 0.5*max(fabs(um),fabs(up))*(phip-phim);
      break;

    case LevelSetSchemeData::ROE :
      lam = fabs(0.5*(um+up));
      tol   = max(fabs(um), fabs(up))*iod_ls.delta;
      if(lam<tol) lam = (lam*lam + tol*tol)/(2*tol);
      flux = 0.5*(phim*um + phip*up) - 0.5*lam*(phip-phim);
      break;

    default:
      print_error("*** Error: Level set flux function not recognized (code = %d).\n", iod_ls.flux);
      exit_mpi();      
  }
  return flux;
} 
//-----------------------------------------------------

void LevelSetOperator::AddSourceTerm(SpaceVariable3D &Phi, SpaceVariable3D &R)
{
  double*** phi = (double***) Phi.GetDataPointer();
  double*** u_x = (double***) dudx.GetDataPointer();
  double*** v_y = (double***) dvdy.GetDataPointer();
  double*** w_z = (double***) dwdz.GetDataPointer();
  double*** res = (double***) R.GetDataPointer();

  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++) 
      for(int i=i0; i<imax; i++)
        res[k][j][i] += phi[k][j][i]*(u_x[k][j][i] + v_y[k][j][i] + w_z[k][j][i]);

  Phi.RestoreDataPointerToLocalVector();
  dudx.RestoreDataPointerToLocalVector();
  dvdy.RestoreDataPointerToLocalVector();
  dwdz.RestoreDataPointerToLocalVector();
  R.RestoreDataPointerAndInsert();
}

//-----------------------------------------------------

void LevelSetOperator::AddSourceTermInBand(SpaceVariable3D &Phi, SpaceVariable3D &R)
{
  double*** phi = (double***) Phi.GetDataPointer();
  double*** u_x = (double***) dudx.GetDataPointer();
  double*** v_y = (double***) dvdy.GetDataPointer();
  double*** w_z = (double***) dwdz.GetDataPointer();
  double*** res = (double***) R.GetDataPointer();

  for(auto it = active_nodes.begin(); it != active_nodes.end(); it++) {
    int i((*it)[0]), j((*it)[1]), k((*it)[2]);
    if(!Phi.IsHere(i,j,k,false)) //not in the interior of the subdomain
      continue; 
    res[k][j][i] += phi[k][j][i]*(u_x[k][j][i] + v_y[k][j][i] + w_z[k][j][i]);
  }

  Phi.RestoreDataPointerToLocalVector();
  dudx.RestoreDataPointerToLocalVector();
  dvdy.RestoreDataPointerToLocalVector();
  dwdz.RestoreDataPointerToLocalVector();
  R.RestoreDataPointerAndInsert();
}


//-----------------------------------------------------

void
LevelSetOperator::ConstructNarrowBandInReinitializer(SpaceVariable3D &Phi)
{
  assert(reinit);
  reinit->ConstructNarrowBand(Phi, Level, UsefulG2, Active, useful_nodes, active_nodes);
}

//-----------------------------------------------------

bool LevelSetOperator::Reinitialize(double time, double dt, int time_step, SpaceVariable3D &Phi,
                                    int special_maxIts, bool must_do)
{
  if(must_do && !reinit) {
    print_error("*** Error: Reinitialization for level set (matid = %d) is requested, but not specified "
                "in the input file.\n", materialid);
    exit_mpi();
  }

  if(!reinit)
    return false; //nothing to do

  if(!isTimeToWrite(time, dt, time_step, iod_ls.reinit.frequency_dt,
                    iod_ls.reinit.frequency, -100.0, must_do))
    return false; //nothing to do (not the right time)

  if(narrow_band) {
//    print("- Reinitializing the level set function (material id: %d), bandwidth = %d.\n", materialid, iod_ls.bandwidth);
    reinit->ReinitializeInBand(Phi, Level, UsefulG2, Active, useful_nodes, active_nodes, special_maxIts);
  } else {
//    print("- Reinitializing the level set function (material id: %d).\n", materialid);
    reinit->ReinitializeFullDomain(Phi, special_maxIts);
  }

  return true;
}

//-----------------------------------------------------
//int debug_counter = 0;
void LevelSetOperator::ReinitializeAfterPhaseTransition(SpaceVariable3D &Phi, vector<Int3> &new_nodes)
{
  if(!reinit) {
    print_error("*** Error: Reinitialization for level set (matid = %d) is requested, but not specified "
                "in the input file.\n", materialid);
    exit_mpi();
  }

  if(narrow_band) {
//    print("- Reinitializing the level set function (id: %d) after phase transition, bandwidth = %d.\n", materialid, iod_ls.bandwidth);

    for(auto it = new_nodes.begin(); it != new_nodes.end(); it++)
      if(std::find(useful_nodes.begin(), useful_nodes.end(), *it) == useful_nodes.end())
        useful_nodes.push_back(*it);

    reinit->ReinitializeInBand(Phi, Level, UsefulG2, Active, useful_nodes, active_nodes);
  } else {
//    print("- Reinitializing the level set function (id: %d) after phase transition.\n", materialid);
    reinit->ReinitializeFullDomain(Phi);
  }

/*
  print_error("Done with reinit. after phase transition!\n");
  std::string str = "Phi"+std::to_string(++debug_counter)+".vtr";
  Phi.WriteToVTRFile(str.c_str());
*/
}

//-----------------------------------------------------

void
LevelSetOperator::PrescribeVelocityFieldForTesting([[maybe_unused]] SpaceVariable3D &V,
                                                   [[maybe_unused]] SpaceVariable3D &Phi,
                                                   [[maybe_unused]] double time)
{
 
#if LEVELSET_TEST == 1
  // rotation of a sloted disk (Problem 4.2.1 of Hartmann, Meinke, Schroder 2008
  Vec5D*** v = (Vec5D***) V.GetDataPointer();
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++) {
 
        double pi = 2.0*acos(0);
        double v1 = fabs(coords[k][j][i][1]-50.0)<1.0e-10 ? 1.0e-10 : 50.0-coords[k][j][i][1];
        double v2 = fabs(coords[k][j][i][0]-50.0)<1.0e-10 ? 1.0e-10 : coords[k][j][i][0]-50.0;
        v[k][j][i][1] = pi/314.0*v1;
        v[k][j][i][2] = pi/314.0*v2;
        v[k][j][i][3] = 0.0;

      }
  coordinates.RestoreDataPointerToLocalVector();
  V.RestoreDataPointerAndInsert();

#elif LEVELSET_TEST == 2
  // vortex deformation of a circle (section 5.2 of Hartmann et al. 2010)
  Vec5D*** v = (Vec5D***) V.GetDataPointer();
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++) {

        double T = 8.0;
        double x = coords[k][j][i][0], y = coords[k][j][i][1];
        double pi = 2.0*acos(0);
        v[k][j][i][1] =  2.0*pow(sin(pi*x),2)*sin(pi*y)*cos(pi*y)*cos(pi*time/T);
        v[k][j][i][2] = -2.0*sin(pi*x)*cos(pi*x)*pow(sin(pi*y),2)*cos(pi*time/T);
        v[k][j][i][3] = 0.0; 

      }
  coordinates.RestoreDataPointerToLocalVector();
  V.RestoreDataPointerAndInsert();

#elif LEVELSET_TEST == 3

  // merging and separation of two circles
  double aa = 1.75, rr = 1.5;
  double s = time<1.0 ? 1.0 : -1.0;
  SpaceVariable3D NPhi(comm, &(dms.ghosted1_3dof)); //inefficient, but OK just for testing level set 

  if(!narrow_band) {
    ComputeNormalDirectionCentralDifferencing_FullDomain(Phi, NPhi);
    Vec5D*** v = (Vec5D***) V.GetDataPointer();
    Vec3D*** normal = (Vec3D***) NPhi.GetDataPointer();
    for(int k=kk0; k<kkmax; k++)
      for(int j=jj0; j<jjmax; j++)
        for(int i=ii0; i<iimax; i++) {
          if(!coordinates.IsHere(i,j,k,false)) {
            v[k][j][i][1] = v[k][j][i][2] = v[k][j][i][3] = 0.0;
            continue;
          }
          v[k][j][i][1] = s*normal[k][j][i][0];
          v[k][j][i][2] = s*normal[k][j][i][1];
          v[k][j][i][3] = 0.0;
        }
    V.RestoreDataPointerAndInsert();
    NPhi.RestoreDataPointerToLocalVector();
  } 
  else {
    ComputeNormalDirectionCentralDifferencing_NarrowBand(Phi, NPhi);
    Vec5D*** v = (Vec5D***) V.GetDataPointer();
    Vec3D*** normal = (Vec3D***) NPhi.GetDataPointer();
    for(auto it = useful_nodes.begin(); it != useful_nodes.end(); it++) {
      int i((*it)[0]), j((*it)[1]), k((*it)[2]);
      if(!coordinates.IsHere(i,j,k,false)) {
        v[k][j][i][1] = v[k][j][i][2] = v[k][j][i][3] = 0.0;
        continue;
      }
      v[k][j][i][1] = s*normal[k][j][i][0];
      v[k][j][i][2] = s*normal[k][j][i][1];
      v[k][j][i][3] = 0.0;
    }
    V.RestoreDataPointerAndInsert();
    NPhi.RestoreDataPointerToLocalVector();
  }
        
/*
  NPhi.WriteToVTRFile("NPhi.vtr");
  exit_mpi();
*/
  NPhi.Destroy();

#endif


}

//-----------------------------------------------------

void
LevelSetOperator::ComputeNormalDirection(SpaceVariable3D &Phi, SpaceVariable3D &NPhi)
{
  if(!narrow_band)
    ComputeNormalDirectionCentralDifferencing_FullDomain(Phi,NPhi);
  else
    ComputeNormalDirectionCentralDifferencing_NarrowBand(Phi,NPhi);
}

//-----------------------------------------------------

void
LevelSetOperator::ComputeNormalDirectionCentralDifferencing_FullDomain(SpaceVariable3D &Phi, SpaceVariable3D &NPhi)
{
  double*** phi   = Phi.GetDataPointer();
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
  Vec3D*** normal = (Vec3D***)NPhi.GetDataPointer();

  double mynorm;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        normal[k][j][i][0] = CentralDifferenceLocal(phi[k][j][i-1],       phi[k][j][i],       phi[k][j][i+1],
                                                 coords[k][j][i-1][0], coords[k][j][i][0], coords[k][j][i+1][0]);
        normal[k][j][i][1] = CentralDifferenceLocal(phi[k][j-1][i],       phi[k][j][i],       phi[k][j+1][i],
                                                 coords[k][j-1][i][1], coords[k][j][i][1], coords[k][j+1][i][1]);
        normal[k][j][i][2] = CentralDifferenceLocal(phi[k-1][j][i],       phi[k][j][i],       phi[k+1][j][i],
                                                 coords[k-1][j][i][2], coords[k][j][i][2], coords[k+1][j][i][2]);

        mynorm = normal[k][j][i].norm();
        if(mynorm!=0.0)
          normal[k][j][i] /= mynorm; 
/*
        if(i==548 && j==6 && k==0) {
          fprintf(stdout,"(%d,%d,%d): normal = %e %e %e. Phi values...\n", i,j,k, normal[k][j][i][0], normal[k][j][i][1],
                  normal[k][j][i][2]);
          fprintf(stdout,"%e  %e  %e\n", phi[k][j+1][i-1], phi[k][j+1][i], phi[k][j+1][i+1]);
          fprintf(stdout,"%e  %e  %e\n", phi[k][j][i-1], phi[k][j][i], phi[k][j][i+1]);
          fprintf(stdout,"%e  %e  %e\n", phi[k][j-1][i-1], phi[k][j-1][i], phi[k][j-1][i+1]);
          fprintf(stdout,"%e  %e  %e\n", phi[k][j-2][i-1], phi[k][j-2][i], phi[k][j-2][i+1]);
        }
*/
      }

  Phi.RestoreDataPointerToLocalVector();
  coordinates.RestoreDataPointerToLocalVector();
  NPhi.RestoreDataPointerAndInsert();
}

//-----------------------------------------------------

void
LevelSetOperator::ComputeNormalDirectionCentralDifferencing_NarrowBand(SpaceVariable3D &Phi, SpaceVariable3D &NPhi)
{
  double*** phi   = Phi.GetDataPointer();
  double*** useful= UsefulG2.GetDataPointer();
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
  Vec3D*** normal = (Vec3D***)NPhi.GetDataPointer();

  double mynorm;
  for(auto it = useful_nodes.begin(); it != useful_nodes.end(); it++) {
    int i((*it)[0]), j((*it)[1]), k((*it)[2]);
    if(!coordinates.IsHere(i,j,k,false))
      continue; //only visit the interior of the subdomain

    // dphi/dx
    if(useful[k][j][i+1] && useful[k][j][i-1])
      normal[k][j][i][0] = CentralDifferenceLocal(phi[k][j][i-1],       phi[k][j][i],       phi[k][j][i+1],
                                               coords[k][j][i-1][0], coords[k][j][i][0], coords[k][j][i+1][0]);
    else if(useful[k][j][i+1])
      normal[k][j][i][0] = (phi[k][j][i+1] - phi[k][j][i])/(coords[k][j][i+1][0] - coords[k][j][i][0]);
    else if(useful[k][j][i-1])
      normal[k][j][i][0] = (phi[k][j][i] - phi[k][j][i-1])/(coords[k][j][i][0] - coords[k][j][i-1][0]);
    else
      normal[k][j][i][0] = 0.0;

    // dphi/dy
    if(useful[k][j+1][i] && useful[k][j+1][i])
      normal[k][j][i][1] = CentralDifferenceLocal(phi[k][j-1][i],       phi[k][j][i],       phi[k][j+1][i],
                                               coords[k][j-1][i][1], coords[k][j][i][1], coords[k][j+1][i][1]);
    else if(useful[k][j+1][i])
      normal[k][j][i][1] = (phi[k][j+1][i] - phi[k][j][i])/(coords[k][j+1][i][1] - coords[k][j][i][1]);
    else if(useful[k][j-1][i])
      normal[k][j][i][1] = (phi[k][j][i] - phi[k][j-1][i])/(coords[k][j][i][1] - coords[k][j-1][i][1]);
    else
      normal[k][j][i][1] = 0.0;

    // dphi/dz
    if(useful[k+1][j][i] && useful[k-1][j][i])
      normal[k][j][i][2] = CentralDifferenceLocal(phi[k-1][j][i],       phi[k][j][i],       phi[k+1][j][i],
                                               coords[k-1][j][i][2], coords[k][j][i][2], coords[k+1][j][i][2]);
    else if(useful[k+1][j][i])
      normal[k][j][i][2] = (phi[k+1][j][i] - phi[k][j][i])/(coords[k+1][j][i][2] - coords[k][j][i][2]);
    else if(useful[k-1][j][i])
      normal[k][j][i][2] = (phi[k][j][i] - phi[k-1][j][i])/(coords[k][j][i][2] - coords[k-1][j][i][2]);
    else
      normal[k][j][i][2] = 0.0;

    mynorm = normal[k][j][i].norm();
    if(mynorm!=0.0)
      normal[k][j][i] /= mynorm; 


/*
    if(i==259 && j==0 && k==0) {
      fprintf(stdout,"(%d,%d,%d): normal = %e %e %e. Phi values...\n", i,j,k, normal[k][j][i][0], normal[k][j][i][1],
              normal[k][j][i][2]);
      fprintf(stdout,"%e  %e  %e\n", phi[k][j+1][i-1], phi[k][j+1][i], phi[k][j+1][i+1]);
      fprintf(stdout,"%e  %e  %e\n", phi[k][j][i-1], phi[k][j][i], phi[k][j][i+1]);
      fprintf(stdout,"%e  %e  %e\n", phi[k][j-1][i-1], phi[k][j-1][i], phi[k][j-1][i+1]);
    }
*/

  }



  Phi.RestoreDataPointerToLocalVector();
  UsefulG2.RestoreDataPointerToLocalVector();
  coordinates.RestoreDataPointerToLocalVector();
  NPhi.RestoreDataPointerAndInsert();
}

//-----------------------------------------------------
void LevelSetOperator::ComputeNormal(SpaceVariable3D &Phi, SpaceVariable3D &NPhi) { //currently only implemented the full domain version
  double*** phi   = Phi.GetDataPointer();
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
  Vec3D*** normal = (Vec3D***)NPhi.GetDataPointer();

  // Compute and store first direvatives, needed for the computation of second order mixed partial derivatives
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        normal[k][j][i][0] = CentralDifferenceLocal(phi[k][j][i-1],       phi[k][j][i],       phi[k][j][i+1],
                                                 coords[k][j][i-1][0], coords[k][j][i][0], coords[k][j][i+1][0]);
        normal[k][j][i][1] = CentralDifferenceLocal(phi[k][j-1][i],       phi[k][j][i],       phi[k][j+1][i],
                                                 coords[k][j-1][i][1], coords[k][j][i][1], coords[k][j+1][i][1]);
        normal[k][j][i][2] = CentralDifferenceLocal(phi[k-1][j][i],       phi[k][j][i],       phi[k+1][j][i],
                                                 coords[k-1][j][i][2], coords[k][j][i][2], coords[k+1][j][i][2]);
      }

  Phi.RestoreDataPointerToLocalVector();
  coordinates.RestoreDataPointerToLocalVector();
  NPhi.RestoreDataPointerAndInsert();
}  

//-----------------------------------------------------

void LevelSetOperator::ComputeUnitNormalAndCurvature(SpaceVariable3D &Phi, SpaceVariable3D &NPhi, SpaceVariable3D &KappaPhi) { // currently only implemented the full domain version

  double*** phi   = Phi.GetDataPointer();
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
  double*** curvature = (double***)KappaPhi.GetDataPointer();
  Vec3D*** normal = (Vec3D***)NPhi.GetDataPointer();

  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
	double phi_xx = SecondOrderDifference(phi[k][j][i-1],       phi[k][j][i],       phi[k][j][i+1],
	                                        coords[k][j][i-1][0], coords[k][j][i][0], coords[k][j][i+1][0]);
	double phi_yy = SecondOrderDifference(phi[k][j-1][i],       phi[k][j][i],       phi[k][j+1][i],
	                                        coords[k][j-1][i][1], coords[k][j][i][1], coords[k][j+1][i][1]);
	double phi_zz = SecondOrderDifference(phi[k-1][j][i],       phi[k][j][i],       phi[k+1][j][i],
	                                        coords[k-1][j][i][2], coords[k][j][i][2], coords[k+1][j][i][2]);
        double phi_xy = CentralDifferenceLocal(/* y-direvative of dphi/dx */
                                               normal[k][j-1][i][0], normal[k][j][i][0], normal[k][j+1][i][0],
                                               coords[k][j-1][i][1], coords[k][j][i][1], coords[k][j+1][i][1]);

        double phi_xz = CentralDifferenceLocal(/* z-direvative of dphi/dx */
                                               normal[k-1][j][i][0], normal[k][j][i][0], normal[k+1][j][i][0], 
                                               coords[k-1][j][i][2], coords[k][j][i][2], coords[k+1][j][i][2]);
        double phi_yz = CentralDifferenceLocal(/* z-direvative of dphi/dy */
                                               normal[k-1][j][i][1], normal[k][j][i][1], normal[k+1][j][i][1],
                                               coords[k-1][j][i][2], coords[k][j][i][2], coords[k+1][j][i][2]);

        double phi_x = normal[k][j][i][0];
        double phi_y = normal[k][j][i][1];
        double phi_z = normal[k][j][i][2];   

        curvature[k][j][i] = phi_x*phi_x*phi_yy - 2.*phi_x*phi_y*phi_xy + phi_y*phi_y*phi_xx
                            +phi_x*phi_x*phi_zz - 2.*phi_x*phi_z*phi_xz + phi_z*phi_z*phi_xx
                            +phi_y*phi_y*phi_zz - 2.*phi_y*phi_z*phi_yz + phi_z*phi_z*phi_yy;  
      }


  double mynorm;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        mynorm = normal[k][j][i].norm();
        if(mynorm!=0.0) {
          curvature[k][j][i] /= (mynorm*mynorm*mynorm); 
          normal[k][j][i] /= mynorm;
        }
      }

  Phi.RestoreDataPointerToLocalVector();
  coordinates.RestoreDataPointerToLocalVector();
  NPhi.RestoreDataPointerAndInsert();
  KappaPhi.RestoreDataPointerAndInsert();
}

//-----------------------------------------------------


