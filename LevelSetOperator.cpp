#include <cfloat> //DBL_MAX
#include <Utils.h>
#include <LevelSetOperator.h>
#include <SpaceOperator.h>
#include <Vector5D.h>
#include <GeoTools.h>
#include <map>

using std::min;
using std::max;
//-----------------------------------------------------

LevelSetOperator::LevelSetOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_,
                                   LevelSetSchemeData &iod_ls_, SpaceOperator &spo)
  : comm(comm_), dm_all(dm_all_), iod(iod_), iod_ls(iod_ls_),
    coordinates(spo.GetMeshCoordinates()),
    delta_xyz(spo.GetMeshDeltaXYZ()),
    volume(spo.GetMeshCellVolumes()),
    rec(comm_, dm_all_, iod_ls_.rec, coordinates, delta_xyz),
    scalar(comm_, &(dm_all_.ghosted1_1dof)),
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
    Phif(comm_, &(dm_all_.ghosted1_1dof))
{

  materialid = iod_ls.materialid;

  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);

  rec.Setup(); //this function requires mesh info (dxyz)

}

//-----------------------------------------------------

LevelSetOperator::~LevelSetOperator()
{

}

//-----------------------------------------------------

void LevelSetOperator::Destroy()
{
  rec.Destroy();
  scalar.Destroy();
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
}

//-----------------------------------------------------

void LevelSetOperator::SetInitialCondition(SpaceVariable3D &Phi)
{
  //get info
  Vec3D***  coords = (Vec3D***)coordinates.GetDataPointer();
  double*** phi    = (double***)Phi.GetDataPointer();

  //first, set phi to be +inf everywhere
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++)
        phi[k][j][i] = DBL_MAX;

  //initialize phi based on plane, cylinder-cones, and sphere.
  MultiInitialConditionsData &ic(iod.ic.multiInitialConditions);

  // planes
  for(auto it=ic.planeMap.dataMap.begin(); it!=ic.planeMap.dataMap.end(); it++) {

    if(it->second->initialConditions.materialid != materialid)
      continue; //not the right one

    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
    Vec3D dir(it->second->nx, it->second->ny, it->second->nz);
    dir /= dir.norm();
    double dist;

    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++) {
          dist = (coords[k][j][i]-x0)*dir;
          if(fabs(dist)<fabs(phi[k][j][i]))
            phi[k][j][i] = -dist; //phi is negative inside the material subdomain, positive outside
        }
  }

  // cylinder-cone
  for(auto it=ic.cylinderconeMap.dataMap.begin(); it!=ic.cylinderconeMap.dataMap.end(); it++) {

    if(it->second->initialConditions.materialid != materialid)
      continue; //not the right one

    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
    Vec3D dir(it->second->nx, it->second->ny, it->second->nz);
    dir /= dir.norm();

    double L = it->second->L; //cylinder height
    double R = it->second->r; //cylinder radius
    double tan_alpha = tan(it->second->opening_angle_degrees/180.0*acos(-1.0));//opening angle
    double Hmax = R/tan_alpha;
    double H = min(it->second->cone_height, Hmax); //cone's height
 
    // define the geometry
    vector<pair<Vec3D, Vec3D> > lineSegments; //boundary of the cylinder-cone
    Vec3D p0(0.0, 0.0, 0.0); //x0
    Vec3D p1(0.0, R, 0.0);
    Vec3D p2(L, R, 0.0);
    Vec3D p3(L+H, (Hmax-H)*tan_alpha, 0.0);
    Vec3D p4(L+H, 0.0, 0.0); 
    lineSegments.push_back(std::make_pair(p0,p1));
    lineSegments.push_back(std::make_pair(p1,p2));
    lineSegments.push_back(std::make_pair(p2,p3));
    lineSegments.push_back(std::make_pair(p3,p4)); //may be a degenerate edge (point) but fine

    double x, r, toto, dist = DBL_MAX;
    Vec3D xp;

    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++) {
    
          x = (coords[k][j][i]-x0)*dir;
          r = (coords[k][j][i] - x0 - x*dir).norm();
          xp = Vec3D(x,r,0.0);

          //calculate unsigned distance from node to the boundary of the cylinder-cone
          for(int l=0; l<lineSegments.size(); l++)
            dist = min(dist, GeoTools::GetShortestDistanceFromPointToLineSegment(xp, 
                                 lineSegments[l].first, lineSegments[l].second, toto));

          //figure out the sign, and update phi
          if(dist<fabs(phi[k][j][i])) {
            if( (x>0 && x<L && r<R) || (x>=L && x<L+H && r<(L+H-x)*tan_alpha) ) //inside
              phi[k][j][i] = -dist;
            else 
              phi[k][j][i] = dist; 
          }
        }
  }

  // spheres 
  for(auto it=ic.sphereMap.dataMap.begin(); it!=ic.sphereMap.dataMap.end(); it++) {

    if(it->second->initialConditions.materialid != materialid)
      continue; //not the right one

    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
    double dist;

    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++) {
          dist = (coords[k][j][i]-x0).norm() - it->second->radius;
          if(fabs(dist)<fabs(phi[k][j][i]))
            phi[k][j][i] = dist; //>0 outside the sphere
        }
  }

  coordinates.RestoreDataPointerToLocalVector();
  Phi.RestoreDataPointerAndInsert();

  ApplyBoundaryConditions(Phi);
}

//-----------------------------------------------------

// Apply boundary conditions by populating ghost cells of Phi
void LevelSetOperator::ApplyBoundaryConditions(SpaceVariable3D &Phi)
{
  double*** phi = (double***) Phi.GetDataPointer();

  int NX, NY, NZ;
  Phi.GetGlobalSize(&NX, &NY, &NZ);

  //! Left boundary
  if(ii0==-1) { 
    switch (iod.mesh.bc_x0) {
      case MeshData::INLET :
      case MeshData::OUTLET :
      case MeshData::WALL :
      case MeshData::SYMMETRY :
        for(int k=k0; k<kmax; k++)
          for(int j=j0; j<jmax; j++)
            phi[k][j][ii0] = phi[k][j][ii0+1];
        break;
      default :
        print_error("Error: Level set boundary condition at x=x0 cannot be specified!\n");
        exit_mpi();
    }
  }

  //! Right boundary
  if(iimax==NX+1) { 
    switch (iod.mesh.bc_xmax) {
      case MeshData::INLET :
      case MeshData::OUTLET :
      case MeshData::WALL :
      case MeshData::SYMMETRY :
        for(int k=k0; k<kmax; k++)
          for(int j=j0; j<jmax; j++)
            phi[k][j][iimax-1] = phi[k][j][iimax-2];
        break;
      default :
        print_error("Error: Level set boundary condition at x=xmax cannot be specified!\n");
        exit_mpi();
    }
  }

  //! Bottom boundary
  if(jj0==-1) { 
    switch (iod.mesh.bc_y0) {
      case MeshData::INLET :
      case MeshData::OUTLET :
      case MeshData::WALL :
      case MeshData::SYMMETRY :
        for(int k=k0; k<kmax; k++)
          for(int i=i0; i<imax; i++)
            phi[k][jj0][i] = phi[k][jj0+1][i];
        break;
      default :
        print_error("Error: Level set boundary condition at y=y0 cannot be specified!\n");
        exit_mpi();
    }
  }

  //! Top boundary
  if(jjmax==NY+1) { 
    switch (iod.mesh.bc_ymax) {
      case MeshData::INLET :
      case MeshData::OUTLET :
      case MeshData::WALL :
      case MeshData::SYMMETRY :
        for(int k=k0; k<kmax; k++)
          for(int i=i0; i<imax; i++)
            phi[k][jjmax-1][i] = phi[k][jjmax-2][i]; 
        break;
      default :
        print_error("Error: Level set boundary condition at y=ymax cannot be specified!\n");
        exit_mpi();
    }
  }

  //! Back boundary (z min)
  if(kk0==-1) { 
    switch (iod.mesh.bc_z0) {
      case MeshData::INLET :
      case MeshData::OUTLET :
      case MeshData::WALL :
      case MeshData::SYMMETRY :
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++)
            phi[kk0][j][i] = phi[kk0+1][j][i]; 
        break;
      default :
        print_error("Error: Level set boundary condition at z=z0 cannot be specified!\n");
        exit_mpi();
    }
  }

  //! Front boundary (z max)
  if(kkmax==NZ+1) { 
    switch (iod.mesh.bc_zmax) {
      case MeshData::INLET :
      case MeshData::OUTLET :
      case MeshData::WALL :
      case MeshData::SYMMETRY :
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++)
            phi[kkmax-1][j][i] = phi[kkmax-2][j][i];
        break;
      default :
        print_error("Error: Boundary condition at z=zmax cannot be specified!\n");
        exit_mpi();
    }
  }

  Phi.RestoreDataPointerAndInsert();
}

//-----------------------------------------------------

void LevelSetOperator::ComputeResidual(SpaceVariable3D &V, SpaceVariable3D &Phi, SpaceVariable3D &R)
{
  Reconstruct(V, Phi); // => ul, ur, dudx, vb, vt, dvdy, wk, wf, dwdz, Phil, Phir, Phib, Phit, Phik, Phif

  ComputeAdvectionFlux(R);  //Advection flux, on the right-hand-side

  AddSourceTerm(Phi, R);
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
  rec.ReconstructIn1D(0/*x-dir*/, scalar, ul, ur, &dudx);

  // Reconstruction: y-velocity
  s = (double***) scalar.GetDataPointer();
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++)
        s[k][j][i] = v[k][j][i][2]; //y-velocity
  scalar.RestoreDataPointerToLocalVector(); //no need to communicate w/ neighbors
  rec.ReconstructIn1D(1/*y-dir*/, scalar, vb, vt, &dvdy);

  // Reconstruction: z-velocity
  s = (double***) scalar.GetDataPointer();
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++)
        s[k][j][i] = v[k][j][i][3]; //z-velocity
  scalar.RestoreDataPointerToLocalVector(); //no need to communicate w/ neighbors
  rec.ReconstructIn1D(2/*z-dir*/, scalar, wk, wf, &dwdz);

  // Reconstruction: Phi
  rec.Reconstruct(Phi, Phil, Phir, Phib, Phit, Phik, Phif);

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

double LevelSetOperator::ComputeLocalAdvectionFlux(double phim, double phip, double um, double up)
{
  double flux = 0.5*(phim*um + phip*up);
  double lam, tol;

  switch (iod_ls.flux) {
    case LevelSetSchemeData::LOCAL_LAX_FRIEDRICHS :
      flux -= 0.5*max(fabs(um),fabs(up))*(phip-phim);
      break;
    case LevelSetSchemeData::ROE :
      lam = abs(0.5*(um+up));
      tol   = max(fabs(um), fabs(up))*iod_ls.delta;
      if(lam<tol) lam = (lam*lam + tol*tol)/(2*tol);
      flux -= 0.5*lam*(phip-phim);
      break;
    default:
      print_error("Error: Level set flux function not recognized (code = %d).\n", iod_ls.flux);
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

void LevelSetOperator::Reinitialize(SpaceVariable3D &Phi)
{
  print_error("Error: Level set reinitialization has not been implemented yet.\n");
  exit_mpi();
}

//-----------------------------------------------------




