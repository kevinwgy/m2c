#include <cfloat> //DBL_MAX
#include <Utils.h>
#include <LevelSetOperator.h>
#include <SpaceOperator.h>
#include <Vector5D.h>
#include <GeoTools.h>
#include <DistancePointToSpheroid.h>
#include <map>

#ifdef LEVELSET_TEST
  #include <Vector2D.h>
#endif

using std::min;
using std::max;
//-----------------------------------------------------

LevelSetOperator::LevelSetOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_,
                                   LevelSetSchemeData &iod_ls_, SpaceOperator &spo)
  : comm(comm_), iod(iod_), iod_ls(iod_ls_),
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

  rec.Setup(spo.GetPointerToInnerGhostNodes(),
            spo.GetPointerToOuterGhostNodes()); //this function requires mesh info (dxyz)

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


  //TODO: Add this feature later
  if(iod.ic.specified[IcData::LEVELSET] || iod.ic.specified[IcData::MATERIALID]) {
    print_error("*** Error: Currently, cannot read level set and material id from user-specified data file.\n");
    exit_mpi();
  }

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

    double x, r, toto, dist;
    Vec3D xp;

    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++) {
    
          x = (coords[k][j][i]-x0)*dir;
          r = (coords[k][j][i] - x0 - x*dir).norm();
          xp = Vec3D(x,r,0.0);

          //calculate unsigned distance from node to the boundary of the cylinder-cone
          dist = DBL_MAX;
          for(int l=0; l<lineSegments.size(); l++)
            dist = min(dist, GeoTools::GetShortestDistanceFromPointToLineSegment(xp, 
                                 lineSegments[l].first, lineSegments[l].second, toto));

          //figure out the sign, and update phi
          if(dist<fabs(phi[k][j][i])) {
            if( (x>0 && x<L && r<R) || (x>=L && x<L+H && r<(L+Hmax-x)*tan_alpha) ) //inside
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


  // spheroids
  for(auto it=ic.spheroidMap.dataMap.begin(); it!=ic.spheroidMap.dataMap.end(); it++) {

    if(it->second->initialConditions.materialid != materialid)
      continue; //not the right one

    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
    Vec3D axis(it->second->axis_x, it->second->axis_y, it->second->axis_z);

    GeoTools::DistanceFromPointToSpheroid distCal(x0, axis, it->second->length, it->second->diameter);

    double dist;

    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++) {
          dist = distCal.Calculate(coords[k][j][i]); //>0 outside the spheroid
          if(fabs(dist)<fabs(phi[k][j][i]))
            phi[k][j][i] = dist; //>0 outside the spheroid
        }
  }


#if LEVELSET_TEST == 1
  // rotation of a sloted disk (Problem 4.2.1 of Hartmann, Meinke, Schroder 2008
  Vec2D x0(50.0,75.0); //center of disk
  double r = 15.0; //radius of disk
  double wid = 5.0, len = 25.0; //slot

  // derivated params
  Vec2D x1(x0[0]-0.5*wid, x0[1]-sqrt(r*r - pow(0.5*wid ,2))); //lowerleft corner of slot
  Vec2D x2(x0[0]-0.5*wid, x0[1]+(len-r)); //upperleft corner of slot
  Vec2D x3(x0[0]+0.5*wid, x0[1]+(len-r)); //upperright corner of slot
  Vec2D x4(x0[0]+0.5*wid, x0[1]-sqrt(r*r - pow(0.5*wid ,2))); //lowerright corner of slot
  double alpha = atan((x4[0]-x0[0])/(x0[1]-x4[1])); //half of the `opening angle'
  //fprintf(stderr,"alpha = %e.\n",alpha);

  double dist;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        Vec2D x(coords[k][j][i][0], coords[k][j][i][1]);

        if(x[1]<=x1[1] && atan(fabs((x[0]-x0[0])/(x[1]-x0[1])))<=alpha) { //Zone 1: Outside-bottom
          if(x[0]<=x0[0])
            phi[k][j][i] = (x-x1).norm();
          else
            phi[k][j][i] = (x-x4).norm();
        }

        else if((x-x0).norm()>=r) { //Zone 2: Outside-rest
          phi[k][j][i] = (x-x0).norm() - r;
        }

        else if(x[0]<=x4[0] && x[0]>=x1[0] && x[1]<=x2[1]) { //Zone 3: Slot
          phi[k][j][i] = min(x4[0]-x[0], min(x[0]-x1[0], x2[1]-x[1]));
        }

        else if(x[0]<=x2[0] && x[1]<=x2[1]) { //Zone 4: Inside-left
          phi[k][j][i] = -min(r-(x-x0).norm(), x2[0]-x[0]); 
        }

        else if(x[0]>=x3[0] && x[1]<=x3[1]) { //Zone 5: Inside-right
          phi[k][j][i] = -min(r-(x-x0).norm(), x[0]-x3[0]);
        }

        else if(x[0]<=x2[0]) { //Zone 6: Inside-topleft
          phi[k][j][i] = -min(r-(x-x0).norm(), (x-x2).norm());
        }

        else if(x[0]<=x3[0]) { //Zone 7: Inside-top
          phi[k][j][i] = -min(r-(x-x0).norm(), x[1]-x2[1]);
        }

        else { //Zone 8: Inside-topright
          phi[k][j][i] = -min(r-(x-x0).norm(), (x-x3).norm());
        }

      }

#endif


  coordinates.RestoreDataPointerToLocalVector();
  Phi.RestoreDataPointerAndInsert();

  ApplyBoundaryConditions(Phi);
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
  //! Left boundary
  if(ii0==-1) { 
    switch (iod_ls.bc_x0) {

      case LevelSetSchemeData::ZERO_NEUMANN :
        for(int k=k0; k<kmax; k++)
          for(int j=j0; j<jmax; j++)
            phi[k][j][ii0] = phi[k][j][ii0+1];
        break;

      case LevelSetSchemeData::LINEAR_EXTRAPOLATION :
        for(int k=k0; k<kmax; k++)
          for(int j=j0; j<jmax; j++) {
            if(ii0+2<NX) { //ii0+2 is within physical domain
              r  = coords[k][j][ii0][0];
              r1 = coords[k][j][ii0+1][0];  f1 = phi[k][j][ii0+1];
              r2 = coords[k][j][ii0+2][0];  f2 = phi[k][j][ii0+2];
              phi[k][j][ii0] = f1 + (f2-f1)/(r2-r1)*(r-r1);
            } else
              phi[k][j][ii0] = phi[k][j][ii0+1];
          }
        break;

      case LevelSetSchemeData::NONE :
        break;

      default :
        print_error("*** Error: Level set boundary condition at x=x0 cannot be specified!\n");
        exit_mpi();
    }
  }

  //! Right boundary
  if(iimax==NX+1) { 
    switch (iod_ls.bc_xmax) {

      case LevelSetSchemeData::ZERO_NEUMANN :
        for(int k=k0; k<kmax; k++)
          for(int j=j0; j<jmax; j++)
            phi[k][j][iimax-1] = phi[k][j][iimax-2];
        break;

      case LevelSetSchemeData::LINEAR_EXTRAPOLATION :
        for(int k=k0; k<kmax; k++)
          for(int j=j0; j<jmax; j++) {
            if(iimax-3>=0) { //iimax-3 is within physical domain
              r  = coords[k][j][iimax-1][0];
              r1 = coords[k][j][iimax-2][0];  f1 = phi[k][j][iimax-2];
              r2 = coords[k][j][iimax-3][0];  f2 = phi[k][j][iimax-3];
              phi[k][j][iimax-1] = f1 + (f2-f1)/(r2-r1)*(r-r1);
            } else
              phi[k][j][iimax-1] = phi[k][j][iimax-2];
          }
        break;
 
      case LevelSetSchemeData::NONE :
        break;

      default :
        print_error("*** Error: Level set boundary condition at x=xmax cannot be specified!\n");
        exit_mpi();
    }
  }

  //! Bottom boundary
  if(jj0==-1) { 
    switch (iod_ls.bc_y0) {

      case LevelSetSchemeData::ZERO_NEUMANN :
        for(int k=k0; k<kmax; k++)
          for(int i=i0; i<imax; i++)
            phi[k][jj0][i] = phi[k][jj0+1][i];
        break;

      case LevelSetSchemeData::LINEAR_EXTRAPOLATION :
        for(int k=k0; k<kmax; k++)
          for(int i=i0; i<imax; i++) {
            if(jj0+2<NY) { //jj0+2 is within physical domain
              r  = coords[k][jj0][i][1];
              r1 = coords[k][jj0+1][i][1];  f1 = phi[k][jj0+1][i];
              r2 = coords[k][jj0+2][i][1];  f2 = phi[k][jj0+2][i];
              phi[k][jj0][i] = f1 + (f2-f1)/(r2-r1)*(r-r1);
            } else
              phi[k][jj0][i] = phi[k][jj0+1][i];
          }
        break;

      case LevelSetSchemeData::NONE :
        break;

      default :
        print_error("*** Error: Level set boundary condition at y=y0 cannot be specified!\n");
        exit_mpi();
    }
  }

  //! Top boundary
  if(jjmax==NY+1) { 
    switch (iod_ls.bc_ymax) {

      case LevelSetSchemeData::ZERO_NEUMANN :
        for(int k=k0; k<kmax; k++)
          for(int i=i0; i<imax; i++)
            phi[k][jjmax-1][i] = phi[k][jjmax-2][i]; 
        break;

      case LevelSetSchemeData::LINEAR_EXTRAPOLATION :
        for(int k=k0; k<kmax; k++)
          for(int i=i0; i<imax; i++) {
            if(jjmax-3>=0) { //jjmax-3 is within physical domain
              r  = coords[k][jjmax-1][i][1];
              r1 = coords[k][jjmax-2][i][1];  f1 = phi[k][jjmax-2][i];
              r2 = coords[k][jjmax-3][i][1];  f2 = phi[k][jjmax-3][i];
              phi[k][jjmax-1][i] = f1 + (f2-f1)/(r2-r1)*(r-r1);
            } else
              phi[k][jjmax-1][i] = phi[k][jjmax-2][i];
          }
        break;
 
      case LevelSetSchemeData::NONE :
        break;

      default :
        print_error("*** Error: Level set boundary condition at y=ymax cannot be specified!\n");
        exit_mpi();
    }
  }

  //! Back boundary (z min)
  if(kk0==-1) { 
    switch (iod_ls.bc_z0) {

      case LevelSetSchemeData::ZERO_NEUMANN :
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++)
            phi[kk0][j][i] = phi[kk0+1][j][i]; 
        break;

      case LevelSetSchemeData::LINEAR_EXTRAPOLATION :
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {
            if(kk0+2<NZ) { //kk0+2 is within physical domain
              r  = coords[kk0][j][i][2];
              r1 = coords[kk0+1][j][i][2];  f1 = phi[kk0+1][j][i];
              r2 = coords[kk0+2][j][i][2];  f2 = phi[kk0+2][j][i];
              phi[kk0][j][i] = f1 + (f2-f1)/(r2-r1)*(r-r1);
            } else
              phi[kk0][j][i] = phi[kk0+1][j][i]; 
          }
        break;

      case LevelSetSchemeData::NONE :
        break;

      default :
        print_error("*** Error: Level set boundary condition at z=z0 cannot be specified!\n");
        exit_mpi();
    }
  }

  //! Front boundary (z max)
  if(kkmax==NZ+1) { 
    switch (iod_ls.bc_zmax) {

      case LevelSetSchemeData::ZERO_NEUMANN :
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++)
            phi[kkmax-1][j][i] = phi[kkmax-2][j][i];
        break;

      case LevelSetSchemeData::LINEAR_EXTRAPOLATION :
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {
            if(kkmax-3>=0) { //kkmax-3 is within physical domain
              r  = coords[kkmax-1][j][i][2];
              r1 = coords[kkmax-2][j][i][2];  f1 = phi[kkmax-2][j][i];
              r2 = coords[kkmax-3][j][i][2];  f2 = phi[kkmax-3][j][i];
              phi[kkmax-1][j][i] = f1 + (f2-f1)/(r2-r1)*(r-r1);
            } else
              phi[kkmax-1][j][i] = phi[kkmax-2][j][i];
          }
        break;

      case LevelSetSchemeData::NONE :
        break;

      default :
        print_error("*** Error: Boundary condition at z=zmax cannot be specified!\n");
        exit_mpi();
    }
  }

  Phi.RestoreDataPointerAndInsert();

  coordinates.RestoreDataPointerToLocalVector();
}

//-----------------------------------------------------

void LevelSetOperator::ComputeResidual(SpaceVariable3D &V, SpaceVariable3D &Phi, SpaceVariable3D &R)
{

#ifdef LEVELSET_TEST
  PrescribeVelocityFieldForTesting(V);
#endif

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
      lam = abs(0.5*(um+up));
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

void LevelSetOperator::Reinitialize(SpaceVariable3D &Phi)
{
  print_error("*** Error: Level set reinitialization has not been implemented yet.\n");
  exit_mpi();
}

//-----------------------------------------------------

void LevelSetOperator::PrescribeVelocityFieldForTesting(SpaceVariable3D &V)
{
  Vec5D*** v = (Vec5D***) V.GetDataPointer();
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();

  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++) {
 
#if LEVELSET_TEST == 1
        // rotation of a sloted disk (Problem 4.2.1 of Hartmann, Meinke, Schroder 2008
        double pi = 2.0*acos(0);
        double v1 = fabs(coords[k][j][i][1]-50.0)<1.0e-10 ? 1.0e-10 : 50.0-coords[k][j][i][1];
        double v2 = fabs(coords[k][j][i][0]-50.0)<1.0e-10 ? 1.0e-10 : coords[k][j][i][0]-50.0;
        v[k][j][i][1] = pi/314.0*v1;
        v[k][j][i][2] = pi/314.0*v2;
        v[k][j][i][3] = 0.0;
#endif

      }

  coordinates.RestoreDataPointerToLocalVector();

  //no need to communicate as the same formula is used in all the subdomains
  V.RestoreDataPointerToLocalVector(); 
}

//-----------------------------------------------------





