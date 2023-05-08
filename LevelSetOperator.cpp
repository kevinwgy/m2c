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

void LevelSetOperator::SetInitialCondition(SpaceVariable3D &Phi, 
                                           unique_ptr<vector<unique_ptr<EmbeddedBoundaryDataSet> > > EBDS,
                                           vector<pair<int,int> > *surf_and_color)
{
  //get info
  Vec3D***  coords = (Vec3D***)coordinates.GetDataPointer();
  double*** phi    = (double***)Phi.GetDataPointer();

  //first, set phi to be a large number everywhere
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++)
        phi[k][j][i] = domain_diagonal;


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

    if(it->second->inclusion != PlaneData::OVERRIDE) {
      print_error("*** Error: Level set initialization only supports Inclusion = Override at the moment.\n");
      exit_mpi();
    }

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

    if(it->second->side != CylinderConeData::INTERIOR) {
      print_error("*** Error: Level set initialization only supports Side = Interior at the moment.\n");
      exit_mpi();
    }

    if(it->second->inclusion != CylinderConeData::OVERRIDE) {
      print_error("*** Error: Level set initialization only supports Inclusion = Override at the moment.\n");
      exit_mpi();
    }

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
          for(int l=0; l<(int)lineSegments.size(); l++)
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


  // cylinder with spherical cap(s)
  for(auto it=ic.cylindersphereMap.dataMap.begin(); it!=ic.cylindersphereMap.dataMap.end(); it++) {

    if(it->second->initialConditions.materialid != materialid)
      continue; //not the right one

    if(it->second->side != CylinderSphereData::INTERIOR) {
      print_error("*** Error: Level set initialization only supports Side = Interior at the moment.\n");
      exit_mpi();
    }

    if(it->second->inclusion != CylinderSphereData::OVERRIDE) {
      print_error("*** Error: Level set initialization only supports Inclusion = Override at the moment.\n");
      exit_mpi();
    }


    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
    Vec3D dir(it->second->nx, it->second->ny, it->second->nz);
    dir /= dir.norm();

    double L = it->second->L; //cylinder height
    double Lhalf = L/2.0;
    double R = it->second->r; //cylinder radius
    bool front_cap = (it->second->front_cap == CylinderSphereData::On);
    bool back_cap = (it->second->back_cap == CylinderSphereData::On);
 
    // define the geometry
    vector<pair<Vec3D, Vec3D> > lineSegments; //boundary of the cylinder
    Vec3D p0(-Lhalf, 0.0, 0.0); //x0
    Vec3D p1(-Lhalf, R, 0.0);
    Vec3D p2(Lhalf, R, 0.0);
    Vec3D p3(Lhalf, 0.0, 0.0);
    lineSegments.push_back(std::make_pair(p1,p2));
    if(!back_cap)
      lineSegments.push_back(std::make_pair(p0,p1));
    if(!front_cap)
      lineSegments.push_back(std::make_pair(p2,p3));

    double x, r, rf(0.0), rb(0.0), toto, dist;
    Vec3D xf = x0 + Lhalf*dir; //center of the front face
    Vec3D xb = x0 - Lhalf*dir; //center of the back face
    Vec3D xp;

    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++) {
    
          x = (coords[k][j][i]-x0)*dir;
          r = (coords[k][j][i] - x0 - x*dir).norm();
          xp = Vec3D(x,r,0.0);

          //calculate unsigned distance from node to the boundary of the cylinder
          if(back_cap && x<=-Lhalf) {
            rb = (coords[k][j][i] - xb).norm();
            dist = fabs(rb - R);
          } 
          else if(front_cap && x>= Lhalf) {
            rf = (coords[k][j][i] - xf).norm();
            dist = fabs(rf - R);
          }
          else {
            dist = DBL_MAX;
            for(int l=0; l<(int)lineSegments.size(); l++)
              dist = min(dist, GeoTools::GetShortestDistanceFromPointToLineSegment(xp, 
                                   lineSegments[l].first, lineSegments[l].second, toto));
          } 

          //figure out the sign, and update phi
          if(dist<fabs(phi[k][j][i])) {
            if( (x>-Lhalf && x<Lhalf && r<R) || (back_cap && x<=-Lhalf && rb<R) ||
                (front_cap && x>=Lhalf && rf<R) ) //inside
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

    if(it->second->side != SphereData::INTERIOR) {
      print_error("*** Error: Level set initialization only supports Side = Interior at the moment.\n");
      exit_mpi();
    }

    if(it->second->inclusion != SphereData::OVERRIDE) {
      print_error("*** Error: Level set initialization only supports Inclusion = Override at the moment.\n");
      exit_mpi();
    }

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


  // parallelepipeds 
  for(auto it=ic.parallelepipedMap.dataMap.begin(); it!=ic.parallelepipedMap.dataMap.end(); it++) {

    if(it->second->initialConditions.materialid != materialid)
      continue; //not the right one

    if(it->second->side != ParallelepipedData::INTERIOR) {
      print_error("*** Error: Level set initialization only supports Side = Interior at the moment.\n");
      exit_mpi();
    }

    if(it->second->inclusion != ParallelepipedData::OVERRIDE) {
      print_error("*** Error: Level set initialization only supports Inclusion = Override at the moment.\n");
      exit_mpi();
    }

    Vec3D x0(it->second->x0, it->second->y0, it->second->z0);
    Vec3D oa(it->second->ax, it->second->ay, it->second->az); oa -= x0;
    Vec3D ob(it->second->bx, it->second->by, it->second->bz); ob -= x0;
    Vec3D oc(it->second->cx, it->second->cy, it->second->cz); oc -= x0;

    if(oa.norm()==0 || ob.norm()==0 || oc.norm()==0 || (oa^ob)*oc<=0.0) {
      print_error("*** Error: Detected error in a user-specified parallelepiped. "
                  "Overlapping vertices or violation of right-hand rule.\n");
      exit_mpi();
    }

    GeoTools::DistanceFromPointToParallelepiped distCal(x0, oa, ob, oc);

    double dist;

    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++) {
          dist = distCal.Calculate(coords[k][j][i]); //>0 outside the spheroid
          if(fabs(dist)<fabs(phi[k][j][i]))
            phi[k][j][i] = dist; //>0 outside the spheroid
        }
  }


  // spheroids
  for(auto it=ic.spheroidMap.dataMap.begin(); it!=ic.spheroidMap.dataMap.end(); it++) {

    if(it->second->initialConditions.materialid != materialid)
      continue; //not the right one

    if(it->second->side != SpheroidData::INTERIOR) {
      print_error("*** Error: Level set initialization only supports Side = Interior at the moment.\n");
      exit_mpi();
    }

    if(it->second->inclusion != SpheroidData::OVERRIDE) {
      print_error("*** Error: Level set initialization only supports Inclusion = Override at the moment.\n");
      exit_mpi();
    }

    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
    Vec3D axis(it->second->axis_x, it->second->axis_y, it->second->axis_z);

    GeoTools::DistanceFromPointToSpheroid distCal(x0, axis, it->second->semi_length, it->second->radius);

    double dist;

    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++) {
          dist = distCal.Calculate(coords[k][j][i]); //>0 outside the spheroid
          if(fabs(dist)<fabs(phi[k][j][i]))
            phi[k][j][i] = dist; //>0 outside the spheroid
        }
  }

  // arbitrary enclosures
  bool active_arbitrary_enclosure = false;
  for(auto&& enclosure : ic.enclosureMap.dataMap) {

    if(enclosure.second->initialConditions.materialid != materialid)
      continue; //not the right one

    if(enclosure.second->inclusion != UserSpecifiedEnclosureData::OVERRIDE) {
      print_error("*** Error: Level set initialization only supports Inclusion = Override at the moment.\n");
      exit_mpi();
    }

    Phi.RestoreDataPointerAndInsert();
    coordinates.RestoreDataPointerToLocalVector();
    bool active = ApplyInitialConditionWithinEnclosure(*enclosure.second, Phi);
    phi = Phi.GetDataPointer();
    coords = (Vec3D***)coordinates.GetDataPointer();

    if(active)
      active_arbitrary_enclosure = true;
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
  //fprintf(stdout,"alpha = %e.\n",alpha);

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

#elif LEVELSET_TEST == 2

  double r = 0.15;
  Vec2D x0(0.5, 0.75);
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        Vec2D x(coords[k][j][i][0], coords[k][j][i][1]);
 
        phi[k][j][i] = (x-x0).norm() - r;

      }

#elif LEVELSET_TEST == 3

  // merging and separation of two circles (4.2.3 of Hartmann et al. 2008)
   
  double a = 1.75, r = 1.5;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        Vec2D x(coords[k][j][i][0], coords[k][j][i][1]);
        double phi1 = sqrt((x[0]-a)*(x[0]-a) + x[1]*x[1]) - r;
        double phi2 = sqrt((x[0]+a)*(x[0]+a) + x[1]*x[1]) - r;

        if(phi1<=0)
          phi[k][j][i] = phi1;
        else if(phi2<=0)
          phi[k][j][i] = phi2;
        else
          phi[k][j][i] = min(phi1, phi2);

      }

#endif


  int ebds_counter = 0;

  if(EBDS && surf_and_color) {
    if(!reinit) {
      print_error("*** Error: Material %d is tracked by both a level set function and an embedded boundary. In this case, \n"
                  "           level set reinitialization must be turned on.\n", materialid);
      exit_mpi();
    }    

    if(iod_ls.reinit.firstLayerTreatment != LevelSetReinitializationData::FIXED)
      print_warning("Warning: Material %d is tracked by both a level set function and an embeddeed boundary. In this case,\n"
                    "         it may be better to 'fix' the first layer nodes when reinitializing the level set.\n",
                    materialid);


    for(auto&& mypair : *surf_and_color) {

      int surf = mypair.first;
      int mycolor = mypair.second;

      int nLayer = (*EBDS)[surf]->Phi_nLayer;
      if(nLayer<2) 
        print_warning("Warning: Material %d is tracked by both a level set function and an embeddeed boundary.\n"
                      "         The intersector should compute unsigned distance ('phi') for at least two layers\n"
                      "         of nodes. Currently, this value is set to %d.\n", (*EBDS)[surf]->Phi_nLayer); 

      double*** color = (*EBDS)[surf]->Color_ptr->GetDataPointer();
      double*** psi   = (*EBDS)[surf]->Phi_ptr->GetDataPointer();
      double*** cpi   = (*EBDS)[surf]->ClosestPointIndex_ptr->GetDataPointer();
     
      for(int k=k0; k<kmax; k++)
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {
    
/*
            //Check if this node is within nLayer of the interface
            bool near = false;
            for(int k1 = k-nLayer; k1 <= k+nLayer; k1++)
              for(int j1 = j-nLayer; j1 <= j+nLayer; j1++)
                for(int i1 = i-nLayer; i1 <= i+nLayer; i1++) {
                  if(!coordinates.IsHereOrInternalGhost(i1,j1,k1))
                    continue;
                  if(color[k1][j1][i1] == mycolor) {
                    near = true;
                    goto DONE;  
                  }
                }
DONE:

            if(!near) 
              continue; //nothing to do
*/
            if(cpi[k][j][i]<0 && color[k][j][i] != mycolor)
              continue;

            ebds_counter++;

            // determine phi at [k][j][i]
            if(color[k][j][i] == mycolor) {//"inside"
              if(phi[k][j][i]>=0) 
                phi[k][j][i] = -psi[k][j][i];
              else
                phi[k][j][i] = std::max(phi[k][j][i], -psi[k][j][i]);
            } else if(color[k][j][i] == 0) {//occluded
              phi[k][j][i] = 0.0; //note that psi may not be 0, because surface has finite thickness
            } else { //"outside"
              if(phi[k][j][i]>=0)
                phi[k][j][i] = std::min(phi[k][j][i], psi[k][j][i]);
              else //this is really bad... the user should avoid this...
                phi[k][j][i] = std::max(phi[k][j][i], -psi[k][j][i]);
            }
          }


      (*EBDS)[surf]->Color_ptr->RestoreDataPointerToLocalVector();
      (*EBDS)[surf]->Phi_ptr->RestoreDataPointerToLocalVector();
      (*EBDS)[surf]->ClosestPointIndex_ptr->RestoreDataPointerToLocalVector();
    }
  }

  coordinates.RestoreDataPointerToLocalVector();
  Phi.RestoreDataPointerAndInsert();

  ApplyBoundaryConditions(Phi);

  if(iod_ls.bandwidth < INT_MAX) {//narrow-band level set method
    assert(reinit);
    reinit->ConstructNarrowBand(Phi, Level, UsefulG2, Active, useful_nodes, active_nodes);
  }

  MPI_Allreduce(MPI_IN_PLACE, &ebds_counter, 1, MPI_INT, MPI_SUM, comm);

  if(ebds_counter>0) { 
    //reinit must have been created (not NULL)
    print("- Updated phi (material id: %d) at %d nodes based on embedded boundary. Going to reinitialize phi.\n",
          materialid, ebds_counter);
    Reinitialize(0.0, 1.0, 0.0, Phi, 600, true/*must do*/); //the first 3 inputs are irrelevant because of "must do"
  }
  else if(active_arbitrary_enclosure && reinit) {
    //if reinit is NULL, should have reinitialized phi (full-domain) in ApplyInitialConditionWithinEnclosure
    print("- Updated phi (material id: %d) based on surface mesh(es). Going to reinitialize phi.\n", materialid);
    Reinitialize(0.0, 1.0, 0.0, Phi, 600, true/*must do*/); //the first 3 inputs are irrelevant because of "must do"
  }

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

void LevelSetOperator::ComputeResidual(SpaceVariable3D &V, SpaceVariable3D &Phi, SpaceVariable3D &R,
                                       [[maybe_unused]] double time, [[maybe_unused]] double dt)
{

#ifdef LEVELSET_TEST
  PrescribeVelocityFieldForTesting(V, Phi, time, dt);
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

void LevelSetOperator::PrescribeVelocityFieldForTesting([[maybe_unused]] SpaceVariable3D &V, [[maybe_unused]] SpaceVariable3D &Phi,
                                                        [[maybe_unused]] double time, [[maybe_unused]] double dt)
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
        v[k][j][i][1] =  2.0*pow(sin(pi*x),2)*sin(pi*y)*cos(pi*y)*cos(pi*(time-0.5*dt)/T);
        v[k][j][i][2] = -2.0*sin(pi*x)*cos(pi*x)*pow(sin(pi*y),2)*cos(pi*(time-0.5*dt)/T);
        v[k][j][i][3] = 0.0; 

      }
  coordinates.RestoreDataPointerToLocalVector();
  V.RestoreDataPointerAndInsert();

#elif LEVELSET_TEST == 3

  // merging and separation of two circles
  double aa = 1.75, rr = 1.5;
  double s = (time-0.5*dt)<1.0 ? 1.0 : -1.0;
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


