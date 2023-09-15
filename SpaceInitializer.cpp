/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <SpaceInitializer.h>
#include <SpaceOperator.h>
#include <LevelSetOperator.h>
#include <EmbeddedBoundaryDataSet.h>
#include <EmbeddedBoundaryOperator.h>
#include <GravityHandler.h>
#include <DistancePointToParallelepiped.h>
#include <DistancePointToSpheroid.h>
#include <DistancePointToSphere.h>
#include <DistancePointToCylinderCone.h>
#include <DistancePointToCylinderSphere.h>
#include <DistancePointToPlane.h>
#include <KDTree.h>
#include <rbf_interp.hpp>

#ifdef LEVELSET_TEST
  #include <Vector2D.h>
#endif

using std::vector;
using std::multimap;
using std::set;
using std::pair;
using std::make_pair;
using std::unique_ptr;

extern double domain_diagonal;
extern int verbose;
extern int INACTIVE_MATERIAL_ID;
//---------------------------------------------------------------------

SpaceInitializer::SpaceInitializer(MPI_Comm &comm_, DataManagers3D &dms_, IoData &iod_,
                                   GlobalMeshInfo &global_mesh_, SpaceVariable3D &coordinates)
                : comm(comm_), dms(dms_), iod(iod_), global_mesh(global_mesh_)
{
  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);
  coordinates.GetGlobalSize(&NX, &NY, &NZ);
}

//---------------------------------------------------------------------

SpaceInitializer::~SpaceInitializer()
{ }

//---------------------------------------------------------------------

void
SpaceInitializer::Destroy()
{ }

//---------------------------------------------------------------------

multimap<int, pair<int,int> >
SpaceInitializer::SetInitialCondition(SpaceVariable3D &V, SpaceVariable3D &ID, vector<SpaceVariable3D*> &Phi,
                                      SpaceOperator &spo, vector<LevelSetOperator*> &lso,
                                      GravityHandler* ghand,
                                      unique_ptr<vector<unique_ptr<EmbeddedBoundaryDataSet> > > EBDS)
{

  // Step 1: Setup the operation order for user-specified geometries.
  //-------------------------------------
  // Internal Geometry ID:
  // 0: Plane: 0,
  // 1: CylinderCone
  // 2: CylinderSphere
  // 3: Sphere
  // 4: Parallelepiped
  // 5: Spheroid
  // 6: Custom-Geometry ("enclosure")
  // 7: Point
  //-------------------------------------
  vector<pair<int,int> > order;  //<geom type, geom dataMap index>
  [[maybe_unused]] int nGeom = OrderUserSpecifiedGeometries(order);

  // Step 2: Apply initial conditions to V, ID (B.C. applied to V, but not ID!)
  multimap<int, pair<int,int> > id2closure = InitializeVandID(V, ID, spo, EBDS.get(), order);

  // Step 3:: Apply initial conditions to Phi
  InitializePhi(ID, spo, Phi, lso, EBDS.get(), order, id2closure);

  // Step 4:: Apply flooding if specified by user
  if(ghand)
    ghand->UpdateInitialConditionByFlooding(V, ID, lso, Phi, EBDS.get());

  return id2closure;
}

//---------------------------------------------------------------------

int
SpaceInitializer::OrderUserSpecifiedGeometries(vector<std::pair<int,int> > &order)
{
  //-------------------------------------
  // Internal Geometry ID:
  // 0: Plane: 0,
  // 1: CylinderCone
  // 2: CylinderSphere
  // 3: Sphere
  // 4: Parallelepiped
  // 5: Spheroid
  // 6: Custom-Geometry ("enclosure")
  // 7: Point
  //-------------------------------------

  MultiInitialConditionsData &ic(iod.ic.multiInitialConditions);
  int nGeom = ic.planeMap.dataMap.size() //0
            + ic.cylinderconeMap.dataMap.size() //1
            + ic.cylindersphereMap.dataMap.size() //2
            + ic.sphereMap.dataMap.size() //3
            + ic.parallelepipedMap.dataMap.size() //4
            + ic.spheroidMap.dataMap.size() //5
            + ic.enclosureMap.dataMap.size() //6
            + ic.pointMap.dataMap.size(); //7
  order.assign(nGeom, make_pair(-1,-1));
  vector<int> user_specified_order(nGeom, -1); //for verification of user's input

  for(auto&& obj : ic.planeMap.dataMap)
    AddGeomToVector(obj.second->order, 0, obj.first, "Plane", order, user_specified_order);
  for(auto&& obj : ic.cylinderconeMap.dataMap)
    AddGeomToVector(obj.second->order, 1, obj.first, "CylinderCone", order, user_specified_order);
  for(auto&& obj : ic.cylindersphereMap.dataMap)
    AddGeomToVector(obj.second->order, 2, obj.first, "CylinderSphere", order, user_specified_order);
  for(auto&& obj : ic.sphereMap.dataMap)
    AddGeomToVector(obj.second->order, 3, obj.first, "Sphere", order, user_specified_order);
  for(auto&& obj : ic.parallelepipedMap.dataMap)
    AddGeomToVector(obj.second->order, 4, obj.first, "Parallelepiped", order, user_specified_order);
  for(auto&& obj : ic.spheroidMap.dataMap)
    AddGeomToVector(obj.second->order, 5, obj.first, "Spheroid", order, user_specified_order);
  for(auto&& obj : ic.enclosureMap.dataMap)
    AddGeomToVector(obj.second->order, 6, obj.first, "ArbitraryEnclosure", order, user_specified_order);
  for(auto&& obj : ic.pointMap.dataMap)
    AddGeomToVector(obj.second->order, 7, obj.first, "Point", order, user_specified_order);

  //verification
  for(int i=0; i<(int)user_specified_order.size()-1; i++)
    assert(user_specified_order[i]<=user_specified_order[i+1]);

  return nGeom;
}

//---------------------------------------------------------------------

void
SpaceInitializer::AddGeomToVector(int o, int type, int ind, string name, vector<std::pair<int,int> > &order,
                                  vector<int> &user_specified_order)
{
  int nGeom = order.size();
  if(o<0 || o>=nGeom) {
    print_error("*** Error: Detected incorrect %s order (%d). Range: [0, %d)\n",
                name.c_str(), o, nGeom);
    exit_mpi();
  }

  int i;
  for(i=o; i<nGeom; i++) {
    if(order[i].first == -1) { //open spot, take it.
      order[i].first = type;
      order[i].second = ind;
      user_specified_order[i] = o;
      return;
    } 
    else { //make sure there is not a conflict
      if(user_specified_order[i]>o) {
        print_error("*** Error: Detected conflicting order in user-specified geometries.\n");
        exit_mpi();
      }
    }
  }

  if(i==nGeom) {
    print_error("*** Error: Unable to specify a valid order for user-specified geometries. "
                "(Rank/Order starts at 0).\n");
    exit_mpi();
  }
}

//---------------------------------------------------------------------

multimap<int, pair<int,int> >
SpaceInitializer::InitializeVandID(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceOperator &spo,
                                   vector<unique_ptr<EmbeddedBoundaryDataSet> >* EBDS,
                                   vector<pair<int,int> > &order)
{

  multimap<int, pair<int,int> > id2closure; //only relevant when there are embedded boundaries

  //----------------------------------------------------------------
  // Step 1: Get the needed data
  //----------------------------------------------------------------
  SpaceVariable3D &coordinates(spo.GetMeshCoordinates());

  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
  Vec5D*** v = (Vec5D***) V.GetDataPointer();
  double*** id = (double***) ID.GetDataPointer();

  print("\n- Initializing the state variables (V) and material id. (ID).\n");

  //----------------------------------------------------------------
  // Step 2: Apply the default initial condition
  //----------------------------------------------------------------
  //! apply the default initial condition (usually, the farfield/inlet b.c.)
  if(iod.ic.default_ic.materialid != 0)
    print_warning("Warning: Material ID of the default initial state is %d. In most cases, it should be 0.\n",
                  iod.ic.default_ic.materialid);

  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++) {
        v[k][j][i][0] = iod.ic.default_ic.density;
        v[k][j][i][1] = iod.ic.default_ic.velocity_x;
        v[k][j][i][2] = iod.ic.default_ic.velocity_y;
        v[k][j][i][3] = iod.ic.default_ic.velocity_z;
        v[k][j][i][4] = iod.ic.default_ic.pressure;
        id[k][j][i]   = iod.ic.default_ic.materialid;
      }

  
  //----------------------------------------------------------------
  // Step 2: Apply initial condition from data file (If specified here)
  //----------------------------------------------------------------
  if(iod.ic.apply_user_file_before_geometries==IcData::YES)
    ApplyUserSpecifiedInitialConditionFile(coords, v, id);


  //----------------------------------------------------------------
  // Step 3: Apply dummy i.c. within embedded closed surfaces
  //----------------------------------------------------------------
  vector<double***> color;
  if(EBDS != nullptr) {
    color.assign(EBDS->size(), NULL);
    for(int i=0; i<(int)EBDS->size(); i++) {
      assert((*EBDS)[i]->Color_ptr);
      color[i] = (*EBDS)[i]->Color_ptr->GetDataPointer();
    }

    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++) {

          for(int surf=0; surf<(int)color.size(); surf++) {
            if(color[surf][k][j][i] <= 0) {//occluded or enclosed
              v[k][j][i][0] = iod.eqs.dummy_state.density;
              v[k][j][i][1] = iod.eqs.dummy_state.velocity_x;
              v[k][j][i][2] = iod.eqs.dummy_state.velocity_y;
              v[k][j][i][3] = iod.eqs.dummy_state.velocity_z;
              v[k][j][i][4] = iod.eqs.dummy_state.pressure;
              id[k][j][i] = INACTIVE_MATERIAL_ID;
              break;
            }
          }
        }
  }  


  //----------------------------------------------------------------
  // Step 4: Apply initial condition based on geometric entities
  //----------------------------------------------------------------
  MultiInitialConditionsData &ic(iod.ic.multiInitialConditions);

  if(!ic.pointMap.dataMap.empty() && EBDS == nullptr) {
    print_error("*** Error: Unable to specify point-based initial conditions w/o embedded boundaries.\n"
                "           Background state can be defined using the 'Farfield' or 'Inlet' structure,\n"
                "           even if the domain does not contain a far-field or inlet boundary.\n");
    exit_mpi();
  }

  for(auto&& obj : order) {
    //-------------------------------------
    // Internal Geometry ID:
    // 0: Plane: 0,
    // 1: CylinderCone
    // 2: CylinderSphere
    // 3: Sphere
    // 4: Parallelepiped
    // 5: Spheroid
    // 6: Custom-Geometry ("enclosure")
    // 7: Point
    //-------------------------------------
    if(obj.first == 0) {//plane
      auto it = ic.planeMap.dataMap.find(obj.second);
      assert(it != ic.planeMap.dataMap.end());

      if(it->second->inclusion != PlaneData::OVERRIDE) {
        print_error("*** Error: Initialization only supports Inclusion = Override at the moment.\n");
        exit_mpi();
      }

      Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
      Vec3D dir(it->second->nx, it->second->ny, it->second->nz);
      GeoTools::DistanceFromPointToPlane distCal(x0, dir); //dir will be normalized by the constructor

      print("  o Found Plane[%d]: (%e %e %e) | normal: (%e %e %e). material id: %d.\n",
            it->first, x0[0], x0[1], x0[2], dir[0], dir[1], dir[2],
            it->second->initialConditions.materialid);

      for(int k=k0; k<kmax; k++)
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {
            if(distCal.Calculate(coords[k][j][i])>0) {
              v[k][j][i][0] = it->second->initialConditions.density;
              v[k][j][i][1] = it->second->initialConditions.velocity_x;
              v[k][j][i][2] = it->second->initialConditions.velocity_y;
              v[k][j][i][3] = it->second->initialConditions.velocity_z;
              v[k][j][i][4] = it->second->initialConditions.pressure;
              id[k][j][i]   = it->second->initialConditions.materialid;
            }
          }
    }
    else if(obj.first == 1) {//cylinder-cone
      auto it = ic.cylinderconeMap.dataMap.find(obj.second);
      assert(it != ic.cylinderconeMap.dataMap.end());

      bool interior = (it->second->side == CylinderConeData::INTERIOR);

      if(it->second->inclusion != CylinderConeData::OVERRIDE) {
        print_error("*** Error: Initialization only supports Inclusion = Override at the moment.\n");
        exit_mpi();
      }

      Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
      Vec3D dir(it->second->nx, it->second->ny, it->second->nz);
      double L = it->second->L; //cylinder height
      double R = it->second->r; //cylinder radius
      double tan_alpha = tan(it->second->opening_angle_degrees/180.0*acos(-1.0));//opening angle
      double Hmax = R/tan_alpha;
      double H = std::min(it->second->cone_height, Hmax); //cone's height
      GeoTools::DistanceFromPointToCylinderCone distCal(x0,dir,L,R,tan_alpha,H);

      print("  o Found Cylinder-Cone[%d]: Base center: (%e %e %e), Axis: (%e %e %e), Mat. ID: %d.\n",
            it->first, x0[0], x0[1], x0[2], dir[0], dir[1], dir[2],
            it->second->initialConditions.materialid);

      for(int k=k0; k<kmax; k++)
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {
            if((interior  && distCal.Calculate(coords[k][j][i])<0) ||
               (!interior && distCal.Calculate(coords[k][j][i])>0)) {
              v[k][j][i][0] = it->second->initialConditions.density;
              v[k][j][i][1] = it->second->initialConditions.velocity_x;
              v[k][j][i][2] = it->second->initialConditions.velocity_y;
              v[k][j][i][3] = it->second->initialConditions.velocity_z;
              v[k][j][i][4] = it->second->initialConditions.pressure;
              id[k][j][i]   = it->second->initialConditions.materialid;
            }
          }
    }
    else if(obj.first == 2) { //cylinder-sphere (i.e., possibly with spherical cap(s))
      auto it = ic.cylindersphereMap.dataMap.find(obj.second);
      assert(it != ic.cylindersphereMap.dataMap.end());

      bool interior = (it->second->side == CylinderSphereData::INTERIOR);

      if(it->second->inclusion != CylinderSphereData::OVERRIDE) {
        print_error("*** Error: Initialization only supports Inclusion = Override at the moment.\n");
        exit_mpi();
      }

      Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z); //center of base
      Vec3D dir(it->second->nx, it->second->ny, it->second->nz);
      double L = it->second->L; //cylinder height
      double R = it->second->r; //cylinder radius
      bool front_cap = (it->second->front_cap == CylinderSphereData::On);
      bool back_cap = (it->second->back_cap == CylinderSphereData::On);
      GeoTools::DistanceFromPointToCylinderSphere distCal(x0, dir, L, R, front_cap, back_cap);

      print("  o Found Cylinder-Sphere[%d]: Base center: (%e %e %e), Axis: (%e %e %e). Mat. ID: \n",
            it->first, x0[0], x0[1], x0[2], dir[0], dir[1], dir[2],
            it->second->initialConditions.materialid);

      for(int k=k0; k<kmax; k++)
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {
            if((interior  && distCal.Calculate(coords[k][j][i])<0) ||
               (!interior && distCal.Calculate(coords[k][j][i])>0)) {
              v[k][j][i][0] = it->second->initialConditions.density;
              v[k][j][i][1] = it->second->initialConditions.velocity_x;
              v[k][j][i][2] = it->second->initialConditions.velocity_y;
              v[k][j][i][3] = it->second->initialConditions.velocity_z;
              v[k][j][i][4] = it->second->initialConditions.pressure;
              id[k][j][i]   = it->second->initialConditions.materialid;
            }
          }
    }
    else if(obj.first == 3) { //sphere
      auto it = ic.sphereMap.dataMap.find(obj.second);
      assert(it != ic.sphereMap.dataMap.end());

      bool interior = (it->second->side == SphereData::INTERIOR);

      if(it->second->inclusion != SphereData::OVERRIDE) {
        print_error("*** Error: Initialization only supports Inclusion = Override at the moment.\n");
        exit_mpi();
      }

      Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
      double R = it->second->radius;
      GeoTools::DistanceFromPointToSphere distCal(x0, R);

      print("  o Found Sphere[%d]: Center: (%e %e %e), Radius: %e. Mat. ID: %d\n",
            it->first, x0[0], x0[1], x0[2], R, it->second->initialConditions.materialid);

      for(int k=k0; k<kmax; k++)
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {
            if((interior  && distCal.Calculate(coords[k][j][i])<0) ||
               (!interior && distCal.Calculate(coords[k][j][i])>0)) {
              v[k][j][i][0] = it->second->initialConditions.density;
              v[k][j][i][1] = it->second->initialConditions.velocity_x;
              v[k][j][i][2] = it->second->initialConditions.velocity_y;
              v[k][j][i][3] = it->second->initialConditions.velocity_z;
              v[k][j][i][4] = it->second->initialConditions.pressure;
              id[k][j][i]   = it->second->initialConditions.materialid;
            }
          }
    }
    else if(obj.first == 4) { //parallelepiped
      auto it = ic.parallelepipedMap.dataMap.find(obj.second);
      assert(it != ic.parallelepipedMap.dataMap.end());

      bool interior = (it->second->side == ParallelepipedData::INTERIOR);

      if(it->second->inclusion != ParallelepipedData::OVERRIDE) {
        print_error("*** Error: Initialization only supports Inclusion = Override at the moment.\n");
        exit_mpi();
      }

      Vec3D x0(it->second->x0, it->second->y0, it->second->z0);
      Vec3D oa(it->second->ax, it->second->ay, it->second->az);  oa -= x0;
      Vec3D ob(it->second->bx, it->second->by, it->second->bz);  ob -= x0;
      Vec3D oc(it->second->cx, it->second->cy, it->second->cz);  oc -= x0;

      if(oa.norm()==0 || ob.norm()==0 || oc.norm()==0 || (oa^ob)*oc<=0.0) {
        print_error("*** Error: Detected error in a user-specified parallelepiped. "
                    "Overlapping vertices or violation of right-hand rule.\n");
        exit_mpi();
      }

      GeoTools::DistanceFromPointToParallelepiped distCal(x0, oa, ob, oc);
      print("  o Found Parallelepiped[%d]: X(%e %e %e). Mat. ID: %d\n", it->first,
            x0[0], x0[1], x0[2], it->second->initialConditions.materialid);

      for(int k=k0; k<kmax; k++)
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {
            if((interior  && distCal.Calculate(coords[k][j][i])<0) ||
               (!interior && distCal.Calculate(coords[k][j][i])>0)) {
              v[k][j][i][0] = it->second->initialConditions.density;
              v[k][j][i][1] = it->second->initialConditions.velocity_x;
              v[k][j][i][2] = it->second->initialConditions.velocity_y;
              v[k][j][i][3] = it->second->initialConditions.velocity_z;
              v[k][j][i][4] = it->second->initialConditions.pressure;
              id[k][j][i]   = it->second->initialConditions.materialid;
            }
          }
    }
    else if(obj.first == 5) { //spheroid
      auto it = ic.spheroidMap.dataMap.find(obj.second);
      assert(it != ic.spheroidMap.dataMap.end());

      bool interior = (it->second->side == SpheroidData::INTERIOR);

      if(it->second->inclusion != SpheroidData::OVERRIDE) {
        print_error("*** Error: Initialization only supports Inclusion = Override at the moment.\n");
        exit_mpi();
      }

      Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
      Vec3D axis(it->second->axis_x, it->second->axis_y, it->second->axis_z);
      GeoTools::DistanceFromPointToSpheroid distCal(x0, axis, it->second->semi_length, it->second->radius);

      print("  o Found Spheroid[%d]: Center: (%e %e %e), Axis: (%e %e %e). Mat. ID: %d\n", it->first,
            x0[0], x0[1], x0[2], axis[0], axis[1], axis[2], it->second->initialConditions.materialid);

      for(int k=k0; k<kmax; k++)
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {
            if((interior  && distCal.Calculate(coords[k][j][i])<0) ||
               (!interior && distCal.Calculate(coords[k][j][i])>0)) {
              v[k][j][i][0] = it->second->initialConditions.density;
              v[k][j][i][1] = it->second->initialConditions.velocity_x;
              v[k][j][i][2] = it->second->initialConditions.velocity_y;
              v[k][j][i][3] = it->second->initialConditions.velocity_z;
              v[k][j][i][4] = it->second->initialConditions.pressure;
              id[k][j][i]   = it->second->initialConditions.materialid;
            }
          }
    }
    else if(obj.first == 6) { //user-specified enclosure
      auto it = ic.enclosureMap.dataMap.find(obj.second);
      assert(it != ic.enclosureMap.dataMap.end());

      if(it->second->inclusion != UserSpecifiedEnclosureData::OVERRIDE) {
        print_error("*** Error: Initialization only supports Inclusion = Override at the moment.\n");
        exit_mpi();
      }

      coordinates.RestoreDataPointerToLocalVector(); //because the next function calls some emb fun.
      InitializeVandIDWithinEnclosure(*it->second, coordinates, *(spo.GetPointerToInnerGhostNodes()),
                                      *(spo.GetPointerToOuterGhostNodes()), v, id);
      coords = (Vec3D***)coordinates.GetDataPointer();
    }
    else if(obj.first == 7) { //point (must be inside an embedded closed surface)
      auto it = ic.pointMap.dataMap.find(obj.second);
      assert(it != ic.pointMap.dataMap.end());

      print("  o Found Point[%d]: (%e %e %e). Mat. ID: %d\n",
            it->first, it->second->x, it->second->y, it->second->z,
            it->second->initialConditions.materialid);

      if(it->second->inclusion != PointData::OVERRIDE) {
        print_error("*** Error: Initialization only supports Inclusion = Override at the moment.\n");
        exit_mpi();
      }

      id2closure.insert(InitializeVandIDByPoint(*it->second, *EBDS, color, v, id));
    }
    else {
      print_error("*** Error: Found unrecognized object id (%d).\n", obj.first);
      exit_mpi();
    }
  }


  // Verification (can be deleted): occluded nodes should have inactive_material_id
  if(EBDS != nullptr) {
    int i,j,k;
    for(int surf=0; surf<(int)EBDS->size(); surf++) {
      set<Int3> *occluded = (*EBDS)[surf]->occluded_ptr;
      set<Int3> *imposed_occluded = (*EBDS)[surf]->imposed_occluded_ptr;
      for(auto it = occluded->begin(); it != occluded->end(); it++) {
        i = (*it)[0];
        j = (*it)[1];
        k = (*it)[2];
        if(!coordinates.IsHere(i,j,k,false))
          continue;
        assert(id[k][j][i] == INACTIVE_MATERIAL_ID);
      }
      for(auto it = imposed_occluded->begin(); it != imposed_occluded->end(); it++) {
        i = (*it)[0];
        j = (*it)[1];
        k = (*it)[2];
        if(!coordinates.IsHere(i,j,k,false))
          continue;
        assert(id[k][j][i] == INACTIVE_MATERIAL_ID);
      }
    }
  }


  //----------------------------------------------------------------
  // Step 5: Apply initial condition from data file (If specified here)
  //----------------------------------------------------------------
  if(iod.ic.apply_user_file_before_geometries==IcData::NO)
    ApplyUserSpecifiedInitialConditionFile(coords, v, id);


  //----------------------------------------------------------------
  // Step 6: Restore Data
  //----------------------------------------------------------------
  if(EBDS != nullptr) {
    for(int i=0; i<(int)EBDS->size(); i++)
      (*EBDS)[i]->Color_ptr->RestoreDataPointerToLocalVector();
  }
  V.RestoreDataPointerAndInsert();
  ID.RestoreDataPointerAndInsert();
  coordinates.RestoreDataPointerToLocalVector(); //!< data was not changed.

  //----------------------------------------------------------------
  // Step 7: Apply boundary conditions to V (ID will be handled in Main by mpo)
  //----------------------------------------------------------------
  spo.ClipDensityAndPressure(V, ID);
  spo.ApplyBoundaryConditions(V);

  return id2closure;

}

//---------------------------------------------------------------------

void
SpaceInitializer::ApplyUserSpecifiedInitialConditionFile(Vec3D*** coords, Vec5D*** v, double*** id)
{

  if(iod.ic.type != IcData::NONE) {

    //! Get coordinates
    Vec3D    x0(iod.ic.x0[0], iod.ic.x0[1], iod.ic.x0[2]); 

    if (iod.ic.type == IcData::PLANAR || iod.ic.type == IcData::CYLINDRICAL) {

      if(iod.ic.type == IcData::PLANAR)
        print("  o Applying the initial condition specified in %s (planar).\n\n", 
              iod.ic.user_specified_ic);
      else
        print("  o Applying the initial condition specified in %s (with cylindrical symmetry).\n\n", 
              iod.ic.user_specified_ic);
 
      Vec3D dir(iod.ic.dir[0], iod.ic.dir[1], iod.ic.dir[2]); 
      dir /= dir.norm();

      double x = 0.0;
      int n = iod.ic.user_data[IcData::COORDINATE].size(); //!< number of data points provided by user (in the axial dir)

      double r = 0.0;
      int nrad = iod.ic.user_data2[IcData::COORDINATE].size(); //!< number of data points in the radial dir.

      int t0, t1;    
      double a0, a1;

      for(int k=k0; k<kmax; k++)
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {

            x = (coords[k][j][i] - x0)*dir; //!< projection onto the 1D axis
//            cout << "coords: " << coords[j][i][0] << ", " << coords[j][i][1] << "; x = " << x << endl;
            if(x<0 || x>iod.ic.user_data[IcData::COORDINATE][n-1])
              continue;

            if(nrad>0) {
              r = (coords[k][j][i] - x0 - x*dir).norm();
              if(r>iod.ic.user_data2[IcData::COORDINATE][nrad-1])
                continue;
            }
 
            //! Find the first 1D coordinate greater than x
            auto upper_it = std::upper_bound(iod.ic.user_data[IcData::COORDINATE].begin(),
                                             iod.ic.user_data[IcData::COORDINATE].end(),
                                             x); 
            t1 = (int)(upper_it - iod.ic.user_data[IcData::COORDINATE].begin());

            if(t1==0) // exactly the first node in 1D
              t1 = 1;

            t0 = t1 - 1;

            //! calculate interpolation weights
            a0 = (iod.ic.user_data[IcData::COORDINATE][t1] - x) /
                 (iod.ic.user_data[IcData::COORDINATE][t1] - iod.ic.user_data[IcData::COORDINATE][t0]);
            a1 = 1.0 - a0;

            //! specify i.c. on node (cell center)
            if(iod.ic.specified[IcData::DENSITY])
              v[k][j][i][0] = a0*iod.ic.user_data[IcData::DENSITY][t0]  
                            + a1*iod.ic.user_data[IcData::DENSITY][t1];
            if(iod.ic.specified[IcData::VELOCITY]) {
              v[k][j][i][1] = (a0*iod.ic.user_data[IcData::VELOCITY][t0] 
                            +  a1*iod.ic.user_data[IcData::VELOCITY][t1])*dir[0];
              v[k][j][i][2] = (a0*iod.ic.user_data[IcData::VELOCITY][t0]
                            +  a1*iod.ic.user_data[IcData::VELOCITY][t1])*dir[1];
              v[k][j][i][3] = (a0*iod.ic.user_data[IcData::VELOCITY][t0] 
                            +  a1*iod.ic.user_data[IcData::VELOCITY][t1])*dir[2];
            }
            if(iod.ic.specified[IcData::PRESSURE]) 
              v[k][j][i][4] = a0*iod.ic.user_data[IcData::PRESSURE][t0]
                            + a1*iod.ic.user_data[IcData::PRESSURE][t1];
            if(iod.ic.specified[IcData::MATERIALID])
              id[k][j][i]   = std::round(a0*iod.ic.user_data[IcData::MATERIALID][t0] 
                                       + a1*iod.ic.user_data[IcData::MATERIALID][t1]);

            //! apply radial variation (if provided by user)
            if(nrad>0) { 
 
              //! Find the first radial coordinate greater than r
              auto upper_it = std::upper_bound(iod.ic.user_data2[IcData::COORDINATE].begin(),
                                               iod.ic.user_data2[IcData::COORDINATE].end(),
                                               r); 
              t1 = (int)(upper_it - iod.ic.user_data2[IcData::COORDINATE].begin());

              if(t1==0) // exactly the first node
                t1 = 1;

              t0 = t1 - 1;

              //! calculate interpolation weights
              a0 = (iod.ic.user_data2[IcData::COORDINATE][t1] - r) /
                   (iod.ic.user_data2[IcData::COORDINATE][t1] - iod.ic.user_data2[IcData::COORDINATE][t0]);
              a1 = 1.0 - a0;

              if(iod.ic.specified[IcData::DENSITY])
                v[k][j][i][0] *= a0*iod.ic.user_data2[IcData::DENSITY][t0] 
                               + a1*iod.ic.user_data2[IcData::DENSITY][t1];
              if(iod.ic.specified[IcData::VELOCITY])
                for(int p=1; p<=3; p++)
                  v[k][j][i][p] *= a0*iod.ic.user_data2[IcData::VELOCITY][t0] 
                                 + a1*iod.ic.user_data2[IcData::VELOCITY][t1];
              if(iod.ic.specified[IcData::PRESSURE])
                v[k][j][i][4] *= a0*iod.ic.user_data2[IcData::PRESSURE][t0]
                               + a1*iod.ic.user_data2[IcData::PRESSURE][t1];
            }
          }

    }

    else if (iod.ic.type == IcData::GENERALCYLINDRICAL) {

      Vec3D dir(iod.ic.dir[0], iod.ic.dir[1], iod.ic.dir[2]); 
      dir /= dir.norm();

      int N = iod.ic.user_data[IcData::COORDINATE].size(); //!< number of data points provided by user

      print("  o Applying the initial condition specified in %s (%d points, cylindrical symmetry).\n\n", 
            iod.ic.user_specified_ic, N);

      //Store sample points in a KDTree (K = 2)
      Vec2D *xy = new Vec2D[N]; 
      PointIn2D *p = new PointIn2D[N];
      for(int i=0; i<N; i++) {
        xy[i] = Vec2D(iod.ic.user_data[IcData::COORDINATE][i], iod.ic.user_data[IcData::RADIALCOORDINATE][i]);
        p[i]  = PointIn2D(i, xy[i]);
      }
      KDTree<PointIn2D,2/*dim*/> tree(N, p); //"tree" uses "p" and changes its ordering. "p" cannot be deleted before "tree"
      print("    * Constructed a KDTree (K=2) to store the data points.\n");

      int numPoints = 15; //this is the number of points for interpolation
      int maxCand = numPoints*20;
      PointIn2D candidates[maxCand*10]; //maxCand may be temporarily increased for fail-safe


      // tentative cutoff distance --- will be adjusted
      double cutoff = sqrt((iod.ic.xmax[0]-iod.ic.xmin[0])*(iod.ic.xmax[1]-iod.ic.xmin[1])/N);

      Vec2D pnode;
      for(int k=k0; k<kmax; k++)
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {

            //3D -> 2D
            pnode[0] = (coords[k][j][i] - x0)*dir; //!< projection onto the 1D axis
            if(pnode[0]<iod.ic.xmin[0] || pnode[0]>iod.ic.xmax[0])
              continue;
            pnode[1] = (coords[k][j][i] - x0 - pnode[0]*dir).norm();
            if(pnode[1]<iod.ic.xmin[1] || pnode[1]>iod.ic.xmax[1])
              continue;
 
            //RBF interpolation w/ KDTree
            int nFound = 0, counter = 0;
            double low_cut = 0.0, high_cut = DBL_MAX;
            while(nFound<numPoints || nFound>maxCand) {
              if(++counter>2000) {
                fprintf(stdout,"*** Error: Cannot find required number of sample points for "
                               "interpolation (by RBF) after %d iterations. "
                               "Coord(3D):%e %e %e  | Coord(2D):%e %e. "
                               "Candidates: %d, cutoff = %e.\n", 
                               counter, coords[k][j][i][0], coords[k][j][i][1], coords[k][j][i][2],
                               pnode[0], pnode[1], nFound, cutoff);
                exit(-1);
              }
              nFound = tree.findCandidatesWithin(pnode, candidates, maxCand, cutoff);

              if(nFound<numPoints) {
                low_cut = std::max(low_cut, cutoff);
                if(high_cut>0.5*DBL_MAX)
                  cutoff *= 4.0;
                else
                  cutoff = 0.5*(low_cut + high_cut);
              }
              else if(nFound>maxCand) {
                high_cut = std::min(high_cut, cutoff);
                cutoff = 0.5*(low_cut + high_cut);
              }

              if((high_cut - low_cut)/high_cut<1e-6) { //fail-safe
                nFound = tree.findCandidatesWithin(pnode, candidates, 10*maxCand, high_cut);
                if(nFound>10*maxCand) {
                  fprintf(stdout,"\033[0;31m*** Error: Cannot find required number of sample points at any"
                                 " cutoff distance.\033[0m\n");
                  exit(-1);
                }

                assert(nFound>=numPoints);
                fprintf(stdout,"\033[0;35mWarning: Unusual behavior. Found %d candidates with cutoff = %e "
                               "(node: %e %e).\033[0m\n", nFound, high_cut, pnode[0], pnode[1]);
                break;
              }
            }

            //figure out the actual points for interpolation (numPoints)
            vector<pair<double,int> > dist2node;
            for(int i=0; i<nFound; i++) 
              dist2node.push_back(make_pair((candidates[i].x-pnode).norm(), candidates[i].id));

            if(nFound>numPoints)
              sort(dist2node.begin(), dist2node.end());         
            dist2node.resize(numPoints);

            //prepare to interpolate
            double xd[2*numPoints];
            for(int i=0; i<numPoints; i++) {
              xd[2*i]   = xy[dist2node[i].second][0];
              xd[2*i+1] = xy[dist2node[i].second][1];
            }
            double r0; //slightly larger than maximum separation
            r0 = dist2node.front().first + dist2node.back().first;
            double fd[numPoints];
            double *rbf_weight = new double[numPoints];
	    double *interp = new double[numPoints];

            //choose a radial basis function
            void (*phi)(int, double[], double, double[]); //a function pointer
            switch (iod.ic.rbf) {
              case IcData::MULTIQUADRIC :
                phi = MathTools::phi1;  break;
              case IcData::INVERSE_MULTIQUADRIC :
                phi = MathTools::phi2;  break;
              case IcData::THIN_PLATE_SPLINE :
                phi = MathTools::phi3;  break;
              case IcData::GAUSSIAN :
                phi = MathTools::phi4;  break;
              default : 
                phi = MathTools::phi1;  break;
            }

            if(iod.ic.specified[IcData::DENSITY]) {
              for(int i=0; i<numPoints; i++)
                fd[i] = iod.ic.user_data[IcData::DENSITY][dist2node[i].second];
              MathTools::rbf_weight(2, numPoints, xd, r0, phi, fd, rbf_weight);
              MathTools::rbf_interp(2, numPoints, xd, r0,
                                    phi, rbf_weight,
                                    1, pnode, interp);
              v[k][j][i][0] = interp[0];
            }

            if(iod.ic.specified[IcData::VELOCITY] || iod.ic.specified[IcData::RADIALVELOCITY])
              v[k][j][i][1] = v[k][j][i][2] = v[k][j][i][3] = 0.0;

            if(iod.ic.specified[IcData::VELOCITY]) {
              for(int i=0; i<numPoints; i++)
                fd[i] = iod.ic.user_data[IcData::VELOCITY][dist2node[i].second];
              MathTools::rbf_weight(2, numPoints, xd, r0, phi, fd, rbf_weight);
              MathTools::rbf_interp(2, numPoints, xd, r0,
                                    phi, rbf_weight,
                                    1, pnode, interp);
              v[k][j][i][1] = interp[0]*dir[0]; 
              v[k][j][i][2] = interp[0]*dir[1]; 
              v[k][j][i][3] = interp[0]*dir[2];
            }

            if(iod.ic.specified[IcData::RADIALVELOCITY]) {
              for(int i=0; i<numPoints; i++)
                fd[i] = iod.ic.user_data[IcData::RADIALVELOCITY][dist2node[i].second];
              MathTools::rbf_weight(2, numPoints, xd, r0, phi, fd, rbf_weight);
              MathTools::rbf_interp(2, numPoints, xd, r0,
                                    phi, rbf_weight,
                                    1, pnode, interp);
              Vec3D dir2 = coords[k][j][i] - x0 - pnode[0]*dir;
              if(dir2.norm()>0)
                dir2 /= dir2.norm();
              v[k][j][i][1] += interp[0]*dir2[0]; 
              v[k][j][i][2] += interp[0]*dir2[1]; 
              v[k][j][i][3] += interp[0]*dir2[2];
            }

            if(iod.ic.specified[IcData::PRESSURE]) {
              for(int i=0; i<numPoints; i++)
                fd[i] = iod.ic.user_data[IcData::PRESSURE][dist2node[i].second];
              MathTools::rbf_weight(2, numPoints, xd, r0, phi, fd, rbf_weight);
              MathTools::rbf_interp(2, numPoints, xd, r0,
                                    phi, rbf_weight,
                                    1, pnode, interp);
              v[k][j][i][4] = interp[0];
            }

            if(iod.ic.specified[IcData::MATERIALID]) {
              for(int i=0; i<numPoints; i++)
                fd[i] = iod.ic.user_data[IcData::MATERIALID][dist2node[i].second];
              MathTools::rbf_weight(2, numPoints, xd, r0, phi, fd, rbf_weight);
              MathTools::rbf_interp(2, numPoints, xd, r0,
                                    phi, rbf_weight,
                                    1, pnode, interp);
              id[k][j][i] = std::round(interp[0]);
            }

            delete [] rbf_weight;
            delete [] interp;
          }

      delete [] xy;
      delete [] p;

    }

    else if (iod.ic.type == IcData::SPHERICAL) {

      print("  o Applying the initial condition specified in %s (with spherical symmetry).\n\n", 
            iod.ic.user_specified_ic);
 
      double x;
      int n = iod.ic.user_data[IcData::COORDINATE].size(); //!< number of data points provided by user
      Vec3D dir;

      int t0, t1;    
      double a0, a1;

      for(int k=k0; k<kmax; k++)
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {

            dir = coords[k][j][i] - x0;
            x = dir.norm();
            dir /= x;
   
//            cout << "coords: " << coords[j][i][0] << ", " << coords[j][i][1] << "; x = " << x << endl;
            if(x>iod.ic.user_data[IcData::COORDINATE][n-1])
              continue;
 
            //! Find the first 1D coordinate greater than x
            auto upper_it = std::upper_bound(iod.ic.user_data[IcData::COORDINATE].begin(),
                                             iod.ic.user_data[IcData::COORDINATE].end(),
                                             x); 
            t1 = (int)(upper_it - iod.ic.user_data[IcData::COORDINATE].begin());

            if(t1==0) // exactly the first node in 1D
              t1 = 1;

            t0 = t1 - 1;

            //! calculate interpolation weights
            a0 = (iod.ic.user_data[IcData::COORDINATE][t1] - x) /
                 (iod.ic.user_data[IcData::COORDINATE][t1] - iod.ic.user_data[IcData::COORDINATE][t0]);
            a1 = 1.0 - a0;

            //! specify i.c. on node (cell center)
            if(iod.ic.specified[IcData::DENSITY])
              v[k][j][i][0] = a0*iod.ic.user_data[IcData::DENSITY][t0]
                            + a1*iod.ic.user_data[IcData::DENSITY][t1];
            if(iod.ic.specified[IcData::VELOCITY]) {
              v[k][j][i][1] = (a0*iod.ic.user_data[IcData::VELOCITY][t0]
                            +  a1*iod.ic.user_data[IcData::VELOCITY][t1])*dir[0];
              v[k][j][i][2] = (a0*iod.ic.user_data[IcData::VELOCITY][t0] 
                            +  a1*iod.ic.user_data[IcData::VELOCITY][t1])*dir[1];
              v[k][j][i][3] = (a0*iod.ic.user_data[IcData::VELOCITY][t0]
                            +  a1*iod.ic.user_data[IcData::VELOCITY][t1])*dir[2];
            }
            if(iod.ic.specified[IcData::PRESSURE])
              v[k][j][i][4] = a0*iod.ic.user_data[IcData::PRESSURE][t0]
                            + a1*iod.ic.user_data[IcData::PRESSURE][t1];
            if(iod.ic.specified[IcData::MATERIALID])
              id[k][j][i]   = std::round(a0*iod.ic.user_data[IcData::MATERIALID][t0] 
                                       + a1*iod.ic.user_data[IcData::MATERIALID][t1]);
          }

    }
  }

}

//---------------------------------------------------------------------
void
SpaceInitializer::InitializeVandIDWithinEnclosure(UserSpecifiedEnclosureData &enclosure,
                                                  SpaceVariable3D &coordinates,
                                                  vector<GhostPoint> &ghost_nodes_inner,
                                                  vector<GhostPoint> &ghost_nodes_outer, 
                                                  Vec5D*** v, double*** id)
{

  //create an ad hoc EmbeddedSurfaceData
  EmbeddedSurfaceData esd;
  esd.provided_by_another_solver = EmbeddedSurfaceData::NO;
  esd.surface_thickness          = enclosure.surface_thickness;
  esd.filename                   = enclosure.surface_filename;

  print("  o Found User-specified enclosure (%s). Mat. ID: %d\n", esd.filename,
        enclosure.initialConditions.materialid);

  //create the embedded operator to track the surface
  EmbeddedBoundaryOperator embed(comm, esd);
  embed.SetCommAndMeshInfo(dms, coordinates, ghost_nodes_inner, ghost_nodes_outer,
                           global_mesh);
  embed.SetupIntersectors();

  //track the surface, and get access to the results 
  embed.TrackSurfaces();
  auto EBDS = embed.GetPointerToEmbeddedBoundaryData(0);
  double*** color = EBDS->Color_ptr->GetDataPointer();

  //impose i.c. inside enclosures --- EXCLUDING occluded nodes
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        if(color[k][j][i] < 0) {//enclosed
          v[k][j][i][0] = enclosure.initialConditions.density;
          v[k][j][i][1] = enclosure.initialConditions.velocity_x;
          v[k][j][i][2] = enclosure.initialConditions.velocity_y;
          v[k][j][i][3] = enclosure.initialConditions.velocity_z;
          v[k][j][i][4] = enclosure.initialConditions.pressure;
          id[k][j][i]   = enclosure.initialConditions.materialid;
        }
      }

  EBDS->Color_ptr->RestoreDataPointerToLocalVector();

  embed.Destroy();

}

//---------------------------------------------------------------------
//! Apply initial condition based on user-specified point and "flood-fill"
pair<int, pair<int,int> > 
SpaceInitializer::InitializeVandIDByPoint(PointData& point,
                                          vector<unique_ptr<EmbeddedBoundaryDataSet> > &EBDS,
                                          vector<double***> &color,
                                          Vec5D*** v, double*** id)
{
  assert(EBDS.size()>0);

  //Step 1. Locate the point within the mesh
  Vec3D this_point(point.x, point.y, point.z);
  Int3 ijk0;
  bool validpoint = global_mesh.FindElementCoveringPoint(this_point, ijk0, NULL, true);
  if(!validpoint) {
    print_error("*** Error: User-specified point (%e %e %e) is outside the computational domain.\n",
                point.x, point.y, point.z);
    exit_mpi();
  }
  vector<Int3> vertices; 
  vertices.reserve(8);
  vector<bool> owner;
  owner.reserve(8);

  int i,j,k;
  for(int dk=0; dk<=1; dk++)
    for(int dj=0; dj<=1; dj++)
      for(int di=0; di<=1; di++) {
        i = ijk0[0] + di;
        j = ijk0[1] + dj;
        k = ijk0[2] + dk;
        vertices.push_back(Int3(i,j,k));
        owner.push_back(i>=i0 && i<imax && j>=j0 && j<jmax && k>=k0 && k<kmax);
      }


  //Step 2. Loop through embedded boundaries and determine the color at user-specified point from each boundary
  vector<int> mycolor(EBDS.size(), INT_MIN);
  for(int surf = 0; surf < (int)EBDS.size(); surf++) {

    vector<int> max_color(vertices.size(), INT_MIN);
    vector<int> min_color(vertices.size(), INT_MAX);
    vector<int> occluded(vertices.size(), 0);

    for(int n=0; n<8; n++) {
      if(owner[n]) {
        i = vertices[n][0];
        j = vertices[n][1];
        k = vertices[n][2];
        max_color[n] = color[surf][k][j][i];
        min_color[n] = color[surf][k][j][i];
        if(color[surf][k][j][i]==0)
          occluded[n] = 1;
      }
    }

    MPI_Allreduce(MPI_IN_PLACE, max_color.data(), 8, MPI_INT, MPI_MAX, comm);
    MPI_Allreduce(MPI_IN_PLACE, min_color.data(), 8, MPI_INT, MPI_MIN, comm);
    MPI_Allreduce(MPI_IN_PLACE, occluded.data(), 8, MPI_INT, MPI_MAX, comm);

    // check if any node in the box is occluded. If yes, exit.
    for(int n=0; n<8; n++) {
      if(occluded[n]>0) {
        print_error("*** Error: User-specified point (%e %e %e) is too close to embedded surface %d.\n",
                    point.x, point.y, point.z, surf);
        exit_mpi();
      }
    }

    // determine mycolor[surf]
    mycolor[surf] = INT_MIN;
    for(int n=0; n<8; n++) {
      if(max_color[n] != INT_MIN) {
        mycolor[surf] = max_color[n];
        break;
      }
    }
    if(mycolor[surf] == INT_MIN) {
      print_error("*** Error: Unable to determine the 'color' of point (%e %e %e) from surface %d.\n",
                  point.x, point.y, point.z, surf);
      exit_mpi();
    }
    for(int n=0; n<8; n++) {
      if(max_color[n] != INT_MIN && max_color[n] != mycolor[surf]) {
        print_error("*** Error: Found different colors (%d and %d) around point (%e %e %e), "
                    "using surface %d.\n",
                    mycolor[surf], max_color[n], point.x, point.y, point.z, surf);
        print_error("           This point may be too close to the embedded surface.\n");
        exit_mpi();
      }
    }
    for(int n=0; n<8; n++) {
      if(min_color[n] != INT_MAX && min_color[n] != mycolor[surf]) {
        print_error("*** Error: Found different colors (%d and %d) around point (%e %e %e), "
                    "using embedded surface %d.\n",
                    mycolor[surf], min_color[n], point.x, point.y, point.z, surf);
        print_error("           This point may be too close to the embedded surface.\n");
        exit_mpi();
      }
    }

  }

  //Step 3: Determine the "ruling surface" and color of this point.
  int ruling_surface = -1;
  int mycolor_final = INT_MIN;
  for(i=0; i<(int)mycolor.size(); i++) {
    assert(mycolor[i] != 0);
    if(mycolor[i]<0) { 
      if(ruling_surface == -1) {
        ruling_surface = i;
        mycolor_final = mycolor[i];
      }
      else { //need to resolve a "conflict"
        if((*EBDS[i]->ColorReachesBoundary_ptr)[-mycolor[i]] == 0) {
          if((*EBDS[ruling_surface]->ColorReachesBoundary_ptr)[-mycolor_final] == 0)
            print_warning("Warning: User-specified point (%e %e %e) is inside isolated regions in "
                          "two surfaces: %d and %d. Enforcing the latter.\n", point.x, point.y, point.z,
                          ruling_surface, i);
          ruling_surface = i;
          mycolor_final = mycolor[i];
        }
        else {
          if((*EBDS[ruling_surface]->ColorReachesBoundary_ptr)[-mycolor_final] == 1) {
            print_warning("Warning: User-specified point (%e %e %e) is inside isolated(*) regions in "
                          "two surfaces: %d and %d. Enforcing the latter.\n", point.x, point.y, point.z,
                          ruling_surface, i);
            ruling_surface = i;
            mycolor_final = mycolor[i];
          }
        }
      }
    }
  }
  if(ruling_surface<0) {
    print_error("*** Error: Unable to impose boundary condition based on point (%e %e %e). This point is connected"
                " to the far-field.\n", point.x, point.y, point.z);
    exit_mpi();
  }


  //Step 4: Update state variables and ID
  double*** ruling_colors = color[ruling_surface];
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        if (ruling_colors[k][j][i] == mycolor_final) { 
            v[k][j][i][0] = point.initialConditions.density;
            v[k][j][i][1] = point.initialConditions.velocity_x;
            v[k][j][i][2] = point.initialConditions.velocity_y;
            v[k][j][i][3] = point.initialConditions.velocity_z;
            v[k][j][i][4] = point.initialConditions.pressure;
            id[k][j][i]   = point.initialConditions.materialid;
          }
        }

  return make_pair(point.initialConditions.materialid, make_pair(ruling_surface, mycolor_final));
}

//---------------------------------------------------------------------

void
SpaceInitializer::InitializePhi(SpaceVariable3D &ID, SpaceOperator &spo,
                                vector<SpaceVariable3D*> &Phi, vector<LevelSetOperator*> &lso,
                                vector<unique_ptr<EmbeddedBoundaryDataSet> >* EBDS,
                                vector<pair<int,int> > &order,
                                multimap<int, pair<int,int> > &id2closure)
{
  for(auto it = iod.schemes.ls.dataMap.begin(); it != iod.schemes.ls.dataMap.end(); it++) {
    int matid = it->second->materialid;

    print("\n- Initializing Phi[%d] for tracking the boundary of material %d.\n", it->first, matid);

    auto closures = id2closure.equal_range(matid);
    if(closures.first != closures.second) {//found this matid
      vector<std::pair<int,int> > surf_and_color;
      for(auto it2 = closures.first; it2 != closures.second; it2++)
        surf_and_color.push_back(it2->second);
      if(it->second->init == LevelSetSchemeData::DISTANCE_CALCULATION)
        InitializePhiByDistance(order, spo.GetMeshCoordinates(), spo.GetMeshDeltaXYZ(),
                                *(spo.GetPointerToInnerGhostNodes()),
                                *(spo.GetPointerToOuterGhostNodes()),
                                *Phi[it->first], *lso[it->first], EBDS,
                                &surf_and_color);
      else //by reinitialization
        InitializePhiByReinitialization(ID, *Phi[it->first], *lso[it->first], &surf_and_color);
    }
    else {
      if(it->second->init == LevelSetSchemeData::DISTANCE_CALCULATION)
        InitializePhiByDistance(order, spo.GetMeshCoordinates(), spo.GetMeshDeltaXYZ(),
                                *(spo.GetPointerToInnerGhostNodes()),
                                *(spo.GetPointerToOuterGhostNodes()),
                                *Phi[it->first], *lso[it->first]);
      else //by reinitialization
        InitializePhiByReinitialization(ID, *Phi[it->first], *lso[it->first]);
    }
  }
}

//---------------------------------------------------------------------

void
SpaceInitializer::InitializePhiByDistance(vector<pair<int,int> > &order,
                                          SpaceVariable3D &coordinates,
                                          SpaceVariable3D &delta_xyz,
                                          vector<GhostPoint> &spo_ghost_nodes_inner,
                                          vector<GhostPoint> &spo_ghost_nodes_outer,
                                          SpaceVariable3D &Phi, LevelSetOperator &lso,
                                          vector<unique_ptr<EmbeddedBoundaryDataSet> >* EBDS,
                                          vector<pair<int,int> > *surf_and_color)
{

  int materialid = lso.GetMaterialID();

  //get info
  Vec3D***  coords = (Vec3D***)coordinates.GetDataPointer();
  double*** phi    = (double***)Phi.GetDataPointer();

  //first, set phi to be a large number everywhere
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++)
        phi[k][j][i] = domain_diagonal;

  //TODO: Add this feature later
  if(iod.ic.specified[IcData::LEVELSET]) {
    print_error("*** Error: Currently, cannot read level set from user-specified data file.\n");
    exit_mpi();
  }


  //initialize phi based on plane, cylinder-cones, and sphere.
  MultiInitialConditionsData &ic(iod.ic.multiInitialConditions);
  bool active_arbitrary_enclosure = false;

  for(auto&& obj : order) {
    //-------------------------------------
    // Internal Geometry ID:
    // 0: Plane: 0,
    // 1: CylinderCone
    // 2: CylinderSphere
    // 3: Sphere
    // 4: Parallelepiped
    // 5: Spheroid
    // 6: Custom-Geometry ("enclosure")
    // 7: Point
    //-------------------------------------
    double dist;
    if(obj.first == 0) {//plane
      auto it = ic.planeMap.dataMap.find(obj.second);
      assert(it != ic.planeMap.dataMap.end()); 
      if(it->second->initialConditions.materialid != materialid)
        continue; //not the right one

      if(it->second->inclusion != PlaneData::OVERRIDE) {
        print_error("*** Error: Level set initialization only supports Inclusion = Override at the moment.\n");
        exit_mpi();
      }

      Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
      Vec3D dir(it->second->nx, it->second->ny, it->second->nz);
      GeoTools::DistanceFromPointToPlane distCal(x0, dir); //dir will be normalized by the constructor

      for(int k=k0; k<kmax; k++)
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {
            dist = distCal.Calculate(coords[k][j][i]);
            if(fabs(dist)<fabs(phi[k][j][i]))
              phi[k][j][i] = -dist; //phi is negative inside the material subdomain, positive outside
          }
    }
    else if(obj.first == 1) {//cylinder-cone
      auto it = ic.cylinderconeMap.dataMap.find(obj.second);
      assert(it != ic.cylinderconeMap.dataMap.end());
      if(it->second->initialConditions.materialid != materialid)
        continue; //not the right one

      if(it->second->side != CylinderConeData::INTERIOR) {
        print_error("*** Error: Level set initialization only supports Side = Interior at the moment.\n");
        exit_mpi();
      }

      bool interior = (it->second->side == CylinderConeData::INTERIOR);

      Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
      Vec3D dir(it->second->nx, it->second->ny, it->second->nz);
      double L = it->second->L; //cylinder height
      double R = it->second->r; //cylinder radius
      double tan_alpha = tan(it->second->opening_angle_degrees/180.0*acos(-1.0));//opening angle
      double Hmax = R/tan_alpha;
      double H = std::min(it->second->cone_height, Hmax); //cone's height
      GeoTools::DistanceFromPointToCylinderCone distCal(x0,dir,L,R,tan_alpha,H);

      for(int k=k0; k<kmax; k++)
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {
            dist = distCal.Calculate(coords[k][j][i]);
            if(fabs(dist)<fabs(phi[k][j][i]))
              phi[k][j][i] = interior ? dist : -dist;
          }
    }
    else if(obj.first == 2) { //cylinder-sphere (i.e., possibly with spherical cap(s))
      auto it = ic.cylindersphereMap.dataMap.find(obj.second);
      assert(it != ic.cylindersphereMap.dataMap.end());
      if(it->second->initialConditions.materialid != materialid)
        continue; //not the right one

      if(it->second->inclusion != CylinderSphereData::OVERRIDE) {
        print_error("*** Error: Level set initialization only supports Inclusion = Override at the moment.\n");
        exit_mpi();
      }

      bool interior = (it->second->side == CylinderSphereData::INTERIOR);

      Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z); //center of base
      Vec3D dir(it->second->nx, it->second->ny, it->second->nz);
      double L = it->second->L; //cylinder height
      double R = it->second->r; //cylinder radius
      bool front_cap = (it->second->front_cap == CylinderSphereData::On);
      bool back_cap = (it->second->back_cap == CylinderSphereData::On);
      GeoTools::DistanceFromPointToCylinderSphere distCal(x0, dir, L, R, front_cap, back_cap);

      for(int k=k0; k<kmax; k++)
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {
            dist = distCal.Calculate(coords[k][j][i]);
            if(fabs(dist)<fabs(phi[k][j][i]))
              phi[k][j][i] = interior ? dist : -dist;
          }
    }
    else if(obj.first == 3) { //sphere
      auto it = ic.sphereMap.dataMap.find(obj.second);
      assert(it != ic.sphereMap.dataMap.end());
      if(it->second->initialConditions.materialid != materialid)
        continue; //not the right one

      if(it->second->inclusion != SphereData::OVERRIDE) {
        print_error("*** Error: Level set initialization only supports Inclusion = Override at the moment.\n");
        exit_mpi();
      }

      bool interior = (it->second->side == SphereData::INTERIOR);

      Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
      double R = it->second->radius;
      GeoTools::DistanceFromPointToSphere distCal(x0, R);

      for(int k=k0; k<kmax; k++)
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {
            dist = distCal.Calculate(coords[k][j][i]);
            if(fabs(dist)<fabs(phi[k][j][i]))
              phi[k][j][i] = interior ? dist : -dist;
          }

    }
    else if(obj.first == 4) { //parallelepiped
      auto it = ic.parallelepipedMap.dataMap.find(obj.second);
      assert(it != ic.parallelepipedMap.dataMap.end());
      if(it->second->initialConditions.materialid != materialid)
        continue; //not the right one

      if(it->second->inclusion != ParallelepipedData::OVERRIDE) {
        print_error("*** Error: Level set initialization only supports Inclusion = Override at the moment.\n");
        exit_mpi();
      }

      bool interior = (it->second->side == ParallelepipedData::INTERIOR);

      Vec3D x0(it->second->x0, it->second->y0, it->second->z0);
      Vec3D oa(it->second->ax, it->second->ay, it->second->az);  oa -= x0;
      Vec3D ob(it->second->bx, it->second->by, it->second->bz);  ob -= x0;
      Vec3D oc(it->second->cx, it->second->cy, it->second->cz);  oc -= x0;

      if(oa.norm()==0 || ob.norm()==0 || oc.norm()==0 || (oa^ob)*oc<=0.0) {
        print_error("*** Error: Detected error in a user-specified parallelepiped. "
                    "Overlapping vertices or violation of right-hand rule.\n");
        exit_mpi();
      }

      GeoTools::DistanceFromPointToParallelepiped distCal(x0, oa, ob, oc);

      for(int k=k0; k<kmax; k++)
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {
            dist = distCal.Calculate(coords[k][j][i]); //>0 outside the spheroid
            if(fabs(dist)<fabs(phi[k][j][i]))
              phi[k][j][i] = interior ? dist : -dist;
          }
    }
    else if(obj.first == 5) { //spheroid
      auto it = ic.spheroidMap.dataMap.find(obj.second);
      assert(it != ic.spheroidMap.dataMap.end());
      if(it->second->initialConditions.materialid != materialid)
        continue; //not the right one

      if(it->second->inclusion != SpheroidData::OVERRIDE) {
        print_error("*** Error: Level set initialization only supports Inclusion = Override at the moment.\n");
        exit_mpi();
      }

      bool interior = (it->second->side == SpheroidData::INTERIOR);

      Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
      Vec3D axis(it->second->axis_x, it->second->axis_y, it->second->axis_z);
      GeoTools::DistanceFromPointToSpheroid distCal(x0, axis, it->second->semi_length, it->second->radius);

      for(int k=k0; k<kmax; k++)
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {
            dist = distCal.Calculate(coords[k][j][i]);
            if(fabs(dist)<fabs(phi[k][j][i]))
              phi[k][j][i] = interior ? dist : -dist;
          }
    }
    else if(obj.first == 6) { //user-specified enclosure
      auto it = ic.enclosureMap.dataMap.find(obj.second);
      assert(it != ic.enclosureMap.dataMap.end());
      if(it->second->initialConditions.materialid != materialid)
        continue; //not the right one

      if(it->second->inclusion != UserSpecifiedEnclosureData::OVERRIDE) {
        print_error("*** Error: Level set initialization only supports Inclusion = Override at the moment.\n");
        exit_mpi();
      }

      Phi.RestoreDataPointerAndInsert();
      coordinates.RestoreDataPointerToLocalVector();
      bool active = InitializePhiWithinEnclosure(*it->second, coordinates, delta_xyz,
                                                 spo_ghost_nodes_inner, spo_ghost_nodes_outer, Phi, lso);
      phi = Phi.GetDataPointer();
      coords = (Vec3D***)coordinates.GetDataPointer();

      if(active)
        active_arbitrary_enclosure = true;
    }
    else if(obj.first == 7) { //point (must be inside an embedded closed surface)
      /* nothing to do here */
    }
    else {
      print_error("*** Error: Level-set initializer found unrecognized object id (%d).\n", obj.first);
      exit_mpi();
    }
  }


//SPECIAL TEST CASES
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


  // Now, take care of embedded surfaces
  int ebds_counter = 0;
  LevelSetSchemeData& iod_ls(lso.GetLevelSetSchemeData());

  if(EBDS && surf_and_color) {
    if(!lso.HasReinitializer()) {
      print_error("*** Error: Material %d is tracked by both a level set function and an embedded boundary.\n"
                  "           In this case, level set reinitialization must be turned on.\n", materialid);
      exit_mpi();
    }

    if(iod_ls.reinit.firstLayerTreatment != LevelSetReinitializationData::FIXED)
      print_warning("Warning: Material %d is tracked by both a level set function and an embedded boundary.\n"
                    "         In this case, it may be better to 'fix' the first layer nodes when "
                              "reinitializing the level set.\n", materialid);

    for(auto&& mypair : *surf_and_color) {

      int surf = mypair.first;
      int mycolor = mypair.second;

      int nLayer = (*EBDS)[surf]->Phi_nLayer;
      if(nLayer<2)
        print_warning("Warning: Material %d is tracked by a level set function AND an embedded boundary.\n"
                      "         Intersector should compute unsigned distance ('phi') for at least two layers\n"
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


  // Finalize
  coordinates.RestoreDataPointerToLocalVector();
  Phi.RestoreDataPointerAndInsert();

  lso.ApplyBoundaryConditions(Phi);

  if(iod_ls.bandwidth < INT_MAX) {//narrow-band level set method
    assert(lso.HasReinitializer());
    lso.ConstructNarrowBandInReinitializer(Phi);
  }

  MPI_Allreduce(MPI_IN_PLACE, &ebds_counter, 1, MPI_INT, MPI_SUM, comm);

/*
  Phi.StoreMeshCoordinates(coordinates);
  Phi.WriteToVTRFile("LS.vtr", "phi");
  MPI_Barrier(comm);
  exit_mpi();
*/

  if(ebds_counter>0) {
    //reinit must have been created (not NULL)
    print("  o Updated Phi (material id: %d) at %d nodes based on embedded boundary. "
          "Going to reinitialize phi.\n", materialid, ebds_counter);
    lso.Reinitialize(0.0, 1.0, 0.0, Phi, 600, true/*must do*/); //first 3 inputs are not used ("must do"!)
  }
  else if(active_arbitrary_enclosure && lso.HasReinitializer()) {
    //if reinit is NULL, should have reinitialized phi (full-domain) in InitializePhiWithinEnclosure
    print("  o Updated Phi (material id: %d) based on surface mesh(es). "
          "Going to reinitialize phi.\n", materialid);
    lso.Reinitialize(0.0, 1.0, 0.0, Phi, 600, true/*must do*/); //first 3 inputs are not used ("must do"!)
  }

}

//---------------------------------------------------------------------

bool 
SpaceInitializer::InitializePhiWithinEnclosure(UserSpecifiedEnclosureData &enclosure,
                                               SpaceVariable3D &coordinates, 
                                               SpaceVariable3D &delta_xyz,
                                               vector<GhostPoint> &spo_ghost_nodes_inner,
                                               vector<GhostPoint> &spo_ghost_nodes_outer,
                                               SpaceVariable3D &Phi, LevelSetOperator &lso)
{
  bool applied_ic = false;
  double*** phi = Phi.GetDataPointer();

  // Create a reinitializer (if not available) for one-time use within this function
  // This is a full-domain reinitializer
  LevelSetReinitializer *tmp_reinit = NULL; 
  if(!lso.HasReinitializer())
    tmp_reinit = new LevelSetReinitializer(comm, dms, lso.GetLevelSetSchemeData(), coordinates, delta_xyz,
                                           *(lso.GetPointerToInnerGhostNodes()),
                                           *(lso.GetPointerToOuterGhostNodes()));//they carry level set b.c.

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
  if(tmp_reinit) {
    lso.ApplyBoundaryConditions(Phi); //need this before reinitialization!
    tmp_reinit->ReinitializeFullDomain(Phi, 600); //one-time use
    tmp_reinit->Destroy(); 
    delete tmp_reinit;
  }

  //clean-up
  EBDS->Color_ptr->RestoreDataPointerToLocalVector();
  EBDS->Phi_ptr->RestoreDataPointerToLocalVector();
  EBDS->ClosestPointIndex_ptr->RestoreDataPointerToLocalVector();

  embed.Destroy();

  return applied_ic;
}


//---------------------------------------------------------------------

void
SpaceInitializer::InitializePhiByReinitialization(SpaceVariable3D &ID,
                                                  SpaceVariable3D &Phi, LevelSetOperator &lso,
                                                  vector<pair<int,int> > *surf_and_color)
{
  int materialid = lso.GetMaterialID();
  double*** id  = ID.GetDataPointer();
  double*** phi = Phi.GetDataPointer();

  // Calculate phi near subdomain boundary
  int myid;
  double dist;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        myid = id[k][j][i];
        dist = global_mesh.GetMinDXYZ(Int3(i,j,k));
        if(myid == materialid) {
          if(i-1>=0 && id[k][j][i-1] != myid)
            dist = std::min(dist, 0.5*(global_mesh.GetX(i)-global_mesh.GetX(i-1)));
          if(i+1<NX && id[k][j][i+1] != myid)
            dist = std::min(dist, 0.5*(global_mesh.GetX(i+1)-global_mesh.GetX(i)));
          if(j-1>=0 && id[k][j-1][i] != myid)
            dist = std::min(dist, 0.5*(global_mesh.GetY(j)-global_mesh.GetY(j-1)));
          if(j+1<NY && id[k][j+1][i] != myid)
            dist = std::min(dist, 0.5*(global_mesh.GetY(j+1)-global_mesh.GetY(j)));
          if(k-1>=0 && id[k-1][j][i] != myid)
            dist = std::min(dist, 0.5*(global_mesh.GetZ(k)-global_mesh.GetZ(k-1)));
          if(k+1<NZ && id[k+1][j][i] != myid)
            dist = std::min(dist, 0.5*(global_mesh.GetZ(k+1)-global_mesh.GetZ(k)));

          phi[k][j][i] = -dist;
        }
        else {
          if(i-1>=0 && id[k][j][i-1] == materialid)
            dist = std::min(dist, 0.5*(global_mesh.GetX(i)-global_mesh.GetX(i-1)));
          if(i+1<NX && id[k][j][i+1] == materialid)
            dist = std::min(dist, 0.5*(global_mesh.GetX(i+1)-global_mesh.GetX(i)));
          if(j-1>=0 && id[k][j-1][i] == materialid)
            dist = std::min(dist, 0.5*(global_mesh.GetY(j)-global_mesh.GetY(j-1)));
          if(j+1<NY && id[k][j+1][i] == materialid)
            dist = std::min(dist, 0.5*(global_mesh.GetY(j+1)-global_mesh.GetY(j)));
          if(k-1>=0 && id[k-1][j][i] == materialid)
            dist = std::min(dist, 0.5*(global_mesh.GetZ(k)-global_mesh.GetZ(k-1)));
          if(k+1<NZ && id[k+1][j][i] == materialid)
            dist = std::min(dist, 0.5*(global_mesh.GetZ(k+1)-global_mesh.GetZ(k)));

          phi[k][j][i] = dist;
        }
      }

  ID.RestoreDataPointerToLocalVector();
  Phi.RestoreDataPointerAndInsert();

  lso.ApplyBoundaryConditions(Phi);

  // Reinitialize Phi
  if(!lso.HasReinitializer()) {
    print_error("*** Error: A reinitializer needs to be specified for material %d's level-set.\n",
                materialid); 
    exit_mpi();
  }

  LevelSetSchemeData& iod_ls(lso.GetLevelSetSchemeData());

  if(surf_and_color && !surf_and_color->empty()) {
    if(iod_ls.reinit.firstLayerTreatment != LevelSetReinitializationData::FIXED)
      print_warning("Warning: Material %d is tracked by both a level set function and an embedded boundary.\n"
                    "         In this case, it may be better to 'fix' the first layer nodes when "
                              "reinitializing the level set.\n", materialid);
  }

  if(iod_ls.bandwidth < INT_MAX) //narrow-band level set method
    lso.ConstructNarrowBandInReinitializer(Phi);

  print("  o Set Phi (material id: %d) near subdomain boundary. Going to reinitialize it.\n", materialid);
  lso.Reinitialize(0.0, 1.0, 0.0, Phi, 600, true/*must do*/);
}

//---------------------------------------------------------------------

