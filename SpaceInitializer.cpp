/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <SpaceInitializer.h>
#include <SpaceOperator.h>
#include <LevelSetOperator.h>
#include <EmbeddedBoundaryDataSet.h>
#include <DistancePointToParallelepiped.h>
#include <DistancePointToSpheroid.h>
#include <DistancePointToSphere.h>
#include <DistancePointToCylinderCone.h>
#include <DistancePointToCylinderSphere.h>
#include <DistancePointToPlane.h>

using std::vector;
using std::multimap;
using std::set;
using std::pair;
using std::make_pair;
using std::unique_ptr;

//---------------------------------------------------------------------

SpaceInitializer::SpaceInitializer(MPI_Comm &comm_, IoData &iod_, GlobalMeshInfo &global_mesh_,
                                   SpaceVariable3D &coordinates)
                : comm(comm_), iod(iod_), global_mesh(global_mesh_)
{
  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);
  coordinates.GetGlobalSize(&NX, &NY, &NZ);
}

//---------------------------------------------------------------------

SpaceInitializer::~SpaceInitialier()
{ }

//---------------------------------------------------------------------

void
SpaceInitializer::Destroy()
{ }

//---------------------------------------------------------------------

multimap<int, pair<int,int> >
SpaceInitializer::SetInitialCondition(SpaceVariable3D &V, SpaceVariable3D &ID, vector<SpaceVariable3D*> &Phi,
                                      SpaceOperator &spo, vector<SpaceVariable3D*> &lso,
                                      unique_ptr<vector<unique_ptr<EmbeddedBoundaryDataSet> > > EBDS)
{

  //----------------------------------------------------------------
  // Step 1: Setup the operation order for user-specified geometries.
  //         Also, create signed distance calculators for later use.
  //----------------------------------------------------------------
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
  int nGeom = OrderUserSpecifiedGeometries(order);

  // Step 2: Apply initial conditions to V, ID (B.C. applied to V, but not ID!)
  multimap<int, pair<int,int> > id2closure = InitializeVandID(V, ID, spo, EBDS, order);

  // Step 3:: Apply initial conditions to Phi
  InitializePhi(ID, Phi, lso, EBDS, order, id2closure);

  // Step 4:: Apply flooding if specified by user
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
SpaceInitializer::AddGeomToVector(int o, int type, int ind, string& name, vector<std::pair<int,int> > &order,
                                  vector<int> &user_specified_order)
{
  int nGeom = order.size();
  if(o<0 || o>=nGeom) {
    print_error("*** Error: Detected incorrect Plane order (%d). Range: [0, %d)\n", o, nGeom);
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
                                   unique_ptr<vector<unique_ptr<EmbeddedBoundaryDataSet> > > EBDS,
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
            it->first, x0[0], x0[1], x0[2], dir[0], dir[1], dir[2]
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
            if(interior  && distCal.Calculate(coords[k][j][i])<0 ||
               !interior && distCal.Calculate(coords[k][j][i])>0) {
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
            if(interior  && distCal.Calculate(coords[k][j][i])<0 ||
               !interior && distCal.Calculate(coords[k][j][i])>0) {
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
            if(interior  && distCal.Calculate(coords[k][j][i])<0 ||
               !interior && distCal.Calculate(coords[k][j][i])>0) {
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
            if(interior  && distCal.Calculate(coords[k][j][i])<0 ||
               !interior && distCal.Calculate(coords[k][j][i])>0) {
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
            if(interior  && distCal.Calculate(coords[k][j][i])<0 ||
               !interior && distCal.Calculate(coords[k][j][i])>0) {
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
      InitializeVandIDWithinEnclosure(*it->second, v, id);
      coords = (Vec3D***)coordinates.GetDataPointer();
    }
    else if(obj.first == 7) { //point (must be inside an embedded closed surface)

      for(auto it=ic.pointMap.dataMap.begin(); it!=ic.pointMap.dataMap.end(); it++) {

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
  embed.SetCommAndMeshInfo(dm_all, coordinates, ghost_nodes_inner, ghost_nodes_outer,
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
        owner.push_back(coordinates.IsHere(i,j,k,false));
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
SpaceInitializer::InitializePhi(SpaceVariable3D &ID, vector<SpaceVariable3D*> &Phi,
                                vector<SpaceVariable3D*> &lso,
                                unique_ptr<vector<unique_ptr<EmbeddedBoundaryDataSet> > > EBDS,
                                vector<pair<int,int> > &order,
                                multimap<int, pair<int,int> > &id2closure)
{
  set<int> ls_tracker;
  for(auto it = iod.schemes.ls.dataMap.begin(); it != iod.schemes.ls.dataMap.end(); it++) {
    int matid = it->second->materialid;

    print("\n- Initializing level set function (%d) for tracking the boundary of material %d.\n",

    if(it->second->init = LevelSetSchemeData::DISTANCE_CALCULATION) {
      auto closures = id2closure.equal_range(matid);
      if(closures.first != closures.second) {//found this matid
        vector<std::pair<int,int> > surf_and_color;
        for(auto it2 = closures.first; it2 != closures.second; it2++)
          surf_and_color.push_back(it2->second);
        InitializePhiByDistanceCalculation(*Phi[it->first],
                                            embed->GetPointerToEmbeddedBoundaryData(),
                                            &surf_and_color);
      } else
        InitializePhiByDistanceCalculation(*Phi[it->first]);
    }
    else { //by reinitialization
      InitializePhiByReinitialization(*Phi[it->first]);
    }
  }
}

//---------------------------------------------------------------------

I AM HERE!

//---------------------------------------------------------------------


//---------------------------------------------------------------------

