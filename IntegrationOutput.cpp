/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/
#include <IntegrationOutput.h>
#include <MathTools.h>
#include <DistancePointToParallelepiped.h>
#include <DistancePointToSpheroid.h>
#include <DistancePointToSphere.h>
#include <DistancePointToCylinderCone.h>
#include <DistancePointToCylinderSphere.h>
#include <DistancePointToPlane.h>
using std::vector;
using std::pair;

//-----------------------------------------------------------------------------------------------------

IntegrationOutput::IntegrationOutput(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod,
                                     LaserAbsorptionSolver* laser_,
                                     std::vector<VarFcnBase*> &vf_, SpaceVariable3D& coordinates_,
                                     SpaceVariable3D& delta_xyz_, SpaceVariable3D& cell_volume_)
                        :comm(comm_), iod_output(iod.output), laser(laser_), vf(vf_), coordinates(coordinates_),
                         delta_xyz(delta_xyz_), cell_volume(cell_volume_), mesh_type(iod.mesh.type),
                         Tag(comm_, &(dm_all_.ghosted1_1dof))
{

  numMaterials = iod.eqs.materials.dataMap.size() + 1; //an extra one for "ghost/inactive"

  int num_ints = iod.output.integrations.dataMap.size();
  if(num_ints==0)
    return; //nothing to do...

  if(num_ints>32) { //very unlikely... may exceed the number of bits that "tag" can handle
    print_error("*** Error: Cannot handle more than 32 integration outputs.\n");
    exit_mpi();
  }

  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);

  Tag.SetConstantValue(0.0, true); //set default tag to 0
  double*** tag = Tag.GetDataPointer();

  frequency.resize(num_ints, -1.0);
  frequency_dt.resize(num_ints, -1.0);
  last_snapshot_time.resize(num_ints, -1.0);
  files.resize(num_ints);
  std::set<int> used_indices;
  for(auto&& integral : iod.output.integrations.dataMap) {
    if(integral.second->frequency<0 && integral.second->frequency_dt<0) {
      print_error("*** Error: Detected integration plot without frequency specification.\n");
      exit_mpi();
    }
    int index = integral.first;
    if(index<0 || index>=num_ints) {
      print_error("*** Error: Integration output indices must start at 0 and have no gaps. Found %d.\n", index);
      exit_mpi();
    }
    frequency[index]    = integral.second->frequency;
    frequency_dt[index] = integral.second->frequency_dt;
    
    SetupIntegrationDomain(tag, index, *integral.second);

    // Setup files
    vector<FILE*> &file(files[index]);
    file.resize(IntegrationData::SIZE, NULL);

    int spn = strlen(iod_output.prefix) + 1;

    if(integral.second->volume[0] != 0) {
      char *filename = new char[spn + strlen(integral.second->volume)];
      snprintf(filename, spn + strlen(integral.second->volume), "%s%s", iod_output.prefix, integral.second->volume);
      file[IntegrationData::VOLUME] = fopen(filename, "w");
      delete [] filename;
    }

    if(integral.second->mass[0] != 0) {
      char *filename = new char[spn + strlen(integral.second->mass)];
      snprintf(filename, spn + strlen(integral.second->mass), "%s%s", iod_output.prefix, integral.second->mass);
      file[IntegrationData::MASS] = fopen(filename, "w");
      delete [] filename;
    }

    if(integral.second->momentum[0] != 0) {
      char *filename = new char[spn + strlen(integral.second->momentum)];
      snprintf(filename, spn + strlen(integral.second->momentum), "%s%s", iod_output.prefix,
               integral.second->momentum);
      file[IntegrationData::MOMENTUM] = fopen(filename, "w");
      delete [] filename;
    }

    if(integral.second->total_energy[0] != 0) {
      char *filename = new char[spn + strlen(integral.second->total_energy)];
      snprintf(filename, spn + strlen(integral.second->total_energy), "%s%s",
               iod_output.prefix, integral.second->total_energy);
      file[IntegrationData::TOTAL_ENERGY] = fopen(filename, "w");
      delete [] filename;
    }

    if(integral.second->total_enthalpy[0] != 0) {
      char *filename = new char[spn + strlen(integral.second->total_enthalpy)];
      snprintf(filename, spn + strlen(integral.second->total_enthalpy), "%s%s",
               iod_output.prefix, integral.second->total_enthalpy);
      file[IntegrationData::TOTAL_ENTHALPY] = fopen(filename, "w");
      delete [] filename;
    }

    if(integral.second->kinetic_energy[0] != 0) {
      char *filename = new char[spn + strlen(integral.second->kinetic_energy)];
      snprintf(filename, spn + strlen(integral.second->kinetic_energy), "%s%s",
               iod_output.prefix, integral.second->kinetic_energy);
      file[IntegrationData::KINETIC_ENERGY] = fopen(filename, "w");
      delete [] filename;
    }

    if(integral.second->internal_energy[0] != 0) {
      char *filename = new char[spn + strlen(integral.second->internal_energy)];
      snprintf(filename, spn + strlen(integral.second->internal_energy), "%s%s",
               iod_output.prefix, integral.second->internal_energy);
      file[IntegrationData::INTERNAL_ENERGY] = fopen(filename, "w");
      delete [] filename;
    }
  
    if(integral.second->potential_energy[0] != 0) {
      char *filename = new char[spn + strlen(integral.second->potential_energy)];
      snprintf(filename, spn + strlen(integral.second->potential_energy), "%s%s",
               iod_output.prefix, integral.second->potential_energy);
      file[IntegrationData::POTENTIAL_ENERGY] = fopen(filename, "w");
      delete [] filename;
    }

    if(integral.second->laser_radiation[0] != 0) {
      if(!laser) {
        print_error("*** Error: Laser radiation integration requested without specifying laser source.\n");
        exit_mpi();
      }
      char *filename = new char[spn + strlen(integral.second->laser_radiation)];
      snprintf(filename, spn + strlen(integral.second->laser_radiation), "%s%s",
               iod_output.prefix, integral.second->laser_radiation);
      file[IntegrationData::LASER_RADIATION] = fopen(filename, "w");
      delete [] filename;
    }
  
    for(int i=0; i<IntegrationData::SIZE; i++)
      if(file[i]) { //write header
        print(file[i], "## Number of materials: %d (0 - %d); ID for ghost/inactive cells: %d.\n",
              numMaterials-1, numMaterials-2, numMaterials-1);
        print(file[i], "## Time step  |  Time  ");
        for(int j=0; j<numMaterials; j++)
          print(file[i], "|  Solution (Mat. %d)  ", j);
        print(file[i],"   |  Sum (including ghost/inactive)\n");
        mpi_barrier();
        print_flush(file[i]);
      }

    used_indices.insert(index);
  }

  Tag.RestoreDataPointerAndInsert();

  if(used_indices.size() != iod_output.integrations.dataMap.size()) {
    print_error("*** Error: Detected error (e.g., duplicates) in integration output indices.\n");
    exit_mpi();
  }

}

//----------------------------------------------------------------------------------------------------

IntegrationOutput::~IntegrationOutput()
{
  for(auto&& file : files)
    for(auto&& f : file)
      if(f) fclose(f); 
}

//------------------------------------------------------------------------------------------------------

void
IntegrationOutput::SetupIntegrationDomain(double*** tag, int index, IntegrationData &integral)
{

  print("\n- Setting up integration domain %d.\n", index);

  vector<pair<int,int> > order;  //<geom type, geom dataMap index>
  int nGeom = OrderUserSpecifiedGeometries(integral, order);

  int mask = 1 << index; //only the index-th bit is 1; 0 eleswhere

  // start with setting the entire domain to be integrated
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) 
        tag[k][j][i] = (int)tag[k][j][i] | mask; //set the index-th bit to 1

  if(nGeom==0)
    return;

  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();

  for(auto&& obj : order) {
    //-------------------------------------
    // Internal Geometry ID:
    // 0: Plane 
    // 1: CylinderCone
    // 2: CylinderSphere
    // 3: Sphere
    // 4: Parallelepiped
    // 5: Spheroid
    //-------------------------------------
    if(obj.first == 0) {//plane
      auto it = integral.planeMap.dataMap.find(obj.second);
      assert(it != integral.planeMap.dataMap.end());

      bool set_intersection = (it->second->inclusion == PlaneData::INTERSECTION);
      if(it->second->inclusion == PlaneData::OVERRIDE) {
        print_warning("  o Warning: 'Override' is not supported. Treating it as 'Intersection'.\n");
        set_intersection = true;
      }

      Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
      Vec3D dir(it->second->nx, it->second->ny, it->second->nz);
      GeoTools::DistanceFromPointToPlane distCal(x0, dir); //dir will be normalized by the constructor

      print("  o Found Plane[%d]: (%e %e %e) | normal: (%e %e %e).\n",
            it->first, x0[0], x0[1], x0[2], dir[0], dir[1], dir[2]);

      for(int k=k0; k<kmax; k++)
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {
            if(distCal.Calculate(coords[k][j][i])>0) {
              if(!set_intersection)
                tag[k][j][i] = (int)tag[k][j][i] | mask;
            } else {
              if(set_intersection)
                tag[k][j][i] = (int)tag[k][j][i] & ~mask; //set bit to 0
            }
          }
    }
    else if(obj.first == 1) {//cylinder-cone
      auto it = integral.cylinderconeMap.dataMap.find(obj.second);
      assert(it != integral.cylinderconeMap.dataMap.end());

      bool interior = (it->second->side == CylinderConeData::INTERIOR);
      bool set_intersection = (it->second->inclusion == CylinderConeData::INTERSECTION);
      if(it->second->inclusion == CylinderConeData::OVERRIDE) {
        print_warning("  o Warning: 'Override' is not supported. Treating it as 'Intersection'.\n");
        set_intersection = true;
      }

      Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
      Vec3D dir(it->second->nx, it->second->ny, it->second->nz);
      double L = it->second->L; //cylinder height
      double R = it->second->r; //cylinder radius
      double tan_alpha = tan(it->second->opening_angle_degrees/180.0*acos(-1.0));//opening angle
      double Hmax = R/tan_alpha;
      double H = std::min(it->second->cone_height, Hmax); //cone's height
      GeoTools::DistanceFromPointToCylinderCone distCal(x0,dir,L,R,tan_alpha,H);

      print("  o Found Cylinder-Cone[%d]: Base center: (%e %e %e), Axis: (%e %e %e).\n",
            it->first, x0[0], x0[1], x0[2], dir[0], dir[1], dir[2]);

      for(int k=k0; k<kmax; k++)
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {
            if((interior  && distCal.Calculate(coords[k][j][i])<0) ||
               (!interior && distCal.Calculate(coords[k][j][i])>0)) {
              if(!set_intersection) 
                tag[k][j][i] = (int)tag[k][j][i] | mask;
            } else {
              if(set_intersection)
                tag[k][j][i] = (int)tag[k][j][i] & ~mask; //set bit to 0
            }
          }
    }
    else if(obj.first == 2) { //cylinder-sphere (i.e., possibly with spherical cap(s))
      auto it = integral.cylindersphereMap.dataMap.find(obj.second);
      assert(it != integral.cylindersphereMap.dataMap.end());

      bool interior = (it->second->side == CylinderSphereData::INTERIOR);
      bool set_intersection = (it->second->inclusion == CylinderSphereData::INTERSECTION);
      if(it->second->inclusion == CylinderSphereData::OVERRIDE) {
        print_warning("  o Warning: 'Override' is not supported. Treating it as 'Intersection'.\n");
        set_intersection = true;
      }

      Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z); //center of base
      Vec3D dir(it->second->nx, it->second->ny, it->second->nz);
      double L = it->second->L; //cylinder height
      double R = it->second->r; //cylinder radius
      bool front_cap = (it->second->front_cap == CylinderSphereData::On);
      bool back_cap = (it->second->back_cap == CylinderSphereData::On);
      GeoTools::DistanceFromPointToCylinderSphere distCal(x0, dir, L, R, front_cap, back_cap);

      print("  o Found Cylinder-Sphere[%d]: Base center: (%e %e %e), Axis: (%e %e %e).\n",
            it->first, x0[0], x0[1], x0[2], dir[0], dir[1], dir[2]);

      for(int k=k0; k<kmax; k++)
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {
            if((interior  && distCal.Calculate(coords[k][j][i])<0) ||
               (!interior && distCal.Calculate(coords[k][j][i])>0)) {
              if(!set_intersection) 
                tag[k][j][i] = (int)tag[k][j][i] | mask;
            } else {
              if(set_intersection)
                tag[k][j][i] = (int)tag[k][j][i] & ~mask; //set bit to 0
            }
          }
    }
    else if(obj.first == 3) { //sphere
      auto it = integral.sphereMap.dataMap.find(obj.second);
      assert(it != integral.sphereMap.dataMap.end());

      bool interior = (it->second->side == SphereData::INTERIOR);
      bool set_intersection = (it->second->inclusion == SphereData::INTERSECTION);
      if(it->second->inclusion == SphereData::OVERRIDE) {
        print_warning("  o Warning: 'Override' is not supported. Treating it as 'Intersection'.\n");
        set_intersection = true;
      }

      Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
      double R = it->second->radius;
      GeoTools::DistanceFromPointToSphere distCal(x0, R);

      print("  o Found Sphere[%d]: Center: (%e %e %e), Radius: %e.\n", it->first, x0[0], x0[1], x0[2], R);

      for(int k=k0; k<kmax; k++)
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {
            if((interior  && distCal.Calculate(coords[k][j][i])<0) ||
               (!interior && distCal.Calculate(coords[k][j][i])>0)) {
              if(!set_intersection) 
                tag[k][j][i] = (int)tag[k][j][i] | mask;
            } else {
              if(set_intersection)
                tag[k][j][i] = (int)tag[k][j][i] & ~mask; //set bit to 0
            }
          }
    }
    else if(obj.first == 4) { //parallelepiped
      auto it = integral.parallelepipedMap.dataMap.find(obj.second);
      assert(it != integral.parallelepipedMap.dataMap.end());

      bool interior = (it->second->side == ParallelepipedData::INTERIOR);
      bool set_intersection = (it->second->inclusion == ParallelepipedData::INTERSECTION);
      if(it->second->inclusion == ParallelepipedData::OVERRIDE) {
        print_warning("  o Warning: 'Override' is not supported. Treating it as 'Intersection'.\n");
        set_intersection = true;
      }

      Vec3D x0(it->second->x0, it->second->y0, it->second->z0);
      Vec3D oa(it->second->ax, it->second->ay, it->second->az);  oa -= x0;
      Vec3D ob(it->second->bx, it->second->by, it->second->bz);  ob -= x0;
      Vec3D oc(it->second->cx, it->second->cy, it->second->cz);  oc -= x0;

      if(oa.norm()==0 || ob.norm()==0 || oc.norm()==0 || (oa^ob)*oc<=0.0) {
        print_error("*** Error: (Integration) Detected error in a user-specified parallelepiped. "
                    "Overlapping vertices or violation of right-hand rule.\n");
        exit_mpi();
      }

      GeoTools::DistanceFromPointToParallelepiped distCal(x0, oa, ob, oc);
      print("  o Found Parallelepiped[%d]: X(%e %e %e).\n", it->first, x0[0], x0[1], x0[2]);

      for(int k=k0; k<kmax; k++)
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {
            if((interior  && distCal.Calculate(coords[k][j][i])<0) ||
               (!interior && distCal.Calculate(coords[k][j][i])>0)) {
              if(!set_intersection) 
                tag[k][j][i] = (int)tag[k][j][i] | mask;
            } else {
              if(set_intersection)
                tag[k][j][i] = (int)tag[k][j][i] & ~mask; //set bit to 0
            }
          }
    }
    else if(obj.first == 5) { //spheroid
      auto it = integral.spheroidMap.dataMap.find(obj.second);
      assert(it != integral.spheroidMap.dataMap.end());

      bool interior = (it->second->side == SpheroidData::INTERIOR);
      bool set_intersection = (it->second->inclusion == SpheroidData::INTERSECTION);
      if(it->second->inclusion == SpheroidData::OVERRIDE) {
        print_warning("  o Warning: 'Override' is not supported. Treating it as 'Intersection'.\n");
        set_intersection = true;
      }

      Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
      Vec3D axis(it->second->axis_x, it->second->axis_y, it->second->axis_z);
      GeoTools::DistanceFromPointToSpheroid distCal(x0, axis, it->second->semi_length, it->second->radius);

      print("  o Found Spheroid[%d]: Center: (%e %e %e), Axis: (%e %e %e).\n", it->first,
            x0[0], x0[1], x0[2], axis[0], axis[1], axis[2]);

      for(int k=k0; k<kmax; k++)
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {
            if((interior  && distCal.Calculate(coords[k][j][i])<0) ||
               (!interior && distCal.Calculate(coords[k][j][i])>0)) {
              if(!set_intersection) 
                tag[k][j][i] = (int)tag[k][j][i] | mask;
            } else {
              if(set_intersection)
                tag[k][j][i] = (int)tag[k][j][i] & ~mask; //set bit to 0
            }
          }
    }
    else {
      print_error("*** Error: Found unrecognized object id (%d).\n", obj.first);
      exit_mpi();
    }

  }

  coordinates.RestoreDataPointerToLocalVector();

}

//------------------------------------------------------------------------------------------------------
// This function is VERY similar to SpaceInitializer::OrderUserSpecifiedGeometries
int
IntegrationOutput::OrderUserSpecifiedGeometries(IntegrationData &integral, vector<std::pair<int,int> > &order)
{
  //-------------------------------------
  // Internal Geometry ID:
  // 0: Plane
  // 1: CylinderCone
  // 2: CylinderSphere
  // 3: Sphere
  // 4: Parallelepiped
  // 5: Spheroid
  //-------------------------------------

  int nGeom = integral.planeMap.dataMap.size() //0
            + integral.cylinderconeMap.dataMap.size() //1
            + integral.cylindersphereMap.dataMap.size() //2
            + integral.sphereMap.dataMap.size() //3
            + integral.parallelepipedMap.dataMap.size() //4
            + integral.spheroidMap.dataMap.size(); //5
  order.assign(nGeom, std::make_pair(-1,-1));
  vector<int> user_specified_order(nGeom, -1); //for verification of user's input

  for(auto&& obj : integral.planeMap.dataMap)
    AddGeomToVector(obj.second->order, 0, obj.first, "Plane", order, user_specified_order);
  for(auto&& obj : integral.cylinderconeMap.dataMap)
    AddGeomToVector(obj.second->order, 1, obj.first, "CylinderCone", order, user_specified_order);
  for(auto&& obj : integral.cylindersphereMap.dataMap)
    AddGeomToVector(obj.second->order, 2, obj.first, "CylinderSphere", order, user_specified_order);
  for(auto&& obj : integral.sphereMap.dataMap)
    AddGeomToVector(obj.second->order, 3, obj.first, "Sphere", order, user_specified_order);
  for(auto&& obj : integral.parallelepipedMap.dataMap)
    AddGeomToVector(obj.second->order, 4, obj.first, "Parallelepiped", order, user_specified_order);
  for(auto&& obj : integral.spheroidMap.dataMap)
    AddGeomToVector(obj.second->order, 5, obj.first, "Spheroid", order, user_specified_order);

  //verification
  for(int i=0; i<(int)user_specified_order.size()-1; i++)
    assert(user_specified_order[i]<=user_specified_order[i+1]);

  return nGeom;
}

//------------------------------------------------------------------------------------------------------
// This function is the same as SpaceInitializer::AddGeomToVector (this duplication allows better modularity)
void
IntegrationOutput::AddGeomToVector(int o, int type, int ind, string name, vector<std::pair<int,int> > &order,
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

//------------------------------------------------------------------------------------------------------

void
IntegrationOutput::WriteIntegrationResults(double time, double dt, int time_step, SpaceVariable3D &V,
                                           SpaceVariable3D &ID, SpaceVariable3D* L, bool force_write)
{
  if(files.empty()) //nothing to write
    return;

  // Get data
  double***   tag = Tag.GetDataPointer();
  Vec5D***      v = (Vec5D***)V.GetDataPointer();
  double***    id = ID.GetDataPointer();
  double***  cell = cell_volume.GetDataPointer();
  Vec3D*** coords = NULL;
  Vec3D***   dxyz = NULL;
  double***     l = NULL;

  if(mesh_type == MeshData::SPHERICAL || mesh_type == MeshData::CYLINDRICAL) {
    coords = (Vec3D***)coordinates.GetDataPointer();
    dxyz   = (Vec3D***)delta_xyz.GetDataPointer();
  }
  for(auto&& file : files) {
    if(file[IntegrationData::LASER_RADIATION]) {
      if(!L || !laser) {
        print_error("*** Error: Requested laser radiation integration, but laser source is not specified.\n");
        exit_mpi();
      }
      l = L->GetDataPointer();
      break;
    }
  }
  
  int index = 0;
  for(auto&& file : files) {
    if(!isTimeToWrite(time, dt, time_step, frequency_dt[index], frequency[index], last_snapshot_time[index],
                      force_write))
      continue;

    if(file[IntegrationData::VOLUME]) {
      print(file[IntegrationData::VOLUME], "%10d    %16.8e    ", time_step, time);
      double volume[numMaterials];
      double sum = 0.0;
      IntegrateVolume(index, tag, coords, dxyz, cell, id, volume);
      for(int i=0; i<numMaterials; i++) {
        print(file[IntegrationData::VOLUME], "%16.8e  ", volume[i]);
        sum += volume[i];
      } 
      print(file[IntegrationData::VOLUME], "%16.8e\n", sum);
      print_flush(file[IntegrationData::VOLUME]);
    }

    if(file[IntegrationData::MASS]) {
      print(file[IntegrationData::MASS], "%10d    %16.8e    ", time_step, time);
      double mass[numMaterials];
      double sum = 0.0;
      IntegrateMass(index, tag, coords, dxyz, cell, v, id, mass);
      for(int i=0; i<numMaterials; i++) {
        print(file[IntegrationData::MASS], "%16.8e  ", mass[i]);
        sum += mass[i];
      }
      print(file[IntegrationData::MASS], "%16.8e\n", sum);
      print_flush(file[IntegrationData::MASS]);
    } 

    if(file[IntegrationData::MOMENTUM]) {
      print(file[IntegrationData::MOMENTUM], "%10d    %16.8e    ", time_step, time);
      Vec3D momentum[numMaterials];
      Vec3D sum = 0.0;
      IntegrateMomentum(index, tag, coords, dxyz, cell, v, id, momentum);
      for(int i=0; i<numMaterials; i++) {
        print(file[IntegrationData::MOMENTUM], "%16.8e  %16.8e  %16.8e  ",
              momentum[i][0], momentum[i][1], momentum[i][2]);
        sum += momentum[i];
      }
      print(file[IntegrationData::MOMENTUM], "%16.8e  %16.8e  %16.8e\n", sum[0], sum[1], sum[2]);
      print_flush(file[IntegrationData::MOMENTUM]);
    } 

    if(file[IntegrationData::TOTAL_ENERGY]) {
      print(file[IntegrationData::TOTAL_ENERGY], "%10d    %16.8e    ", time_step, time);
      double E[numMaterials];
      double sum = 0.0;
      IntegrateTotalEnergy(index, tag, coords, dxyz, cell, v, id, E);
      for(int i=0; i<numMaterials; i++) {
        print(file[IntegrationData::TOTAL_ENERGY], "%16.8e  ", E[i]);
        sum += E[i];
      }
      print(file[IntegrationData::TOTAL_ENERGY], "%16.8e\n", sum);
      print_flush(file[IntegrationData::TOTAL_ENERGY]);
    }

    if(file[IntegrationData::TOTAL_ENTHALPY]) {
      print(file[IntegrationData::TOTAL_ENTHALPY], "%10d    %16.8e    ", time_step, time);
      double H[numMaterials];
      double sum = 0.0;
      IntegrateTotalEnthalpy(index, tag, coords, dxyz, cell, v, id, H);
      for(int i=0; i<numMaterials; i++) {
        print(file[IntegrationData::TOTAL_ENTHALPY], "%16.8e  ", H[i]);
        sum += H[i];
      }
      print(file[IntegrationData::TOTAL_ENTHALPY], "%16.8e\n", sum);
      print_flush(file[IntegrationData::TOTAL_ENTHALPY]);
    }

    if(file[IntegrationData::KINETIC_ENERGY]) {
      print(file[IntegrationData::KINETIC_ENERGY], "%10d    %16.8e    ", time_step, time);
      double kinetic[numMaterials];
      double sum = 0.0;
      IntegrateKineticEnergy(index, tag, coords, dxyz, cell, v, id, kinetic);
      for(int i=0; i<numMaterials; i++) {
        print(file[IntegrationData::KINETIC_ENERGY], "%16.8e  ", kinetic[i]);
        sum += kinetic[i];
      }
      print(file[IntegrationData::KINETIC_ENERGY], "%16.8e\n", sum);
      print_flush(file[IntegrationData::KINETIC_ENERGY]);
    }

    if(file[IntegrationData::INTERNAL_ENERGY]) {
      print(file[IntegrationData::INTERNAL_ENERGY], "%10d    %16.8e    ", time_step, time);
      double internal[numMaterials];
      double sum = 0.0;
      IntegrateInternalEnergy(index, tag, coords, dxyz, cell, v, id, internal);
      for(int i=0; i<numMaterials; i++) {
        print(file[IntegrationData::INTERNAL_ENERGY], "%16.8e  ", internal[i]);
        sum += internal[i];
      }
      print(file[IntegrationData::INTERNAL_ENERGY], "%16.8e\n", sum);
      print_flush(file[IntegrationData::INTERNAL_ENERGY]);
    }

    if(file[IntegrationData::POTENTIAL_ENERGY]) {
      print(file[IntegrationData::POTENTIAL_ENERGY], "%10d    %16.8e    ", time_step, time);
      double potential[numMaterials];
      double sum = 0.0;
      IntegratePotentialEnergy(index, tag, coords, dxyz, cell, v, id, potential);
      for(int i=0; i<numMaterials; i++) {
        print(file[IntegrationData::POTENTIAL_ENERGY], "%16.8e  ", potential[i]);
        sum += potential[i];
      }
      print(file[IntegrationData::POTENTIAL_ENERGY], "%16.8e\n", sum);
      print_flush(file[IntegrationData::POTENTIAL_ENERGY]);
    }

    if(file[IntegrationData::LASER_RADIATION]) {
      print(file[IntegrationData::LASER_RADIATION], "%10d    %16.8e    ", time_step, time);
      assert(l);
      double radiation[numMaterials];
      double sum = 0.0;
      IntegrateLaserRadiation(index, tag, coords, dxyz, cell, v, id, l, radiation);
      for(int i=0; i<numMaterials; i++) {
        print(file[IntegrationData::LASER_RADIATION], "%16.8e  ", radiation[i]);
        sum += radiation[i];
      }
      print(file[IntegrationData::LASER_RADIATION], "%16.8e\n", sum);
      print_flush(file[IntegrationData::LASER_RADIATION]);
    }

    last_snapshot_time[index] = time;

    index++;
  }    

  Tag.RestoreDataPointerToLocalVector();
  V.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();
  cell_volume.RestoreDataPointerToLocalVector();
  if(coords) coordinates.RestoreDataPointerToLocalVector();
  if(dxyz) delta_xyz.RestoreDataPointerToLocalVector();
  if(l) L->RestoreDataPointerToLocalVector();
}   
    
//--------------------------------------------------------------------------------------------------------------

void
IntegrationOutput::IntegrateVolume(int index, double*** tag, Vec3D*** coords, Vec3D*** dxyz, double*** cell,
                                        double*** id, double* volume)
{
  //initialize volumes to 0
  for(int i=0; i<numMaterials; i++)
      volume[i] = 0.0;
 
  // loop through the true domain (excluding the ghost layer)
  double PI = acos(0.0)*2.0;
  [[maybe_unused]] double scalar;
  int myid;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        myid = id[k][j][i];
        if(myid<0 || myid>=numMaterials) {
          fprintf(stderr,"*** Error: Detected an unrecognized material id (%d)\n",
                  myid);
          exit(-1);
        }
        if(MathTools::GetBit((int)tag[k][j][i], index)) {
          if(mesh_type == MeshData::SPHERICAL){
            scalar = PI*4.0*coords[k][j][i][0]*coords[k][j][i][0]/dxyz[k][j][i][2]/dxyz[k][j][i][1];
            volume[myid] += cell[k][j][i]*scalar;
          }
          else if(mesh_type == MeshData::CYLINDRICAL){
            scalar = PI*2.0*coords[k][j][i][1]/dxyz[k][j][i][2];
            volume[myid] += cell[k][j][i]*scalar;
          }
          else
            volume[myid] += cell[k][j][i];
        }
      }

  MPI_Allreduce(MPI_IN_PLACE, volume, numMaterials, MPI_DOUBLE, MPI_SUM, comm);
}

//--------------------------------------------------------------------------------------------------------------

void
IntegrationOutput::IntegrateMass(int index, double*** tag, Vec3D*** coords, Vec3D*** dxyz, double*** cell,
                                 Vec5D*** v, double*** id, double* mass)
{

  for(int i=0; i<numMaterials; i++)
      mass[i] = 0.0;

  double PI = acos(0.0)*2.0;
  [[maybe_unused]] double scalar;
  int myid;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        myid = id[k][j][i];
        if(myid<0 || myid>=numMaterials) {
          fprintf(stderr,"*** Error: Detected an unrecognized material id (%d)\n",
                  myid);
          exit(-1);
        }
        if(MathTools::GetBit((int)tag[k][j][i], index)) {
          if(mesh_type == MeshData::SPHERICAL){
            scalar = PI*4.0*coords[k][j][i][0]*coords[k][j][i][0]/dxyz[k][j][i][2]/dxyz[k][j][i][1];
            mass[myid] += v[k][j][i][0]*cell[k][j][i]*scalar;
          }
          else if(mesh_type == MeshData::CYLINDRICAL){
            scalar = PI*2.0*coords[k][j][i][1]/dxyz[k][j][i][2];
            mass[myid] += v[k][j][i][0]*cell[k][j][i]*scalar;
          }
          else
             mass[myid] += v[k][j][i][0]*cell[k][j][i];
        }
      }

  MPI_Allreduce(MPI_IN_PLACE, mass, numMaterials, MPI_DOUBLE, MPI_SUM, comm);
}

//----------------------------------------------------------------------------------------------------------------
void
IntegrationOutput::IntegrateMomentum(int index, double*** tag, Vec3D*** coords, Vec3D*** dxyz, double*** cell,
                                     Vec5D*** v, double*** id, Vec3D* momentum)
{

  for(int i=0; i<numMaterials; i++)
      momentum[i] = 0.0;

  double PI = acos(0.0)*2.0;
  [[maybe_unused]] double scalar;
  Vec3D velocity;
  int myid;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        myid = id[k][j][i];
        if(myid<0 || myid>=numMaterials) {
          fprintf(stderr,"*** Error: Detected an unrecognized material id (%d)\n",
                  myid);
          exit(-1);
        }
        if(MathTools::GetBit((int)tag[k][j][i], index)) {
          velocity = Vec3D(v[k][j][i][1],v[k][j][i][2],v[k][j][i][3]);
          if(mesh_type == MeshData::SPHERICAL){
            scalar = PI*4.0*coords[k][j][i][0]*coords[k][j][i][0]/dxyz[k][j][i][2]/dxyz[k][j][i][1];
            momentum[myid] += v[k][j][i][0]*cell[k][j][i]*scalar*velocity;
          }
          else if(mesh_type == MeshData::CYLINDRICAL){
            scalar = PI*2.0*coords[k][j][i][1]/dxyz[k][j][i][2];
            momentum[myid] += v[k][j][i][0]*cell[k][j][i]*scalar*velocity;
          }
          else
             momentum[myid] += v[k][j][i][0]*cell[k][j][i]*velocity;
        }
      }

  MPI_Allreduce(MPI_IN_PLACE, (double*)momentum, 3*numMaterials, MPI_DOUBLE, MPI_SUM, comm);

}

//----------------------------------------------------------------------------------------------------------------

void
IntegrationOutput::IntegrateTotalEnergy(int index, double*** tag, Vec3D*** coords, Vec3D*** dxyz,
                                        double*** cell, Vec5D*** v, double*** id, double* E)
{

  for(int i=0; i<numMaterials; i++)
    E[i] = 0.0;

  double PI = acos(0.0)*2.0;
  [[maybe_unused]] double scalar;
  int myid;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        myid = id[k][j][i];
        if(myid<0 || myid>=numMaterials) {
          fprintf(stderr,"*** Error: Detected an unrecognized material id (%d)\n",
                  myid);
          exit(-1);
        }
        if(MathTools::GetBit((int)tag[k][j][i], index)) {
          double vsquare = v[k][j][i][1]*v[k][j][i][1] + v[k][j][i][2]*v[k][j][i][2]
                         + v[k][j][i][3]*v[k][j][i][3];
          double e = vf[myid]->GetInternalEnergyPerUnitMass(v[k][j][i][0], v[k][j][i][4]);
          if(mesh_type == MeshData::SPHERICAL){
            scalar = PI*4.0*coords[k][j][i][0]*coords[k][j][i][0]/dxyz[k][j][i][2]/dxyz[k][j][i][1];
            E[myid] += v[k][j][i][0]*(0.5*vsquare + e)*cell[k][j][i]*scalar;
          }
          else if(mesh_type == MeshData::CYLINDRICAL){
            scalar = PI*2.0*coords[k][j][i][1]/dxyz[k][j][i][2];
            E[myid] += v[k][j][i][0]*(0.5*vsquare + e)*cell[k][j][i]*scalar;
          }
          else
            E[myid] += v[k][j][i][0]*(0.5*vsquare + e)*cell[k][j][i];
        }
      }

  MPI_Allreduce(MPI_IN_PLACE, E, numMaterials, MPI_DOUBLE, MPI_SUM, comm);

}

//-----------------------------------------------------------------------------------------------------------------

void
IntegrationOutput::IntegrateTotalEnthalpy(int index, double*** tag, Vec3D*** coords, Vec3D*** dxyz,
                                          double*** cell, Vec5D*** v, double*** id, double* H)
{

  for(int i=0; i<numMaterials; i++)
    H[i] = 0.0;

  double PI = acos(0.0)*2.0;
  [[maybe_unused]] double scalar;
  int myid;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        myid = id[k][j][i];
        if(myid<0 || myid>=numMaterials) {
          fprintf(stderr,"*** Error: Detected an unrecognized material id (%d)\n",
                  myid);
          exit(-1);
        }
        if(MathTools::GetBit((int)tag[k][j][i], index)) {
          double vsquare = v[k][j][i][1]*v[k][j][i][1] + v[k][j][i][2]*v[k][j][i][2]
                         + v[k][j][i][3]*v[k][j][i][3];
          double e = vf[myid]->GetInternalEnergyPerUnitMass(v[k][j][i][0], v[k][j][i][4]);
          if(mesh_type == MeshData::SPHERICAL){
            scalar = PI*4.0*coords[k][j][i][0]*coords[k][j][i][0]/dxyz[k][j][i][2]/dxyz[k][j][i][1];
            H[myid] += (v[k][j][i][0]*(0.5*vsquare+e) + v[k][j][i][4])*cell[k][j][i]*scalar;
          }
          else if(mesh_type == MeshData::CYLINDRICAL){
            scalar = PI*2.0*coords[k][j][i][1]/dxyz[k][j][i][2];
            H[myid] += (v[k][j][i][0]*(0.5*vsquare+e) + v[k][j][i][4])*cell[k][j][i]*scalar;
          }
          else
            H[myid] += (v[k][j][i][0]*(0.5*vsquare+e) + v[k][j][i][4])*cell[k][j][i];
        }
      }
 
  MPI_Allreduce(MPI_IN_PLACE, H, numMaterials, MPI_DOUBLE, MPI_SUM, comm);
 
}

//-----------------------------------------------------------------------------------------------------------------

void
IntegrationOutput::IntegrateKineticEnergy(int index, double*** tag, Vec3D*** coords, Vec3D*** dxyz,
                                          double*** cell, Vec5D*** v, double*** id, double* kinetic)
{

  for(int i=0; i<numMaterials; i++)
    kinetic[i] = 0.0;

  double PI = acos(0.0)*2.0;
  [[maybe_unused]] double scalar;
  int myid;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        myid = id[k][j][i];
        if(myid<0 || myid>=numMaterials) {
          fprintf(stderr,"*** Error: Detected an unrecognized material id (%d)\n",
                  myid);
          exit(-1);
        }
        if(MathTools::GetBit((int)tag[k][j][i], index)) {
          double vsquare = v[k][j][i][1]*v[k][j][i][1] + v[k][j][i][2]*v[k][j][i][2]
                         + v[k][j][i][3]*v[k][j][i][3];
          if(mesh_type == MeshData::SPHERICAL){
            scalar = PI*4.0*coords[k][j][i][0]*coords[k][j][i][0]/dxyz[k][j][i][2]/dxyz[k][j][i][1];
            kinetic[myid] += 0.5*v[k][j][i][0]*vsquare*cell[k][j][i]*scalar;
          }
          else if(mesh_type == MeshData::CYLINDRICAL){
            scalar = PI*2.0*coords[k][j][i][1]/dxyz[k][j][i][2];
            kinetic[myid] += 0.5*v[k][j][i][0]*vsquare*cell[k][j][i]*scalar;
          }
          else
            kinetic[myid] += 0.5*v[k][j][i][0]*vsquare*cell[k][j][i];
        }
      }

  MPI_Allreduce(MPI_IN_PLACE, kinetic, numMaterials, MPI_DOUBLE, MPI_SUM, comm);

}

//---------------------------------------------------------------------------------------------------------

void
IntegrationOutput::IntegrateInternalEnergy(int index, double*** tag, Vec3D*** coords, Vec3D*** dxyz,
                                           double*** cell, Vec5D*** v, double*** id, double* internal)
{

  for(int i=0; i<numMaterials; i++)
    internal[i] = 0.0;

  double PI = acos(0.0)*2.0;
  [[maybe_unused]] double scalar;
  int myid;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        myid = id[k][j][i];
        if(myid<0 || myid>=numMaterials) {
          fprintf(stderr,"*** Error: Detected an unrecognized material id (%d)\n",
                  myid);
          exit(-1);
        }
        if(MathTools::GetBit((int)tag[k][j][i], index)) {
          double e = vf[myid]->GetInternalEnergyPerUnitMass(v[k][j][i][0], v[k][j][i][4]);
          if(mesh_type == MeshData::SPHERICAL){
            scalar = PI*4.0*coords[k][j][i][0]*coords[k][j][i][0]/dxyz[k][j][i][2]/dxyz[k][j][i][1];
            internal[myid] += v[k][j][i][0]*e*cell[k][j][i]*scalar;
          }
          else if(mesh_type == MeshData::CYLINDRICAL){
            scalar = PI*2.0*coords[k][j][i][1]/dxyz[k][j][i][2];
            internal[myid] += v[k][j][i][0]*e*cell[k][j][i]*scalar;
          }
          else
            internal[myid] += v[k][j][i][0]*e*cell[k][j][i];
        }
      }
 
  MPI_Allreduce(MPI_IN_PLACE, internal, numMaterials, MPI_DOUBLE, MPI_SUM, comm);

}
    
//---------------------------------------------------------------------------------------------------------
 
void
IntegrationOutput::IntegratePotentialEnergy(int index, double*** tag, Vec3D*** coords, Vec3D*** dxyz,
                                            double*** cell, Vec5D*** v, double*** id, double* potential)
{

  for(int i=0; i<numMaterials; i++)
    potential[i] = 0.0;

  double PI = acos(0.0)*2.0;
  [[maybe_unused]] double scalar;
  int myid;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        myid = id[k][j][i];
        if(myid<0 || myid>=numMaterials) {
          fprintf(stderr,"*** Error: Detected an unrecognized material id (%d)\n",
                  myid);
          exit(-1);
        }
        if(MathTools::GetBit((int)tag[k][j][i], index)) {
          if(mesh_type == MeshData::SPHERICAL){
            scalar = PI*4.0*coords[k][j][i][0]*coords[k][j][i][0]/dxyz[k][j][i][2]/dxyz[k][j][i][1];
            potential[myid] += v[k][j][i][4]*cell[k][j][i]*scalar;
          }
          else if(mesh_type == MeshData::CYLINDRICAL){
            scalar = PI*2.0*coords[k][j][i][1]/dxyz[k][j][i][2];
            potential[myid] += v[k][j][i][4]*cell[k][j][i]*scalar;
          }
          else
            potential[myid] += v[k][j][i][4]*cell[k][j][i];
        }
      }

  MPI_Allreduce(MPI_IN_PLACE, potential, numMaterials, MPI_DOUBLE, MPI_SUM, comm);

}

//---------------------------------------------------------------------------------------------------------

void
IntegrationOutput::IntegrateLaserRadiation(int index, double*** tag, Vec3D*** coords, Vec3D*** dxyz,
                                           double*** cell, Vec5D*** v, double*** id, double*** l, double* radiation)
{
  assert(l);
  assert(laser);

  for(int i=0; i<numMaterials; i++)
    radiation[i] = 0.0;

  double PI = acos(0.0)*2.0;
  [[maybe_unused]] double scalar;
  int myid;
  double e, myT, eta;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        myid = id[k][j][i];
        if(myid<0 || myid>=numMaterials) {
          fprintf(stderr,"*** Error: Detected an unrecognized material id (%d)\n",
                  myid);
          exit(-1);
        }
        if(MathTools::GetBit((int)tag[k][j][i], index)) {
          e   = vf[myid]->GetInternalEnergyPerUnitMass(v[k][j][i][0], v[k][j][i][4]);
          myT = vf[myid]->GetTemperature(v[k][j][i][0], e);
          eta = laser->GetAbsorptionCoefficient(myT, myid);
          if(mesh_type == MeshData::SPHERICAL){
            scalar = PI*4.0*coords[k][j][i][0]*coords[k][j][i][0]/dxyz[k][j][i][2]/dxyz[k][j][i][1];
            radiation[myid] += eta*l[k][j][i]*cell[k][j][i]*scalar;
          }
          else if(mesh_type == MeshData::CYLINDRICAL){
            scalar = PI*2.0*coords[k][j][i][1]/dxyz[k][j][i][2];
            radiation[myid] += eta*l[k][j][i]*cell[k][j][i]*scalar;
          }
          else
            radiation[myid] += eta*l[k][j][i]*cell[k][j][i];
        }
      }

  MPI_Allreduce(MPI_IN_PLACE, radiation, numMaterials, MPI_DOUBLE, MPI_SUM, comm);

}
