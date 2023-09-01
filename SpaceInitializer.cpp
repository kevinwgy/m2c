/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <SpaceInitializer.h>
#include <SpaceOperator.h>
#include <LevelSetOperator.h>
#include <EmbeddedBoundaryDataSet.h>

using std::vector;

//---------------------------------------------------------------------

SpaceInitializer::SpaceInitializer(MPI_Comm &comm_, IoData &iod_, GlobalMeshInfo &global_mesh_)
                : comm(comm_), iod(iod_), global_mesh(global_mesh_)
{ }

//---------------------------------------------------------------------

SpaceInitializer::~SpaceInitialier()
{ }

//---------------------------------------------------------------------

void
SpaceInitializer::Destroy()
{ }

//---------------------------------------------------------------------

std::multimap<int, std::pair<int,int> >
SpaceInitializer::SetInitialCondition(SpaceVariable3D &V, SpaceVariable3D &ID, vector<SpaceVariable3D*> &Phi,
                                      SpaceOperator &spo, vector<SpaceVariable3D*> &lso,
                                      std::unique_ptr<vector<std::unique_ptr<EmbeddedBoundaryDataSet> > > EBDS)
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
  vector<std::pair<int,int> > order;  //<geom type, geom dataMap index>
  int nGeom = OrderUserSpecifiedGeometries(order);

  // Step 2: Apply initial conditions to V, ID
  I AM HERE

  // Step 3:: Apply initial conditions to Phi
 
}

//---------------------------------------------------------------------

int
SpaceInitializer::OrderUserSpecifiedGeometries(vector<std::pair<int,int> > &order)
{
  MultiInitialConditionsData &ic(iod.ic.multiInitialConditions);
  int nGeom = ic.planeMap.dataMap.size() //0
            + ic.cylinderconeMap.dataMap.size() //1
            + ic.cylindersphereMap.dataMap.size() //2
            + ic.sphereMap.dataMap.size() //3
            + ic.parallelepipedMap.dataMap.size() //4
            + ic.spheroidMap.dataMap.size() //5
            + ic.enclosureMap.dataMap.size() //6
            + ic.pointMap.dataMap.size(); //7
  order.assign(nGeom, std::make_pair(-1,-1));
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


//---------------------------------------------------------------------





















//---------------------------------------------------------------------

