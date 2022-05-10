#include<EmbeddedBoundaryOperator.h>
#include<deque>
#include<list>
#include<utility>
#include<Utils.h>
#include<memory.h> //unique_ptr
#include<dlfcn.h> //dlopen, dlclose

using std::string;
using std::map;
using std::vector;
using std::unique_ptr;

//------------------------------------------------------------------------------------------------

EmbeddedBoundaryOperator::EmbeddedBoundaryOperator(MPI_Comm &comm_, IoData &iod_, bool surface_from_other_solver) 
                        : comm(comm_), iod(iod_), hasSurfFromOtherSolver(surface_from_other_solver),
                          dms_ptr(NULL), coordinates_ptr(NULL), ghost_nodes_inner_ptr(NULL),
                          ghost_nodes_outer_ptr(NULL), x_glob_ptr(NULL), y_glob_ptr(NULL), z_glob_ptr(NULL),
                          dx_glob_ptr(NULL), dy_glob_ptr(NULL), dz_glob_ptr(NULL)
{
  // count surfaces from files
  int counter = iod.ebm.embed_surfaces.surfaces.dataMap.size();

  // get the correct size of surfaces and F
  if(surface_from_other_solver) {
    surfaces.resize(counter+1);
    F.resize(counter+1);
  } else {
    surfaces.resize(counter);
    F.resize(counter);
  }

  surfaces.resize(surfaces.size());
  F_prev.resize(F.size());
 
  // set default boundary type to "None"
  surface_type.resize(surfaces.size(), EmbeddedSurfaceData::None);
  iod_embedded_surfaces.resize(surfaces.size(), NULL);

  // read surfaces from files
  for(auto it = iod.ebm.embed_surfaces.surfaces.dataMap.begin();
           it != iod.ebm.embed_surfaces.surfaces.dataMap.end(); it++) {

    int index = it->first;
    if(index<0 || index>=counter) {
      print_error("*** Error: Detected error in the indices of embedded surfaces (id = %d).", index);
      exit_mpi();
    }

    iod_embedded_surfaces[index] = it->second;

    if(index==0) {
      if(surface_from_other_solver) {
        if(it->second->surface_provided_by_other_solver != EmbeddedSurfaceData::YES) {
          print_error("*** Error: Conflict input about EmbeddedSurface[%d]. Should mesh be provided by another solver?", 
                      index); 
          exit_mpi();
        } 
        continue; //no file to read
      } else {
        if(it->second->surface_provided_by_other_solver != EmbeddedSurfaceData::NO) {
          print_error("*** Error: Conflict input about EmbeddedSurface[%d]. Should mesh be provided by user?", index); 
          exit_mpi();
        }
      }
    } else {
      if(it->second->surface_provided_by_other_solver != EmbeddedSurfaceData::NO) {
        print_error("*** Error: Currently, only one embedded surface (with id 0) can be provided by another solver.");
        exit_mpi();
      }
    }

    surface_type[index] = it->second->type;

    ReadMeshFile(it->second->filename, surface_type[index], surfaces[index].X, surfaces[index].elems);

    surfaces[index].X0 = surfaces[index].X;

    surfaces[index].BuildConnectivities();
    surfaces[index].CalculateNormalsAndAreas();
/*
    bool orientation = surfaces[index].CheckSurfaceOrientation(); 
    bool closed = surfaces[index].CheckSurfaceClosedness();
*/
  }

  print("- Activated the Embedded Boundary Method. Detected %d surface(s) (%d from concurrent program(s)).\n",
        surfaces.size(), surfaces.size() - counter); 

  // set NULL to intersector pointers
  intersector.resize(surfaces.size(), NULL);

  // setup dynamics_calculator
  SetupUserDefinedDynamicsCalculator();
}

//------------------------------------------------------------------------------------------------

EmbeddedBoundaryOperator::~EmbeddedBoundaryOperator()
{
  for(auto it = dynamics_calculator.begin(); it != dynamics_calculator.end(); it++) {
    if(std::get<0>(*it)) {
      assert(std::get<2>(*it)); //this is the destruction function
      (std::get<2>(*it))(std::get<0>(*it)); //use the destruction function to destroy the calculator
      assert(std::get<1>(*it));
      dlclose(std::get<1>(*it));
    }
  }
}

//------------------------------------------------------------------------------------------------

void
EmbeddedBoundaryOperator::Destroy()
{
  for(int i=0; i<intersector.size(); i++)
    if(intersector[i])
      intersector[i]->Destroy();
}

//------------------------------------------------------------------------------------------------

void
EmbeddedBoundaryOperator::SetCommAndMeshInfo(DataManagers3D &dms_, SpaceVariable3D &coordinates_,
                              vector<GhostPoint> &ghost_nodes_inner_, vector<GhostPoint> &ghost_nodes_outer_,
                              vector<double> &x_, vector<double> &y_, vector<double> &z_,
                              vector<double> &dx_, vector<double> &dy_, vector<double> &dz_)
{
  dms_ptr               = &dms_;
  coordinates_ptr       = &coordinates_;
  ghost_nodes_inner_ptr = &ghost_nodes_inner_; 
  ghost_nodes_outer_ptr = &ghost_nodes_outer_; 
  x_glob_ptr            = &x_;
  y_glob_ptr            = &y_;
  z_glob_ptr            = &z_;
  dx_glob_ptr           = &dx_;
  dy_glob_ptr           = &dy_;
  dz_glob_ptr           = &dz_;
}

//------------------------------------------------------------------------------------------------

void
EmbeddedBoundaryOperator::SetupIntersectors()
{
  for(int i=0; i<intersector.size(); i++) {
    intersector[i] = new Intersector(comm, *dms_ptr, *iod_embedded_surfaces[i], surfaces[i],
                                     *coordinates_ptr, *ghost_nodes_inner_ptr, *ghost_nodes_outer_ptr,
                                     *x_glob_ptr, *y_glob_ptr, *z_glob_ptr,
                                     *dx_glob_ptr, *dy_glob_ptr, *dz_glob_ptr);
  }
}

//------------------------------------------------------------------------------------------------

void
EmbeddedBoundaryOperator::StoreID2Closure(std::map<int, std::pair<int,int> > &id2closure_)
{
  id2color = id2closure_;

  inactive_colors.clear();
  for(int i=0; i<surfaces.size(); i++) {
    unique_ptr<EmbeddedBoundaryDataSet> EBDS = GetPointerToIntersectoResultsOnSurface(surf);
    int nRegions = EBDS->nRegions; //this is the number of *closures*. Colors -1, -2, ...
    for(int color = -1; color>=-nRegions; color--) {
      bool found = false;
      for(auto it = id2color.begin(); it != id2color.end(); it++) {
        if(it->second->first == i && it->second->second == color) {
          found = true;
          break;
        }
      }
      if(!found)
        inactive_colors.push_back(std::make_pair(i, color));
    } 
  }
}

//------------------------------------------------------------------------------------------------

void
EmbeddedBoundaryOperator::ReadMeshFile(const char *filename, EmbeddedSurfaceData::Type& surface_type, 
                                       vector<Vec3D> &Xs, vector<Int3> &Es)
{

  // read data from the surface input file.
  FILE *topFile;
  topFile = fopen(filename, "r");
  if(topFile == NULL) {
    print_error("*** Error: embedded structure surface mesh doesn't exist (%s).\n", filename);
    exit_mpi();
  }
 
  int MAXLINE = 500;
  char line[MAXLINE], key1[MAXLINE], key2[MAXLINE], copyForType[MAXLINE];

  int num0 = 0;
  int num1 = 0;
  double x1, x2, x3;
  int node1, node2, node3;
  int type_reading = 0; //1 means reading node set; 2 means reading element set
  std::deque<std::pair<int, Vec3D>> nodeList;
  std::deque<std::array<int, 4>> elemList; // element ID + three node IDs
  int maxNode = 0, maxElem = 0;
  bool found_nodes = false;
  bool found_elems = false;


  // --------------------
  // Read the file
  // --------------------
  while(fgets(line, MAXLINE, topFile) != 0) {

    sscanf(line, "%s", key1);
    string key1_string(key1);

    if(key1[0] == '#') {
      //Do nothing. This is user's comment
    }
    else if(same_strings_insensitive(key1_string,"Nodes")){
      if(found_nodes) {//already found nodes... This is a conflict
        print_error("*** Error: Found multiple sets of nodes (keyword 'Nodes') in %s.\n", filename);
        exit_mpi();
      }
      sscanf(line, "%*s %s", key2);
      type_reading = 1;
      found_nodes = true;
    }
    else if(same_strings_insensitive(key1_string, "Elements")) {

      if(found_elems) {//already found elements... This is a conflict
        print_error("*** Error: Found multiple sets of elements (keyword 'Elements') in %s.\n", filename);
        exit_mpi();
      }
      sscanf(line, "%*s %s", key2);
      type_reading = 2;
      found_elems = true;
/*
      // now we look for keywords for the type of structure
      strcpy(copyForType, key2);
      int l = 0;
      while((copyForType[l] != '\0') && (l < MAXLINE)) {
        copyForType[l] = (char)std::tolower(static_cast<unsigned char>(copyForType[l]));
        l++;
      }
      // read the name of the file and detects keyword for type
      if(strstr(copyForType, "symmetry") != NULL)
        surface_type = EmbeddedSurfaceData::Symmetry;
      else if(strstr(copyForType, "porouswall") != NULL)
        surface_type = EmbeddedSurfaceData::PorousWall;
      else if(strstr(copyForType, "wall") != NULL) //this check has to be placed after "porouswall"
        surface_type = EmbeddedSurfaceData::Wall;
      else if(strstr(copyForType, "none") != NULL) 
        surface_type = EmbeddedSurfaceData::None;
      else {
        print_error("*** Error: Detected unsupported surface type in %s (%s).\n", filename, key2);
        exit_mpi();
      } 
*/
    }
    else if(type_reading == 1) { //reading a node (following "Nodes Blabla")
      int count = sscanf(line, "%d %lf %lf %lf", &num1, &x1, &x2, &x3);
      if(count != 4) {
        print_error("*** Error: Cannot interpret line %s (in %s). Expecting a node.\n", line, filename);
        exit_mpi();
      }
      if(num1 < 1) {
        print_error("*** Error: detected a node with index %d in embedded surface file %s.\n", num1, filename);
        exit_mpi();
      }
      if(num1 > maxNode)
        maxNode = num1;

      nodeList.push_back({num1, {x1, x2, x3}});
    }
    else if(type_reading == 2) { // we are reading an element --- HAS TO BE A TRIANGLE!
      int count = sscanf(line, "%d %d %d %d %d", &num0, &num1, &node1, &node2, &node3);
      if(count != 5) {
        print_error("*** Error: Cannot interpret line %s (in %s). Expecting a triangular element.\n", line, filename);
        exit_mpi();
      }
      if(num0 < 1) {
        print_error("*** Error: detected an element with index %d in embedded surface file %s.\n", num0, filename);
        exit_mpi();
      }
      if(num0 > maxElem)
        maxElem = num0;

      elemList.push_back({num0, node1, node2, node3});
    }
    else { // found something I cannot understand...
      print_error("*** Error: Unable to interpret line %s (in %s).\n", line, filename);
      exit_mpi();
    }

  }

  fclose(topFile);

  if(!found_nodes) {
    print_error("*** Error: Unable to find node set in %s.\n", filename);
    exit_mpi();
  }
  if(!found_elems) {
    print_error("*** Error: Unable to find element set in %s.\n", filename);
    exit_mpi();
  }

  // ----------------------------
  // Now, check and store nodes
  // ----------------------------
  int nNodes = nodeList.size();
  map<int,int> old2new;
  Xs.resize(nNodes);
  int id(-1);
  if(nNodes != maxNode) { // need to renumber nodes, i.e. create "old2new"
    print_warning("Warning: The node indices of an embedded surface may have a gap: "
                  "max index = %d, number of nodes = %d. Renumbering nodes. (%s)\n",
                  maxNode, nNodes, filename);
//    assert(nNodes < maxNode);

    int current_id = 0; 
    vector<bool> nodecheck(maxNode+1, false);
    for(auto it1 = nodeList.begin(); it1 != nodeList.end(); it1++) {
      id = it1->first;
      if(nodecheck[id]) {
        print_error("*** Error: Found duplicate node (id: %d) in embedded surface file %s.\n", id, filename);
        exit(-1);
      }
      nodecheck[id] = true;
      Xs[current_id] = it1->second; 
      old2new[id] = current_id;
      current_id++;
    }
    assert(current_id==Xs.size());
  } 
  else { //in good shape
    vector<bool> nodecheck(nNodes, false);
    for(auto it1 = nodeList.begin(); it1 != nodeList.end(); it1++) {
      id = it1->first - 1; 
      if(nodecheck[id]) {
        print_error("*** Error: Found duplicate node (id: %d) in embedded surface file %s.\n", id+1, filename);
        exit(-1);
      }
      nodecheck[id] = true;
      Xs[it1->first - 1] = it1->second;
    }
  }


  // ------------------------------
  // check nodes used by elements
  // ------------------------------
  for(auto it = elemList.begin(); it != elemList.end(); it++) {

    id = (*it)[0];
    node1 = (*it)[1];
    node2 = (*it)[2];
    node3 = (*it)[3];
      
    if(old2new.empty()) {//node set is original order

      if(node1<=0 || node1 > nNodes) {
        print_error("*** Error: Detected unknown node number (%d) in element %d (%s).\n", node1, id, filename);
        exit_mpi();
      }

      if(node2<=0 || node2 > nNodes) {
        print_error("*** Error: Detected unknown node number (%d) in element %d (%s).\n", node2, id, filename);
        exit_mpi();
      }

      if(node3<=0 || node3 > nNodes) {
        print_error("*** Error: Detected unknown node number (%d) in element %d (%s).\n", node3, id, filename);
        exit_mpi();
      }
    }
    else {// nodes are renumbered

      auto p1 = old2new.find(node1);
      if(p1 == old2new.end()) { 
        print_error("*** Error: Detected unknown node number (%d) in element %d (%s).\n", node1, id, filename);
        exit_mpi();
      }

      auto p2 = old2new.find(node2);
      if(p2 == old2new.end()) { 
        print_error("*** Error: Detected unknown node number (%d) in element %d (%s).\n", node2, id, filename);
        exit_mpi();
      }

      auto p3 = old2new.find(node3);
      if(p3 == old2new.end()) { 
        print_error("*** Error: Detected unknown node number (%d) in element %d (%s).\n", node3, id, filename);
        exit_mpi();
      }
    }
  }


  // ----------------------------
  // check and store elements
  // ----------------------------
  int nElems = elemList.size();
  Es.resize(nElems);
  if(nElems != maxElem) { // need to renumber elements.
    print_warning("Warning: The element indices of an embedded surface may have a gap: "
                  "max index = %d, number of elements = %d. Renumbering elements. (%s)\n",
                  maxElem, nElems, filename);
//    assert(nElems < maxElem);
    
    int current_id = 0; 
    vector<bool> elemcheck(maxElem+1, false);
    for(auto it = elemList.begin(); it != elemList.end(); it++) {
      id = (*it)[0];
      if(elemcheck[id]) {
        print_error("*** Error: Found duplicate element (id: %d) in embedded surface file %s.\n", id, filename);
        exit_mpi();
      }
      elemcheck[id] = true;

      node1 = (*it)[1];
      node2 = (*it)[2];
      node3 = (*it)[3];
      
      if(old2new.empty()) //node set is original order
        Es[current_id] = Int3(node1-1, node2-1, node3-1);
      else {// nodes are renumbered
        auto p1 = old2new.find(node1);
        auto p2 = old2new.find(node2);
        auto p3 = old2new.find(node3);
        Es[current_id] = Int3(p1->second, p2->second, p3->second);
      }      
      current_id++;
    }
  } 
  else { //element numbers in good shape

    vector<bool> elemcheck(nElems, false);
    for(auto it = elemList.begin(); it != elemList.end(); it++) {
      id = (*it)[0] - 1;
      if(elemcheck[id]) {
        print_error("*** Error: Found duplicate element (id: %d) in embedded surface file %s.\n", id, filename);
        exit_mpi();
      }
      elemcheck[id] = true;

      node1 = (*it)[1];
      node2 = (*it)[2];
      node3 = (*it)[3];
      
      if(old2new.empty()) //node set is original order
        Es[id] = Int3(node1-1, node2-1, node3-1);
      else {// nodes are renumbered
        auto p1 = old2new.find(node1);
        auto p2 = old2new.find(node2);
        auto p3 = old2new.find(node3);
        Es[id] = Int3(p1->second, p2->second, p3->second);
      }
    }
  }

}

//------------------------------------------------------------------------------------------------

void
EmbeddedBoundaryOperator::UpdateSurfacesPrevAndFPrev(bool partial_copy)
{
  //partial_copy means only the nodal coordinates (not topology) gets copied to surfaces_prev
  //
  assert(partial_copy);

  assert(F.size() == surfaces.size());
  assert(surfaces.size() == F.size());
  assert(F.size() == F_prev.size());

  // loop through all the surfaces
  for(int i = 0; i < F.size(); i++) {
    // copy force
    assert(F[i].size() == F_prev[i].size());
    for(int j=0; j<F[i].size(); j++)
      F_prev[i][j] = F[i][j];

    // copy nodal coords
    assert(surfaces[i].X.size() == surfaces_prev[i].X.size());
    for(int j=0; j<surfaces[i].X.size(); j++)
      surfaces_prev[i].X[j] = surfaces[i].X[j];
  }
}

//------------------------------------------------------------------------------------------------

void
EmbeddedBoundaryOperator::ComputeForces(SpaceVariable3D &V, SpaceVariable3D &ID)
{

  Vec5D***  v  = (Vec5D***) V.GetDataPointer();
  double*** id = ID.GetDataPointer();
  
  // loop through all the embedded surfaces
  for(int surf=0; surf<surfaces.size(); surf++) {

    unique_ptr<EmbeddedBoundaryDataSet> EBDS = GetPointerToIntersectoResultsOnSurface(surf);

    I AM HERE



  }

  //TODO: Must assemble force

  V.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();
}

//------------------------------------------------------------------------------------------------

void
EmbeddedBoundaryOperator::TrackSurfaces()
{
  vector<bool> hasInlet(intersector.size(), false);
  vector<bool> hasOutlet(intersector.size(), false);
  vector<bool> hasOccluded(intersector.size(), false);
  vector<int> numRegions(intersector.size(), 0);
  int phi_layers = 3;
  for(int i = 0; i < intersector.size(); i++) {
    bool a, b, c;
    int d;
    intersector[i]->TrackSurfaceFullCourse(a, b, c, d, phi_layers);
    hasInlet[i] = a;
    hasOutlet[i] = b;
    hasOccluded[i] = c;
    numRegions[i] = d;
  }
}

//------------------------------------------------------------------------------------------------

void
EmbeddedBoundaryOperator::TrackUpdatedSurfaces()
{
  int phi_layers = 3;
  for(int i=0; i<intersector.size(); i++) {
    if(iod_embedded_surfaces[i]->surface_provided_by_other_solver == EmbeddedSurfaceData::NO &&
       strcmp(iod_embedded_surfaces[i]->dynamics_calculator, "") == 0)
      continue; //this surface is fixed...
    intersector[i]->RecomputeFullCourse(surfaces_prev[i].X, phi_layers);
  }
}

//------------------------------------------------------------------------------------------------

void
EmbeddedBoundaryOperator::ApplyUserDefinedSurfaceDynamics(double t, double dt)
{
  for(int i=0; i<surfaces.size(); i++) {
    if(strcmp(iod_embedded_surfaces[i]->dynamics_calculator, "") == 0)
      continue; //not specified for this surface
    UserDefinedDynamics *calculator(std::get<0>(dynamics_calculator[i]));
    assert(calculator);

    vector<Vec3D> &Xs(surfaces[i].X);
    vector<Vec3D> &X0(surfaces[i].X0);
    vector<Vec3D> disp(Xs.size(), 0.0);
    calculator->GetUserDefinedDynamics(t, Xs.size(), (double*)Xs.data(), 
                                       (double*)disp.data(), (double*)surfaces[i].Udot.data());
    for(int i=0; i<Xs.size(); i++) {
      for(int j=0; j<3; j++)
        Xs[i][j] = X0[i][j] + disp[i][j];
    }
  }
}

//------------------------------------------------------------------------------------------------

unique_ptr<vector<unique_ptr<EmbeddedBoundaryDataSet> > >
EmbeddedBoundaryOperator::GetPointerToIntersectorResults()
{
  unique_ptr<vector<unique_ptr<EmbeddedBoundaryDataSet> > > p(new vector<unique_ptr<EmbeddedBoundaryDataSet> >());

  for(int i=0; i<intersector.size(); i++) {
    assert(intersector[i]); //should not call this function if surfaces are not "tracked"
    p->push_back(intersector[i]->GetPointerToResults());
  }

  return p;
}

//------------------------------------------------------------------------------------------------

unique_ptr<EmbeddedBoundaryDataSet>
EmbeddedBoundaryOperator::GetPointerToIntersectoResultsOnSurface(int i)
{
  assert(i>=0 && i<intersector.size());
  assert(intersector[i]);
  return intersector[i]->GetPointerToResults();
}

//------------------------------------------------------------------------------------------------

void
EmbeddedBoundaryOperator::SetupUserDefinedDynamicsCalculator()
{
  dynamics_calculator.resize(surfaces.size(), std::make_tuple(nullptr,nullptr,nullptr));
  for(int i=0; i<surfaces.size(); i++) {
    if(strcmp(iod_embedded_surfaces[i]->dynamics_calculator, "") == 0)
      continue; //not specified for this surface
    if(iod_embedded_surfaces[i]->surface_provided_by_other_solver == EmbeddedSurfaceData::YES) {
      print_error("*** Error: Unable to apply user-defined dynamics for Surface %d, which is owned by "
                  "another solver.\n", i);
      exit_mpi();
    }

    // setup the calculator: dynamically load the object
    std::tuple<UserDefinedDynamics*, void*, DestroyUDD*> &me(dynamics_calculator[i]);
    std::get<1>(me) = dlopen(iod_embedded_surfaces[i]->dynamics_calculator, RTLD_NOW);
    if(!std::get<1>(me)) {
      print_error("*** Error: Unable to load object %s.\n", iod_embedded_surfaces[i]->dynamics_calculator);
      exit_mpi();
    }
    dlerror();
    CreateUDD* create = (CreateUDD*) dlsym(std::get<1>(me), "Create");
    const char* dlsym_error = dlerror();
    if(dlsym_error) {
      print_error("*** Error: Unable to find function Create in %s.\n",
                  iod_embedded_surfaces[i]->dynamics_calculator);
      exit_mpi();
    }
    std::get<2>(me) = (DestroyUDD*) dlsym(std::get<1>(me), "Destroy");
    dlsym_error = dlerror();
    if(dlsym_error) {
      print_error("*** Error: Unable to find function Destroy in %s.\n",
                  iod_embedded_surfaces[i]->dynamics_calculator);
      exit_mpi();
    }
     
    //This is the actual calculator
    std::get<0>(me) = create();

    print("- Loaded user-defined dynamics calculator for surface %d from %s.\n",
          i, iod_embedded_surfaces[i]->dynamics_calculator);
  }

}

//------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------









