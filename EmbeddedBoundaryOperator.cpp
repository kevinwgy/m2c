#include<EmbeddedBoundaryOperator.h>
#include<deque>
#include<list>
#include<utility>
#include<Utils.h>
#include<Vector5D.h>
#include<GeoTools.h>
#include<trilinear_interpolation.h>
#include<gauss_quadratures.h>
#include<memory.h> //unique_ptr
#include<dlfcn.h> //dlopen, dlclose

using std::string;
using std::map;
using std::vector;
using std::unique_ptr;

extern int INACTIVE_MATERIAL_ID;

//------------------------------------------------------------------------------------------------

EmbeddedBoundaryOperator::EmbeddedBoundaryOperator(MPI_Comm &comm_, IoData &iod_, bool surface_from_other_solver) 
                        : comm(comm_), iod(iod_), hasSurfFromOtherSolver(surface_from_other_solver),
                          dms_ptr(NULL), coordinates_ptr(NULL), ghost_nodes_inner_ptr(NULL),
                          ghost_nodes_outer_ptr(NULL), global_mesh_ptr(NULL)
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

  surfaces_prev.resize(surfaces.size());
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
        if(it->second->provided_by_another_solver != EmbeddedSurfaceData::YES) {
          print_error("*** Error: Conflict input about EmbeddedSurface[%d]. Should mesh be provided by another solver?", 
                      index); 
          exit_mpi();
        } 
        continue; //no file to read
      } else {
        if(it->second->provided_by_another_solver != EmbeddedSurfaceData::NO) {
          print_error("*** Error: Conflict input about EmbeddedSurface[%d]. Should mesh be provided by user?", index); 
          exit_mpi();
        }
      }
    } else {
      if(it->second->provided_by_another_solver != EmbeddedSurfaceData::NO) {
        print_error("*** Error: Currently, only one embedded surface (with id 0) can be provided by another solver.");
        exit_mpi();
      }
    }

    surface_type[index] = it->second->type;

    ReadMeshFile(it->second->filename, surface_type[index], surfaces[index].X, surfaces[index].elems);

    surfaces[index].X0 = surfaces[index].X;

    // initialize the velocity vector (to 0)
    surfaces[index].Udot.resize(surfaces[index].X.size(), 0.0);

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

  // setup output
  for(int i=0; i<surfaces.size(); i++)
    lagout.push_back(LagrangianOutput(comm, iod_embedded_surfaces[i]->output));

  // setup dynamics_calculator
  SetupUserDefinedDynamicsCalculator();
}

//------------------------------------------------------------------------------------------------

EmbeddedBoundaryOperator::~EmbeddedBoundaryOperator()
{

  for(int i=0; i<intersector.size(); i++)
    if(intersector[i])
      delete intersector[i];

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
                              GlobalMeshInfo &global_mesh_)
{
  dms_ptr               = &dms_;
  coordinates_ptr       = &coordinates_;
  ghost_nodes_inner_ptr = &ghost_nodes_inner_; 
  ghost_nodes_outer_ptr = &ghost_nodes_outer_; 
  global_mesh_ptr       = &global_mesh_;
}

//------------------------------------------------------------------------------------------------

void
EmbeddedBoundaryOperator::SetupIntersectors()
{
  for(int i=0; i<intersector.size(); i++) {
    intersector[i] = new Intersector(comm, *dms_ptr, *iod_embedded_surfaces[i], surfaces[i],
                                     *coordinates_ptr, *ghost_nodes_inner_ptr, *ghost_nodes_outer_ptr,
                                     *global_mesh_ptr);
  }
}

//------------------------------------------------------------------------------------------------

void
EmbeddedBoundaryOperator::FindSolidBodies(std::multimap<int, std::pair<int,int> > &id2closure)
{

  // Part 1: Find inactive colors. Warning: When multiple surfaces have inactive regions that are close
  //         to each other or overlapping, the information collected here is invalid.
  inactive_colors.clear();
  for(int i=0; i<surfaces.size(); i++) {
    unique_ptr<EmbeddedBoundaryDataSet> EBDS = GetPointerToEmbeddedBoundaryData(i);
    int nRegions = EBDS->nRegions; //this is the number of *closures*. Colors -1, -2, ...
    for(int color = -1; color>=-nRegions; color--) {
      bool found = false;
      for(auto it = id2closure.begin(); it != id2closure.end(); it++) {
        if(it->second.first == i && it->second.second == color) {
          found = true;
          break;
        }
      }
      if(!found)
        inactive_colors.insert(std::make_pair(i, color));
    } 
  }


  // Part 2: Find inactive_elem_status. Needed for force computation
  inactive_elem_status.resize(surfaces.size());
  for(int surf=0; surf<surfaces.size(); surf++)
    inactive_elem_status[surf].resize(surfaces[surf].elems.size(), 0);

  vector<bool> touched(surfaces.size(), false);
  for(auto it = inactive_colors.begin(); it != inactive_colors.end(); it++) {
    int surf = it->first;
    int this_color = it->second;
    assert(intersector[surf]);
    vector<int> &status(inactive_elem_status[surf]);
    if(touched[surf]) { //be careful! Create a new vector for FindColorBoundary, then merge w/ existing one
      vector<int> tmp;
      intersector[surf]->FindColorBoundary(this_color, tmp);
      assert(tmp.size() == status.size());
      for(int i=0; i<tmp.size(); i++) {
        if(tmp[i]==1) {
          if(status[i]==0)
            status[i] = 1;
          else if(status[i]==2)
            status[i] = 3;
        }
        else if(tmp[i]==2) {
          if(status[i]==0)
            status[i] = 2;
          else if(status[i]==1)
            status[i] = 3;
        }
        else if(tmp[i]==3)
          status[i] = 3;
      }
    }
    else {
      intersector[surf]->FindColorBoundary(this_color, status);
      touched[surf] = true;
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

    if(F_prev[i].size()>0) //not the first time
      assert(F[i].size() == F_prev[i].size());
    else
      F_prev[i].resize(F[i].size(),Vec3D(0.0));

    // copy force
    for(int j=0; j<F[i].size(); j++)
      F_prev[i][j] = F[i][j];

    if(surfaces_prev[i].X.size()>0) //not the first time
      assert(surfaces[i].X.size() == surfaces_prev[i].X.size());
    else
      surfaces_prev[i].X.resize(surfaces[i].X.size(),Vec3D(0.0));
      
    // copy nodal coords
    for(int j=0; j<surfaces[i].X.size(); j++)
      surfaces_prev[i].X[j] = surfaces[i].X[j];
  }
}

//------------------------------------------------------------------------------------------------

void
EmbeddedBoundaryOperator::ComputeForces(SpaceVariable3D &V, SpaceVariable3D &ID)
{

  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);

  Vec5D***  v  = (Vec5D***) V.GetDataPointer();
  double*** id = ID.GetDataPointer();
  
  // loop through all the embedded surfaces
  for(int surf=0; surf<surfaces.size(); surf++) {

    // Clear force vector
    vector<Vec3D>&  Fs(F[surf]); //Nodal loads (TO BE COMPUTED)
    Fs.resize(surfaces[surf].X.size(), 0.0);

    // Get quadrature info
    int np = 0; //number of Gauss points
    switch (iod_embedded_surfaces[surf]->quadrature) {
      case EmbeddedSurfaceData::NONE :
        np = 0; //force is always 0 --> one-way coupling
        break;
      case EmbeddedSurfaceData::ONE_POINT :
        np = 1;
        break;
      case EmbeddedSurfaceData::THREE_POINT :
        np = 3;
        break;
      case EmbeddedSurfaceData::FOUR_POINT :
        np = 4;
        break;
      case EmbeddedSurfaceData::SIX_POINT :
        np = 6;
        break;
      default :
        print_error("*** Error: Detected unknown Gauss quadrature rule (%d).\n", 
                    (int)iod_embedded_surfaces[surf]->quadrature);
        exit_mpi();
    }

    // ------------------------------------
    if(np==0) //one-way coupling
      continue;
    // ------------------------------------

    // Collect info about the surface and intersection results
    vector<Vec3D>&  Xs(surfaces[surf].X);
    vector<Int3>&   Es(surfaces[surf].elems);
    vector<Vec3D>&  Ns(surfaces[surf].elemNorm);
    vector<double>& As(surfaces[surf].elemArea);
    vector<int>&    status(inactive_elem_status[surf]);

    vector<double> gweight(np, 0.0);
    vector<Vec3D>  gbary(np, 0.0); //barycentric coords of Gauss points (symmetric) 
    MathTools::GaussQuadraturesTriangle::GetParameters(np, gweight.data(), gbary.data());
 
    
    vector<int> scope;
    intersector[surf]->GetElementsInScope(scope);

    //Note that different subdomain scopes overlap. We need to avoid repetition!
    for(auto it = scope.begin(); it != scope.end(); it++) {

      int tid = *it; //triangle id
      Int3 n(Es[tid][0], Es[tid][1], Es[tid][2]);

      vector<Vec3D> tg(np, 0.0); //traction at each Gauss point, (-pI + tau)n --> a Vec3D

      assert(fabs(Ns[tid].norm()-1.0)<1.0e-12); //normal must be valid!

      for(int side=0; side<2; side++) { //loop through the two sides

        Vec3D normal = Ns[tid];
        if(side==1)
          normal *= -1.0;

        for(int p=0; p<np; p++) { //loop through Gauss points

          Vec3D xg = gbary[p][0]*Xs[n[0]] + gbary[p][1]*Xs[n[1]] + gbary[p][2]*Xs[n[2]];

          // Lofting (Multiple processors may process the same point (xg). Make sure they produce the same result
          double loft = CalculateLoftingHeight(xg, iod_embedded_surfaces[surf]->gauss_points_lofting);
          xg += loft*normal;
          //fprintf(stderr,"side = %d, p = %d, loft = %e, xg = %e %e %e.\n", side, p, loft,xg[0], xg[1], xg[2]);
          
          // Check if this Gauss point is in this subdomain.
          Int3 ijk;
          bool foundit = global_mesh_ptr->FindCellCoveringPoint(xg, ijk, false);
          if(!foundit || !coordinates_ptr->IsHere(ijk[0],ijk[1],ijk[2],false))
            continue;

          // Calculate traction at Gauss point on this "side"
          if(status[tid]==3 || status[tid]==side+1) //this side faces the interior of a solid body 
            tg[p] += -1.0*iod_embedded_surfaces[surf]->internal_pressure*normal;
          else 
            tg[p] += CalculateTractionAtPoint(xg, side, As[tid], normal, n, Xs, v, id); 
        }

      }

      // Now, tg carries traction from both sides of the triangle

      // Integrate (See KW's notes for the formula)
      for(int p=0; p<np; p++) {
        tg[p] *= As[tid];
        // each node of the triangle gets some load from this Gauss point
        for(int node=0; node<3; node++)
          Fs[n[node]] += gweight[p]*gbary[p][node]*tg[p];
      }
    }

    // Processor 0 assembles the loads on the entire surface
    if(mpi_rank==0)
      MPI_Reduce(MPI_IN_PLACE, (double*)Fs.data(), 3*Fs.size(), MPI_DOUBLE, MPI_SUM, 0, comm);
    else
      MPI_Reduce((double*)Fs.data(), NULL, 3*Fs.size(), MPI_DOUBLE, MPI_SUM, 0, comm);

/*
    Vec3D sum(0.0);
    for(auto&& f : Fs) {
      sum+= f;
      print("%e %e %e\n", f[0], f[1], f[2]);
    }
    print("sum = %e %e %e.\n", sum[0], sum[1], sum[2]);
*/
  }
  
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
    if(iod_embedded_surfaces[i]->provided_by_another_solver == EmbeddedSurfaceData::NO &&
       strcmp(iod_embedded_surfaces[i]->dynamics_calculator, "") == 0)
      continue; //this surface is fixed...
    intersector[i]->RecomputeFullCourse(surfaces_prev[i].X, phi_layers);
  }
}

//------------------------------------------------------------------------------------------------

void
EmbeddedBoundaryOperator::ApplyUserDefinedSurfaceDynamics(double t, double dt)
{
  for(int surf=0; surf<surfaces.size(); surf++) {
    if(strcmp(iod_embedded_surfaces[surf]->dynamics_calculator, "") == 0)
      continue; //not specified for this surface
    UserDefinedDynamics *calculator(std::get<0>(dynamics_calculator[surf]));
    assert(calculator);

    vector<Vec3D> &Xs(surfaces[surf].X);
    vector<Vec3D> &X0(surfaces[surf].X0);
    vector<Vec3D> disp(Xs.size(), 0.0);
    fprintf(stderr,"Hello... udot size is %d\n", (int)surfaces[surf].Udot.size());
    calculator->GetUserDefinedDynamics(t, disp.size(), (double*)X0.data(), (double*)Xs.data(),
                                       (double*)disp.data(), (double*)surfaces[surf].Udot.data());
    fprintf(stderr,"Hello again...\n");
    for(int i=0; i<Xs.size(); i++) {
      for(int j=0; j<3; j++)
        Xs[i][j] = X0[i][j] + disp[i][j];
    }
  }
}

//------------------------------------------------------------------------------------------------

unique_ptr<vector<unique_ptr<EmbeddedBoundaryDataSet> > >
EmbeddedBoundaryOperator::GetPointerToEmbeddedBoundaryData()
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
EmbeddedBoundaryOperator::GetPointerToEmbeddedBoundaryData(int i)
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
    if(iod_embedded_surfaces[i]->provided_by_another_solver == EmbeddedSurfaceData::YES) {
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

void
EmbeddedBoundaryOperator::OutputSurfaces()
{
  for(int i=0; i<lagout.size(); i++)
    lagout[i].OutputTriangulatedMesh(surfaces[i].X0, surfaces[i].elems);
}

//------------------------------------------------------------------------------------------------

void
EmbeddedBoundaryOperator::OutputResults(double time, double dt, int time_step, bool force_write)
{
  for(int i=0; i<lagout.size(); i++)
    lagout[i].OutputResults(time, dt, time_step, surfaces[i].X0, surfaces[i].X, F[i], force_write);
}

//------------------------------------------------------------------------------------------------

double
EmbeddedBoundaryOperator::CalculateLoftingHeight(Vec3D &p, double factor)
{
  if(factor==0)
    return 0.0;

  assert(factor>0);

  Int3 ijk;
  bool foundit = global_mesh_ptr->FindCellCoveringPoint(p, ijk, true);
  if(!foundit)
    return 0.0; //no lofting applied to triangles outside. They will not (and should not) get a force.

  double size = std::min(global_mesh_ptr->dx_glob[ijk[0]],
                         std::min(global_mesh_ptr->dy_glob[ijk[1]], global_mesh_ptr->dz_glob[ijk[2]]));

  return factor*size;
}

//------------------------------------------------------------------------------------------------

Vec3D
EmbeddedBoundaryOperator::CalculateTractionAtPoint(Vec3D &p, int side, double area,
                                                   Vec3D &normal, Int3 &tnodes, vector<Vec3D> &Xs, 
                                                   Vec5D*** v, double*** id)
{

  Int3 ijk0(INT_MAX);
  Vec3D xi;
  global_mesh_ptr->FindElementCoveringPoint(p, ijk0, &xi, true);

  int i,j,k;

  // find which nodes in the element are on the correct side given by "normal"
  bool sameside[2][2][2];
  for(int dk=0; dk<=1; dk++)
    for(int dj=0; dj<=1; dj++)
      for(int di=0; di<=1; di++) {
        i = ijk0[0] + di;
        j = ijk0[1] + dj;
        k = ijk0[2] + dk;

        assert(coordinates_ptr->IsHere(i,j,k,true)); //Not in the ghosted subdomain? Something is wrong

        if(coordinates_ptr->OutsidePhysicalDomainAndUnpopulated(i,j,k)) { //state variable unavailable here.
          sameside[dk][dj][di] = false;
          continue;
        }

        if(id[k][j][i] == INACTIVE_MATERIAL_ID) {
          sameside[dk][dj][di] = false;
          continue;
        }


        Vec3D x(global_mesh_ptr->GetX(i), global_mesh_ptr->GetY(j), global_mesh_ptr->GetZ(k));
        
        double tmp[3]; 
        double signed_distance = (side==0) ? GeoTools::ProjectPointToPlane(x, Xs[tnodes[0]], Xs[tnodes[1]],
                                                           Xs[tnodes[2]], tmp, &area, &normal)
                                           : GeoTools::ProjectPointToPlane(x, Xs[tnodes[1]], Xs[tnodes[0]],
                                                           Xs[tnodes[2]], tmp, &area, &normal);

        sameside[dk][dj][di] = (signed_distance>0);
      }


  // interpolate pressure at the point
  // We populate opposite side and inactive nodes by average of same side / active nodes.
  // TODO: This can be done more carefully.
  double pressure[2][2][2];
  double total_pressure = 0.0;
  int n_pressure = 0;
  for(int dk=0; dk<=1; dk++)
    for(int dj=0; dj<=1; dj++)
      for(int di=0; di<=1; di++) {

        if(!sameside[dk][dj][di])
          continue;

        i = ijk0[0] + di;
        j = ijk0[1] + dj;
        k = ijk0[2] + dk;

        pressure[dk][dj][di] = v[k][j][i][4]; //get pressure
        total_pressure += pressure[dk][dj][di];
        n_pressure++;
      }

  if(n_pressure==0) {
    fprintf(stderr,"\033[0;31m*** Error: No valid active nodes for interpolating pressure at "
                   "Gauss point (%e, %e, %e). Try reducing surface thickness.\n\033[0m", p[0], p[1], p[2]);
    exit(-1);
  }

  double avg_pressure = total_pressure/n_pressure;
  for(int dk=0; dk<=1; dk++)
    for(int dj=0; dj<=1; dj++)
      for(int di=0; di<=1; di++) {

        if(sameside[dk][dj][di])
          continue;

        i = ijk0[0] + di;
        j = ijk0[1] + dj;
        k = ijk0[2] + dk;

        pressure[dk][dj][di] = avg_pressure;
      }
  
  //Now, perform trilinear interpolation to get p at the point
  double my_pressure = MathTools::trilinear_interpolation(pressure[0][0][0], pressure[0][0][1],
                                      pressure[0][1][0], pressure[0][1][1], pressure[1][0][0], 
                                      pressure[1][0][1], pressure[1][1][0], pressure[1][1][1], (double*)xi);

  //TODO: Add viscous force later!

  return -1.0*my_pressure*normal;  

}

//------------------------------------------------------------------------------------------------









