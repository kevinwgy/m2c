/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<EmbeddedBoundaryOperator.h>
#include<Vector5D.h>
#include<CommunicationTools.h>
#include<GeoTools.h>
#include<trilinear_interpolation.h>
#include<gauss_quadratures.h>
#include<rbf_interp.hpp>
#include<deque>
#include<list>
#include<utility>
#include<fstream>
#include<string>
#include<iomanip> //std::setw
#include<memory.h> //unique_ptr
#include<dlfcn.h> //dlopen, dlclose

using std::string;
using std::map;
using std::vector;
using std::unique_ptr;

extern int INACTIVE_MATERIAL_ID;
extern int verbose;

//------------------------------------------------------------------------------------------------

EmbeddedBoundaryOperator::EmbeddedBoundaryOperator(MPI_Comm &comm_, IoData &iod_, bool surface_from_other_solver) 
                        : comm(comm_), iod(iod_), hasSurfFromOtherSolver(surface_from_other_solver),
                          dms_ptr(NULL), coordinates_ptr(NULL), ghost_nodes_inner_ptr(NULL),
                          ghost_nodes_outer_ptr(NULL), global_mesh_ptr(NULL),
                          cylindrical_symmetry(false) 
{
  // count surfaces from files
  int counter = iod_.ebm.embed_surfaces.surfaces.dataMap.size();

  // get the correct size of surfaces and F
  surfaces.resize(counter);
  F.resize(counter);
  F_over_A.resize(counter);
  Anodal.resize(counter);

  surfaces_prev.resize(surfaces.size());
  F_prev.resize(F.size());
  F_over_A_prev.resize(F_over_A.size());
 
  // set default boundary type to "None"
  surface_type.assign(surfaces.size(), EmbeddedSurfaceData::None);
  iod_embedded_surfaces.assign(surfaces.size(), NULL);

  // read surfaces from files
  int nConcurrent = 0;
  for(auto it = iod_.ebm.embed_surfaces.surfaces.dataMap.begin();
           it != iod_.ebm.embed_surfaces.surfaces.dataMap.end(); it++) {

    int index = it->first;
    if(index<0 || index>=counter) {
      print_error("*** Error: Detected error in the indices of embedded surfaces (id = %d).", index);
      exit_mpi();
    }

    iod_embedded_surfaces[index] = it->second;

    if(index==0) {
      if(surface_from_other_solver) {
        if(it->second->provided_by_another_solver != EmbeddedSurfaceData::YES) {
          print_error("*** Error: Conflict input in EmbeddedSurface[%d]. Should mesh be provided by another solver?", 
                      index); 
          exit_mpi();
        } 
        nConcurrent++;
        continue; //no file to read
      } else {
        if(it->second->provided_by_another_solver != EmbeddedSurfaceData::NO) {
          print_error("*** Error: Conflict input in EmbeddedSurface[%d]. Should mesh be provided by user?", index); 
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

    ReadMeshFile(it->second->filename, surfaces[index].X, surfaces[index].elems);

    surfaces[index].X0 = surfaces[index].X;

    // initialize the velocity vector (to 0)
    surfaces[index].Udot.assign(surfaces[index].X.size(), 0.0);

    surfaces[index].BuildConnectivities();
    surfaces[index].CalculateNormalsAndAreas();
/*
    bool orientation = surfaces[index].CheckSurfaceOrientation(); 
    bool closed = surfaces[index].CheckSurfaceClosedness();
*/
  }

  print("- Activated the Embedded Boundary Method. Detected %d surface(s) (%d from concurrent program(s)).\n\n",
        surfaces.size(), nConcurrent); 

  // set NULL to intersector pointers
  intersector.assign(surfaces.size(), NULL);

  // set NULL to multi-surface intersector pointers (after santity check)
  std::set<std::pair<int,int> > surface_pairs;
  for(auto&& my_pair : iod_.ebm.embed_surfaces.surface_intersections.dataMap) {
    SurfaceIntersectionData* inter(my_pair.second);
    surface_pairs.insert(std::make_pair(std::min(inter->surface1_id, inter->surface2_id),
                                        std::max(inter->surface1_id, inter->surface2_id)));
  }
  if(surface_pairs.size() != iod_.ebm.embed_surfaces.surface_intersections.dataMap.size()) {
    print_error("*** Error: Detected identical SurfaceIntersection pairs.\n");
    exit_mpi();
  }
  multi_intersector.assign(surface_pairs.size(), NULL);


  // setup output
  for(int i=0; i<(int)surfaces.size(); i++)
    lagout.push_back(LagrangianOutput(comm, iod_embedded_surfaces[i]->output));

  // setup dynamics_calculator
  SetupUserDefinedDynamicsCalculator();

  // set cylindrical_symmetry
  cylindrical_symmetry = (iod_.mesh.type == MeshData::CYLINDRICAL);

  // initialize 2d_to_3d to false
  twoD_to_threeD.assign(counter, false);

}

//------------------------------------------------------------------------------------------------
// A constructor for tracking a single embedded surface provided using a mesh file
// The surface may contain multiple enclosures
// NOTE: surface intersections IGNORED in this constructor
EmbeddedBoundaryOperator::EmbeddedBoundaryOperator(MPI_Comm &comm_, IoData &iod_, EmbeddedSurfaceData &iod_surface)
                        : comm(comm_), iod(iod_), hasSurfFromOtherSolver(false),
                          dms_ptr(NULL), coordinates_ptr(NULL), ghost_nodes_inner_ptr(NULL),
                          ghost_nodes_outer_ptr(NULL), global_mesh_ptr(NULL)
{

  // get the correct size of surfaces and F
  surfaces.resize(1);
  F.resize(1);
  F_over_A.resize(1);
  Anodal.resize(1);

  surfaces_prev.resize(1);
  F_prev.resize(1);
  F_over_A_prev.resize(1);
 
  // set default boundary type to "None"
  iod_embedded_surfaces.assign(1, &iod_surface);
  surface_type.assign(1, iod_surface.type);

  ReadMeshFile(iod_surface.filename, surfaces[0].X, surfaces[0].elems);

  surfaces[0].X0 = surfaces[0].X;
  surfaces[0].Udot.assign(surfaces[0].X.size(), 0.0);

  surfaces[0].BuildConnectivities();
  surfaces[0].CalculateNormalsAndAreas();
/*
  bool orientation = surfaces[index].CheckSurfaceOrientation(); 
  bool closed = surfaces[index].CheckSurfaceClosedness();
*/

  print("- Activated the Embedded Boundary Method to track the surface provided in %s\n\n",
        iod_surface.filename);

  // set NULL to intersector pointers
  intersector.assign(1, NULL);

  // setup output
  lagout.push_back(LagrangianOutput(comm, iod_embedded_surfaces[0]->output));

  // setup dynamics_calculator
  SetupUserDefinedDynamicsCalculator();
}

//------------------------------------------------------------------------------------------------

EmbeddedBoundaryOperator::~EmbeddedBoundaryOperator()
{

  for(int i=0; i<(int)intersector.size(); i++)
    if(intersector[i])
      delete intersector[i];

  for(int i=0; i<(int)multi_intersector.size(); i++)
    if(multi_intersector[i])
      delete multi_intersector[i];

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
  for(int i=0; i<(int)intersector.size(); i++)
    if(intersector[i])
      intersector[i]->Destroy();

  for(int i=0; i<(int)multi_intersector.size(); i++)
    if(multi_intersector[i])
      multi_intersector[i]->Destroy();
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


  // determine twoD_to_threeD (in the constructor, it is set to false by default)
  assert(global_mesh_ptr);
  if(global_mesh_ptr->IsMesh2D()) {
    bool involves_2D_3D_mapping = false;
    for(int surf=0; surf<(int)twoD_to_threeD.size(); surf++) {
      twoD_to_threeD[surf] = IsEmbeddedSurfaceIn3D(surf);
      if(twoD_to_threeD[surf])
        involves_2D_3D_mapping = true;
    }

    if(involves_2D_3D_mapping) {
      print("- Detected embedded surface(s) in 3D, while M2C domain is in 2D.\n");
      if(!global_mesh_ptr->two_dimensional_xy) {
        print_error("*** Error: 2D Mesh must be in x-y plane (for 2D->3D mapping).\n");
        exit_mpi();
      }
      else
        print("  o Activated 2D->3D mapping in load calculation. Lofting may be increased internally.\n");
    }
  }

}

//------------------------------------------------------------------------------------------------

bool
EmbeddedBoundaryOperator::IsEmbeddedSurfaceIn3D(int surf)
{
  assert(surf<(int)surfaces.size());
  vector<Vec3D> &X(surfaces[surf].X);
  assert(X.size()>0);

  std::set<double> xyz[3]; //using *set* to eliminate duplicates
  for(auto&& x : X)
    for(int d=0; d<3; d++)
      xyz[d].insert(x[d]);

  for(int d=0; d<3; d++) {
    if(xyz[d].size()<3)
      return false; //this dimension is not discretized. Cannot be 3D.

    vector<double> coords(xyz[d].begin(), xyz[d].end()); //build a vector to store the coords
    std::sort(coords.begin(), coords.end()); //increasing order
    double hmax = 0.0;
    for(int i=0; i<(int)coords.size()-1; i++)
      hmax = std::max(hmax, coords[i+1]-coords[i]);

    if(hmax>0.99*(coords.back() - coords.front())) //direction d has not been discretized
      return false; //not 3D
  }

  return true;
} 

//------------------------------------------------------------------------------------------------

void
EmbeddedBoundaryOperator::SetupIntersectors()
{
  for(int i=0; i<(int)intersector.size(); i++) {
    intersector[i] = new Intersector(comm, *dms_ptr, *iod_embedded_surfaces[i], surfaces[i],
                                     *coordinates_ptr, *ghost_nodes_inner_ptr, *ghost_nodes_outer_ptr,
                                     *global_mesh_ptr);

  }

  assert(multi_intersector.size() == iod.ebm.embed_surfaces.surface_intersections.dataMap.size());
  int iter=0;
  for(auto&& my_pair : iod.ebm.embed_surfaces.surface_intersections.dataMap) {
    multi_intersector[iter] = new MultiSurfaceIntersector(comm, *dms_ptr, *my_pair.second, *coordinates_ptr,
                                                          surfaces, intersector, *ghost_nodes_inner_ptr,
                                                          *ghost_nodes_outer_ptr, *global_mesh_ptr);
    iter++;
  }
 
}

//------------------------------------------------------------------------------------------------

void
EmbeddedBoundaryOperator::FindSolidBodies(std::multimap<int, std::pair<int,int> > &id2closure)
{

  // Part 1: Find inactive colors. Warning: When multiple surfaces have inactive regions that are close
  //         to each other or overlapping, the information collected here is invalid. (TODO)
  inactive_colors.clear();
  for(int i=0; i<(int)surfaces.size(); i++) {
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
  for(int surf=0; surf<(int)surfaces.size(); surf++)
    inactive_elem_status[surf].assign(surfaces[surf].elems.size(), 0);

  vector<bool> touched(surfaces.size(), false);
  for(auto it = inactive_colors.begin(); it != inactive_colors.end(); it++) {
    int surf = it->first;
    int this_color = it->second;
    assert(intersector[surf]);
    vector<int> &status(inactive_elem_status[surf]);
    if(touched[surf]) { //be careful! Create a new vector for FindColorBoundary, then merge w/ existing one
      vector<int> tmp;
      if(twoD_to_threeD[surf])
        intersector[surf]->FindColorBoundary3Dto2D(this_color, cylindrical_symmetry, tmp);
      else
        intersector[surf]->FindColorBoundary(this_color, tmp);

      assert(tmp.size() == status.size());
      for(int i=0; i<(int)tmp.size(); i++) {
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
      if(twoD_to_threeD[surf])
        intersector[surf]->FindColorBoundary3Dto2D(this_color, cylindrical_symmetry, status);
      else
        intersector[surf]->FindColorBoundary(this_color, status);
      touched[surf] = true;
    }
  }


  // output the wetted sides (i.e. active_elem_status)
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);
  for(int surf=0; surf<(int)surfaces.size(); surf++) {

    if(mpi_rank != 0)
      continue;

    if(strcmp(iod_embedded_surfaces[surf]->output.wetting_output_filename,"")) {

      char outname[512];
      snprintf(outname, 512, "%s%s", iod_embedded_surfaces[surf]->output.prefix,
              iod_embedded_surfaces[surf]->output.wetting_output_filename);

      std::fstream out(outname, std::fstream::out);
      if(!out.is_open()) {
        fprintf(stdout,"\033[0;31m*** Error: Cannot write file %s.\n\033[0m", outname);
        exit(-1);
      }

      vector<Vec3D>&  Xs(surfaces[surf].X);
      vector<Int3>&   Es(surfaces[surf].elems);
      vector<Vec3D>&  Ns(surfaces[surf].elemNorm);

      // find the median element "size" --> length of markers
      vector<double> tmp = surfaces[surf].elemArea; //make a copy
      auto mid = tmp.begin() + tmp.size()/2;
      std::nth_element(tmp.begin(), mid, tmp.end()); //partial sort to find median
      double midarea = tmp[tmp.size()/2];
      assert(midarea>=0);
      double amplification_factor = 2.0;
      double marker_length = amplification_factor*sqrt(midarea*2.0);

      vector<int> &status(inactive_elem_status[surf]);

      // Write nodes
      out << "Nodes WettedSurfacePoints" << std::endl;
      Vec3D p,q;
      for(int i=0; i<(int)Es.size(); i++) {
        Int3 &nod(Es[i]);
        switch (status[i]) {
          case 0 : //both sides are wetted
            p = (Xs[nod[0]]+Xs[nod[1]]+Xs[nod[2]])/3.0 - marker_length*Ns[i];
            q = p + 2.0*marker_length*Ns[i];
            break;
          case 1 : //negative side is wetted
            p = (Xs[nod[0]]+Xs[nod[1]]+Xs[nod[2]])/3.0;
            q = p - marker_length*Ns[i];
            break;
          case 2 : //positive side is wetted
            p = (Xs[nod[0]]+Xs[nod[1]]+Xs[nod[2]])/3.0;
            q = p + marker_length*Ns[i];
            break;
          case 3 : //neither side is wetted
            p = (Xs[nod[0]]+Xs[nod[1]]+Xs[nod[2]])/3.0;
            q = p;
            break;
        }
        out << std::setw(10) << 2*i+1
            << std::setw(14) << std::scientific << p[0]
            << std::setw(14) << std::scientific << p[1]
            << std::setw(14) << std::scientific << p[2] << "\n";
        out << std::setw(10) << 2*i+2
            << std::setw(14) << std::scientific << q[0]
            << std::setw(14) << std::scientific << q[1]
            << std::setw(14) << std::scientific << q[2] << "\n";
      }

      // Write line segments / "markers"
      out << "Elements Markers using WettedSurfacePoints" << std::endl;
      for(int i=0; i<(int)Es.size(); i++) {
        out << std::setw(10) << i+1 << "  1  "  //"1" for line segment 
            << std::setw(10) << 2*i+1
            << std::setw(10) << 2*i+2 << std::endl;
      }

      out.flush();
      out.close();
    }
  }

  MPI_Barrier(comm);

}

//------------------------------------------------------------------------------------------------

void
EmbeddedBoundaryOperator::ReadMeshFile(const char *filename, vector<Vec3D> &Xs, vector<Int3> &Es)
{
  string fname(filename);
  assert(fname.size()>0);
  auto loc = fname.find_last_of(".");
  if(loc>=fname.size()-1) {//assume the default format (top) if file extension not detected
    ReadMeshFileInTopFormat(filename, Xs, Es); 
    return;
  }

  string ext = fname.substr(loc+1);
  if(ext == "obj" || ext == "Obj" || ext == "OBJ")
    ReadMeshFileInOBJFormat(filename, Xs, Es);
  else if(ext == "stl" || ext == "Stl" || ext == "STL")
    ReadMeshFileInSTLFormat(filename, Xs, Es);
  else
    ReadMeshFileInTopFormat(filename, Xs, Es); 
}

//------------------------------------------------------------------------------------------------

void
EmbeddedBoundaryOperator::ReadMeshFileInTopFormat(const char *filename, vector<Vec3D> &Xs, vector<Int3> &Es)
{

  // read data from the surface input file.
  FILE *topFile;
  topFile = fopen(filename, "r");
  if(topFile == NULL) {
    print_error("*** Error: Cannot open embedded surface mesh file (%s).\n", filename);
    exit_mpi();
  }
 
  int MAXLINE = 500;
  char line[MAXLINE], key1[MAXLINE], key2[MAXLINE]; //, copyForType[MAXLINE];

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
      type_reading = 2;
      found_elems = true;

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
    assert(current_id==(int)Xs.size());
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
EmbeddedBoundaryOperator::ReadMeshFileInOBJFormat(const char *filename, vector<Vec3D> &Xs, vector<Int3> &Es)
{
  // This function is adapted from toys::obj2top. But unlike the "toy", it does not separate different groups

  Xs.clear();
  Es.clear();

  std::ifstream input(filename, std::ifstream::in);
  if(input.fail()) {
    print_error("*** Error: Cannot open embedded surface mesh file (%s).\n", filename);
    exit_mpi();
  }

  string line, word;

  double area_tol = 1e-12;

  Vec3D xyz;
  int nVerts, maxVerts = 1024, n1, n2, n3;
  vector<int> ind(maxVerts); //should not have polygons with more than 1024 vertices!
  vector<std::pair<string, vector<Int3> > > elem_groups;
  vector<Int3> *elems = NULL; 

  std::set<string> ignored_keywords;

  while(getline(input, line)) {

    auto first_nonspace_id = line.find_first_not_of(" ");
    if((unsigned)first_nonspace_id<line.size() && line[first_nonspace_id] == '#')
      continue;  //this line is comment

    std::stringstream linestream(line);
    if(!(linestream >> word)) //the first word in the line
      continue;

    if(word == "v") { //vertex
      linestream >> xyz[0] >> xyz[1] >> xyz[2]; //skipping vertex properties (if present)
      Xs.push_back(xyz);
    }
    else if(word == "g") { //element group
      string group_name;
      linestream >> group_name; //the second word in the line -> group name  
      while(linestream >> word) //in stl, group name may have multiple words...
        group_name = group_name + "_" + word;

      bool found(false);
      for(auto&& eg : elem_groups) {
        if(eg.first == group_name) {
          elems = &eg.second;
          found = true;
          break;
        }
      }
      if(!found) {
        elem_groups.push_back(std::make_pair(group_name, vector<Int3>()));
        elems = &elem_groups.back().second;
      }
    }
    else if(word == "f") { //face element

      if(elem_groups.empty()) { //user did not specify group name
        assert(elems==NULL);
        elem_groups.push_back(std::make_pair("Default", vector<Int3>()));
        elems = &elem_groups.back().second;
      }

      nVerts = 0;
      while(linestream >> word) { 
        ind[nVerts++] = std::stoi(word.substr(0,word.find("/")));
        if(nVerts>maxVerts) {
          print_error("*** Error: Found a face element in %s with more than %d nodes.\n",
                      filename, maxVerts); 
          exit_mpi(); 
        }
      }
      if(nVerts<3) {
        print_error("*** Error: Found a face element in %s with only %d nodes.\n",
                    filename, nVerts); 
        exit_mpi(); 
      }
      for(int i=1; i<nVerts-1; i++)
        elems->push_back(Int3(ind[0], ind[i], ind[i+1]));
    }
    else
      ignored_keywords.insert(word);
  }


  for(auto&& key : ignored_keywords) {
    if(!key.empty())
      print_warning("Warning: Ignored lines in %s starting with %s.\n",
                    filename, key.c_str());
  }

  if(verbose>=1) {
    print("- Found %d nodes in %s.\n", (int)Xs.size(), filename);
    for(auto&& eg : elem_groups)
      print("- Obtained %d triangular elements in group %s from %s.\n",
            (int)eg.second.size(), eg.first.c_str(), filename);

    if(elem_groups.size()>=2)
      print_warning("Warning: Merging multiple (%d) groups into one.\n", (int)elem_groups.size()); 
  }


  for(auto&& eg : elem_groups) {
    for(auto&& e : eg.second) {
      n1 = e[0]; n2 = e[1]; n3 = e[2];
      if(n1<=0 || n1>(int)Xs.size() || n2<=0 || n2>(int)Xs.size() ||
         n3<=0 || n3>(int)Xs.size()) {
        print_error("*** Error: Found element (%d %d %d) in %s with unknown node(s).\n",
                n1, n2, n3, filename);
        exit_mpi();
      }
      Vec3D cr = (Xs[n2-1] - Xs[n1-1])^(Xs[n3-1] - Xs[n1-1]);
      if(cr.norm() < area_tol) {
        print_warning("Warning: Detected a degenerate triangle with area %e --- dropped from the list.\n",
                      cr.norm());
        continue;
      }
      Es.push_back(Int3(n1-1,n2-1,n3-1)); //node id starts at 0
    }
  }

  input.close();
}

//------------------------------------------------------------------------------------------------

void
EmbeddedBoundaryOperator::ReadMeshFileInSTLFormat(const char *filename, vector<Vec3D> &Xs, vector<Int3> &Es)
{
  // This function is adapted from toys::stl2top.

  Xs.clear();
  Es.clear();
  
  std::ifstream input(filename, std::ifstream::in);

  string line;
  string word;
  getline(input, line);

  double area_tol = 1e-12;

  int nodeid(0);
  double x,y,z;
 
  while(true) {

    getline(input, line);

    if(line.compare(0,8,"endsolid") == 0)
      break;

    getline(input, line);
    input >> word >> x >> y >> z;
    Vec3D X1(x,y,z);
    input >> word >> x >> y >> z;
    Vec3D X2(x,y,z);
    input >> word >> x >> y >> z;
    Vec3D X3(x,y,z);
    Vec3D cr = (X2 - X1)^(X3 - X1);
    getline(input, line); //get the end of line
    getline(input, line);
    getline(input, line);

    if(cr.norm() < area_tol) {
      print_warning("Warning: Detected a degenerate triangle with area %e --- dropped from the list.\n",
                    cr.norm());
      continue;
    }

    Xs.push_back(X1);
    Xs.push_back(X2);
    Xs.push_back(X3);
    Es.push_back(Int3(nodeid, nodeid+1, nodeid+2));
    nodeid += 3;

  }

  if(verbose>=1)
    print("Found %d triangular elements in %s.\n", (int)Es.size(), filename);

  input.close();
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
  assert(F_over_A.size() == F_over_A_prev.size());

  // loop through all the surfaces
  for(int i = 0; i < (int)F.size(); i++) {

    if(F_prev[i].size()>0) {//not the first time
      assert(F[i].size() == F_prev[i].size());
      assert(F_over_A[i].size() == F_over_A_prev[i].size());
    } else {
      F_prev[i].assign(F[i].size(),Vec3D(0.0));
      F_over_A_prev[i].assign(F_over_A[i].size(),Vec3D(0.0));
    }

    // copy force & nodal area
    for(int j=0; j<(int)F[i].size(); j++) {
      F_prev[i][j] = F[i][j];
      F_over_A_prev[i][j] = F_over_A[i][j];
    }

    if(surfaces_prev[i].X.size()>0) //not the first time
      assert(surfaces[i].X.size() == surfaces_prev[i].X.size());
    else
      surfaces_prev[i].X.assign(surfaces[i].X.size(),Vec3D(0.0));
      
    // copy nodal coords
    for(int j=0; j<(int)surfaces[i].X.size(); j++)
      surfaces_prev[i].X[j] = surfaces[i].X[j];
  }
}

//------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// NOTE: ONLY PROC 0 HAS THE CORRECT FORCES FOR THE ENTIRE SURFACE
//----------------------------------------------------------------

void
EmbeddedBoundaryOperator::ComputeForces(SpaceVariable3D &V, SpaceVariable3D &ID)
{

/*
  V.StoreMeshCoordinates(*coordinates_ptr);
  V.WriteToVTRFile("V.vtr","sol");
  ID.StoreMeshCoordinates(*coordinates_ptr);
  ID.WriteToVTRFile("ID.vtr","ID");
*/

  Vec5D***  v  = (Vec5D***) V.GetDataPointer();
  double*** id = ID.GetDataPointer();
  
  // loop through all the embedded surfaces
  for(int surf=0; surf<(int)surfaces.size(); surf++) {

    // Get force vector
    vector<Vec3D>&  Fs(F[surf]); //Nodal loads (TO BE COMPUTED)
    vector<Vec3D>&  FAs(F_over_A[surf]); //Nodal loads divided by nodal area (TO BE COMPUTED)

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

    if(!twoD_to_threeD[surf]) 
      ComputeForcesOnSurfaceDirectly(surf, np, v, id, Fs, FAs);
    else
      ComputeForcesOnSurface2DTo3D(surf, np, v, id, Fs, FAs);

/*
    Vec3D sum(0.0);
    for(int i=0; i<Fs.size(); i++) {
      sum+= Fs[i];
      print("Tri %d: %e %e %e (%e)\n", i, Fs[i][0], Fs[i][1], Fs[i][2], Fs[i].norm());
    }
    print("sum = %e %e %e.\n", sum[0], sum[1], sum[2]);
*/
  }
  
  V.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();

}

//------------------------------------------------------------------------------------------------

double
EmbeddedBoundaryOperator::TrackSurfaces(int phi_layers)
{
  assert(phi_layers>0);

  double max_dist = -DBL_MAX;

  for(auto&& surf : surfaces)
    surf.CalculateNormalsAndAreas();

  for(int i = 0; i < (int)intersector.size(); i++) {
    intersector[i]->BuildKDTreeAndFindIntersections();
    intersector[i]->FloodFillColors();
  }

  //check & handle potential surface intersections (if specified by user)
  vector<bool> modified(intersector.size(), false);
  for(auto&& multiX : multi_intersector) {

    if(!multiX->CheckSurfaceIntersections())
      continue;

    //Now, we know these two surfaces intersect.
    
    int numNew = multiX->FindNewEnclosures();
    if(numNew>0) {
      if(verbose>=1)
        print("    o Found %d new enclosures due to embedded surface intersections.\n", numNew);
      int xid = multiX->UpdateIntersectors(); // modifies the involved single-surface intersectors
      if(xid>=0) //indeed, an interesector was modified
        modified[xid] = true;

      //NOTE: COLOR IN Intersectors ARE NOT UPDATED!
    }
  }

  // compute unsigned shortest distance
  for(int i = 0; i < (int)intersector.size(); i++) {
    double max_dist0 = intersector[i]->ComputeUnsignedDistance(phi_layers, modified[i] ? 1 : -1);
    if(max_dist0>max_dist)
      max_dist = max_dist0;
  }

  return max_dist;
}

//------------------------------------------------------------------------------------------------

void
EmbeddedBoundaryOperator::TrackUpdatedSurfaces()
{
  double max_dist = -DBL_MAX;

  int phi_layers = 3;
  for(int i=0; i<(int)intersector.size(); i++) {

    //avoid unncessary/repetitive work
    if(iod_embedded_surfaces[i]->provided_by_another_solver == EmbeddedSurfaceData::NO &&
       strcmp(iod_embedded_surfaces[i]->dynamics_calculator, "") == 0) {
       //this surface is static

      bool involved = false;
      for(auto&& multiX : multi_intersector) {
        for(int surf=0; surf<(int)multiX->GetNumSurfaces(); surf++)
          if(multiX->GetSurfaceID(surf) == i) {
            involved = true;
            break;
          }
        if(involved)
          break;
      }
      if(!involved)
        continue; //skip current intersector/surface
    }

    // actual computations
    surfaces[i].CalculateNormalsAndAreas();
    intersector[i]->BuildKDTreeAndFindIntersections();
    intersector[i]->FindSweptNodes(surfaces_prev[i].X);
    intersector[i]->RefillAfterSurfaceUpdate();
  }

  //check & handle potential surface intersections (if specified by user)
  vector<bool> modified(intersector.size(), false);
  for(auto&& multiX : multi_intersector) {

    bool already_found = multiX->FoundMultiSurfIntersection();
    if(!multiX->CheckSurfaceIntersections())
      continue;

    //Now, we know these two surfaces intersect.
    //
    int numNew = already_found ? multiX->FindNewEnclosuresAfterSurfaceUpdate() //true->1
                               : multiX->FindNewEnclosures();
    if(numNew>0) {
      if(verbose>=1)
        print("    o Found %d new enclosures due to embedded surface intersections.\n", numNew);
      int xid = multiX->UpdateIntersectors(); // modifies the involved single-surface intersectors
                                              // NOTE: "color" is NOT updated
      if(xid>=0) //indeed, an interesector was modified
        modified[xid] = true;

      //NOTE: COLOR IS NOT UPDATED!
    }

  }

  for(int i = 0; i < (int)intersector.size(); i++) {
    if(iod_embedded_surfaces[i]->provided_by_another_solver == EmbeddedSurfaceData::NO &&
       strcmp(iod_embedded_surfaces[i]->dynamics_calculator, "") == 0 &&
       modified[i] == false) {

      if(!multi_intersector.empty())
        //in this case, no need to re-compute distance, but may need to rebuild scopeN
        intersector[i]->BuildLayerNSubdomainScopeAndKDTree(phi_layers);

      continue;
    }
    //"ComputeUnsignedDistance": first build scopN and KDTree.
    double max_dist0 = intersector[i]->ComputeUnsignedDistance(phi_layers, modified[i] ? 1 : -1);
    if(max_dist0>max_dist)
      max_dist = max_dist0;
  }


  //return max_dist; //WARNING: IF ALL SURFACES ARE FIXED, MAX_DIST gets garbage
  return;
}

//------------------------------------------------------------------------------------------------

void
EmbeddedBoundaryOperator::ApplyUserDefinedSurfaceDynamics(double t, double dt)
{
  for(int surf=0; surf<(int)surfaces.size(); surf++) {
    if(strcmp(iod_embedded_surfaces[surf]->dynamics_calculator, "") == 0)
      continue; //not specified for this surface
    UserDefinedDynamics *calculator(std::get<0>(dynamics_calculator[surf]));
    assert(calculator);

    vector<Vec3D> &Xs(surfaces[surf].X);
    vector<Vec3D> &X0(surfaces[surf].X0);
    vector<Vec3D> disp(Xs.size(), 0.0);
    calculator->GetUserDefinedDynamics(t, dt, disp.size(), (double*)X0.data(), (double*)Xs.data(),
                                       (double*)disp.data(), (double*)surfaces[surf].Udot.data());
    for(int i=0; i<(int)Xs.size(); i++) {
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

  for(int i=0; i<(int)intersector.size(); i++) {
    assert(intersector[i]); //should not call this function if surfaces are not "tracked"
    p->push_back(intersector[i]->GetPointerToResults());
  }

  return p;
}

//------------------------------------------------------------------------------------------------

unique_ptr<EmbeddedBoundaryDataSet>
EmbeddedBoundaryOperator::GetPointerToEmbeddedBoundaryData(int i)
{
  assert(i>=0 && i<(int)intersector.size());
  assert(intersector[i]);
  return intersector[i]->GetPointerToResults();
}

//------------------------------------------------------------------------------------------------

void
EmbeddedBoundaryOperator::SetupUserDefinedDynamicsCalculator()
{
  dynamics_calculator.assign(surfaces.size(), std::make_tuple(nullptr,nullptr,nullptr));
  for(int i=0; i<(int)surfaces.size(); i++) {
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
  for(int i=0; i<(int)lagout.size(); i++)
    lagout[i].OutputTriangulatedMesh(surfaces[i].X0, surfaces[i].elems);
}

//------------------------------------------------------------------------------------------------

void
EmbeddedBoundaryOperator::OutputResults(double time, double dt, int time_step, bool force_write)
{
  for(int i=0; i<(int)lagout.size(); i++)
    lagout[i].OutputResults(time, dt, time_step, surfaces[i].X0, surfaces[i].X, F[i], &F_over_A[i], force_write);
}

//------------------------------------------------------------------------------------------------

void
EmbeddedBoundaryOperator::ComputeForcesOnSurfaceDirectly(int surf, int np, Vec5D*** v, double*** id,
                                                         vector<Vec3D> &Fs, vector<Vec3D> &FAs)
{

  Fs.assign(surfaces[surf].X.size(), 0.0);

  if(FAs.size()!=surfaces[surf].X.size())
    FAs.assign(surfaces[surf].X.size(), 0.0);

  vector<double>& An(Anodal[surf]);
  An.assign(surfaces[surf].X.size(), 0.0);

  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);

  // Collect info about the surface and intersection results
  vector<Vec3D>&  Xs(surfaces[surf].X);
  vector<Int3>&   Es(surfaces[surf].elems);
  vector<Vec3D>&  Ns(surfaces[surf].elemNorm);
  vector<double>& As(surfaces[surf].elemArea);
  vector<int>&    status(inactive_elem_status[surf]);

  vector<double> gweight(np, 0.0);
  vector<Vec3D>  gbary(np, 0.0); //barycentric coords of Gauss points (symmetric) 
  MathTools::GaussQuadraturesTriangle::GetParameters(np, gweight.data(), gbary.data());
 
  vector<Vec3D> xgs(np, 0.0); //internal var.
    
  vector<int> scope;
  intersector[surf]->GetElementsInScope1(scope);

  //Note that different subdomain scopes overlap. We need to avoid repetition!
  for(auto it = scope.begin(); it != scope.end(); it++) {

    int tid = *it; //triangle id

    if(status[tid]==3)
      continue; //both sides facing interactive regions (interior of solid body)

    if(!multi_intersector.empty() && //there can be dropped elements...
       (int)surfaces[surf].elemtag.size()>tid && surfaces[surf].elemtag[tid] == 1) //"1" means dropped
      continue;

    Int3 n(Es[tid][0], Es[tid][1], Es[tid][2]);

    // Get Gauss points (before lofting)
    for(int p=0; p<np; p++) 
      xgs[p] = gbary[p][0]*Xs[n[0]] + gbary[p][1]*Xs[n[1]] + gbary[p][2]*Xs[n[2]];

    // traction at each Gauss point, (-pI + tau)n --> a Vec3D
    vector<Vec3D> tg(np, 0.0); 

    // whether the fraction of area "controlled" by this Gauss point should be calculated by this cpu core.
    // Only affects An and FAs, not Fs
    vector<bool> calculate_An(np, false);

    assert(fabs(Ns[tid].norm()-1.0)<1.0e-12); //normal must be valid!

    for(int side=0; side<2; side++) { //loop through the two sides

      Vec3D normal = Ns[tid];
      if(side==1)
        normal *= -1.0;

      for(int p=0; p<np; p++) { //loop through Gauss points

        Vec3D xg = xgs[p];

        Int3 ijk0;
        if(!global_mesh_ptr->FindCellCoveringPoint(xg, ijk0, false))
          continue; //before lofting, we need to make sure the original Gauss point is in the domain.
                    //otherwise, there can be complexities.

        if(coordinates_ptr->IsHere(ijk0[0],ijk0[1],ijk0[2],false))
          calculate_An[p] = true;  //this is to avoid double-counting area
        
        // Lofting (Multiple processors may process the same point (xg). Make sure they produce the same result
        double loft = CalculateLoftingHeight(xg, iod_embedded_surfaces[surf]->gauss_points_lofting);
        xg += loft*normal;
          
        // Check if this Gauss point is in this subdomain.
        Int3 ijk;
        bool foundit = global_mesh_ptr->FindCellCoveringPoint(xg, ijk, false);
        if(!foundit) { //pull it back to the domain (not necessarily the current subdomain)
          Vec3D xg0 = xg - loft*normal;
          double pull_back = 0.5*loft;
          for(int step=0; step<5; step++) {
            xg -= pull_back*normal;
            foundit = global_mesh_ptr->FindCellCoveringPoint(xg, ijk, false);
            if(foundit)
              break;
            pull_back /= 2.0;
          }
          if(!foundit) {
            xg = xg0;
            ijk = ijk0;
            foundit = true;
          }
        }

        if(!coordinates_ptr->IsHere(ijk[0],ijk[1],ijk[2],false))
          continue;

        // Calculate traction at Gauss point on this "side"
        if(status[tid]==side+1) //this side faces the interior of a solid body 
          tg[p] += -1.0*iod_embedded_surfaces[surf]->internal_pressure*normal;
        else 
          tg[p] += CalculateTractionAtPoint(xg, normal, v, id); 

      }

    }

    // Now, tg carries traction from both sides of the triangle

    // Integrate (See KW's notes for the formula)
    for(int p=0; p<np; p++) {
      tg[p] *= As[tid];
      // each node of the triangle gets some load from this Gauss point
      for(int node=0; node<3; node++) {
        double coeff = gweight[p]*gbary[p][node];
        Fs[n[node]] += coeff*tg[p];
        if(calculate_An[p])
          An[n[node]] += coeff*As[tid];
      } 
    }
  }

  // Processor 0 assembles the loads on the entire surface
  if(mpi_rank==0) {
    MPI_Reduce(MPI_IN_PLACE, (double*)Fs.data(), 3*Fs.size(), MPI_DOUBLE, MPI_SUM, 0, comm);
    MPI_Reduce(MPI_IN_PLACE, An.data(), An.size(), MPI_DOUBLE, MPI_SUM, 0, comm);
  } else {
    MPI_Reduce((double*)Fs.data(), NULL, 3*Fs.size(), MPI_DOUBLE, MPI_SUM, 0, comm);
    MPI_Reduce(An.data(), NULL, An.size(), MPI_DOUBLE, MPI_SUM, 0, comm);
  }

  if(mpi_rank==0) {
    for(int i=0; i<(int)Fs.size(); i++)
      FAs[i] = An[i]==0.0 ? 0.0 : Fs[i]/An[i];
  }

  MPI_Barrier(comm);
}

//------------------------------------------------------------------------------------------------
// Similar to ComputeForcesOnSurfaceDirectly. Project each Gauss point onto x-y plane.
// Note: For robustness, larger lofting (relative to triangle size) is imposed, regardless of
//       user input. Otherwise, it could happen that after converting to x-y plane, the point lands
//       on the wrong (i.e., opposite) side of the embedded surface.
void
EmbeddedBoundaryOperator::ComputeForcesOnSurface2DTo3D(int surf, int np, Vec5D*** v, double*** id,
                                                       vector<Vec3D> &Fs, vector<Vec3D> &FAs)
{

  assert(global_mesh_ptr->two_dimensional_xy);
  assert(twoD_to_threeD[surf]);

  Fs.assign(surfaces[surf].X.size(), 0.0);

  if(FAs.size()!=surfaces[surf].X.size())
    FAs.assign(surfaces[surf].X.size(), 0.0);

  vector<double>& An(Anodal[surf]);
  An.assign(surfaces[surf].X.size(), 0.0);

  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);

  // Collect info about the surface and intersection results
  vector<Vec3D>&  Xs(surfaces[surf].X);
  vector<Int3>&   Es(surfaces[surf].elems);
  vector<Vec3D>&  Ns(surfaces[surf].elemNorm);
  vector<double>& As(surfaces[surf].elemArea);
  vector<int>&    status(inactive_elem_status[surf]);

  vector<double> gweight(np, 0.0);
  vector<Vec3D>  gbary(np, 0.0); //barycentric coords of Gauss points (symmetric) 
  MathTools::GaussQuadraturesTriangle::GetParameters(np, gweight.data(), gbary.data());
 
  vector<Vec3D> xgs(np, 0.0); //internal var.

  int worker_proc = 0; //mpi_rank 0 collects the data and send it to the structural solver 
  set<int> my_indices_to_send;
   

  // -----------------------------------------------------------
  // Step 1: Each processor loops through all embedded boundary elements.
  //         Find equiv point on x-y plane & compute traction there.
  // -----------------------------------------------------------
  double z0 = global_mesh_ptr->GetZ(0); //It is assumed that the mesh is on x-y plane
  for(size_t tid=0; tid<Es.size(); tid++) {

    if(status[tid]==3)
      continue; //both sides facing interactive regions (interior of solid body)

    if(!multi_intersector.empty() && //there can be dropped elements...
       surfaces[surf].elemtag.size()>tid && surfaces[surf].elemtag[tid] == 1) //"1" means dropped
      continue;

    // Get the three vertices of the triangle
    Int3 n(Es[tid][0], Es[tid][1], Es[tid][2]);

    // Get Gauss points (before lofting)
    for(int p=0; p<np; p++)
      xgs[p] = gbary[p][0]*Xs[n[0]] + gbary[p][1]*Xs[n[1]] + gbary[p][2]*Xs[n[2]];

    // Get element normal and convert it to x-y plane
    assert(fabs(Ns[tid].norm()-1.0)<1.0e-12); //normal must be valid!
    Vec3D normal = Ns[tid];
    Vec3D normal_xy(0.0); //all three components are 0
    if(cylindrical_symmetry) {
      normal_xy[0] = normal[0];
      if(xgs[0][1]*normal[1]+xgs[0][2]*normal[2]>=0) //any point on the triangle is fine
        normal_xy[1] =  sqrt(normal[1]*normal[1] + normal[2]*normal[2]);
      else
        normal_xy[1] = -sqrt(normal[1]*normal[1] + normal[2]*normal[2]);
    } else {//truly 2D, normal[2] should be 0 or extremely small.
      double denom = sqrt(normal[0]*normal[0] + normal[1]*normal[1]);
      normal_xy[0] = normal[0]/denom;
      normal_xy[1] = normal[1]/denom;
    }

    // traction at each Gauss point, (-pI + tau)n --> a Vec3D
    vector<Vec3D> tg(np, 0.0);

    // whether the fraction of area "controlled" by this Gauss point should be calculated by this cpu core.
    // Only affects An and FAs, not Fs
    vector<bool> calculate_An(np, false);

    double circumradius = (Xs[n[0]]-Xs[n[1]]).norm()*(Xs[n[1]]-Xs[n[2]]).norm()*(Xs[n[2]]-Xs[n[0]]).norm()
                        / (4.0*As[tid]);  //a*b*c/(4*area)
    bool in_subdomain = false; //at least one side/Gauss point in subdomain ==> needs to be assembled
    for(int side=0; side<2; side++) { //loop through the two sides

      if(side==1) { //directly modify 'normal' and 'normal_xy'!
        normal    *= -1.0;
        normal_xy *= -1.0;
      }

      for(int p=0; p<np; p++) { //loop through Gauss points

        Vec3D xg = xgs[p];
        if(cylindrical_symmetry)
          xg[1] = sqrt(xg[1]*xg[1]+xg[2]*xg[2]);
        xg[2] = z0; //xg should be within the fluid mesh

        Int3 ijk0;
        if(!global_mesh_ptr->FindCellCoveringPoint(xg, ijk0, false))
            continue; //before lofting, we need to make sure the un-lofted Gauss point is in the domain.
                      //otherwise, no force calculated at this Gauss point

        if(coordinates_ptr->IsHere(ijk0[0],ijk0[1],ijk0[2],false))
          calculate_An[p] = true;  //this is to avoid double-counting area

        // Lofting (Multiple processors may process the same point. They must produce the same result)
        double loft = iod_embedded_surfaces[surf]->gauss_points_lofting * global_mesh_ptr->GetMinDXYZ(ijk0);
        loft = std::max(loft, 0.2*circumradius); //NOTE: avoid landing on wrong side when converted to 2D
        xg += loft*normal_xy;

        // Check if this lofted Gauss point is in current subdomain.
        Int3 ijk;
        bool foundit = global_mesh_ptr->FindCellCoveringPoint(xg, ijk, false);
        if(!foundit) { //pull it back to the domain (not necessarily the current subdomain)
          Vec3D xg0 = xg - loft*normal_xy; //backup
          double pull_back = 0.5*loft;
          for(int step=0; step<5; step++) {
            xg -= pull_back*normal_xy;
            foundit = global_mesh_ptr->FindCellCoveringPoint(xg, ijk, false);
            if(foundit)
              break;
            pull_back /= 2.0;
          }
          if(!foundit) {
            xg = xg0;
            ijk = ijk0;
            foundit = true;
          }
        }

        if(!coordinates_ptr->IsHere(ijk[0],ijk[1],ijk[2],false))
          continue;

        if(!in_subdomain)
          in_subdomain = true;

        if(status[tid]==side+1) //this side faces the interior of a solid body
          tg[p] += -1.0*iod_embedded_surfaces[surf]->internal_pressure*normal;
        else {
          Vec3D tg_xy = CalculateTractionAtPoint(xg, normal_xy, v, id);
          tg[p][0] += tg_xy[0];
          if(cylindrical_symmetry) {
            double norm_yz = sqrt(xgs[p][1]*xgs[p][1] + xgs[p][2]*xgs[p][2]);
            if(norm_yz != 0.0) {
              tg[p][1] += tg_xy[1]*xgs[p][1]/norm_yz;
              tg[p][2] += tg_xy[1]*xgs[p][2]/norm_yz;
            } else
              tg[p][1] += tg[p][2] = 0.0; 
          } else { //truly 2D
            tg[p][1] += tg_xy[1];
            tg[p][2] += 0.0;
          }
        }
      }
    }

    if(!in_subdomain)
      continue;

    // Now, tg carries traction from both sides of the triangle

    // Integrate (See KW's notes for the formula)
    for(int p=0; p<np; p++) {
      tg[p] *= As[tid];
      // each node of the triangle gets some load from this Gauss point
      for(int node=0; node<3; node++) {
        double coeff = gweight[p]*gbary[p][node];
        Fs[n[node]] += coeff*tg[p];
        if(calculate_An[p])
          An[n[node]] += coeff*As[tid];
        my_indices_to_send.insert(n[node]); //avoid duplication (set)
      } 
    }
  }

  // Prepare data to gather
  vector<Vec3D> my_loads; //non-workers send there loads and indices to worker_proc
  vector<double> my_areas;
  vector<int> my_indices;
  if(mpi_rank!=worker_proc) { //worker_proc should NOT send!
    for(auto&& node : my_indices_to_send) {
      my_indices.push_back(node);
      my_areas.push_back(An[node]); //could be 0 if calculate_An is false --- ok.
      my_loads.push_back(Fs[node]);
    }
  }

  // Processor 0 gets data from others and add them to Fs and An
  vector<Vec3D> other_loads;
  vector<double> other_areas;
  vector<int> other_indices;
  CommunicationTools::GatherVector<Vec3D>(comm, worker_proc, my_loads, &other_loads);
  CommunicationTools::GatherVector<double>(comm, worker_proc, my_areas, &other_areas);
  CommunicationTools::GatherVector<int>(comm, worker_proc, my_indices, &other_indices);

  if(mpi_rank==worker_proc) {
    assert(other_loads.size() == other_areas.size() &&
           other_loads.size() == other_indices.size());
    int node;
    for(size_t i=0; i<other_loads.size(); i++) {
      node = other_indices[i];
      Fs[node] += other_loads[i];
      An[node] += other_areas[i];
    } 
    for(size_t i=0; i<Fs.size(); i++)
      FAs[i] = An[i]==0.0 ? 0.0 : Fs[i]/An[i];
  }


  MPI_Barrier(comm);
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

  double size = global_mesh_ptr->GetMinDXYZ(ijk);

  return factor*size;
}

//------------------------------------------------------------------------------------------------
// Calculates the one-sided traction from the side indicated by "normal"
Vec3D
EmbeddedBoundaryOperator::CalculateTractionAtPoint(Vec3D &p, Vec3D &normal, Vec5D*** v, double*** id)
{
  //int mpi_rank;
  //MPI_Comm_rank(comm, &mpi_rank);


  Int3 ijk0(INT_MAX);
  Vec3D xi;
  global_mesh_ptr->FindElementCoveringPoint(p, ijk0, &xi, true);

  int i,j,k;

  // find which nodes in the element are on the correct side given by "normal"
  // first, we move the Gauss point (p) along the normal direction outside "wall_thickness"
  double loft = 0.0;
  for(auto&& xter : intersector)
    loft = std::max(loft, xter->GetSurfaceHalfThickness());
  loft *= 2.0; //moving out of half_thickness

  int iter, max_iter = 10;
  bool sameside[2][2][2];
  bool found_sameside;
  for(iter=0; iter<max_iter; iter++) {//gradually increase "loft", if necessary

    Vec3D ref_point = p + loft*normal;

    found_sameside = false;
    for(int dk=0; dk<=1; dk++)
      for(int dj=0; dj<=1; dj++)
        for(int di=0; di<=1; di++) {
          i = ijk0[0] + di;
          j = ijk0[1] + dj;
          k = ijk0[2] + dk;

          assert(coordinates_ptr->IsHere(i,j,k,true)); //Not in the ghosted subdomain? Something is wrong

          if(coordinates_ptr->OutsidePhysicalDomain(i,j,k)) { //state variable unavailable here.
            sameside[dk][dj][di] = false;
            continue;
          } //NOTE: Must not check (external) ghost nodes, even if the state variable is correct (w/ b.c.).
            //      This is because the surface may not extend into the ghost layer. Hence, the intersection
            //      function may not work as expected!

          if(id[k][j][i] == INACTIVE_MATERIAL_ID) {
            sameside[dk][dj][di] = false;
            continue;
          }


          Vec3D x(global_mesh_ptr->GetX(i), global_mesh_ptr->GetY(j), global_mesh_ptr->GetZ(k));
          sameside[dk][dj][di] = true; //start w/ true
          for(auto&& xter : intersector) {
            if(xter->Intersects(x, ref_point, NULL, false, 1)) { //ignore dropped elems due to multi_intersector
              sameside[dk][dj][di] = false;
              break;
            }
          }

          if(sameside[dk][dj][di])
            found_sameside = true; //at least found one on the "same side"
        }

    if(found_sameside)
      break;
    else
      loft *= 2.0; //increase lofting distance
  }

  if(iter>=5 && verbose>=1) {
    if(found_sameside)
      fprintf(stdout,"\033[0;35mWarning: Applied a lofting height of %e (iter=%d) to find valid nodes for interpolating \n"
                               "         pressure at Gauss point (%e, %e, %e).\033[0m\n",
              loft, iter, p[0], p[1], p[2]);
    //if found_sameside == false, will trigger another warning message below.
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

        //fprintf(stdout,"[%d] side = %d: (%d,%d,%d), p = %e.\n", mpi_rank, side, i,j,k, v[k][j][i][4]);
        pressure[dk][dj][di] = v[k][j][i][4]; //get pressure
        total_pressure += pressure[dk][dj][di];
        n_pressure++;
      }

  double avg_pressure;
  if(n_pressure==0) {
    fprintf(stdout,"\033[0;35mWarning: No valid active nodes for interpolating pressure at "
                   "Gauss point (%e, %e, %e). Try increasing lofting height.\033[0m\n",
                   p[0], p[1], p[2]);
    avg_pressure = 0.0;
  } else
    avg_pressure = total_pressure/n_pressure;

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

int
EmbeddedBoundaryOperator::CombineSharedGaussPointData(vector<double>& data, vector<double>& shared_data)
{
  if(shared_data.empty())
    return 0;

  if(shared_data.size()%3!=0 || data.size()%5!=0) {
    fprintf(stdout,"\033[0;31m*** Error: Encountered issue in load transfer (2D->3D). "
                   "Lofting may be too large.\n\033[0m");
    exit(-1);
  }

  int Nshared = shared_data.size()/3;

  double tolerance = std::max(1.0e-12, 1.0e-13*std::min(fabs(shared_data[0]), fabs(shared_data[1])));

  int counter = 0;

  for(int i=0; i<Nshared; i++) {

    bool already_checked = false;
    for(int j=0; j<i; j++) {
      if(std::max(fabs(shared_data[3*i]   - shared_data[3*j]  ),
                  fabs(shared_data[3*i+1] - shared_data[3*j+1])) < tolerance
         && fabs(shared_data[3*i+2] - shared_data[3*j+2]) < tolerance) {
        already_checked = true;
        break;
      }
    }
    if(already_checked)
      continue;


    assert(data.size()%5==0);
    int N = data.size()/5;
    int first_index = -1;
    int dup_index = -1;
    for(int j=0; j<N; j++) {
      if(std::max(fabs(shared_data[3*i]   - data[5*j]  ),
                  fabs(shared_data[3*i+1] - data[5*j+1])) < tolerance
         && fabs(shared_data[3*i+2] - data[5*j+2]) < tolerance) {
        if(first_index==-1)
          first_index = j;
        else {
          data[5*first_index+3] += data[5*j+3];
          data[5*first_index+4] += data[5*j+4];
          dup_index = j;
          break;
        }
      }
    }

    if(first_index<0 || dup_index<0) {
      fprintf(stdout,"\033[0;31m*** Error: Encountered issue in load transfer (2D->3D). "
                     "Lofting may be too large.\n\033[0m");
      exit(-1);
    }

    data.erase(data.begin()+5*dup_index, data.begin()+5*dup_index+5); //erase 5 elements
    
    counter++;
  }

  assert(counter*2==Nshared);

  return counter;

}

//------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------





//------------------------------------------------------------------------------------------------




//------------------------------------------------------------------------------------------------



