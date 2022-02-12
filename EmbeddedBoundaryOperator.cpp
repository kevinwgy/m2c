#include<EmbeddedBoundaryOperator.h>
#include<deque>
#include<list>
#include<utility>

using std::string;
using std::map;
using std::vector;

//------------------------------------------------------------------------------------------------

EmbeddedBoundaryOperator::EmbeddedBoundaryOperator(IoData &iod_, bool surface_from_other_solver) 
                        : iod(iod_), hasSurfFromOtherSolver(surface_from_other_solver)
{
  // count surfaces from files
  int counter = iod.embed_surfaces.surfaces.dataMap.size();

  // get the correct size of surfaces and F
  if(surface_from_other_solver) {
    surfaces.resize(counter+1);
    F.resize(counter+1);
  } else {
    surfaces.resize(counter);
    F.resize(counter);
  }

  // read surfaces from files
  for(auto it = iod.embed_surfaces.surfaces.dataMap.begin();
           it != iod.embed_surfaces.surfaces.dataMap.end(); it++) {

    int index = it->first;
    if(index<0 || index>=counter) {
      print_error("*** Error: Detected error in the indices of embedded surfaces (id = %d).\n", index);
      exit_mpi();
    }
    if(surface_from_other_solver)
      index++;

    // currently, these parameters are not used in M2C.
    map<int, EmbeddedSurfaceData::Type> &boundaryConditionsMap;
    string nodeSetName; 
    vector<int> faceIDs;
    int *surfaceID = NULL;
    int *faceID = NULL;
    int numStNodes(0), numStElems(0), maxFaceID(0); 
    
    ReadMeshFile(it->second->filename, nodeSetName, boundaryConditionsMap, faceIDs, numStNodes,
                 numStElems, surface[index].X, surfaceID, surface[index].elems, faceID, maxFaceID);

    surface[index].X0 = surface[index].X;

    // delete
    if(surfaceID) delete [] surfaceID;
    if(faceID)    delete [] faceID; 

  }

  print("- Activted the Embedded Boundary Method. Detected %d surfaces (%d from concurrent programs).\n",
        surfaces.size(), surfaces.size() - counter); 

}

//------------------------------------------------------------------------------------------------

EmbeddedBoundaryOperator::~EmbeddedBoundaryOperator()
{ }

//------------------------------------------------------------------------------------------------

void
EmbeddedBoundaryOperator::Destroy()
{ }

//------------------------------------------------------------------------------------------------

// copied from AERO-F. Originally written by KW, extended by other students
void
EmbeddedBoundaryOperator::ReadMeshFile(char *filename, string& nodeSetName,
                                map<int, BoundaryData::Type>& boundaryConditionsMap, vector<int>& faceIDs,
                                int& numStNodes, int& numStElems, vector<Vec3D> &Xs, int *&surfaceID, vector<Int3> &stElem,
                                int *&faceID, int& maxFaceID) 
{
  // read data from the solid surface input file.
  FILE *topFile;
  topFile = fopen(filename, "r");
  if(topFile == NULL) {
    print_error("*** Error: embedded structure surface mesh doesn't exist (%s).\n", filename);
    exit_mpi();
  }
 
  int MAXLINE = 500;
  char line[MAXLINE], key1[MAXLINE], key2[MAXLINE], copyForType[MAXLINE];

  // load the nodes and initialize all node-based variables.
  // load solid nodes at t=0
  int num0 = 0;
  int num1 = 0;
  double x1, x2, x3;
  int node1, node2, node3;
  int type_read = 0;
  int surfaceid = 0;
  std::deque<std::pair<int, Vec3D>> nodeList;
  std::deque<std::array<int, 4>> elemList; // surface ID + three node IDs
  int maxIndex = 0, maxElem = -1;

  while(fgets(line, MAXLINE, topFile) != 0) {

    sscanf(line, "%s", key1);
    bool skip = false;
    if(strcmp(key1, "Nodes") == 0) {
      sscanf(line, "%*s %s", key2);
      nodeSetName = std::string(key2);
      skip = true;
      type_read = 1;
    }
    else if(strcmp(key1, "Elements") == 0) {
      sscanf(line, "%*s %s", key2);
      skip = true;
      type_read = 2;
      int underscore_pos = -1;
      int k = 0;
      while((key2[k] != '\0') && (k < MAXLINE)) {
        if(key2[k] == '_') {
          underscore_pos = k;
        }
        k++;
      }
      if(underscore_pos > -1) {
        sscanf(key2 + (underscore_pos + 1), "%d", &surfaceid);
        // now we look for keywords for the type of structure
        strcpy(copyForType, key2);
        int l = 0;
        while((copyForType[l] != '\0') && (l < MAXLINE)) {
          copyForType[l] = (char)std::tolower(static_cast<unsigned char>(copyForType[l]));
          l++;
        }
        // read the name of the file and detects keyword for type
        if(strstr(copyForType, "symmetry") != NULL) {
          if(boundaryConditionsMap.count(surfaceid) != 0) {
            print_error("*** Error: two embedded surfaces have the same id (%d)\n", surfaceid);
            exit_mpi();
          }
          boundaryConditionsMap[surfaceid] = EmbeddedSurfaceData::Symmetry;
        }
        else if(strstr(copyForType, "porouswall") != NULL) {
          if(boundaryConditionsMap.count(surfaceid) != 0) {
            print_error("*** Error: two embedded surfaces have the same id (%d)\n", surfaceid);
            exit_mpi();
          }
          boundaryConditionsMap[surfaceid] = BoundaryData::PorousWall;
        }
        else {
          faceIDs.push_back(surfaceid);
        }
      }
    }
    if(!skip) { // we are reading a node or an element (not a header)
      if(type_read == 1) {
        sscanf(line, "%d %lf %lf %lf", &num1, &x1, &x2, &x3);
        if(num1 < 1) {
          print_error("*** Error: detected a node with index %d in the embedded surface file!\n", num1);
          exit_mpi();
        }
        if(num1 > maxIndex) {
          maxIndex = num1;
        }
        nodeList.push_back({num1, {x1, x2, x3}});
      }
      if(type_read == 2) { // we are reading an element
        sscanf(line, "%d %d %d %d %d", &num0, &num1, &node1, &node2, &node3);
        elemList.push_back({surfaceid, node1 - 1, node2 - 1, node3 - 1});
      }
    }
  }
  numStNodes = int(nodeList.size());
  if(numStNodes != maxIndex) {
    print_error("*** Error: the node set of the embedded surface has gaps: max index = %d, number of nodes = %d.",
                 maxIndex, numStNodes);
    exit_mpi();
  }
  numStElems = int(elemList.size());
  // feed data to Xs
  Xs.resize(numStNodes);
  surfaceID = new int[numStNodes];
  for(auto it1 = nodeList.begin(); it1 != nodeList.end(); it1++) {
    Xs[it1->first - 1][0] = it1->second[0];
    Xs[it1->first - 1][1] = it1->second[1];
    Xs[it1->first - 1][2] = it1->second[2];
    surfaceID[it1->first - 1] = 0;
  }
  nodeList.clear();
  stElem.resize(numStElems);
  faceID = new int[numStElems];
  maxFaceID = 0;
  auto it2 = elemList.begin();
  for(int i = 0; i < numStElems; i++) {
    stElem[i][0] = (*it2)[1];
    stElem[i][1] = (*it2)[2];
    stElem[i][2] = (*it2)[3];
    faceID[i] = (*it2)[0]; // give the face the ID of the element
    maxFaceID = std::max(maxFaceID, faceID[i]);
    surfaceID[(*it2)[1]] = (*it2)[0]; // gives every node the ID of the element
    surfaceID[(*it2)[2]] = (*it2)[0];
    surfaceID[(*it2)[3]] = (*it2)[0];
    it2++;
  }
  elemList.clear();
  fclose(topFile);
}

//------------------------------------------------------------------------------------------------

void
EmbeddedBoundaryOperator::ComputeForces(SpaceVariable3D &V, SpaceVariable3D &ID)
{
  print_warning("Warning: EmbeddedBoundaryOperator::ComputeForces has not been implemented.\n");
}

//------------------------------------------------------------------------------------------------

void
EmbeddedBoundaryOperator::TrackSurfaces()
{
  print_warning("Warning: EmbeddedBoundaryOperator::TrackSurfaces has not been implemented.\n");
}

//------------------------------------------------------------------------------------------------

void
EmbeddedBoundaryOperator::TrackUpdatedSurfaceFromOtherSolver()
{
  print_warning("Warning: EmbeddedBoundaryOperator::TrackUpdatedSurfaceFromOtherSolver has not been implemented.\n");
}

//------------------------------------------------------------------------------------------------




//------------------------------------------------------------------------------------------------









