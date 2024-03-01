/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef CRACKINGSURFACE_H_
#define CRACKINGSURFACE_H_

#include <map>
#include <set>
#include <fstream>
#include <vector>

struct Int3;

/*****************************************************************************************
 * Classes for handling FSI w/ fracture (copied from AERO-F, written by KW in grad school)
 ****************************************************************************************/

struct PhantomElement {
  int nNodes;
  double *phi;
  int *nodes;

  // constructors
  PhantomElement(): nNodes(-1), phi(0), nodes(0) {}
  PhantomElement(int, int *, double *);
  PhantomElement(int, int, int, int, double, double, double, double);
  // destructor
  ~PhantomElement() {
    if(phi) {
      delete[] phi;
    }
    if(nodes) {
      delete[] nodes;
    }
  }

  // update nodes
  void update(int *, double *); // KW: Doesn't have to update both nodes and phi. Use "NULL".

  void writeCrackingData(std::ofstream&) const;
  static PhantomElement *readCrackingData(std::ifstream&);
};

//------------------------------------------------------------------------------

struct LatestCracking {
  std::set<int> phantomQuads;
  std::map<int, int> phantomNodes; // Note: "phantomNodes" are NOT equivalent to "nodes of phantomQuads"!
};

//------------------------------------------------------------------------------

class CrackingSurface {
  const int elemType; // currently only support quadrangles.
  int nTotalQuads, nUsedQuads;
  int nTotalTrias, nUsedTrias;
  int nTotalNodes, nUsedNodes;

  std::map<int, PhantomElement *> phantoms; // size: number of cracked (quad) elements
  LatestCracking latest;
  bool gotNewCracking;
  int (*tria2quad)[2]; // size: nTotalTrias
  int (*quad2tria)[2]; // size: nTotalQuads
  bool *cracked; // size: nTotalQuads
  bool *deleted; // size: nTotalQuads, in the case of Element Deletion, a "cracked" element is "deleted".

  // For cracking simulations, this contains a map from a set of purely undeleted
  // triangles to the full list of triangles.
  int *triangle_id_map;
  int numRealTriangles;

  void constructTriangleMap();

 public:
  CrackingSurface(int, int, int, int, int);
  ~CrackingSurface();

  // called by EmbeddedStructure only!
  int splitQuads(int *, int, int *);
  int updateCracking(int, int, int *, double *, int *, std::vector<Int3> &, int, int *, int);

  int numCrackedElements() {
    return int(phantoms.size());
  }
  bool hasCracked(int);
  double getPhi(int, double, double, bool * = 0, bool = false);

  double getPhiPhysBAM(int, double, double, bool * = 0, bool = false);

  bool purelyPhantom(int);
  bool purelyPhantomPhysBAM(int);

  bool getNewCrackingFlag() const {
    return gotNewCracking;
  }
  void setNewCrackingFlag(bool flag) {
    gotNewCracking = flag;
  }

  int totNodes() const {
    return nTotalNodes;
  }
  int usedNodes() const {
    return nUsedNodes;
  }
  int totTrias() const {
    return nTotalTrias;
  }
  int usedTrias() const {
    return nUsedTrias;
  }
  std::set<int> getLatestPhantomQuads() const {
    return latest.phantomQuads;
  }
  std::map<int, int> getLatestPhantomNodes() const {
    return latest.phantomNodes;
  }
  void getQuad2Tria(int quad, int& trId1, int& trId2) {
    trId1 = quad2tria[quad][0];
    trId2 = quad2tria[quad][1];
  }

  // for debug
  void printInfo(char *);

  void writeCrackingData(std::ofstream&) const;
  void readCrackingData(std::ifstream&);

  // see triangle_id_map
  int mapTriangleID(int);

  int numberRealTriangles() {
    if(triangle_id_map) {
      return numRealTriangles;
    }
    else {
      return nUsedTrias;
    }
  }
};

//------------------------------------------------------------------------------
#endif
