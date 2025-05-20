/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <TriangulatedSurface.h>
#include <map>
#include <cassert>
using std::set;
using std::vector;
using std::map;
using std::pair;

//----------------------------------------------------

TriangulatedSurface&
TriangulatedSurface::operator=(const TriangulatedSurface &S2)
{
  degenerate = S2.degenerate;

  X0    = S2.X0;
  X     = S2.X;
  Udot  = S2.Udot;
  elems = S2.elems;
  
  active_nodes = S2.active_nodes;
  active_elems = S2.active_elems;

  nodetag = S2.nodetag;
  elemtag = S2.elemtag;

  elemNorm = S2.elemNorm;
  elemArea = S2.elemArea;

  node2node = S2.node2node;
  node2elem = S2.node2elem;
  elem2elem = S2.elem2elem;

  return *this;
}

//----------------------------------------------------

void
TriangulatedSurface::Append(const TriangulatedSurface &S2)
{
  // sanity checks
  assert(degenerate == S2.degenerate);
  assert(X0.empty() == S2.X0.empty());
  assert(X.empty() == S2.X.empty());
  assert(Udot.empty() == S2.Udot.empty());
  assert(elems.empty() == S2.elems.empty());
  assert(elemNorm.empty() == S2.elemNorm.empty()); //either both are empty, or both filled
  assert(elemArea.empty() == S2.elemArea.empty());
  assert(node2node.empty() == S2.node2node.empty());
  assert(node2elem.empty() == S2.node2elem.empty());
  assert(elem2elem.empty() == S2.elem2elem.empty());

  int nNodes = X.size();
  int nElems = elems.size(); //before appending

  X0.insert(X0.end(), S2.X0.begin(), S2.X0.end());
  X.insert(X.end(), S2.X.begin(), S2.X.end());
  Udot.insert(Udot.end(), S2.Udot.begin(), S2.Udot.end());

  for(auto&& e : S2.elems)
    elems.push_back(Int3(e[0]+nNodes, e[1]+nNodes, e[2]+nNodes)); //shift index

  active_nodes += S2.active_nodes;
  active_elems += S2.active_elems;

  nodetag.insert(nodetag.end(), S2.nodetag.begin(), S2.nodetag.end());
  elemtag.insert(elemtag.end(), S2.elemtag.begin(), S2.elemtag.end());

  elemNorm.insert(elemNorm.end(), S2.elemNorm.begin(), S2.elemNorm.end());
  elemArea.insert(elemArea.end(), S2.elemArea.begin(), S2.elemArea.end());

  for(auto&& n2n : S2.node2node) {
    node2node.push_back(set<int>());
    for(auto&& nid : n2n)
      node2node.back().insert(nid + nNodes); //shift index
  }
  for(auto&& n2e : S2.node2elem) {
    node2elem.push_back(set<int>());
    for(auto&& eid : n2e)
      node2elem.back().insert(eid + nElems); //shift index
  }
  for(auto&& e2e : S2.elem2elem) {
    elem2elem.push_back(set<int>());
    for(auto&& eid : e2e)
      elem2elem.back().insert(eid + nElems); //shift index
  }
  
}

//----------------------------------------------------

void TriangulatedSurface::BuildConnectivities()
{
  int nNodes = X.size();
  int nElems = elems.size();

  node2node.clear();
  node2node.resize(nNodes);
  node2elem.clear();
  node2elem.resize(nNodes);
  
  for(int i=0; i<nElems; i++)
    for(int j=0; j<3; j++) {
      node2node[elems[i][j]].insert(elems[i][(j+1)%3]);
      node2node[elems[i][j]].insert(elems[i][(j+2)%3]);
      node2elem[elems[i][j]].insert(i);
    }

  elem2elem.clear();
  elem2elem.resize(nElems);
  for(int i=0; i<nElems; i++)
    for(int j=0; j<3; j++)
      for(auto it = node2elem[elems[i][j]].begin(); it != node2elem[elems[i][j]].end(); it++)
        elem2elem[i].insert(*it); //no risk of duplicates, as elem2elem[i] is a "set"

}

//----------------------------------------------------

void TriangulatedSurface::CalculateNormalsAndAreas()
{
  int nElems = elems.size();

  elemNorm.clear();
  elemNorm.resize(nElems);

  elemArea.clear();
  elemArea.resize(nElems);

  int n1, n2, n3;
  Vec3D dx2, dx3;

  // Also look to determine a point inside the solid but away from the structure.
  for (int i=0; i<nElems; i++) {
    n1 = elems[i][0];
    n2 = elems[i][1];
    n3 = elems[i][2];

    dx2 = X[n2] - X[n1];
    if(n2==n3) {//assuming this is a set of line segments in 2D (x-y)
      degenerate = true;
      elemNorm[i][0] = -dx2[1];
      elemNorm[i][1] = dx2[0];
      elemNorm[i][2] = 0.0;
    } else {
      dx3 = X[n3] - X[n1];
      elemNorm[i] = dx2^dx3; // cross product ==> normal dir.
    }

    double nrm = elemNorm[i].norm();
    // normalize the normal.
    if(nrm > 0.0) {
      elemNorm[i] /= nrm;
      elemArea[i] = 0.5*nrm;
    } else {
      fprintf(stdout, "*** Error: area (length) of triangle (edge) %d is %e.\n", i+1, nrm);
      exit(-1);
    }
  }

  //verify the degenerate case
  if(degenerate) {
    for(int i=0; i<nElems; i++)
      if(n2!=n3) {
        fprintf(stdout, "*** Error: Cannot handle a mix of triangles and line segments.\n");
        exit(-1);
      }
    for(int i=0; i<(int)X.size(); i++)
      if(X[i][2] != 0) {
        fprintf(stdout, "*** Error: A degenerated triangulated surface must be on the x-y plane. Found z = %e.\n", 
                     X[i][2]);
        exit(-1);
      }
  }
}

//----------------------------------------------------
// TODO: This function does not support fracture (i.e. active_nodes...)
typedef pair<int, int> iipair;
typedef pair<int, bool> ibpair;
typedef pair<iipair, ibpair> EdgePair;
EdgePair makeEdgePair(int node1, int node2, int triangleNumber) {
  if(node1 < node2)
    return EdgePair(iipair(node1, node2), ibpair(triangleNumber, true));
  else
    return EdgePair(iipair(node2, node1), ibpair(triangleNumber, false));
}

// Check whether all the elements have consistent orientation.
bool
TriangulatedSurface::CheckSurfaceOrientation()
{
  assert(!degenerate);

  map<iipair, ibpair> edgeMap;
  for (int iTriangle=0; iTriangle<(int)elems.size(); iTriangle++) {
    int from1, to1, from2, to2, from3, to3;
    from1 = elems[iTriangle][0];  to1 = elems[iTriangle][1];
    from2 = elems[iTriangle][1];  to2 = elems[iTriangle][2];
    from3 = elems[iTriangle][2];  to3 = elems[iTriangle][0];

    EdgePair ep[3];
    ep[0] = makeEdgePair(from1, to1, iTriangle);
    ep[1] = makeEdgePair(from2, to2, iTriangle);
    ep[2] = makeEdgePair(from3, to3, iTriangle);

    for(int i=0; i < 3; ++i) {
      map<iipair, ibpair>::iterator it = edgeMap.find(ep[i].first);
      if(it != edgeMap.end()) { // we found this edge
         if(it->second.second == ep[i].second.second)
           return false;
      } else
        edgeMap[ep[i].first] = ep[i].second;
    }
  }
  return true;
}

//----------------------------------------------------
// Checks whether the elements have consistent orientation,
// AND the surface is closed.
bool
TriangulatedSurface::CheckSurfaceOrientationAndClosedness()
{
  assert(!degenerate);

  if(!CheckSurfaceOrientation())
    return false;

  std::multimap<int,int> edges;
  int p1, p2;
  for(int i=0; i<(int)elems.size(); i++) {
    for(int j=0; j<3; j++) {//loop through the three edges
      p1 = elems[i][j];
      p2 = elems[i][(j+1)%3];
      // is p2->p1 in "edges"? if yes, cancel it. otherwise, add p1->p2 to "edges"
      auto res = edges.equal_range(p2);
      bool found = false;
      for(auto it = res.first; it != res.second; it++) {
        if(it->second == p1) {
          edges.erase(it);
          found = true;
          break; 
        }
      } 
      if(!found)
        edges.insert(pair<int,int>(p1,p2));
    }
  }
  return edges.empty();
}

//----------------------------------------------------
// Checks whether the surface is closed.
bool
TriangulatedSurface::CheckSurfaceClosedness()
{
  assert(!degenerate);

  std::multimap<int,int> edges;
  int p1, p2;
  for(int i=0; i<(int)elems.size(); i++) {
    for(int j=0; j<3; j++) {//loop through the three edges
      p1 = elems[i][j];
      p2 = elems[i][(j+1)%3];
      // is p2->p1 or p1->p2 in "edges"? if yes, cancel it. otherwise, add p1->p2 to "edges"
      auto res = edges.equal_range(p2);
      bool found = false;
      for(auto it = res.first; it != res.second; it++) {
        if(it->second == p1) {
          edges.erase(it);
          found = true;
          break; 
        }
      } 
      if(!found) {
        auto res2 = edges.equal_range(p1);
        for(auto it = res2.first; it != res2.second; it++) {
          if(it->second == p2) {
            edges.erase(it);
            found = true;
            break; 
          }
        }
        if(!found)
          edges.insert(pair<int,int>(p1,p2));
      }
    }
  }

  return edges.empty();
}

//----------------------------------------------------

