/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <CrackingSurface.h>
#include <Vector3D.h>

#include <cstdlib>
#include <fstream>
#include <iostream>

using std::map;
using std::set;

//------------------------------------------------------------------------------

PhantomElement::PhantomElement(int n, int *nod, double *ph): nNodes(n) {
  phi = new double[n];
  nodes = new int[n];
  for(int i = 0; i < n; i++) {
    phi[i] = ph[i];
    nodes[i] = nod[i];
  }
}

//------------------------------------------------------------------------------

PhantomElement::PhantomElement(int a, int b, int c, int d,
                               double phia, double phib, double phic, double phid): nNodes(4) {
  phi = new double[4];
  nodes = new int[4];
  phi[0] = phia;
  phi[1] = phib;
  phi[2] = phic;
  phi[3] = phid;
  nodes[0] = a;
  nodes[1] = b;
  nodes[2] = c;
  nodes[3] = d;
}

//------------------------------------------------------------------------------

void PhantomElement::update(int *nod, double *ph) {
  if(nod)
    for(int i = 0; i < nNodes; i++) {
      nodes[i] = nod[i];
    }
  if(ph)
    for(int i = 0; i < nNodes; i++) {
      phi[i] = ph[i];
    }
}

//------------------------------------------------------------------------------

void PhantomElement::writeCrackingData(std::ofstream& restart_file) const {
  restart_file.write(reinterpret_cast<const char *>(&nNodes), sizeof(int));
  restart_file.write(reinterpret_cast<const char *>(phi), sizeof(double)*nNodes);
  restart_file.write(reinterpret_cast<const char *>(nodes), sizeof(int)*nNodes);
}

//------------------------------------------------------------------------------

PhantomElement *PhantomElement::readCrackingData(std::ifstream& restart_file) {
  PhantomElement *E = new PhantomElement;
  restart_file.read(reinterpret_cast<char *>(&E->nNodes), sizeof(int));
  E->phi = new double[E->nNodes];
  E->nodes = new int[E->nNodes];
  restart_file.read(reinterpret_cast<char *>(E->phi), sizeof(double)*E->nNodes);
  restart_file.read(reinterpret_cast<char *>(E->nodes), sizeof(int)*E->nNodes);
  return E;
}

//------------------------------------------------------------------------------

CrackingSurface::CrackingSurface(int eType, int nUsed, int nTotal, int nUsedNd, int nTotNodes): elemType(eType) {
  if(eType != 4) {
    fprintf(stdout, "*** Error: ElemType(%d) is not supported for CrackingSurface!\n", eType);
    exit(1);
  }
  nTotalNodes = nTotNodes;
  nUsedNodes = nUsedNd;
  nTotalQuads = nTotal;
  nUsedQuads = nUsed;
  nTotalTrias = nTotal * 2,  nUsedTrias = 0; // surface not split yet!
  tria2quad = new int[nTotalTrias][2];
  quad2tria = new int[nTotalQuads][2];
  cracked = new bool[nTotalQuads];
  deleted = new bool[nTotalQuads];
  for(int i = 0; i < nTotalTrias; i++) {
    tria2quad[i][0] = tria2quad[i][1] = -1;
  }
  for(int i = 0; i < nTotalQuads; i++) {
    cracked[i] = false;
    deleted[i] = false;
    quad2tria[i][0] = quad2tria[i][1] = -1;
  }
  gotNewCracking = false;
  triangle_id_map = NULL;
  numRealTriangles = 0;
}

//------------------------------------------------------------------------------

CrackingSurface::~CrackingSurface() {
  for(map<int, PhantomElement *>::iterator it = phantoms.begin(); it != phantoms.end(); it++) {
    delete it->second;
  }
  if(tria2quad) {
    delete[] tria2quad;
  }
  if(cracked) {
    delete[] cracked;
  }
  if(deleted) {
    delete[] deleted;
  }
  if(triangle_id_map) {
    delete [] triangle_id_map;
  }
}

//------------------------------------------------------------------------------

int CrackingSurface::splitQuads(int *quadTopo, int nQuads, int *triaTopo) {
  if(nQuads != nUsedQuads) {
    fprintf(stdout, "Software bug in CrackingSurface::splitQuads!\n");
    exit(-1);
  }
  int count = 0;
  for(int i = 0; i < nQuads; i++) {
    triaTopo[count*3]   = quadTopo[i * 4];
    triaTopo[count*3+1] = quadTopo[i * 4 + 1];
    triaTopo[count*3+2] = quadTopo[i * 4 + 2];
    tria2quad[count][0] = i;
    tria2quad[count][1] = 0;
    quad2tria[i][0] = count;
    count++;
    if(quadTopo[i * 4 + 2] == quadTopo[i * 4 + 3]) { // this is a degenerated triangle
      quad2tria[i][1] = -1;
      continue;
    }
    triaTopo[count*3]   = quadTopo[i * 4];
    triaTopo[count*3+1] = quadTopo[i * 4 + 2];
    triaTopo[count*3+2] = quadTopo[i * 4 + 3];
    tria2quad[count][0] = i;
    tria2quad[count][1] = 1;
    quad2tria[i][1] = count;
    count++;
  }
  nUsedTrias = count;
  for(int i = nQuads; i < nTotalQuads; i++) {
    tria2quad[count][0] = i;
    tria2quad[count][1] = 0;
    quad2tria[i][0] = count;
    count++;
    tria2quad[count][0] = i;
    tria2quad[count][1] = 1;
    quad2tria[i][1] = count;
    count++;
  }
  for(int i = nUsedTrias; i < nTotalTrias; i++) {
    triaTopo[i*3] = triaTopo[i*3+1] = triaTopo[i*3+2] = -1;
  }
  return nUsedTrias;
}

//------------------------------------------------------------------------------

int CrackingSurface::updateCracking(int numConnUpdate, int numLSUpdate, int *connUpdate, double *phi,
                                    int *phiIndex, std::vector<Int3> &triaTopo, int nUsedNd, int *new2old, int numNewNodes) {
  if(numConnUpdate == 0) {
    return 0;
  }
  if(gotNewCracking) {
    fprintf(stdout, "*** Warning: last cracking update hasn't been pushed to intersector!\n");
  }
  gotNewCracking = true;
  latest.phantomQuads.clear();
  latest.phantomNodes.clear();
  int quadId, nNew = 0, maxtrId = nUsedTrias - 1;
  int maxQuad = nUsedQuads - 1, nNewQuad = 0;
  // build a reverse map of phiIndex. The cost is trivial since phiIndex is very short.
  map<int, int> temp;
  map<int, int>::iterator itor;
  for(int i = 0; i < numLSUpdate; i++) {
    temp[phiIndex[i]] = i;
  }
  for(int i = 0; i < numConnUpdate; i++) {
    quadId = connUpdate[5 * i];
    if(quadId > maxQuad) {
      maxQuad = quadId;
    }
    bool already_cracked = cracked[quadId];
    cracked[quadId] = true;
    latest.phantomQuads.insert(quadId);
    // insert a new phantom element or modify an existing one
    itor = temp.find(quadId);
    if(already_cracked) {
      if(itor != temp.end()) {
        fprintf(stdout, "*** Warning: Shouldn't modify the levelset of a previously-cracked element!\n");
      }
      phantoms[quadId]->update(&(connUpdate[5 * i + 1]), itor == temp.end() ? NULL : & (phi[itor->second]));
    }
    else {
      if(itor == temp.end()) {
        fprintf(stdout, "*** Error: Need the level-set for a new cracked element!\n");
        exit(-1);
      }
      int ind = itor->second;
      if(phi[4 * ind] < 0.0 && phi[4 * ind + 1] < 0.0 && phi[4 * ind + 2] < 0.0 && phi[4 * ind + 3] < 0.0) {
        deleted[quadId] = true;
      }
      phantoms[quadId] = new PhantomElement(connUpdate[5 * i + 1], connUpdate[5 * i + 2], connUpdate[5 * i + 3], connUpdate[5 * i + 4],
                                            phi[4 * ind], phi[4 * ind + 1], phi[4 * ind + 2], phi[4 * ind + 3]);
    }
    // modify the triangle mesh connectivity
    int trId1, trId2;
    if(quadId >= nUsedQuads) { // this is a new quad
      nNewQuad++;
      if(quadId > nUsedQuads + numConnUpdate / 2) {
        fprintf(stdout, "*** Error: nUsed = %d, newConn/2 = %d, currentId = %d!\n", nUsedQuads, numConnUpdate / 2, quadId + 1);
        exit(-1);
      }
      nNew += 2;
      trId1 = nUsedTrias + 2 * (quadId - nUsedQuads);
      trId2 = trId1 + 1;
      if(maxtrId < trId2) {
        maxtrId = trId2;
      }
    }
    else { // this is an existing quad
      trId1 = quad2tria[quadId][0];
      trId2 = quad2tria[quadId][1];
    }
    triaTopo[trId1][0] = connUpdate[5 * i + 1];
    triaTopo[trId1][1] = connUpdate[5 * i + 2];
    triaTopo[trId1][2] = connUpdate[5 * i + 3];
    // Check to see if the second quad exists.  If the element deletion is done with triangles,
    // then the quad is a degenerated triangle and the second triangle has id -1
    if(trId2 >= 0) {
      triaTopo[trId2][0] = connUpdate[5 * i + 1];
      triaTopo[trId2][1] = connUpdate[5 * i + 3];
      triaTopo[trId2][2] = connUpdate[5 * i + 4];
    }
  }
  // construct latest.phantomNodes
  for(int i = 0; i < numNewNodes; i++) {
    latest.phantomNodes[new2old[2 * i]] = new2old[2 * i + 1];
  }
  if(maxQuad + 1 != nUsedQuads + nNewQuad) {
    fprintf(stdout, "*** Error: Inconsistency in the number of structure quad elements (%d %d %d)! (Could be a software bug.)\n", maxQuad, nUsedQuads,
            nNewQuad);
    exit(-1);
  }
  nUsedQuads += nNewQuad;
  if(maxtrId + 1 != nUsedTrias + nNew) {
    fprintf(stdout, "SOFTWARE BUG: Violated the ordering of new elements (%d v.s. %d)\n", maxtrId + 1, nUsedTrias + 2 * nNew);
    exit(-1);
  }
  nUsedTrias += nNew;
  nUsedNodes = nUsedNd;
  constructTriangleMap();
  return nNew;
}

//------------------------------------------------------------------------------

bool CrackingSurface::hasCracked(int trId) {
  if(trId >= nUsedTrias) {
    fprintf(stdout, "*** Error: Unable to access Triangle %d of the embedded surface(%d)!\n", trId + 1, nUsedTrias);
    exit(-1);
  }
  if(cracked[tria2quad[trId][0]]) {
    return true;
  }
  else {
    return false;
  }
}

//------------------------------------------------------------------------------

double CrackingSurface::getPhi(int trId, double xi1, double xi2, bool *hasCracked, bool debug) {
  if(trId >= nUsedTrias) {
    fprintf(stdout, "*** Error: Unable to access Triangle %d of the embedded surface(%d)!\n", trId + 1, nUsedTrias);
    exit(-1);
  }
  if(!cracked[tria2quad[trId][0]]) { // no cracking
    if(hasCracked) {
      *hasCracked = false;
    }
    if(debug)
      fprintf(stdout, "--- Debug info: trId = %d, tria2quad = %d/%d, cracked = %d, phi = 1.0.\n",
              trId, tria2quad[trId][0] + 1, tria2quad[trId][1], *hasCracked);
    return 1.0;
  }
  // This element really cracked.
  if(hasCracked) {
    *hasCracked = true;
  }
  if(phantoms.find(tria2quad[trId][0]) == phantoms.end()) {
    fprintf(stdout, "*** Error:Triangle %d (in Quad %d) contains no cracking!\n", trId, tria2quad[trId][0]);
    exit(-1);
  }
  double *phi = phantoms[tria2quad[trId][0]]->phi;
  double xi3;
  switch(tria2quad[trId][1]) {
    case 0: // This triangle is ABC
      xi3 = 1.0 - xi1 - xi2;
      if(debug)
        fprintf(stdout, "--- Now in getPhi! input:(%d,%e,%e), Quad %d/%d, phi = (%e %e %e %e), phix = %e\n",
                trId + 1, xi1, xi2, tria2quad[trId][0] + 1, tria2quad[trId][1], phi[0], phi[1], phi[2], phi[3],
                phi[0]*xi1 * (1.0 - xi3) + phi[1] * (1.0 - xi1) * (1.0 - xi3) + phi[2] * (1.0 - xi1)*xi3 + phi[3]*xi1 * xi3);
      return phi[0] * xi1 * (1.0 - xi3) + phi[1] * (1.0 - xi1) * (1.0 - xi3) + phi[2] * (1.0 - xi1) * xi3 + phi[3] * xi1 * xi3;
    case 1: // This triangle is ACD
      if(debug)
        fprintf(stdout, "--- Now in getPhi! input:(%d,%e,%e), Quad %d/%d, phi = (%e %e %e %e), phix = %e\n",
                trId + 1, xi1, xi2, tria2quad[trId][0] + 1, tria2quad[trId][1], phi[0], phi[1], phi[2], phi[3],
                phi[0]*xi1 * (1.0 - xi2) + phi[1]*xi1 * xi2 + phi[2] * (1.0 - xi1)*xi2 + phi[3] * (1.0 - xi1) * (1.0 - xi2));
      return phi[0] * xi1 * (1.0 - xi2) + phi[1] * xi1 * xi2 + phi[2] * (1.0 - xi1) * xi2 + phi[3] * (1.0 - xi1) * (1.0 - xi2);
    default:
      fprintf(stdout, "Software bug in the cracking surface...\n");
      exit(-1);
  }
}

//----------------------------------------------------------------------------

bool CrackingSurface::purelyPhantom(int trId) {
  if(trId >= nUsedTrias) {
    fprintf(stdout, "*** Error: Unable to access Triangle %d of the embedded surface(%d)!\n", trId + 1, nUsedTrias);
    exit(-1);
  }
  return deleted[tria2quad[trId][0]];
}

//----------------------------------------------------------------------------

void CrackingSurface::printInfo(char *filename) {
  FILE *myout = fopen(filename, "w");
  fprintf(myout, "...Information about the cracking surface...\n");
  fprintf(myout, "  elemType = %d, n{Total,Used}Quads = %d/%d, n{Total,Used}Trias = %d/%d, n{Total,Used}Nodes = %d/%d.\n",
          elemType, nTotalQuads, nUsedQuads, nTotalTrias, nUsedTrias, nTotalNodes, nUsedNodes);
  fprintf(myout, "phantoms(%d)...\n", int(phantoms.size()));
  for(map<int, PhantomElement *>::iterator it = phantoms.begin(); it != phantoms.end(); it++) {
    fprintf(myout, "  %d --> nodes(%d): ", it->first + 1, it->second->nNodes);
    int *nod = it->second->nodes;
    for(int i = 0; i < it->second->nNodes; i++) {
      fprintf(myout, "%d ", nod[i] + 1);
    }
    fprintf(myout, "\n");
    fprintf(myout, "         phi: ");
    double *ph = it->second->phi;
    for(int i = 0; i < it->second->nNodes; i++) {
      fprintf(myout, "%e ", ph[i]);
    }
    fprintf(myout, "\n");
  }
  fprintf(myout, "latest...\n");
  fprintf(myout, "  phantomQuads(%d): ", int(latest.phantomQuads.size()));
  for(set<int>::iterator it = latest.phantomQuads.begin(); it != latest.phantomQuads.end(); it++) {
    fprintf(myout, "%d ", (*it) + 1);
  }
  fprintf(myout, "\n");
  fprintf(myout, "  phantomNodes(%d): ", int(latest.phantomNodes.size()));
  for(map<int, int>::iterator it = latest.phantomNodes.begin(); it != latest.phantomNodes.end(); it++) {
    fprintf(myout, "%d(%d)  ", it->first + 1, it->second + 1);
  }
  fprintf(myout, "\n");
  fprintf(myout, "tria2quad(%d)...\n", nTotalTrias);
  for(int i = 0; i < nTotalTrias; i++) {
    fprintf(myout, "  Tria %d <-> Quad %d.%d\n", i + 1, tria2quad[i][0] + 1, tria2quad[i][1] + 1);
  }
  fprintf(myout, "quad2tria(%d)...\n", nTotalQuads);
  for(int i = 0; i < nTotalQuads; i++) {
    fprintf(myout, "  Quad %d <-> Tria %d and %d\n", i + 1, quad2tria[i][0] + 1, quad2tria[i][1] + 1);
  }
  fprintf(myout, "cracked(%d)...\n", nTotalQuads);
  for(int i = 0; i < nTotalQuads; i++) {
    fprintf(myout, "  Quad %d -> %d\n", i + 1, cracked[i]);
  }
  fprintf(myout, "deleted(%d)...\n", nTotalQuads);
  for(int i = 0; i < nTotalQuads; i++) {
    fprintf(myout, "  Quad %d -> %d\n", i + 1, deleted[i]);
  }
  fclose(myout);
}

//----------------------------------------------------------------------------

void CrackingSurface::writeCrackingData(std::ofstream& restart_file) const {
  restart_file.write(reinterpret_cast<const char *>(&nTotalQuads), sizeof(nTotalQuads));
  restart_file.write(reinterpret_cast<const char *>(&nUsedQuads), sizeof(nUsedQuads));
  restart_file.write(reinterpret_cast<const char *>(&nTotalTrias), sizeof(nTotalTrias));
  restart_file.write(reinterpret_cast<const char *>(&nUsedTrias), sizeof(nUsedTrias));
  restart_file.write(reinterpret_cast<const char *>(&nTotalNodes), sizeof(nTotalNodes));
  restart_file.write(reinterpret_cast<const char *>(&nUsedNodes), sizeof(nUsedNodes));
  int sz = int(phantoms.size());
  restart_file.write(reinterpret_cast<const char *>(&sz), sizeof(int));
  for(map<int, PhantomElement *>::const_iterator itr = phantoms.begin();
      itr != phantoms.end(); ++itr) {
    restart_file.write(reinterpret_cast<const char *>(&itr->first),
                       sizeof(int));
    (itr->second)->writeCrackingData(restart_file);
  }
  restart_file.write(reinterpret_cast<const char *>(cracked), sizeof(bool)*nTotalQuads);
  restart_file.write(reinterpret_cast<const char *>(deleted), sizeof(bool)*nTotalQuads);
}

//------------------------------------------------------------------------------

void CrackingSurface::readCrackingData(std::ifstream& restart_file) {
  restart_file.read(reinterpret_cast<char *>(&nTotalQuads), sizeof(nTotalQuads));
  restart_file.read(reinterpret_cast<char *>(&nUsedQuads), sizeof(nUsedQuads));
  restart_file.read(reinterpret_cast<char *>(&nTotalTrias), sizeof(nTotalTrias));
  restart_file.read(reinterpret_cast<char *>(&nUsedTrias), sizeof(nUsedTrias));
  restart_file.read(reinterpret_cast<char *>(&nTotalNodes), sizeof(nTotalNodes));
  restart_file.read(reinterpret_cast<char *>(&nUsedNodes), sizeof(nUsedNodes));
  int sz, tmp;
  restart_file.read(reinterpret_cast<char *>(&sz), sizeof(int));
  for(int i = 0; i < sz; ++i) {
    restart_file.read(reinterpret_cast<char *>(&tmp),
                      sizeof(int));
    phantoms[tmp] = PhantomElement::readCrackingData(restart_file);
  }
  restart_file.read(reinterpret_cast<char *>(cracked), sizeof(bool)*nTotalQuads);
  restart_file.read(reinterpret_cast<char *>(deleted), sizeof(bool)*nTotalQuads);
  constructTriangleMap();
}

//------------------------------------------------------------------------------

void CrackingSurface::constructTriangleMap() {
  int real_tri_count = 0;
  for(int i = 0; i < nUsedTrias; i++) {
    if(!purelyPhantom(i)) {
      ++real_tri_count;
    }
  }
  if(triangle_id_map) {
    delete [] triangle_id_map;
  }
  triangle_id_map = new int[real_tri_count];
  int q = 0;
  for(int i = 0; i < nUsedTrias; i++) {
    if(!purelyPhantom(i)) {
      triangle_id_map[q] = i;
      ++q;
    }
  }
  numRealTriangles = real_tri_count;
}

//------------------------------------------------------------------------------

int CrackingSurface::mapTriangleID(int id) {
  if(triangle_id_map) {
    return triangle_id_map[id];
  }
  else {
    return id;
  }
}

//------------------------------------------------------------------------------

double CrackingSurface::getPhiPhysBAM(int trId, double xi1, double xi2, bool *hasCracked, bool debug) {
  if(triangle_id_map) {
    return getPhi(triangle_id_map[trId], xi1, xi2, hasCracked, debug);
  }
  else {
    return getPhi(trId, xi1, xi2, hasCracked, debug);
  }
}

//------------------------------------------------------------------------------

bool CrackingSurface::purelyPhantomPhysBAM(int trId) {
  if(triangle_id_map) {
    return purelyPhantom(triangle_id_map[trId]);
  }
  else {
    return purelyPhantom(trId);
  }
}

//------------------------------------------------------------------------------
