/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<ClosestTriangle.h>
#include<iostream>
using std::set;
using std::map;

//----------------------------------------------------------------------------

ClosestTriangle::ClosestTriangle(int (*nd)[3], Vec3D *sX, Vec3D *sN, set<int> *n2n, set<int> *n2e) {
  fail = false;
  triNodes = nd;
  structX = sX;
  structNorm = sN;
  node2node = n2n;
  node2elem = n2e;
}

//----------------------------------------------------------------------------

void
ClosestTriangle::start(Vec3D xp) {
  x = xp;
  isFirst = true;
  bestTrId = -1;
  n1 = n2 = -1;
  mode = -1;
  nPairs = 0;
  minDist = 1.0e10;
  nd2tri.clear();
  nd2tester.clear();
  dist2othernode = 1.0e10;
}
//----------------------------------------------------------------------------

double ClosestTriangle::edgeProject(Vec3D x0, int n1, int n2, double &alpha) const
{
  Vec3D xA =   structX[n1];
  Vec3D xB =   structX[n2];
  Vec3D AB= xB-xA;
  Vec3D AX = x0-xA;
  alpha = AB*AX/(AB*AB);
  Vec3D P = xA + alpha*AB;
  return (P-x0).norm();
}

//----------------------------------------------------------------------------

double ClosestTriangle::project(Vec3D x0, int tria, double& xi1, double& xi2) const
{
  int iA = triNodes[tria][0];
  int iB = triNodes[tria][1];
  int iC = triNodes[tria][2];
  Vec3D xA = structX[iA];
  Vec3D xB = structX[iB];
  Vec3D xC = structX[iC];

  Vec3D ABC = 0.5*(xB-xA)^(xC-xA);
  double areaABC = ABC.norm();
  Vec3D dir = 1.0/areaABC*ABC;

  //calculate the projection.
  double dist = (x0-xA)*dir;
  Vec3D xp = x0 - dist*dir;

  //calculate barycentric coords.
  double areaPBC = (0.5*(xB-xp)^(xC-xp))*dir;
  double areaPCA = (0.5*(xC-xp)^(xA-xp))*dir;
  xi1 = areaPBC/areaABC;
  xi2 = areaPCA/areaABC;

  return dist;
}

//----------------------------------------------------------------------------

void
ClosestTriangle::checkTriangle(int trId)
{
  double dist;
  int *nd = triNodes[trId];
  double xi1, xi2;
  dist = project(x, trId, xi1, xi2);
  double xi3 = 1-xi1-xi2;
  const double eps = 0;

  if(xi1 >= -eps && xi2 >= -eps && xi3 >= -eps) {
    if(isFirst || std::fabs(minDist) > std::fabs(dist)) {
      isFirst = false;
      minDist = dist;
      bestTrId = trId;
      n1 = n2 = -1;
      mode = 0;
    }
  }else {
    if(xi1 < -eps)
      checkEdge(trId, nd[1], nd[2], nd[0], dist);
    if(xi2 < -eps)
      checkEdge(trId, nd[0], nd[2], nd[1], dist);
    if(xi3 < -eps)
      checkEdge(trId, nd[0], nd[1], nd[2], dist);
  }
}

//----------------------------------------------------------------------------

void ClosestTriangle::checkVertex(int ip1, int trId, double trDist) {

  // If this node is already our best solution
  if(n1 == ip1 && n2 < 0) {
    if(nPairs == maxNPairs) {
      std::cerr << "Warning: On the embedded structure surface, detected too many (>" << maxNPairs <<
                   ") peripheral triangles to a node!" << std::endl;
    }
    periTri[nPairs++] = trId;
    int repeated1, repeated2;
    registerNodes(ip1,trId, repeated1, repeated2);
    bool thisSign = (trDist >= 0);

    if((thisSign != isPositive || !isConsistent)/* && !hasBeenFixed*/) { // need to figure out the correct sign...
      isConsistent = false;
      if(repeated1 != -1) {// this triangle and another traversed triangle share an edge that has this damn node.
        //Step 1. Determine if it is a "hill" or a "dent".
        int repeated = repeated1;
        double dist2node = (structX[triNodes[trId][repeated]]-x).norm();
        if(dist2node<dist2othernode) {
          dist2othernode = dist2node;
          double xi1_temp, xi2_temp;
          double dist2 = project(structX[nd2tester[triNodes[trId][repeated]]], trId, xi1_temp, xi2_temp);
          int type = (dist2>0) ? 1 : 2; //1: dent, 2: hill
          //Step 2. Check the angles between (ip1->x) and the normal of the two triangles
          Vec3D ip1_x = x - structX[ip1];
          double alpha = ip1_x*structNorm[trId];
          double beta  = ip1_x*structNorm[nd2tri[triNodes[trId][repeated]]];
          if(type==1)
            isPositive = (alpha>=0&&beta>=0) ? true : false;
          else
            isPositive = (alpha<0&&beta<0) ? false : true;

//          hasBeenFixed = true; //the sign is determined.
        }
      }
      if(repeated2 != -1) {// this triangle and another traversed triangle share an edge that has this damn node.
        //Step 1. Determine if it is a "hill" or a "dent".
        int repeated = repeated2;
        double dist2node = (structX[triNodes[trId][repeated]]-x).norm();
        if(dist2node<dist2othernode) {
          dist2othernode = dist2node;
          double xi1_temp, xi2_temp;
          double dist2 = project(structX[nd2tester[triNodes[trId][repeated]]], trId, xi1_temp, xi2_temp);
          int type = (dist2>0) ? 1 : 2; //1: dent, 2: hill
          //Step 2. Check the angles between (ip1->x) and the normal of the two triangles
          Vec3D ip1_x = x - structX[ip1];
          double alpha = ip1_x*structNorm[trId];
          double beta  = ip1_x*structNorm[nd2tri[triNodes[trId][repeated]]];
          if(type==1)
            isPositive = (alpha>=0&&beta>=0) ? true : false;
          else
            isPositive = (alpha<0&&beta<0) ? false : true;

//          hasBeenFixed = true; //the sign is determined.
        }
      }
    }
    return;
  }

  // Compute the distance to the node. Determining the sign of the distance
  // is delayed until the very end
  double dist = (x-structX[ip1]).norm();
  if(isFirst || dist < std::fabs(minDist)) {
    isFirst = false;
    n1 = ip1;
    n2 = -1;
    nPairs = 1;
    minDist = dist;
    dist2othernode = 1.0e10; //distance to the other node on the edge that is shared with another triangle
    periTri[0] = trId;
    bestTrId = trId;
    nd2tri.clear();
    nd2tester.clear();
    int repeated1, repeated2;
    registerNodes(ip1,trId, repeated1, repeated2);
    mode = 2;
    hasBeenFixed = false;
    isConsistent = true;
    isPositive = trDist >= 0;
  }
}

//----------------------------------------------------------------------------

int ClosestTriangle::registerNodes(int ip1, int trId, int& repeated1, int& repeated2)
{
    map<int,int>::iterator it;

    unsigned int repeated_node = 0;
    repeated1 = repeated2 = -1;

    if(ip1==triNodes[trId][0]) {
        it = nd2tri.find(triNodes[trId][1]);
        if(it==nd2tri.end()) {
          nd2tri[triNodes[trId][1]] = trId;
          nd2tester[triNodes[trId][1]] = triNodes[trId][2];
        } else if(it->second!=trId) {
          repeated_node |= 1 << 1;
          repeated1 = 1;
        }
        it = nd2tri.find(triNodes[trId][2]);
        if(it==nd2tri.end()) {
           nd2tri[triNodes[trId][2]] = trId;
           nd2tester[triNodes[trId][2]] = triNodes[trId][1];
        } else if(it->second!=trId) {
          repeated_node |= 1 << 2;
          repeated2 = 2;
        }
    }
    else if(ip1==triNodes[trId][1]) {
        it = nd2tri.find(triNodes[trId][2]);
        if(it==nd2tri.end()) {
          nd2tri[triNodes[trId][2]] = trId;
          nd2tester[triNodes[trId][2]] = triNodes[trId][0];
        } else if(it->second!=trId) {
          repeated_node |= 1 << 2;
          repeated1 = 2;
        }
        it = nd2tri.find(triNodes[trId][0]);
        if(it==nd2tri.end()) {
           nd2tri[triNodes[trId][0]] = trId;
           nd2tester[triNodes[trId][0]] = triNodes[trId][2];
        } else if(it->second!=trId) {
          repeated_node |= 1 << 0;
          repeated2 = 0;
        }
    }
    else if(ip1==triNodes[trId][2]) {
        it = nd2tri.find(triNodes[trId][0]);
        if(it==nd2tri.end()) {
          nd2tri[triNodes[trId][0]] = trId;
          nd2tester[triNodes[trId][0]] = triNodes[trId][1];
        } else if(it->second!=trId) {
          repeated_node |= 1 << 0;
          repeated1 = 0;
        }
        it = nd2tri.find(triNodes[trId][1]);
        if(it==nd2tri.end()) {
           nd2tri[triNodes[trId][1]] = trId;
           nd2tester[triNodes[trId][1]] = triNodes[trId][0];
        } else if(it->second!=trId) {
          repeated_node |= 1 << 1;
          repeated2 = 1;
        }
    } else {
      fprintf(stdout,"*** Error (software bug): node %d doesn't belong to triangle %d!\n", ip1+1, trId+1);
      exit(-1);
    }

    return repeated_node;
}

//----------------------------------------------------------------------------

bool
ClosestTriangle::checkEdge(int trId, int ip1, int ip2, int p3, double trDist) {
  int p1, p2;
  if(ip1 < ip2) {
    p1 = ip1;
    p2 = ip2;
  } else {
    p1 = ip2;
    p2 = ip1;
  }

  if(n1 == p1 && n2 == p2) { //<! If we've already looked at this edge before
    // If we have already examined this edge and the sign opinions disagree, we need to make the right decision
    // Check if the p3 is in the positive or negative side of the other triangle
    // When signs disagree, the true distance has the opposite sign. because the triangle surface is the
    // end of point with the same sign as p3;
    if(trDist*minDist >= 0)
      return true;
    double xi1, xi2;
    double d2 = project(structX[p3], bestTrId, xi1, xi2);
    if(d2*minDist>0)
      minDist = -minDist;
    return true;
  } else {
    double dist, alpha;
    double sign = trDist >= 0 ? 1 : -1;
    dist = sign*edgeProject(x, p1, p2, alpha);

    int cn1 = (alpha > 1) ? p2 : p1;
    int cn2 = p2;
    if(alpha < 0 || alpha > 1) {
      checkVertex(cn1, trId, trDist);
      return true;
    }

    if(isFirst || std::fabs(dist) < std::fabs(minDist)) {
      isFirst = false;
      bestTrId = trId;
      n = structNorm[bestTrId];
      minDist = dist;
      n1 = cn1;
      n2 = cn2;
      mode = 1;
      return true;
    }
  }
  return false;
}

//----------------------------------------------------------------------------

double ClosestTriangle::getSignedVertexDistance() const {

  return isPositive ? minDist : -minDist;

/*  if(isConsistent)
    return isPositive ? minDist : -minDist;
  else if(hasBeenFixed)
    return isPositive ? minDist : -minDist;
  return minDist;
*/
}

//----------------------------------------------------------------------------

double ClosestTriangle::findSignedVertexDistance()
{
  set<int> &vertices = node2node[n1]; //vertices in the direct neighborhood of n1
  Vec3D    &xp       = structX[n1];   //coordinate of n1

  //step 1: find a test direction.
  Vec3D dir(0,0,0); //the test direction (normalized).
  double rr = 0.0; //distance in the test direction.

  double ry = 0.0;
  Vec3D npx = x - xp;
  npx = 1.0/npx.norm()*npx;
  for(set<int>::iterator it=vertices.begin(); it!=vertices.end(); it++) {
    Vec3D xr = structX[*it];

    double d_xp_xr = (xp-xr).norm();
    if(d_xp_xr>=ry)
      ry = d_xp_xr;

    if(npx*(xp-xr)<1.0e-14) continue;
    double t = npx*(x-xr)/(npx*(xp-xr));
    Vec3D xt = xr + t*(xp-xr);
    double d_x_xt = (xt - x).norm();
    if(d_x_xt>=rr) {
      dir = xt - x;
      rr  = d_x_xt;
    }
  }
  if(rr<1.0e-14) {fprintf(stdout,"*** Error: (in ClosestTriangle) distance = %e.\n",rr);exit(-1);}
  dir *= 1.0/rr; //normalize dir
  rr = std::min(1.5*rr, std::max((x-xp).norm(), ry));
  rr *= 0.05;

  int nTrial = 50;
  for(int iTrial=0; iTrial<nTrial; iTrial++) {
    Vec3D x_trial;
    if(iTrial%2==1)
      x_trial = x - (double(iTrial+1.0)*rr)*dir;
    else
      x_trial = x + (double(iTrial+1.0)*rr)*dir;

    int sign = findTesterStatus(x_trial);
    if(sign!=0) {
      minDist = (double)sign*std::fabs(minDist);
//      fprintf(stdout,"NOTE: x = (%e, %e, %e), node = %d, x_trial = (%e, %e, %e). sign = %d, minDist = %e\n", x[0], x[1], x[2], n1+1, x_trial[0], x_trial[1], x_trial[2], sign, minDist);
      return minDist;
    }// else
 //     fprintf(stdout,"NOTE: x = (%e, %e, %e), node = %d, x_trial = (%e, %e, %e). I = %d, sign = %d\n", x[0], x[1], x[2], n1+1, x_trial[0], x_trial[1], x_trial[2], iTrial, sign);
  }

  fprintf(stdout,"*** Error: (in ClosestTriangle) failed in determining node status! nTrial = %d.\n", nTrial);
  fail = true;
//  exit(-1);

  return 0.0;
}

//----------------------------------------------------------------------------

int ClosestTriangle::findTesterStatus(Vec3D xt) const
{
//  Vec3D xdebug(-4.348610e+00, -5.030928e+00, -6.636450e-02);
//  bool debug = (xt-xdebug).norm()<1.0e-5 ? true : false;

  set<int> &elements = node2elem[n1]; //elements in the direct neighborhood of n1
  double mindist = 1.0e14;
  const double eps = 1.0e-14;
  int nn1,nn2;
  int myMode = -1;
  int bestTriangle = -1;
  nn1 = nn2 = -1;

  for(set<int>::iterator it=elements.begin(); it!=elements.end(); it++) {
    double xi[3];
    double dist = project(xt, *it, xi[0], xi[1]);
    xi[2] = 1.0 - xi[0] - xi[1]; // project onto the plane determined by this triangle

    if((n1==triNodes[*it][0] || xi[0] >= -eps) &&
       (n1==triNodes[*it][1] || xi[1] >= -eps) &&
       (n1==triNodes[*it][2] || xi[2] >= -eps)) {
      if(std::fabs(mindist) >= std::fabs(dist)) {
        mindist = dist;
        myMode = 0;
        nn1 = nn2 = -1;
        bestTriangle = *it;
      }
    }
    else {
      if(!(n1==triNodes[*it][0] || xi[0] >= -eps))
        checkEdgeForTester(xt, *it, triNodes[*it][1], triNodes[*it][2], triNodes[*it][0], dist, nn1, nn2, mindist, myMode, bestTriangle, eps);
      if(!(n1==triNodes[*it][1] || xi[1] >= -eps))
        checkEdgeForTester(xt, *it, triNodes[*it][2], triNodes[*it][0], triNodes[*it][1], dist, nn1, nn2, mindist, myMode, bestTriangle, eps);
      if(!(n1==triNodes[*it][2] || xi[2] >= -eps))
        checkEdgeForTester(xt, *it, triNodes[*it][0], triNodes[*it][1], triNodes[*it][2], dist, nn1, nn2, mindist, myMode, bestTriangle, eps);
    }
//    if(debug)
//      fprintf(stdout,"NOTE (debug): node = %d, trId = %d. | mode = %d, nn1 = %d, nn2 = %d, bestTriangle = %d, mindist = %e.\n",
//              n1+1, *it+1, myMode, nn1+1, nn2+1, bestTriangle+1, mindist);
  }
//  if(debug)
//    fprintf(stdout,"NOTE (final): node = %d, mode = %d, nn1 = %d, nn2 = %d, bestTriangle = %d, mindist = %e.\n",
//            n1+1, myMode, nn1+1, nn2+1, bestTriangle+1, mindist);


  int sign;
  if(myMode<0) {
//    fprintf(stdout,"WARNING: (in IntersectorFRG) Mode = %d!\n", myMode);
//    fprintf(stdout,"-- x = (%e %e %e), tester = (%e, %e, %e).\n", x[0], x[1], x[2], xt[0], xt[1], xt[2]);
    sign = 0;
  } else
    sign = mindist>=0.0 ? 1 : -1;

  return sign;
}

//----------------------------------------------------------------------------

void ClosestTriangle::checkEdgeForTester(Vec3D xt, int trId, int ip1, int ip2, int p3, double trDist, int &nn1, int &nn2,
                                         double &mindist, int &myMode, int &bestTriangle, const double eps) const
{
  int p1, p2;
  if(ip1 < ip2) {
    p1 = ip1;  p2 = ip2;
  } else {
    p1 = ip2;  p2 = ip1;
  }

  if(nn1 == p1 && nn2 == p2) { //this edge has been traversed.
    if(trDist*mindist >= 0)
      return;
    double xi1,xi2;
    double d2 = project(structX[p3], bestTriangle, xi1, xi2);
    if(d2*mindist>0)
      mindist = -mindist;
    return;
  } else {
    double dist, alpha;
    double sign = trDist >= 0 ? 1 : -1;
    dist = sign*edgeProject(xt, p1, p2, alpha);
    if((alpha<-eps && p1==n1) || (alpha>1.0+eps && p2==n1))
      return;

    if(std::fabs(dist) < std::fabs(mindist)) {
      bestTriangle = trId;
      nn1 = p1;
      nn2 = p2;
      myMode = 1;
      mindist = dist;
    }
  }
}

//----------------------------------------------------------------------------

