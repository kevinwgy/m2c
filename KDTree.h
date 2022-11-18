#ifndef _KD_TREE_
#define _KD_TREE_

/**********************************************************************************
 * Copyright © Michel Lesoinne, 2008
 **********************************************************************************/

template <class Obj> class DirComp {
    int dir;
  public:
    DirComp(int d) : dir(d) {}
    bool operator()(const Obj &o1, const Obj &o2) const
      { return o1.val(dir) < o2.val(dir); }
};

/** This templated class implements a kD-tree
 * The Obj class can be any type
 * CompType must provide a comparison object for each direction.
 * Objects in the tree have a span and the tree is built based on their
 * lowest-coordinates. They need to return what their width is in each direction*/
template <class Obj, int dim=3, class CompType=DirComp<Obj> >
class KDTree {
   protected:
     int nObj;
     int dir;
     double splitVal; // Value detemining which branch to take
     double w[dim];
     Obj *obj; // Non zero if this is a terminal leaf.
     KDTree<Obj, dim, CompType> *leftTree, *rightTree;

     int getBestSplit(int nobj, Obj *allObjs, int dir);
   public:
     KDTree(int nobj, Obj *allObjs, int depth = 0);
     ~KDTree();
     // Give the maximum width of object in this tree.
     // note this is not the same as the width of the bounding box
     // of the tree.
     double maxWidth(int dir)
       { return w[dir]; }
     int maxdepth() {
       if(obj) return 1;
       else return std::max(leftTree->maxdepth(), rightTree->maxdepth())+1;
     }
     int maxLeaf() {
       if(obj) return nObj;
       else return std::max(leftTree->maxLeaf(), rightTree->maxLeaf());
     }
     int findCandidates(double x[dim], Obj *o, int maxNObj, int depth=0);
     int findCloseCandidates(double x[dim], Obj *o, int maxNObj, double &dist);
     int findCandidatesWithin(double x[dim], Obj *o, int maxNObj, double dist);
     int findCandidatesInBox(double xmin[dim], double xmax[dim],  Obj *o, int maxNObj);
};

#include <algorithm>


template <class Obj, int dim, class CompType>
KDTree<Obj, dim, CompType>::KDTree(int nobj, Obj *allObjs, int depth) {
  nObj = nobj;
  dir = depth % dim;
  int firstDir = dir;
  int split;
  int bestDir = 0;
  double quality[dim], bestQuality = 0.0;

  while(1) {
    if(nObj <= 4) {
      split = nObj;
      quality[dir] = 0;
      bestDir = dir;
      break;
    }
    split = getBestSplit(nObj, allObjs, dir);
    // compute the quality
    quality[dir] = std::min(split, nObj-split);
    // is it the best we got?
    if(dir == firstDir || quality[dir] > bestQuality) {
      bestQuality = quality[dir];
      bestDir = dir;
    }
    // See if we should look for better.
    if(quality[dir] > 0.25*nObj)
      break;
    else {
       int newdir = (dir+1)%dim;
       if(newdir == firstDir) break;
       dir = newdir;
    }
  }

  if(quality[bestDir] < 0.25*nObj) {
    obj = allObjs;
    int i, d;
    for(d = 0; d < dim; ++d)
      w[d] = 0.0;
    for(i = 0; i < nObj; ++i)
      for(d = 0; d < dim; ++d)
        w[d] = std::max(w[d], obj[i].width(d));
    return;
  }
  if(dir != bestDir)
     split = getBestSplit(nobj, allObjs, bestDir);
  dir = bestDir;
  splitVal = ( allObjs[split].val(dir) + allObjs[split-1].val(dir) )/2.0;
  obj = 0; // Indicates we have subtrees
  leftTree = new KDTree<Obj, dim, CompType>(split, allObjs,dir+1);
  rightTree = new KDTree<Obj, dim, CompType>(nObj-split, allObjs+split, dir+1);
  for(int d = 0; d < dim; ++d)
    w[d] = std::max(leftTree->w[d], rightTree->w[d]);
}

template <class Obj, int dim, class CompType>
KDTree<Obj, dim, CompType>::~KDTree() {
  if(obj == 0) {
    delete leftTree;
    delete rightTree;
  }
}


template <class Obj, int dim, class CompType>
int
KDTree<Obj, dim, CompType>::getBestSplit(int nobj, Obj *allObjs, int dir) {
  if(nobj <= 4) {
     return nobj;
  }

  CompType compare(dir);
  std::sort(allObjs, allObjs+nobj, compare);
  int leftSplit = nobj/2;
  int rightSplit = leftSplit;
  // Find the best median. Make sure that the split point
  // is such that the object at the median and the one before it
  // are not collocated.
  while(leftSplit > 0
         && allObjs[leftSplit].val(dir) == allObjs[leftSplit-1].val(dir))
    leftSplit--;
  // Note that since nobj>4 we are certain that rightSplit cannot be 0
  while(rightSplit < nobj
         && allObjs[rightSplit].val(dir) == allObjs[rightSplit-1].val(dir))
    rightSplit++;

  int split = (nobj-leftSplit < rightSplit) ? leftSplit : rightSplit;

  return split;
};

template <class Obj, int dim, class CompType>
int
KDTree<Obj, dim, CompType>::findCandidates(double x[dim], Obj *o, int maxNObj, int depth) {
  if(obj) {
    int nPot = 0;
    for(int i = 0; i < nObj; ++i) {
      // check that we fall within the box.
      int k;

      for(k = 0; k < dim && x[k] >= obj[i].val(k) && x[k]  <= obj[i].val(k)+obj[i].width(k)
               ; )
	   ++k;
      if(k !=  dim) // if one test was failed
        continue;
      if(nPot >= maxNObj)
        return maxNObj+1; // we do not have enough space
      o[nPot++] = obj[i];
    }
    return nPot;
  }
  int nFound = 0;
  if(x[dir] <= splitVal+leftTree->maxWidth(dir))
    nFound = leftTree->findCandidates(x, o, maxNObj,depth+1);

  if(x[dir] >= splitVal)
    nFound += rightTree->findCandidates(x, o+nFound, maxNObj-nFound, depth+1);

  return nFound;
}

template <class Obj, int dim, class CompType>
int
KDTree<Obj, dim, CompType>::findCloseCandidates(double x[dim],
              Obj *o, int maxNObj, double &dist) {
  if(obj) {
    int nPot = 0;
    for(int i = 0; i < nObj; ++i) {
      int k=0;
      // compute the pseudo-distance of x to the bounding box
      double locDist =0;
      for(k = 0; k < dim; ++k)
        locDist = std::max(locDist,
	            std::max(obj[i].val(k)-x[k], x[k]-obj[i].val(k)-obj[i].width(k)));
      if(dist >= 0 && locDist > dist)
        continue;
      if(dist < 0 || locDist < dist) {
        nPot = 0;
	dist = locDist;
      }
      if(nPot >= maxNObj)
        return maxNObj+1; // we do not have enough space

      o[nPot++] = obj[i];
    }
    return nPot;
  }
  int nFound = 0;
  // Examine the two subtrees in an order dependent of x to minimize the number of leafs
  // to be examined.
  if(x[dir] <= splitVal) {
    nFound = leftTree->findCloseCandidates(x, o, maxNObj, dist);
    double potDist = dist;
    int nRightFound=0;
    if(x[dir] >= splitVal-dist) {
      nRightFound = rightTree->findCloseCandidates(x, o+nFound, maxNObj-nFound, potDist);
      if(nRightFound > maxNObj)
        return nRightFound;
      if(potDist < dist) {
        dist = potDist;
        for(int i = 0; i < nRightFound; ++i)
          o[i] = o[nFound+i];
        nFound = nRightFound;
      } else
        nFound += nRightFound;
    }
  } else {
    nFound = rightTree->findCloseCandidates(x, o, maxNObj, dist);
    double potDist = dist;
    int nLeftFound = 0;
    if(x[dir] <= splitVal+leftTree->maxWidth(dir)+dist) {
      nLeftFound = leftTree->findCloseCandidates(x, o+nFound, maxNObj-nFound, potDist);
      if(potDist < dist) {
        dist = potDist;
        for(int i = 0; i < nLeftFound; ++i)
          o[i] = o[nFound+i];
        nFound = nLeftFound;
      } else
        nFound += nLeftFound;
    }
  }
  return nFound;
}

template <class Obj, int dim, class CompType>
int
KDTree<Obj, dim, CompType>::findCandidatesWithin(double x[dim],
              Obj *o, int maxNObj, double dist) {
  if(obj) {
    int nPot = 0;
    for(int i = 0; i < nObj; ++i) {
      int k=0;
      // compute the pseudo-distance of x to the bounding box
      double locDist =0;
      for(k = 0; k < dim; ++k)
        locDist = std::max(locDist,
	            std::max(obj[i].val(k)-x[k], x[k]-obj[i].val(k)-obj[i].width(k)));
      if(locDist > dist)
        continue;
      if(nPot >= maxNObj)
        return maxNObj+1; // we do not have enough space

      o[nPot++] = obj[i];
    }
    return nPot;
  }
  int nFound = 0;

  if(x[dir] <= splitVal+leftTree->maxWidth(dir)+dist)
    nFound = leftTree->findCandidatesWithin(x, o, maxNObj, dist);

  if(x[dir]+dist >= splitVal)
    nFound += rightTree->findCandidatesWithin(x, o+nFound, maxNObj-nFound, dist);

  return nFound;
}

template <class Obj, int dim, class CompType>
int KDTree<Obj, dim, CompType>::findCandidatesInBox(double xmin[dim], double xmax[dim],  Obj *o, int maxNObj) {
  if(obj) {
    int nPot = 0;
    for(int i = 0; i < nObj; ++i) {
      bool isOutside = false;
      for(int k = 0; k < dim; ++k)
        if(xmax[k] < obj[i].val(k) ||
           xmin[k] > obj[i].val(k)+obj[i].width(k))
          isOutside = true;
      if(isOutside)
        continue;
      if(nPot < maxNObj)
        o[nPot] = obj[i];
      nPot++;
    }
    return nPot;
  }
  int nFound = 0;

  if(xmin[dir] <= splitVal+leftTree->maxWidth(dir))
    nFound = leftTree->findCandidatesInBox(xmin, xmax, o, maxNObj);

  if(xmax[dir] >= splitVal)
    nFound += rightTree->findCandidatesInBox(xmin, xmax, o+nFound, maxNObj-nFound);

  return nFound;
}

#endif
