#ifndef _GHOST_FLUID_OPERATOR_H_
#define _GHOST_FLUID_OPERATOR_H_

#include<CustomCommunicator.h>
#include<GhostPoint.h>
#include<EmbeddedBoundaryFormula.h>

class EmbeddedBoundaryDataSet;

/**************************************************************
* Class GhostFluidOperator is responsible for populating ghost
* nodes in ``inactive'' regions next to an embedded boundary.
* It requires and works with Intersector and 
* EmbeddedBoundaryFormula. Multiple population formulas are
* supported.
**************************************************************/

class GhostFluidOperator {

  MPI_Comm &comm;
  

  std::vector<std::pair<Int3, EmbeddedBoundaryFormula> > ghostNodes1; //ghost nodes whose image is inside this subdomain

  std::vector<std::vector<Int3> > ghostNodes2; //ghost nodes whose image is inside another subdomain
  std::vector<int> ghostNodes2_sender;

  std::vector<std::vector<std::pair<Int3, EmbeddedBoundaryFormula> > > friendsGhostNodes;
  std::vector<int> friendsGhostNodes_receiver;

  // laser radiance at ghost nodes
  //   std::vector<double> l1; //ghostNodes1
  //     std::vector<std::vector<double> > l2; //ghostNodes2


public:

  GhostFluidOperator(











};


#endif
