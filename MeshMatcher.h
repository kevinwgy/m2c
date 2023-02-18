/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _MESH_MATCHER_H_
#define _MESH_MATCHER_H_

#include<Vector3D.h>
#include<SpaceVariable.h>
#include<map>

//------------------------------------------------------------
// Class MeshMatcher is responsible for transferring data from 
// one (partitioned) mesh to another. 
//------------------------------------------------------------

class MeshMatcher {

  //! whether exact node-to-node matching is expected
  bool exact_match;

  //! communicator (combination of mesh1 and mesh2)
  MPI_Comm& comm; //all the cores
  int rank, size;

  //! whether the current proc cores owns part of mesh1 or mesh2
  bool mesh1_owner;
  bool mesh2_owner;

  //! matched_nodes
  std::map<int, std::vector<int> >    send_i;
  std::map<int, std::vector<int> >    send_j;
  std::map<int, std::vector<int> >    send_k;
  std::map<int, std::vector<double> > send_buffer_dim1;
  std::map<int, std::vector<double> > send_buffer_dim2;
  std::map<int, std::vector<double> > send_buffer_dim3;
  std::map<int, std::vector<int> >    recv_i;
  std::map<int, std::vector<int> >    recv_j;
  std::map<int, std::vector<int> >    recv_k;
  std::map<int, std::vector<double> > recv_buffer_dim1;
  std::map<int, std::vector<double> > recv_buffer_dim2;
  std::map<int, std::vector<double> > recv_buffer_dim3;

public:

  MeshMatcher(MPI_Comm& comm_, SpaceVariable3D* coordinates1, SpaceVariable3D* coordinates2);
  ~MeshMatcher();
  
  void Transfer(SpaceVariable3D* V1, SpaceVariable3D* V2);

private:

  void FindGlobalMeshInfo(SpaceVariable3D* coordinates, std::vector<double> &x, std::vector<double> &y, 
                          std::vector<double> &z, std::vector<double> &x_minus, std::vector<double> &x_plus, 
                          std::vector<double> &y_minus, std::vector<double> &y_plus, 
                          std::vector<double> &z_minus, std::vector<double> &z_plus);

  void SetupTransferExactMatch(std::vector<double>& x1, std::vector<double>& y1, std::vector<double>& z1,
                               std::vector<double>& x1_minus, std::vector<double>& x1_plus, //coords of ghosts
                               std::vector<double>& y1_minus, std::vector<double>& y1_plus, //coords of ghosts
                               std::vector<double>& z1_minus, std::vector<double>& z1_plus, //coords of ghosts
                               std::vector<double>& x2, std::vector<double>& y2, std::vector<double>& z2,
                               std::vector<double>& x2_minus, std::vector<double>& x2_plus, //coords of ghosts
                               std::vector<double>& y2_minus, std::vector<double>& y2_plus, //coords of ghosts
                               std::vector<double>& z2_minus, std::vector<double>& z2_plus, //coords of ghosts
                               bool mesh1_owner, int ii0_1, int jj0_1, int kk0_1,
                               int iimax_1, int jjmax_1, int kkmax_1,
                               bool mesh2_owner, int ii0_2, int jj0_2, int kk0_2,
                               int iimax_2, int jjmax_2, int kkmax_2);

};



#endif
