#ifndef _MESH_MATCHER_H_
#define _MESH_MATCHER_H_

#include<Vector3D.h>
#include<SpaceVariable.h>

//------------------------------------------------------------
// Class MeshMatcher is responsible for transferring data from 
// one (partitioned) mesh to another. 
//------------------------------------------------------------

class MeshMatcher {

  //! whether exact node-to-node matching is expected
  bool exact_match;

  //! communicators (may be NULL for some processor cores)
  MPI_Comm* comm; //all the cores
  MPI_Comm* comm1; //!< mesh 1
  MPI_Comm* comm2; //!< mesh 2
  int rank, size, rank1, size1, rank2, size2;

  //! mesh coordinates
  SpaceVariable3D* coordinates1;
  SpaceVariable3D* coordinates2;

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

  MeshMatcher(MPI_Comm* comm_, MPI_Comm* comm1_, MPI_Comm* comm2_, 
              SpaceVariable3D* coords1_, SpaceVariable3D* coords2_,
              vector<double>& x1, vector<double>& y1, vector<double>& z1,
              vector<double>& x2, vector<double>& y2, vector<double>& z2);

  ~MeshMatcher();
  
  void SendData(SpaceVariable3D* V1, SpaceVariable3D* V2);



};



#endif
