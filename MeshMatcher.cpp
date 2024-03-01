/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<MeshMatcher.h>
#include<SpaceVariable.h>
#include<Utils.h>
#include<cfloat> //DBL_MAX
#include<cassert>
using std::vector;
using std::pair;
using std::map;

//-------------------------------------------------------------------------------

MeshMatcher::MeshMatcher(MPI_Comm& comm_, SpaceVariable3D* coordinates1, SpaceVariable3D* coordinates2)
           : comm(comm_)
{

  exact_match = true; //can be extended in future if needed

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  mesh1_owner = (coordinates1!=NULL);
  mesh2_owner = (coordinates2!=NULL);

//  int i0_1(-INT_MAX), j0_1(-INT_MAX), k0_1(-INT_MAX), imax_1(-INT_MAX), jmax_1(-INT_MAX), kmax_1(-INT_MAX);    
  int ii0_1(-INT_MAX), jj0_1(-INT_MAX), kk0_1(-INT_MAX), iimax_1(-INT_MAX), jjmax_1(-INT_MAX), kkmax_1(-INT_MAX);    
  if(coordinates1) {
//    coordinates1->GetCornerIndices(&i0_1, &j0_1, &k0_1, &imax_1, &jmax_1, &kmax_1);
    coordinates1->GetGhostedCornerIndices(&ii0_1, &jj0_1, &kk0_1, &iimax_1, &jjmax_1, &kkmax_1);
  }

//  int i0_2(-INT_MAX), j0_2(-INT_MAX), k0_2(-INT_MAX), imax_2(-INT_MAX), jmax_2(-INT_MAX), kmax_2(-INT_MAX);    
  int ii0_2(-INT_MAX), jj0_2(-INT_MAX), kk0_2(-INT_MAX), iimax_2(-INT_MAX), jjmax_2(-INT_MAX), kkmax_2(-INT_MAX);    
  if(coordinates2) {
//    coordinates2->GetCornerIndices(&i0_2, &j0_2, &k0_2, &imax_2, &jmax_2, &kmax_2);
    coordinates2->GetGhostedCornerIndices(&ii0_2, &jj0_2, &kk0_2, &iimax_2, &jjmax_2, &kkmax_2);
  }

  vector<double> x1, y1, z1, x1_minus, x1_plus, y1_minus, y1_plus, z1_minus, z1_plus;
  FindGlobalMeshInfo(coordinates1, x1, y1, z1, x1_minus, x1_plus, y1_minus, y1_plus, z1_minus, z1_plus);

  vector<double> x2, y2, z2, x2_minus, x2_plus, y2_minus, y2_plus, z2_minus, z2_plus;
  FindGlobalMeshInfo(coordinates2, x2, y2, z2, x2_minus, x2_plus, y2_minus, y2_plus, z2_minus, z2_plus);

  SetupTransferExactMatch(x1,y1,z1,x1_minus,x1_plus,y1_minus,y1_plus,z1_minus,z1_plus,
                          x2,y2,z2,x2_minus,x2_plus,y2_minus,y2_plus,z2_minus,z2_plus,
                          mesh1_owner, ii0_1, jj0_1, kk0_1, iimax_1, jjmax_1, kkmax_1,
                          mesh2_owner, ii0_2, jj0_2, kk0_2, iimax_2, jjmax_2, kkmax_2);

}

//-------------------------------------------------------------------------------

MeshMatcher::~MeshMatcher()
{ }

//-------------------------------------------------------------------------------

void
MeshMatcher::FindGlobalMeshInfo(SpaceVariable3D* coordinates, std::vector<double> &x, std::vector<double> &y,
                                std::vector<double> &z, std::vector<double> &x_minus, std::vector<double> &x_plus, 
                                std::vector<double> &y_minus, std::vector<double> &y_plus, 
                                std::vector<double> &z_minus, std::vector<double> &z_plus)
{

  vector<int> info(4, -INT_MAX);
  if(coordinates) {//owner
    coordinates->GetGlobalSize(&(info[0]), &(info[1]), &(info[2]));
    info[3] = coordinates->NumGhostLayers();
  }
  MPI_Allreduce(MPI_IN_PLACE, info.data(), info.size(), MPI_INT, MPI_MAX, comm);

  int i0(-INT_MAX), j0(-INT_MAX), k0(-INT_MAX), imax(-INT_MAX), jmax(-INT_MAX), kmax(-INT_MAX);    
  int ii0(-INT_MAX), jj0(-INT_MAX), kk0(-INT_MAX), iimax(-INT_MAX), jjmax(-INT_MAX), kkmax(-INT_MAX);    
  if(coordinates) {
    coordinates->GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
    coordinates->GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);
  }

  Vec3D*** coords = coordinates ? (Vec3D***)coordinates->GetDataPointer() : NULL;

  //collect x,y,z and the ghost coords
  x.assign(info[0], -INT_MAX);
  y.assign(info[1], -INT_MAX);
  z.assign(info[2], -INT_MAX);
  x_minus.assign(info[3], -INT_MAX);  x_plus.assign(info[3], -INT_MAX);
  y_minus.assign(info[3], -INT_MAX);  y_plus.assign(info[3], -INT_MAX);
  z_minus.assign(info[3], -INT_MAX);  z_plus.assign(info[3], -INT_MAX);
  if(coordinates) {
    for(int i=i0; i<imax; i++)  x[i] = coords[k0][j0][i][0];    
    for(int j=j0; j<jmax; j++)  y[j] = coords[k0][j][i0][1];
    for(int k=k0; k<kmax; k++)  z[k] = coords[k][j0][i0][2];
    for(int i=ii0; i<0; i++)    x_minus[i-ii0] = coords[k0][j0][i][0];
    for(int j=jj0; j<0; j++)    y_minus[j-jj0] = coords[k0][j][i0][1];
    for(int k=kk0; k<0; k++)    z_minus[k-kk0] = coords[k][j0][i0][2];
    for(int i=iimax-1; i>=info[0]; i--)  x_plus[i-info[0]] = coords[k0][j0][i][0];
    for(int j=jjmax-1; j>=info[1]; j--)  y_plus[j-info[1]] = coords[k0][j][i0][1];
    for(int k=kkmax-1; k>=info[2]; k--)  z_plus[k-info[2]] = coords[k][j0][i0][2];
  }

  MPI_Allreduce(MPI_IN_PLACE, x.data(), x.size(), MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(MPI_IN_PLACE, y.data(), y.size(), MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(MPI_IN_PLACE, z.data(), z.size(), MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(MPI_IN_PLACE, x_minus.data(), x_minus.size(), MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(MPI_IN_PLACE, y_minus.data(), y_minus.size(), MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(MPI_IN_PLACE, z_minus.data(), z_minus.size(), MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(MPI_IN_PLACE, x_plus.data(), x_plus.size(), MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(MPI_IN_PLACE, y_plus.data(), y_plus.size(), MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(MPI_IN_PLACE, z_plus.data(), z_plus.size(), MPI_DOUBLE, MPI_MAX, comm);

  if(coordinates)
    coordinates->RestoreDataPointerToLocalVector();

}

//-------------------------------------------------------------------------------
// Setup data transfer for two meshes that match exactly --- that is, the overlapped region
// has no gaps.
void 
MeshMatcher::SetupTransferExactMatch(vector<double>& x1, vector<double>& y1, vector<double>& z1,
                                     vector<double>& x1_minus, vector<double>& x1_plus, //coords of ghosts
                                     vector<double>& y1_minus, vector<double>& y1_plus, //coords of ghosts
                                     vector<double>& z1_minus, vector<double>& z1_plus, //coords of ghosts
                                     vector<double>& x2, vector<double>& y2, vector<double>& z2,
                                     vector<double>& x2_minus, vector<double>& x2_plus, //coords of ghosts
                                     vector<double>& y2_minus, vector<double>& y2_plus, //coords of ghosts
                                     vector<double>& z2_minus, vector<double>& z2_plus, //coords of ghosts
                                     bool mesh1_owner, int ii0_1, int jj0_1, int kk0_1,
                                     int iimax_1, int jjmax_1, int kkmax_1,
                                     bool mesh2_owner, int ii0_2, int jj0_2, int kk0_2,
                                     int iimax_2, int jjmax_2, int kkmax_2)
{
  int inactive = -999;

  int ghost_width_1 = x1_minus.size();
  int ghost_width_2 = x2_minus.size();


  //----------------------------------------------------------
  // Step 0. Find a mesh-specific length tolerance
  //----------------------------------------------------------
  double eps = DBL_MAX;
  for(int i=0; i<(int)x1.size()-1; i++)
    eps = std::min(eps, x1[i+1]-x1[i]);
  for(int i=0; i<(int)x2.size()-1; i++)
    eps = std::min(eps, x2[i+1]-x2[i]);
  eps *= 1.0e-6;


  //----------------------------------------------------------
  // Step 1. Gather information about the ownership of Mesh1
  //----------------------------------------------------------
  int* mesh1_ownership = new int[size]; 
  mesh1_ownership[rank] = mesh1_owner ? 1 : 0;
  MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, mesh1_ownership, 1, MPI_INT, comm);


  //----------------------------------------------------------
  // Step 2. Gather information about the ownership of Mesh2
  //----------------------------------------------------------
  int* mesh2_indices = new int[6*size]; //(ii0,jj0,kk0,iimax,jjmax,kkmax)
  if(mesh2_owner){
    mesh2_indices[6*rank]   = ii0_2;
    mesh2_indices[6*rank+1] = jj0_2;
    mesh2_indices[6*rank+2] = kk0_2;
    mesh2_indices[6*rank+3] = iimax_2;
    mesh2_indices[6*rank+4] = jjmax_2;
    mesh2_indices[6*rank+5] = kkmax_2;
  } else
    for(int p=0; p<6; p++)
      mesh2_indices[6*rank+p] = inactive;
  MPI_Allgather(MPI_IN_PLACE, 6, MPI_INT, mesh2_indices, 6, MPI_INT, comm);


  //----------------------------------------------------------
  // Step 3. Each sender (Mesh1 owner) prepares packages for
  //         the receivers (Mesh2 owners)
  //----------------------------------------------------------
  //internal variables to inform the receiver of the mapped nodes
  map<int, vector<int> > dest_i, dest_j, dest_k;
  
  if(mesh1_owner) {

    Int3 myid0(ii0_1, jj0_1, kk0_1);
    Int3 myidmax(iimax_1, jjmax_1, kkmax_1);
    Vec3D myxyz0, myxyzmax;
    myxyz0[0] = myid0[0]<0 ? x1_minus[ghost_width_1+myid0[0]] : x1[myid0[0]];
    myxyz0[1] = myid0[1]<0 ? y1_minus[ghost_width_1+myid0[1]] : y1[myid0[1]];
    myxyz0[2] = myid0[2]<0 ? z1_minus[ghost_width_1+myid0[2]] : z1[myid0[2]];
    myxyzmax[0] = (myidmax[0]-1)>=(int)x1.size() ? x1_plus[myidmax[0]-1-x1.size()] : x1[myidmax[0]-1];
    myxyzmax[1] = (myidmax[1]-1)>=(int)y1.size() ? y1_plus[myidmax[1]-1-y1.size()] : y1[myidmax[1]-1];
    myxyzmax[2] = (myidmax[2]-1)>=(int)z1.size() ? z1_plus[myidmax[2]-1-z1.size()] : z1[myidmax[2]-1];

//    fprintf(stdout,"[%d] Owner of mesh 1 [%d,%d,%d] (%e,%e,%e) -> [%d,%d,%d] (%e,%e,%e)\n", rank, ii0_1, jj0_1, kk0_1,
//            myxyz0[0], myxyz0[1], myxyz0[2], iimax_1, jjmax_1, kkmax_1, myxyzmax[0], myxyzmax[1], myxyzmax[2]);

    for(int proc=0; proc<size; proc++) {
      if(mesh2_indices[6*proc] == inactive)
        continue;

      Int3 id0(mesh2_indices[6*proc], mesh2_indices[6*proc+1], mesh2_indices[6*proc+2]);
      Int3 idmax(mesh2_indices[6*proc+3], mesh2_indices[6*proc+4], mesh2_indices[6*proc+5]);
      Vec3D xyz0, xyzmax;

      xyz0[0] = id0[0]<0 ? x2_minus[ghost_width_2+id0[0]] : x2[id0[0]];
      xyz0[1] = id0[1]<0 ? y2_minus[ghost_width_2+id0[1]] : y2[id0[1]];
      xyz0[2] = id0[2]<0 ? z2_minus[ghost_width_2+id0[2]] : z2[id0[2]];
      xyzmax[0] = (idmax[0]-1)>=(int)x2.size() ? x2_plus[idmax[0]-1-x2.size()] : x2[idmax[0]-1];
      xyzmax[1] = (idmax[1]-1)>=(int)y2.size() ? y2_plus[idmax[1]-1-y2.size()] : y2[idmax[1]-1];
      xyzmax[2] = (idmax[2]-1)>=(int)z2.size() ? z2_plus[idmax[2]-1-z2.size()] : z2[idmax[2]-1];

//      fprintf(stdout,"[%d] is looking at [%d], owner of mesh 2: [%d,%d,%d] (%e,%e,%e) -> [%d,%d,%d] (%e,%e,%e)\n", rank, proc, id0[0], 
//              id0[1], id0[2], xyz0[0], xyz0[1], xyz0[2], idmax[0], idmax[1], idmax[2], xyzmax[0], 
//              xyzmax[1], xyzmax[2]);
      bool overlap = true;
      for(int p=0; p<3; p++)
        if(myxyzmax[p] < xyz0[p]-eps || myxyz0[p] > xyzmax[p]+eps) {
          overlap = false;
          break;
        }

      if(!overlap)
        continue;


      //------------------------------------------------
      // Find the overlapping between the current proc (Mesh1) and "proc" (Mesh2)
      //------------------------------------------------
      send_i[proc] = vector<int>();
      send_j[proc] = vector<int>();
      send_k[proc] = vector<int>();
      dest_i[proc] = vector<int>();
      dest_j[proc] = vector<int>();
      dest_k[proc] = vector<int>();
      for(int dir=0; dir<3; dir++) {

        vector<int>&           send(dir==0 ? send_i[proc] : (dir==1 ? send_j[proc] : send_k[proc]));
        vector<int>&           dest(dir==0 ? dest_i[proc] : (dir==1 ? dest_j[proc] : dest_k[proc]));
        vector<double>&          s1(dir==0 ? x1           : (dir==1 ?           y1 : z1));
        vector<double>&    s1_minus(dir==0 ? x1_minus     : (dir==1 ?     y1_minus : z1_minus));
        vector<double>&     s1_plus(dir==0 ? x1_plus      : (dir==1 ?      y1_plus : z1_plus));
        vector<double>&          s2(dir==0 ? x2           : (dir==1 ?           y2 : z2));
        vector<double>&    s2_minus(dir==0 ? x2_minus     : (dir==1 ?     y2_minus : z2_minus));
        vector<double>&     s2_plus(dir==0 ? x2_plus      : (dir==1 ?      y2_plus : z2_plus));

        double s, sprime;
        for(int i=id0[dir]; i<idmax[dir]; i++) { 
          if(i<0)                s = s2_minus[ghost_width_2+i];
          else if(i>=(int)s2.size())  s = s2_plus[i-s2.size()];
          else                   s = s2[i];

          for(int i1=myid0[dir]; i1<myidmax[dir]; i1++) {
            if(i1<0)                sprime = s1_minus[ghost_width_1+i1];
            else if(i1>=(int)s1.size())  sprime = s1_plus[i1-s1.size()];
            else                    sprime = s1[i1]; 

            if(fabs(sprime-s)<eps) { //found an exact match
              send.push_back(i1);
              dest.push_back(i);
              break;
            }
          }
        }

//        if(dir==0) fprintf(stdout,"[%d] send_i to [%d]: size = %d.\n", rank, proc, (int)send.size());
//        else if(dir==1) fprintf(stdout,"[%d] send_j to [%d]: size = %d.\n", rank, proc, (int)send.size());
//        else fprintf(stdout,"[%d] send_k to [%d]: size = %d.\n", rank, proc, (int)send.size());
        
        assert(!send.empty());

      }
    }   

  }
  

  //----------------------------------------------------------
  // Step 4. Set up communication channels.
  //----------------------------------------------------------
  int* package_size = new int[size];

  // loop through all the processor cores as senders
  for(int sender = 0; sender < size; sender++) {

    if(!mesh1_ownership[sender])
      continue; //nothing to send

    //send number of matched nodes to everyone (including itself)
    for(int dir = 0; dir < 3; dir++) {//dimension-by-dimension
      int num_nodes_to_receive = 0;
      MPI_Request recv_request;
      MPI_Irecv(&num_nodes_to_receive, 1, MPI_INT, sender, sender, comm, &recv_request);
      if(rank == sender) { 
        map<int, vector<int> >& dest(dir==0 ? dest_i : (dir==1 ? dest_j : dest_k));
        for(int receiver = 0; receiver < size; receiver++) {
          auto it = dest.find(receiver);
          package_size[receiver] = (it==dest.end()) ? 0 : it->second.size();
          //fprintf(stdout,"[%d] dir = %d, I will send %d nodes to [%d].\n", rank, dir, package_size[receiver], receiver);
          MPI_Send(&(package_size[receiver]), 1, MPI_INT, receiver, sender, comm);
        }
      }  
      MPI_Wait(&recv_request, MPI_STATUSES_IGNORE);
      //fprintf(stdout,"[%d] dir = %d, I will receive %d nodes from [%d].\n", rank, dir, num_nodes_to_receive, sender);

      //send node-to-node mappings to everyone (including itself)
      if(num_nodes_to_receive>0) {
        map<int, vector<int> >& recv(dir==0 ? recv_i : (dir==1 ? recv_j : recv_k));
        recv[sender] = vector<int>(num_nodes_to_receive, 0);
        MPI_Irecv(recv[sender].data(), num_nodes_to_receive, MPI_INT, sender, sender, comm, &recv_request);
      }
      if(rank == sender) { 
        map<int, vector<int> >& dest(dir==0 ? dest_i : (dir==1 ? dest_j : dest_k));
        for(auto it = dest.begin(); it != dest.end(); it++) {
          int receiver = it->first;
          MPI_Send(it->second.data(), it->second.size(), MPI_INT, receiver, sender, comm);
        }
      }  
      if(num_nodes_to_receive>0)
        MPI_Wait(&recv_request, MPI_STATUSES_IGNORE);
    }

    //"Sender" creates send buffers
    if(rank == sender) {
      for(int receiver = 0; receiver < size; receiver++) {
        auto it = send_i.find(receiver);
        if(it == send_i.end())
          continue;

        int N = send_i[receiver].size()*send_j[receiver].size()*send_k[receiver].size();
        send_buffer_dim1[receiver] = vector<double>(N, 0.0);
        send_buffer_dim2[receiver] = vector<double>(2*N, 0.0);
        send_buffer_dim3[receiver] = vector<double>(3*N, 0.0);
      }
    }

    //each receiver should create a receive buffer corresponding to this "Sender"
    auto it = recv_i.find(sender); 
    if(it != recv_i.end()) {
      int N = recv_i[sender].size()*recv_j[sender].size()*recv_k[sender].size();
      recv_buffer_dim1[sender] = vector<double>(N, 0.0);
      recv_buffer_dim2[sender] = vector<double>(2*N, 0.0);
      recv_buffer_dim3[sender] = vector<double>(3*N, 0.0);
    }

  }

  
  //----------------------------------------------------------
  // Step 5. Verification
  //----------------------------------------------------------
/*  print("Mesh 1: x-coords:\n");
  for(int i=0; i<x1.size(); i++)
    print("%d    %e.\n", i, x1[i]);
  print("Mesh 1: y-coords:\n");
  for(int i=0; i<y1.size(); i++)
    print("%d    %e.\n", i, y1[i]);
  print("Mesh 1: z-coords:\n");
  for(int i=0; i<z1.size(); i++)
    print("%d    %e.\n", i, z1[i]);
  print("Mesh 2: x-coords:\n");
  for(int i=0; i<x2.size(); i++)
    print("%d    %e.\n", i, x2[i]);
  print("Mesh 2: y-coords:\n");
  for(int i=0; i<y2.size(); i++)
    print("%d    %e.\n", i, y2[i]);
  print("Mesh 2: z-coords:\n");
  for(int i=0; i<z2.size(); i++)
    print("%d    %e.\n", i, z2[i]);
*/
  assert(send_i.size()==send_j.size() && send_j.size()==send_k.size() &&
         send_k.size()==send_buffer_dim1.size() && 
         send_buffer_dim1.size()==send_buffer_dim2.size() &&
         send_buffer_dim2.size()==send_buffer_dim3.size());
  assert(recv_i.size()==recv_j.size() && send_j.size()==send_k.size() &&
         recv_k.size()==recv_buffer_dim1.size() &&
         recv_buffer_dim1.size()==recv_buffer_dim2.size() &&
         recv_buffer_dim2.size()==recv_buffer_dim3.size());

  auto it_i = send_i.begin();
  auto it_j = send_j.begin();
  auto it_k = send_k.begin();
  while(it_i != send_i.end()) {
    int receiver = it_i->first;
    assert(receiver == it_j->first && receiver == it_k->first);
    assert(send_buffer_dim1[receiver].size() == it_i->second.size()*it_j->second.size()*it_k->second.size());

//    fprintf(stdout,"[%d]: Sending %d points to [%d]. i = [%d, %d], j = [%d, %d], k = [%d, %d].\n", 
//            rank, (int)send_buffer_dim1[receiver].size(), receiver,
//            it_i->second.front(), it_i->second.back(),
//            it_j->second.front(), it_j->second.back(),
//            it_k->second.front(), it_k->second.back());

    it_i++;  it_j++;  it_k++;
  }

  it_i = recv_i.begin();
  it_j = recv_j.begin();
  it_k = recv_k.begin();
  while(it_i != recv_i.end()) {
    int sender = it_i->first;
    assert(sender == it_j->first && sender == it_k->first);
    assert(recv_buffer_dim1[sender].size() == it_i->second.size()*it_j->second.size()*it_k->second.size());

//    fprintf(stdout,"[%d]: Receiving %d points from [%d]. i = [%d, %d], j = [%d, %d], k = [%d, %d].\n", 
//            rank, (int)recv_buffer_dim1[sender].size(), sender,
//            it_i->second.front(), it_i->second.back(),
//            it_j->second.front(), it_j->second.back(),
//            it_k->second.front(), it_k->second.back());

    it_i++;  it_j++;  it_k++;
  }


  //----------------------------------------------------------
  // Step 6. Clean-up
  //----------------------------------------------------------
  delete [] mesh1_ownership;
  delete [] mesh2_indices;
  delete [] package_size;

}

//-------------------------------------------------------------------------------

void
MeshMatcher::Transfer(SpaceVariable3D* V1, SpaceVariable3D* V2)
{
  if(mesh1_owner) assert(V1); 
  if(mesh2_owner) assert(V2); 
    
  int dim = 1;
  if(mesh1_owner) dim = V1->NumDOF();
  else if(mesh2_owner) dim = V2->NumDOF();
  // if neither mesh1_owner nor mesh2_owner, dim does not matter (set to 1)
  
  if(dim<1 || dim>3) 
    fprintf(stdout,"*** Error: Cannot transfer %d-dimensional data between two meshes.\n", dim);

  map<int,vector<double> >& send_buffer(dim==1 ? send_buffer_dim1 : (dim==2 ? send_buffer_dim2 : send_buffer_dim3));
  map<int,vector<double> >& recv_buffer(dim==1 ? recv_buffer_dim1 : (dim==2 ? recv_buffer_dim2 : recv_buffer_dim3));
   
  // non-blocking send
  double*** v1 = V1 ? V1->GetDataPointer() : NULL;
  vector<MPI_Request> send_requests;
  for(auto it = send_buffer.begin(); it != send_buffer.end(); it++) {
    int receiver = it->first;
    vector<int>& sendi(send_i[receiver]);
    vector<int>& sendj(send_j[receiver]);
    vector<int>& sendk(send_k[receiver]);
    vector<double>& buffer = it->second;

    // fill buffer
    int counter = 0;
    for(auto it_k = sendk.begin(); it_k != sendk.end(); it_k++) 
      for(auto it_j = sendj.begin(); it_j != sendj.end(); it_j++) 
        for(auto it_i = sendi.begin(); it_i != sendi.end(); it_i++)
          for(int p=0; p<dim; p++) {
            buffer[counter] = v1[*it_k][*it_j][*it_i*dim+p];
            counter++; 
          }
    assert(counter == (int)buffer.size());

    // send package to "receiver" (if not self)
    if(receiver != rank) {
      send_requests.push_back(MPI_Request());
      MPI_Isend(buffer.data(), buffer.size(), MPI_DOUBLE, receiver, rank, comm, &(send_requests.back()));
    }
  }
  if(V1) 
    V1->RestoreDataPointerToLocalVector();


  // non-blocking recv
  vector<MPI_Request> recv_requests;
  for(auto it = recv_buffer.begin(); it != recv_buffer.end(); it++) {
    int sender = it->first;
    vector<double>& buffer = it->second;
    if(sender != rank) {
      recv_requests.push_back(MPI_Request());
      MPI_Irecv(buffer.data(), buffer.size(), MPI_DOUBLE, sender, sender, comm, &(recv_requests.back()));
    } else { //sender is myself ==> direct copy
      vector<double>& my_send_buffer(send_buffer[rank]);
      assert(my_send_buffer.size() == buffer.size());
      for(int i=0; i<(int)my_send_buffer.size(); i++)
        buffer[i] = my_send_buffer[i];
    }
  }

  // wait
  MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE); //not necessary
  MPI_Waitall(recv_requests.size(), recv_requests.data(), MPI_STATUSES_IGNORE);

  // fill V2 w/ received data
  double*** v2 = V2 ? V2->GetDataPointer() : NULL;
  for(auto it = recv_buffer.begin(); it != recv_buffer.end(); it++) {
    int sender = it->first;
    vector<double>& buffer = it->second;

    vector<int>& recvi(recv_i[sender]);
    vector<int>& recvj(recv_j[sender]);
    vector<int>& recvk(recv_k[sender]);
    int counter = 0;
    for(auto it_k = recvk.begin(); it_k != recvk.end(); it_k++) 
      for(auto it_j = recvj.begin(); it_j != recvj.end(); it_j++) 
        for(auto it_i = recvi.begin(); it_i != recvi.end(); it_i++)
          for(int p=0; p<dim; p++) {
            v2[*it_k][*it_j][*it_i*dim+p] = buffer[counter];
            counter++; 
          }
    assert(counter == (int)buffer.size());
  }
  if(V2) V2->RestoreDataPointerAndInsert();
 
}

//-------------------------------------------------------------------------------


