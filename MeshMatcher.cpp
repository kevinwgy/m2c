#include<MeshMatcher.h>
using std::vector;
using std::pair;
using std::map;

//-------------------------------------------------------------------------------

MeshMatcher::MeshMatcher(MPI_Comm* comm_, MPI_Comm* comm1_, MPI_Comm* comm2_, 
                         SpaceVariable3D* coords1_, SpaceVariable3D* coords2_,
                         vector<double>& x1, vector<double>& y1, vector<double>& z1,
                         vector<double>& x2, vector<double>& y2, vector<double>& z2)
           : comm(comm_), comm1(comm1_), comm2(comm2_), 
             coordinates1(coords1_), coordinates2(coords2_)
{

  exact_match = true; //can be extended in future if needed

  assert(comm_); //this cannot be NULL!
  MPI_Comm_rank(*comm, &rank);
  MPI_Comm_size(*comm, &size);

  if(comm1) {
    MPI_Comm_rank(*comm1, &rank1);
    MPI_Comm_size(*comm1, &size1);
  } else 
    rank1 = size1 = -1;

  if(comm2) {
    MPI_Comm_rank(*comm2, &rank2);
    MPI_Comm_size(*comm2, &size2);
  } else 
    rank2 = size2 = -1;

  SetupMesh1ToMesh2Transfer(x1,y1,z1,x2,y2,z2);


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
                                     int iimax_2, int jjmax_2, int kkmax_2,
                                     vector<vector<Int3> > &send_nodes,
                                     vector<vector<Int3> > &receiver_packages,
{
  int inactive = -999;

  int ghost_width_1 = x1_minus.size();
  int ghost_width_2 = x2_minus.size();


  //----------------------------------------------------------
  // Step 0. Find a mesh-specific length tolerance
  //----------------------------------------------------------
  double eps = DBL_MAX;
  for(int i=0; i<x1.size()-1; i++)
    eps = std::min(eps, x1[i+1]-x1[i]);
  for(int i=0; i<x2.size()-1; i++)
    eps = std::min(eps, x2[i+1]-x2[i]);
  eps /= 1.0e-6;


  //----------------------------------------------------------
  // Step 1. Gather information about the ownership of Mesh1
  //----------------------------------------------------------
  int* mesh1_indices = new int[6*size]; //(ii0,jj0,kk0,iimax,jjmax,kkmax)
  if(mesh1_owner){
    mesh1_indices[6*rank]   = ii0_1;
    mesh1_indices[6*rank+1] = jj0_1;
    mesh1_indices[6*rank+2] = kk0_1;
    mesh1_indices[6*rank+3] = iimax_1;
    mesh1_indices[6*rank+4] = jjmax_1;
    mesh1_indices[6*rank+5] = kkmax_1;
  } else
    for(int p=0; p<6; p++)
      mesh1_indices[6*rank+p] = inactive;
  MPI_Allgather(MPI_IN_PLACE, 6, MPI_INT, mesh1_indices, 6, MPI_INT, *comm);


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
  MPI_Allgather(MPI_IN_PLACE, 6, MPI_INT, mesh2_indices, 6, MPI_INT, *comm);


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
    myxyzmax[0] = myidmax[0]>=x1.size() ? x1_plus[myidmax[0]-x1.size()] : x1[myidmax[0]];
    myxyzmax[1] = myidmax[1]>=y1.size() ? y1_plus[myidmax[1]-y1.size()] : y1[myidmax[1]];
    myxyzmax[2] = myidmax[2]>=z1.size() ? z1_plus[myidmax[2]-z1.size()] : z1[myidmax[2]];

    for(int proc=0; proc<size; proc++) {
      if(mesh2_indices[6*proc] == inactive)
        continue;

      Int3 id0(mesh2_indices[6*proc], mesh2_indices[6*proc+1], mesh2_indices[6*proc+2]);
      Int3 idmax(mesh2_indices[6*proc+3], mesh2_indices[6*proc+4], mesh2_indices[6*proc+5]);
      Vec3D xyz0, xyzmax;

      xyz0[0] = id0[0]<0 ? x2_minus[ghost_width_2+id0[0]] : x2[id0[0]];
      xyz0[1] = id0[1]<0 ? y2_minus[ghost_width_2+id0[1]] : y2[id0[1]];
      xyz0[2] = id0[2]<0 ? z2_minus[ghost_width_2+id0[2]] : z2[id0[2]];
      xyzmax[0] = idmax[0]>=x2.size() ? x2_plus[idmax[0]-x2.size()] : x2[idmax[0]];
      xyzmax[1] = idmax[1]>=y2.size() ? y2_plus[idmax[1]-y2.size()] : y2[idmax[1]];
      xyzmax[2] = idmax[2]>=z2.size() ? z2_plus[idmax[2]-z2.size()] : z2[idmax[2]];

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

        int i1 = myid0[dir];
        double s, s1;
        for(int i=id0[dir]; i<idmax[dir]; i++) { 
          if(i<0)                s = s2_minus[ghost_width_2+i];
          else if(i>=s2.size())  s = s2_plus[i-s2.size()];
          else                   s = s2[i];

          while(i1<myidmax[dir]) {
            if(i1<0)                s1 = s1_minus[ghost_width_1+i1];
            else if(i1>=x1.size())  s1 = s1_plus[i1-s1.size()];
            else                    s1 = s1[i1]; 

            if(fabs(s1-s)<eps) { //found an exact match
              send.push_back(i1);
              dest.push_back(i);
              break;
            }

            i1++;
          }

          if(i1>=myidmax[dir])
            break;
        }
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

    if(mesh1_indices[6*sender] == inactive)
      continue; //nothing to send

    //send number of matched nodes to everyone (including itself)
    for(int dir = 0; dir < 3; dir++) {//dimension-by-dimension
      if(rank == sender) { 
        map<int, vector<int> >& dest(dir==0 ? dest_i : (dir==1 ? dest_j : dest_k));
        for(int receiver = 0; receiver < size; receiver++) {
          auto it = dest.find(receiver);
          package_size[receiver] = (it==dest.end()) ? 0 : it->second.size();
          MPI_Send(&(package_size[receiver]), 1, MPI_INT, receiver, sender, *comm);
        }
      }  
      int num_nodes_to_receive = 0;
      MPI_Recv(&num_nodes_to_receive, 1, MPI_INT, sender, sender, *comm);

      //send node-to-node mappings to everyone (including itself)
      if(rank == sender) { 
        map<int, vector<int> >& dest(dir==0 ? dest_i : (dir==1 ? dest_j : dest_k));
        for(auto it = dest.begin(); it != dest.end(); it++) {
          int receiver = it->first;
          MPI_Send(it->second.data(), it->second.size(), MPI_INT, receiver, sender, *comm);
        }
      }  
      if(num_nodes_to_receive>0) {
        map<int, vector<int> >& recv(dir==0 ? recv_i : (dir==1 ? recv_j : recv_k));
        recv[sender] = vector<int>(num_nodes_to_receive, 0);
        MPI_Recv(recv[sender].data(), num_nodes_to_receive, MPI_INT, sender, sender, *comm);
      }
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
  print("Mesh 1: x-coords:\n");
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

    fprintf(stderr,"[%d]: Sending %d points to [%d]. i = [%d, %d), j = [%d, %d), k = [%d, %d).\n", 
            rank, send_buffer_dim1[receiver].size(), receiver,
            it_i->second->front(), it_i->second->back(),
            it_j->second->front(), it_j->second->back(),
            it_k->second->front(), it_k->second->back());

    it_i++;  it_j++;  it_k++;
  }

  it_i = recv_i.begin();
  it_j = recv_i.begin();
  it_k = recv_i.begin();
  while(it_i != recv_i.end()) {
    int sender = it_i->first;
    assert(sender == it_j->first && sender == it_k->first);
    assert(recv_buffer_dim1[receiver].size() == it_i->second.size()*it_j->second.size()*it_k->second.size());

    fprintf(stderr,"[%d]: Receiving %d points from [%d]. i = [%d, %d), j = [%d, %d), k = [%d, %d).\n", 
            rank, recv_buffer_dim1[sender].size(), sender,
            it_i->second->front(), it_i->second->back(),
            it_j->second->front(), it_j->second->back(),
            it_k->second->front(), it_k->second->back());

    it_i++;  it_j++;  it_k++;
  }


  //----------------------------------------------------------
  // Step 6. Clean-up
  //----------------------------------------------------------
  delete [] mesh1_indices;
  delete [] mesh2_indices;
  delete [] package_size;

}

//-------------------------------------------------------------------------------

void
MeshMatcher::SendData(SpaceVariable3D* V1, SpaceVariable3D* V2)
{
  if(mesh1_owner) assert(V1) else assert(V1==NULL); 
  if(mesh2_owner) assert(V2) else assert(V2==NULL); 
    
  int dim = 1;
  if(mesh1_owner) dim = V1->NumDOF();
  else if(mesh2_owner) dim = V2->NumDOF();
  // if neither mesh1_owner nor mesh2_owner, dim does not matter (set to 1)
  
  if(dim<1 || dim>3) 
    fprintf(stderr,"*** Error: Cannot transfer %d-dimensional data between two meshes.\n", dim);

  map<int,vector<double> > send_buffer(dim==1 ? send_buffer_dim1 : (dim==2 ? send_buffer_dim2 : send_buffer_dim3));
  map<int,vector<double> > recv_buffer(dim==1 ? recv_buffer_dim1 : (dim==2 ? recv_buffer_dim2 : recv_buffer_dim3));
   
  // non-blocking send
  double*** v1 = V1 ? V1->GetDataPointer() : NULL;
  vector<MPI_Request> send_requests;
  for(auto it = send_buffer.begin(); it != send_buffer.end(); it++) {
    int receiver = it->first;
    vector<int> &sendi(send_i[receiver]);
    vector<int> &sendj(send_j[receiver]);
    vector<int> &sendk(send_k[receiver]);
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
    assert(counter == buffer.size());

    // send package to "receiver"
    send_requests.push_back(MPI_Request());
    MPI_Isend(buffer.data(), buffer.size(), MPI_DOUBLE, receiver, rank, *comm, &(send_requests.back());
  }
  if(V1) V1.RestoreDataPointerToLocalVector();


  // non-blocking recv
  vector<MPI_Request> recv_requests;
  for(auto it = recv_buffer.begin(); it != recv_buffer.end(); it++) {
    int sender = it->first;
    vector<double>& buffer = it->second;
    recv_requests.push_back(MPI_Request());
    MPI_Irecv(buffer.data(), buffer.size(), MPI_DOUBLE, sender, sender, *comm, &(recv_requests.bac()));
  }

  // wait
  MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE); //not necessary
  MPI_Waitall(recv_requests.size(), recv_requests.data(), MPI_STATUSES_IGNORE);

  // fill V2 w/ received data
  double*** v2 = V2 ? V2->GetDataPointer() : NULL;
  for(auto it = recv_buffer.begin(); it != recv_buffer.end(); it++) {
    int sender = it->first;
    vector<double>& buffer = it->second;

    vector<int> &recvi(recv_i[sender]);
    vector<int> &recvj(recv_j[sender]);
    vector<int> &recvk(recv_k[sender]);
    int counter = 0;
    for(auto it_k = recvk.begin(); it_k != recvk.end(); it_k++) 
      for(auto it_j = recvj.begin(); it_j != recvj.end(); it_j++) 
        for(auto it_i = recvi.begin(); it_i != recvi.end(); it_i++)
          for(int p=0; p<dim; p++) {
            v2[*it_k][*it_j][*it_i*dim+p] = buffer[counter];
            counter++; 
          }
    assert(counter == buffer.size());
  }
  if(V2) V2->RestoreDataPointerAndInsert();
 
}

//-------------------------------------------------------------------------------


