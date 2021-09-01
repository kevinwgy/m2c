#include<CustomCommunicator.h>
using std::vector;

//------------------------------------------------------------------------------------

CustomCommunicator::CustomCommunicator(MPI_Comm& comm_, SpaceVariable3D &V, 
                                       std::vector<Int3> &ghost_nodes)
                  : comm(comm_), dof(V.NumDOF()), ghost_width(V.NumGhostLayers())
{

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  //-------------------------------------
  // Create the communication channels
  //-------------------------------------
  // This routine does not rely on a
  // certain way of mesh partitioning (e.g.,
  // cutting along x, y, z axes as done by
  // PETSc)
  //-------------------------------------
  
  //----------------------------------------------------------------------------------------
  //Step 0: Verify the ghost_nodes are indeed ghost nodes.
  //----------------------------------------------------------------------------------------
  for(auto it = ghost_nodes.begin(); it != ghost_nodes.end(); it ++) {
    int i((*it)[0]), j((*it)[1]), k((*it)[2]);
    if(V.IsHere(i,j,k,true) && !V.IsHere(i,j,k,false)) {
      continue; //good
    } else {
      fprintf(stderr,"*** Error: [Proc %d] Passing an interior node (%d,%d,%d) to CustomCommunicator.\n",
              rank, i,j,k);
      exit(-1);
    }
  }


  //----------------------------------------------------------------------------------------
  //Step 1: We will use V to build the channels. So here we backup the data in V, which
  //        will be applied back to it at the end.
  //----------------------------------------------------------------------------------------
  int ii0, jj0, kk0, iimax, jjmax, kkmax;
  V.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);
  double *backup = new double[(kkmax-kk0)*(jjmax-jj0)*(iimax-ii0)*dof];
  double*** v = V.GetDataPointer();
  int counter = 0;
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++) 
        for(int p=0; p<dof; p++) 
          backup[counter++] = v[k][j][i*dof+p];
  V.RestoreDataPointerToLocalVector();


  //----------------------------------------------------------------------------------------
  //Step 2: Create channels one by one
  //----------------------------------------------------------------------------------------
  V.SetConstantValue(0.0, true);

  for(int proc = 0; proc < size; proc++) {

    //Creating the send/receive channel to populate (selected) ghost nodes in Subdomain #proc

    //2.1. proc sends ghost nodes to everyone else   
    int local_size = 0;
    if(proc == rank) //tell other processors the number of points
      local_size = ghost_nodes.size();
    MPI_Bcast(&local_size, 1, MPI_INT, proc, comm);

    int* nodes = new int[3*local_size];
    if(proc == rank)
      for(int i = 0; i < ghost_nodes.size(); i++) {
        nodes[i*3]   = ghost_nodes[i][0];
        nodes[i*3+1] = ghost_nodes[i][1];
        nodes[i*3+2] = ghost_nodes[i][2];
      }
    MPI_Bcast(&nodes, 3*local_size, MPI_INT, proc, comm);

    //2.2. non-procs figure out if any of the ghost nodes are owned by them. If yes, create
    //     a "send package", and tag V
    v = V.GetDataPointer();
    int package_id = -1;
    if(proc != rank) {
      for(int n = 0; n < local_size; n++) {
        int i(nodes[n*3]), j(nodes[n*3+1]), k(nodes[n*3+2]);

        if(V.IsHere(i,j,k,false)) {
          if(package_id==-1) {
            send_pack.push_back(Package(SEND, proc));
            package_id = send_pack.size() - 1;
          }
          send_pack[package_id].index.push_back(Int3(i,j,k));

          v[k][j][i] = rank + 1; //have to add 1, so that the first proc appears to be 1 (instead of 0)
        }
      }
    }
    V.RestoreDataPointerAndInsert();
    delete [] nodes;

    //2.3. proc figures out who is the owner of each ghost node; create "receive packages"
    v = V.GetDataPointer();
    if(proc == rank) {

      vector<int> recv_package_id(size, -1);

      for(auto it = ghost_nodes.begin(); it != ghost_nodes.end(); it++) {
        int i((*it)[0]), j((*it)[1]), k((*it)[2]);
        assert(v[k][j][i]>0 && v[k][j][i]<=size); //each ghost node must have one (and only one) owner
        int owner = v[k][j][i] - 1;        
        if(recv_package_id[owner] == -1) {
          recv_pack.push_back(Package(RECEIVE, owner)); 
          recv_package_id[owner] = recv_pack.size() - 1;
        } 
        recv_pack[recv_package_id[owner]].index.push_back(Int3(i,j,k));

        v[k][j][i] = 0; //restore for next proc
      }
    }
    else {//restore V for next proc
      if(package_id != -1)
        for(auto it = send_pack[package_id].index.begin(); it != send_pack[package_id].index.end(); it++){
          int i((*it)[0]), j((*it)[1]), k((*it)[2]);
          v[k][j][i] = 0;
        }
    }

    V.RestoreDataPointerToLocalVector();
  }


  //----------------------------------------------------------------------------------------
  //Step 3: Allocate space for buffer in send_pack and recv_pack
  //----------------------------------------------------------------------------------------
  for(auto it = send_pack.begin(); it != send_pack.end(); it++)
    it->buffer.resize(dof*it->index.size());

  for(auto it = recv_pack.begin(); it != recv_pack.end(); it++)
    it->buffer.resize(dof*it->index.size());


  //----------------------------------------------------------------------------------------
  //Step 4: Restore orginal data in V
  //----------------------------------------------------------------------------------------
  v = V.GetDataPointer();
  counter = 0;
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++) 
        for(int p=0; p<dof; p++) 
          v[k][j][i*dof+p] = backup[counter++];
  V.RestoreDataPointerToLocalVector();
  delete [] backup;

}

//------------------------------------------------------------------------------------

void
CustomCommunicator::ExchangeAndInsert(SpaceVariable3D &V)
{
  // verify compatibility
  assert(V.NumDOF()==dof && V.NumGhostLayers()==ghost_width);

  double*** v = V.GetDataPointer();

  // non-blocking send
  vector<MPI_Request> send_requests;
  for(auto it = send_pack.begin(); it != send_pack.end(); it++) {

    int receiver = it->rank;
    vector<double> &buffer(it->buffer);
    vector<Int3> &index(it->index);

    int i,j,k;
    for(int n = 0; n < index.size(); n++) {
      i = index[n][0];
      j = index[n][1];
      k = index[n][2];
      for(int p = 0; p < dof; p++)
        buffer[n*dof+p] = v[k][j][i*dof+p];
    }

    send_requests.push_back();
    MPI_Isend(buffer.data(), buffer.size(), MPI_DOUBLE, receiver, rank, comm, 
              send_requests[send_requests.size()-1]);
  }

  // non-blocking receive
  vector<MPI_Request> recv_requests;
  for(auto it = recv_pack.begin(); it != recv_pack.end(); it ++) {
    int sender = it->rank;
    recv_requests.push_back();
    MPI_Irecv(it->buffer.data(), it->buffer.size(), MPI_DOUBLE, sender, sender, comm,
              recv_requests[recv_requests.size()-1]);
  }

  // wait
  MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE); //not necessary
  MPI_Waitall(recv_requests.size(), recv_requests.data(), MPI_STATUSES_IGNORE);


  //insert received data to vector V
  for(auto it = recv_pack.begin(); it != recv_pack.end(); it++) {
    vector<double> &buffer(it->buffer);
    vector<Int3> &index(it->index);

    int i,j,k;
    for(int n = 0; n < index.size(); n++) {
      i = index[n][0];
      j = index[n][1];
      k = index[n][2];
      for(int p = 0; p < dof; p++)
        v[k][j][i*dof+p] = buffer[n*dof+p];
    }
  }

  // Restore w/o additional MPI comm.
  V.RestoreDataPointerToLocalVector();

}

//------------------------------------------------------------------------------------


