#include<M2CTwinMessenger.h>
#include<cassert>
#include<climits>
#include<algorithm> //std::find


using std::vector;

extern int verbose;
extern int INACTIVE_MATERIAL_ID;

//---------------------------------------------------------------

M2CTwinMessenger::M2CTwinMessenger(IoData &iod_, MPI_Comm &m2c_comm_, MPI_Comm &joint_comm_,
                                   int status_)
              : iod(iod_), m2c_comm(m2c_comm_), joint_comm(joint_comm_),
                coordinates(NULL), ghost_nodes_inner(NULL), ghost_nodes_outer(NULL),
                global_mesh(NULL), TMP(NULL), floodfiller(NULL)
{

  MPI_Comm_rank(m2c_comm, &m2c_rank);
  MPI_Comm_size(m2c_comm, &m2c_size);

  assert(status_==1 || status_==2);
  if(status_==1)
    twinning_status = LEADER;
  else
    twinning_status = FOLLOWER;
}

//---------------------------------------------------------------

M2CTwinMessenger::~M2CTwinMessenger()
{
  if(TMP)         delete TMP;
  if(floodfiller) delete floodfiller;
}

//---------------------------------------------------------------

void
M2CTwinMessenger::Destroy()
{
  if(TMP)
    TMP->Destroy();
  if(floodfiller)
    floodfiller->Destroy();
}

//---------------------------------------------------------------

void
M2CTwinMessenger::CommunicateBeforeTimeStepping(SpaceVariable3D &coordinates_, DataManagers3D &dms_
                                                vector<GhostPoint> &ghost_nodes_inner_,
                                                vector<GhostPoint> &ghost_nodes_outer_,
                                                GlobalMeshInfo &global_mesh_)
{

  coordinates       = &coordinates_;
  ghost_nodes_inner = &ghost_nodes_outer_;
  ghost_nodes_outer = &ghost_nodes_outer_;
  global_mesh       = &global_mesh_;

  if(twinning_status = LEADER) {


  } else {// FOLLOWER
    TMP = new SpaceVariable3D(m2c_comm, &(dms_.ghosted1_1dof));  
    floodfiller = new FloodFill(m2c_comm, dms_, ghost_nodes_inner_, ghost_nodes_outer_);
  }



  if(twinning_status == LEADER) {

    int numFollowerProcs(0);
    MPI_Comm_remote_size(joint_comm, &numFollowerProcs);

    // -----------------------------------------------
    // Step 1. Collect and send import nodes (ghost boundary)
    // -----------------------------------------------
    vector<Int3> import_all;
    vector<Vec3D> import_all_coords;
    for(auto&& ghost : ghost_nodes_outer)
      if(ghost.bcType == (int)MeshData::OVERSET && 
         ghost.type_projection == GhostPoint::FACE) {
        import_all.push_back(ghost.ijk);
        import_all_coords.push_back(global_mesh->GetXYZ(ghost.ijk));
      }

    int nNodes_all = import_all.size();

    vector<MPI_Request> send_requests;
    for(int proc = 0; proc < numFollowerProcs; proc++) {
      send_requests.push_back(MPI_Request());
      MPI_Isend(&nNodes_all, 1, MPI_INT, proc, m2c_rank/*tag*/, joint_comm, 
                &send_requests.back());
    }
    MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);
    send_requests.clear();

    if(nNodes_all>0) { //MPI allows sending data of 0 size. But here we avoid it for clarity
      for(int proc = 0; proc < numFollowerProcs; proc++) {
        send_requests.push_back(MPI_Request());
        MPI_Isend((double*)import_all_coords.data(), 3*nNodes_all, MPI_DOUBLE, 
                  proc, m2c_rank/*tag*/, joint_comm, &send_requests.back());
      }
      MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);
      send_requests.clear();
    }


    // -----------------------------------------------
    // Step 2. Receive adoption info from followers
    // -----------------------------------------------
    vector<MPI_Request> recv_requests;
    vector<int> nNodes(numFollowerProcs, 0);
    if(nNodes_all>0) {
      for(int proc = 0; proc < numFollowerProcs; proc++) {
        recv_requests.push_back(MPI_Request);
        MPI_Irecv(&nNodes[proc], 1, MPI_INT, proc, proc, joint_comm, &recv_requests.back());
      }
      MPI_Waitall(recv_requests.size(), recv_requests.data(), MPI_STATUSES_IGNORE);
      recv_requests.clear();
    }

    vector<int> found[numFollowerProcs];
    for(int proc = 0; proc < numFollowerProcs; proc++) {
      if(nNodes[proc].size()>0) {
        found[proc].resize(nNodes[proc].size());
        recv_requests.push_back(MPI_Request());
        MPI_Irecv(found[proc].data(), found[proc].size(), proc, proc, joint_comm,
                  &recv_requests.back()); 
      }
    }

    MPI_Waitall(recv_requests.size(), recv_requests.data(), MPI_STATUSES_IGNORE);
    recv_requests.clear();

    // -----------------------------------------------
    // Step 3. Construct import_nodes; avoid duplicates
    // -----------------------------------------------
    import_nodes.resize(numFollowerProcs);

    std::multimap<int,int> duplicates;
    vector<int> owner(import_all.size(),-1);
    for(int proc = 0; proc < numFollowerProcs; proc++) {
      for(auto&& id : found[proc]) {
        assert(id>=0 && id<import_all.size());
        if(owner[id]<0) {
          owner[id] = proc;
          import_nodes[proc].push_back(import_all[id]);  
        } else //duplicates
          duplicates.insert(make_pair(proc, id));
      }
    }

    // make sure there are no orphans
    for(int i=0; i<owner.size(); i++) {
      if(own == -1) {
        fprintf(stderr,"\033[0;31m*** Error: [M2C-M2C] Node (%d,%d,%d) (%e,%e,%e) is not picked"
                       " up by any follower processor.\033[0m\n", import_all[i][0], 
                       import_all[i][1], import_all[i][2], import_all_coords[i][0],
                       import_all_coords[i][1], import_all_coords[i][2]);
        exit_mpi();
      }
    } 

    // notify followers about duplication, which should be erased.
    vector<int> dups[numFollowerProcs];
    send_requests.clear();
    if(nNodes_all>0) { //MPI allows size = 0. But here, we try to do a bit better...
      for(auto&& nod : duplicates) 
        dups[nod.first].push_back(nod.second);
 
      for(int proc = 0; proc < numFollowerProcs; proc++) {
        send_requests.push_back(MPI_Request());
        int ndup = dups[proc].size();
        MPI_Isend(&ndup, 1, MPI_INT, proc, m2c_rank/*tag*/, joint_comm, &send_requests.back());
      }
      MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);
      send_requests.clear();
    }

    for(int proc = 0; proc < numFollowerProcs; proc++) {
      if(dups[proc].size()>0) {
        send_requests.push_back(MPI_Request());
        MPI_Isend(dups[proc].data(), dups[proc].size(), proc, m2c_rank, joint_comm,
                  &send_requests.back()); 
      }
    }
    MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);
    

    // -----------------------------------------------
    // Step 4. Send global mesh to follower proc 0, who will then broadcast the info
    // -----------------------------------------------
    if(m2c_rank==0) {
      int mysize;

      mysize = global_mesh->x_glob.size();
      MPI_Send(&mysize, 1, MPI_INT, 0, 0, joint_comm);
      MPI_Send(global_mesh->x_glob.data(), mysize, MPI_DOUBLE, 0, 0, joint_comm);  
      mysize = global_mesh->dx_glob.size();
      MPI_Send(&mysize, 1, MPI_INT, 0, 0/*tag*/, joint_comm);
      MPI_Send(global_mesh->dx_glob.data(), mysize, MPI_DOUBLE, 0, 0, joint_comm);  

      mysize = global_mesh->y_glob.size();
      MPI_Send(&mysize, 1, MPI_INT, 0, 0, joint_comm);
      MPI_Send(global_mesh->y_glob.data(), mysize, MPI_DOUBLE, 0, 0, joint_comm);  
      mysize = global_mesh->dy_glob.size();
      MPI_Send(&mysize, 1, MPI_INT, 0, 0, joint_comm);
      MPI_Send(global_mesh->dy_glob.data(), mysize, MPI_DOUBLE, 0, 0, joint_comm);  

      mysize = global_mesh->z_glob.size();
      MPI_Send(&mysize, 1, MPI_INT, 0, 0, joint_comm);
      MPI_Send(global_mesh->z_glob.data(), mysize, MPI_DOUBLE, 0, 0, joint_comm);  
      mysize = global_mesh->dz_glob.size();
      MPI_Send(&mysize, 1, MPI_INT, 0, 0, joint_comm);
      MPI_Send(global_mesh->dz_glob.data(), mysize, MPI_DOUBLE, 0, 0, joint_comm);  

      
    } 
  }

  else { // Follower

    int numLeaderProcs = 0;
    MPI_Comm_remote_size(joint_comm, &numLeaderProcs);

    // -----------------------------------------------
    // Step 1. Receive nodes from leader procs 
    // -----------------------------------------------
    vector<MPI_Request> recv_requests;
    vector<int> nNodes(numLeaderProcs, 0);
    for(int proc = 0; proc < numLeaderProcs; proc++) {
      recv_requests.push_back(MPI_Request());
      MPI_Irecv(&(nNodes[proc]), 1, MPI_INT, proc, proc, joint_comm, &recv_requests.back());
    }
    MPI_Waitall(recv_requests.size(), recv_requests.data(), MPI_STATUSES_IGNORE);
    recv_requests.clear();

    vector<Vec3D> export_points_all[numLeaderProcs];
    for(int proc = 0; proc < numLeaderProcs; proc++) {
      export_points_all[proc].resize(nNodes[proc]);
      if(nNodes[proc]>0) {
        recv_requests.push_back(MPI_Request());
        MPI_Irecv((double*)export_points_all[proc].data(), 3*nNodes[proc], MPI_DOUBLE, 
                  proc, proc, joint_comm, &recv_requests.back());
      }
    }

    MPI_Waitall(recv_requests.size(), recv_requests.data(), MPI_STATUSES_IGNORE);
    recv_requests.clear();

    // -----------------------------------------------
    // Step 2. Locate nodes and send a package to each leader proc
    // -----------------------------------------------
    vector<MPI_Request> send_requests;
    export_points.resize(numLeaderProcs);
    vector<int> found[numLeaderProcs];
    for(int proc = 0; proc < numLeaderProcs; proc++) {

      if(export_points_all[proc].empty())
        continue; //no MPI exchange

      Int3 ijk0;
      Vec3D xi;
      for(int i=0; i<export_points_all[proc].size(); i++) { 
        if(global_mesh->FindElementCoveringPoint(export_points_all[proc][i], ijk0, &xi, true)) {
          found[proc].push_back(i);
          export_points[proc].push_back(GhostPoint(ijk0, xi));
        }
      }

      int count = export_points[proc].size();
      send_requests.push_back(MPI_Request());
      MPI_Isend(&count, 1, MPI_INT, proc, m2c_rank/*tag*/, joint_comm, &send_requests.back());
    }

    MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);
    send_requests.clear();

    for(int proc = 0; proc < numLeaderProcs; proc++) {
      // It is legit to have size = 0 in MPI. It still transfers "metadata", which is unnecessary
      if(found[proc].size()>0) {
        send_requests.push_back(MPI_Request());
        MPI_Isend(found[proc].data(), found[proc].size(), proc, m2c_rank/*tag*/, joint_comm,
                  &send_requests.back()); 
      }
    }
    MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);
    send_requests.clear();

    // -----------------------------------------------
    // Step 3. Receives duplicates (if any) and remove them
    // -----------------------------------------------
    vector<int> dups[numLeaderProcs];
    vector<int> ndup(numLeaderProcs,0);
    recv_requests.clear();
    for(int proc = 0; proc < numLeaderProcs; proc++) {
      if(export_points_all[proc].empty())
        continue; //no MPI exchange

      recv_requests.push_back(MPI_Request());
      MPI_Irecv(&ndup[proc], 1, MPI_INT, proc, proc, joint_comm, &recv_requests.back());
    } 
    MPI_Waitall(recv_requests.size(), recv_requests.data(), MPI_STATUSES_IGNORE);
    recv_requests.clear();

    for(int proc = 0; proc < numLeaderProcs; proc++) {
      if(ndup[proc]>0) {
        dups[proc].resize(ndups[proc]);

        recv_requests.push_back(MPI_Request());
        MPI_Irecv(dups[proc].data(), dups[proc].size(), MPI_INT, proc, proc,
                  joint_comm, &recv_requests.back());
      }
    }
    MPI_Waitall(recv_requests.size(), recv_requests.data(), MPI_STATUSES_IGNORE);
    recv_requests.clear();
    
    //find duplicates and remove them
    for(int proc = 0; proc < numLeaderProcs; proc++) {
      for(dup : dups[proc]) {
        auto it = std::find(found[proc].begin(), found[proc].end(), dup);
        assert(it != found[proc].end());
        export_points.erase(export_points.begin() + (int)(it - found[proc].begin()));
      }
    }
    

    // -----------------------------------------------
    // Step 4. Get info about leader's mesh
    // -----------------------------------------------
    if(m2c_rank==0) { //receives data from leader proc #0, then broadcast
      int mysize;

      // x
      MPI_Recv(&mysize, 1, MPI_INT, 0, 0, joint_comm);
      assert(mysize>0);
      global_mesh_twin.x_glob.resize(mysize);
      MPI_Bcast(&mysize, 1, MPI_INT, 0, m2c_comm);

      MPI_Recv(global_mesh_twin->x_glob.data(), mysize, MPI_DOUBLE, 0, 0, joint_comm);
      MPI_Bcast(global_mesh_twin->x_glob.data(), mysize, MPI_DOUBLE, 0, m2c_comm);

      // dx
      MPI_Recv(&mysize, 1, MPI_INT, 0, 0, joint_comm);
      assert(mysize>0);
      global_mesh_twin.dx_glob.resize(mysize);
      MPI_Bcast(&mysize, 1, MPI_INT, 0, m2c_comm);

      MPI_Recv(global_mesh_twin->dx_glob.data(), mysize, MPI_DOUBLE, 0, 0, joint_comm);
      MPI_Bcast(global_mesh_twin->dx_glob.data(), mysize, MPI_DOUBLE, 0, m2c_comm);

      // y
      MPI_Recv(&mysize, 1, MPI_INT, 0, 0, joint_comm);
      assert(mysize>0);
      global_mesh_twin.y_glob.resize(mysize);
      MPI_Bcast(&mysize, 1, MPI_INT, 0, m2c_comm);

      MPI_Recv(global_mesh_twin->y_glob.data(), mysize, MPI_DOUBLE, 0, 0, joint_comm);
      MPI_Bcast(global_mesh_twin->y_glob.data(), mysize, MPI_DOUBLE, 0, m2c_comm);

      // dy
      MPI_Recv(&mysize, 1, MPI_INT, 0, 0, joint_comm);
      assert(mysize>0);
      global_mesh_twin.dy_glob.resize(mysize);
      MPI_Bcast(&mysize, 1, MPI_INT, 0, m2c_comm);

      MPI_Recv(global_mesh_twin->dy_glob.data(), mysize, MPI_DOUBLE, 0, 0, joint_comm);
      MPI_Bcast(global_mesh_twin->dy_glob.data(), mysize, MPI_DOUBLE, 0, m2c_comm);

      // z
      MPI_Recv(&mysize, 1, MPI_INT, 0, 0, joint_comm);
      assert(mysize>0);
      global_mesh_twin.z_glob.resize(mysize);
      MPI_Bcast(&mysize, 1, MPI_INT, 0, m2c_comm);

      MPI_Recv(global_mesh_twin->z_glob.data(), mysize, MPI_DOUBLE, 0, 0, joint_comm);
      MPI_Bcast(global_mesh_twin->z_glob.data(), mysize, MPI_DOUBLE, 0, m2c_comm);

      // dz
      MPI_Recv(&mysize, 1, MPI_INT, 0, 0, joint_comm);
      assert(mysize>0);
      global_mesh_twin.dz_glob.resize(mysize);
      MPI_Bcast(&mysize, 1, MPI_INT, 0, m2c_comm);

      MPI_Recv(global_mesh_twin->dz_glob.data(), mysize, MPI_DOUBLE, 0, 0, joint_comm);
      MPI_Bcast(global_mesh_twin->dz_glob.data(), mysize, MPI_DOUBLE, 0, m2c_comm);

    } else {
      int mysize;

      // x
      MPI_Bcast(&mysize, 1, MPI_INT, 0, m2c_comm);
      assert(mysize>0);
      global_mesh_twin.x_glob.resize(mysize);
      MPI_Bcast(global_mesh_twin->x_glob.data(), mysize, MPI_DOUBLE, 0, m2c_comm);
       
      // dx
      MPI_Bcast(&mysize, 1, MPI_INT, 0, m2c_comm);
      assert(mysize>0);
      global_mesh_twin.dx_glob.resize(mysize);
      MPI_Bcast(global_mesh_twin->dx_glob.data(), mysize, MPI_DOUBLE, 0, m2c_comm);
       
      // y
      MPI_Bcast(&mysize, 1, MPI_INT, 0, m2c_comm);
      assert(mysize>0);
      global_mesh_twin.y_glob.resize(mysize);
      MPI_Bcast(global_mesh_twin->y_glob.data(), mysize, MPI_DOUBLE, 0, m2c_comm);
       
      // dy
      MPI_Bcast(&mysize, 1, MPI_INT, 0, m2c_comm);
      assert(mysize>0);
      global_mesh_twin.dy_glob.resize(mysize);
      MPI_Bcast(global_mesh_twin->dy_glob.data(), mysize, MPI_DOUBLE, 0, m2c_comm);
       
      // z
      MPI_Bcast(&mysize, 1, MPI_INT, 0, m2c_comm);
      assert(mysize>0);
      global_mesh_twin.y_glob.resize(mysize);
      MPI_Bcast(global_mesh_twin->y_glob.data(), mysize, MPI_DOUBLE, 0, m2c_comm);
       
      // dz
      MPI_Bcast(&mysize, 1, MPI_INT, 0, m2c_comm);
      assert(mysize>0);
      global_mesh_twin.dz_glob.resize(mysize);
      MPI_Bcast(global_mesh_twin->dz_glob.data(), mysize, MPI_DOUBLE, 0, m2c_comm);
    }

    // Now, global_mesh_twin has been created...


  }

}

//---------------------------------------------------------------

void
M2CTwinMessenger::FirstExchange()
{

}

//---------------------------------------------------------------

void
M2CTwinMessenger::Exchange()
{

}

//---------------------------------------------------------------

void
M2CTwinMessenger::FinalExchange()
{

}

//---------------------------------------------------------------


//---------------------------------------------------------------

