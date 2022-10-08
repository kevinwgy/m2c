#include<M2CTwinMessenger.h>
#include<cassert>
#include<climits>
#include<algorithm> //std::find

#define M2C_M2C_INITIAL_TAG1 38601
#define M2C_M2C_INITIAL_TAG2 38602

using std::vector;

extern int verbose;

//---------------------------------------------------------------

M2CTwinMessenger::M2CTwinMessenger(IoData &iod_, MPI_Comm &m2c_comm_, MPI_Comm &joint_comm_,
                                   int status_)
              : iod(iod_), m2c_comm(m2c_comm_), joint_comm(joint_comm_),
                coordinates(NULL), ghost_nodes_outer(NULL), global_mesh(NULL)
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
{ }

//---------------------------------------------------------------

void
M2CTwinMessenger::Destroy()
{ }

//---------------------------------------------------------------

void
M2CTwinMessenger::CommunicateBeforeTimeStepping(SpaceVariable3D &coordinates_,
                                                vector<GhostPoint> &ghost_nodes_outer_,
                                                GlobalMeshInfo &global_mesh_)
{

  coordinates       = &coordinates_;
  ghost_nodes_outer = &ghost_nodes_outer_;
  global_mesh       = &global_mesh_;

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
    // Step 4. Find import_nodes
    // -----------------------------------------------
    I AM HERE


  }






    int my_data_count = 3*import_nodes.size();

    vector<int> counts;
    if(m2c_rank == worker) // I am the receiver
      counts.resize(m2c_size, -1);

    // Gather the number of import nodes from each processor
    MPI_Gather(&my_data_count, 1, MPI_INT, counts.data(), 1, MPI_INT, worker, m2c_comm);

    vector<int> displacements;
    int counter = 0;
    if(m2c_rank == worker) {
      displacements.resize(m2c_size,-1);
      for(int i=0; i<displacements.size(); i++) {
        displacements[i] = counter;
        counter += counts[i];
      }
    }

    // Gather the import nodes on "worker"
    vector<Int3> all_data;
    if(m2c_rank == worker) //receiver
      all_data.resize(counter);

    MPI_Gatherv((int*)import_nodes.data(), my_data_count, MPI_INT, (int*)all_data.data(),
                counts.data(), displacements.data(), MPI_INT, worker, m2c_comm);


    // -----------------------------------------------
    // Step 2. Worker passes ALL the import nodes to ALL
    //         processors running the M2C twin
    // -----------------------------------------------
    vector<Vec3D> all_import_nodes_coords;
    if(m2c_rank == worker) {
      all_import_nodes_coords.reserve(all_data.size());
      for(auto&& ijk : all_data)
        all_import_nodes_coords.push_back(global_mesh->GetXYZ(ijk));

      int numFollowerProcs(0);
      MPI_Comm_remote_size(joint_comm, &numFollowerProcs);

      // Send number of nodes
      vector<MPI_Request> send_requests;
      for(int proc = 0; proc < numFollowerProcs; proc++) {
        send_requests.push_back(MPI_Request());
        MPI_Isend(&counter, 1, MPI_INT, proc, M2C_M2C_INITIAL_TAG1, joint_comm, &send_requests.back());
      }
      MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);

      // Send data
      send_requests.clear();
      for(int proc = 0; proc < numFollowerProcs; proc++) {
        send_requests.push_back(MPI_Request());
        MPI_Isend((double*)all_import_nodes_coords, 3*all_import_nodes_coords.size(), MPI_DOUBLE,
                  proc, M2C_M2C_INITIAL_TAG2, joint_comm, &send_requests.back());
      }
      MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);
    } 

  }

  else { //FOLLOWER

    // -----------------------------------------------
    // Step 1. Each processor receives ALL the expected
    //         nodes from the "worker" processor of the
    //         Leader M2C.
    // -----------------------------------------------
    int worker = 0;

    int leader_worker = 0; //must be the same as above
    int numLeaderProcs = 0;
    MPI_Comm_remote_size(joint_comm, &numLeaderProcs);


    MPI_Request recv_request;
    int counter = 0;
    MPI_Irecv(&counter, 1, MPI_INT, leader_worker, M2C_M2C_INITIAL_TAG1, joint_comm, &recv_request);
    MPI_Wait(&recv_request, MPI_STATUSES_IGNORE);

    assert(counter%3 == 0);
    vector<Vec3D> all_nodes_coords(counter/3);
    MPI_Irecv((double*)all_nodes_coords, counter, MPI_DOUBLE, leader_worker, M2C_M2C_INITIAL_TAG2, 
              joint_comm, &recv_request);
    MPI_Wait(&recv_request, MPI_STATUSES_IGNORE);


    // -----------------------------------------------
    // Step 2. Find the location of each ``expected node''
    // -----------------------------------------------
    vector<int> found(all_nodes_coords.size(), m2c_size); //default value: m2c_size
    vector<GhostPoint> export_points;
    Int3 ijk0;
    Vec3D xi;
    for(int i=0; i<all_nodes_coords.size(); i++)
      if(global_mesh->FindElementCoveringPoint(all_nodes_coords[i], ijk0, &xi, true)) {
        found[i] = m2c_rank;
        export_points.push_back(GhostPoint(ijk0, xi));
      }
    } 
    MPI_Allreduce(MPI_IN_PLACE, found.data(), found.size(), MPI_INT, MPI_MIN, m2c_comm);
    for(int i=0; i<found.size(); i++)
      if(found[i] == m2c_size) { //this point is not adopted...
        print_error("*** Error: M2C Follower cannot find point (%e %e %e) in "
                    "the domain.\n", all_nodes_coords[i][0], all_nodes_coords[i][1],
                    all_nodes_coords[i][2]);
        exit_mpi();
      }

    // -----------------------------------------------
    // Step 3. Send ownerships to Leader
    // -----------------------------------------------
    if(m2c_rank == worker) {
      
      vector<int> packages[numLeaderProcs];
      




    }

    // -----------------------------------------------
    // Step 4. Find out the import nodes
    // -----------------------------------------------


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

