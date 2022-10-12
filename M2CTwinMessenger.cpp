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
                global_mesh(NULL), TMP(NULL), TMP3(NULL), Color(NULL)floodfiller(NULL)
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
  if(TMP3)        delete TMP3;
  if(Color)       delete Color;
  if(floodfiller) delete floodfiller;
}

//---------------------------------------------------------------

void
M2CTwinMessenger::Destroy()
{
  if(TMP)
    TMP->Destroy();
  if(TMP3)
    TMP3->Destroy();
  if(Color)
    Color->Destroy();
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
    TMP3 = new SpaceVariable3D(m2c_comm, &(dms_.ghosted1_3dof)); //dof: 3 
    Color = new SpaceVariable3D(m2c_comm, &(dms_.ghosted1_1dof));  
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
        recv_requests.push_back(MPI_Request());
        MPI_Irecv(&nNodes[proc], 1, MPI_INT, proc, proc, joint_comm, &recv_requests.back());
      }
      MPI_Waitall(recv_requests.size(), recv_requests.data(), MPI_STATUSES_IGNORE);
      recv_requests.clear();
    }

    vector<int> found[numFollowerProcs];
    for(int proc = 0; proc < numFollowerProcs; proc++) {
      if(nNodes[proc]>0) {
        found[proc].resize(nNodes[proc]);
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
      if(owner[i] == -1) {
        fprintf(stderr,"\033[0;31m*** Error: [M2C-M2C] Node (%d,%d,%d) (%e,%e,%e) is not picked"
                       " up by any follower processor.\033[0m\n", import_all[i][0], 
                       import_all[i][1], import_all[i][2], import_all_coords[i][0],
                       import_all_coords[i][1], import_all_coords[i][2]);
        exit(-1);
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

      for(int proc = 0; proc < numFollowerProcs; proc++) {
        if(dups[proc].size()>0) {
          send_requests.push_back(MPI_Request());
          MPI_Isend(dups[proc].data(), dups[proc].size(), proc, m2c_rank, joint_comm,
                    &send_requests.back()); 
        }
      }
      MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);
      send_requests.clear(); 
    }


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

      // receive the location(s) of the overset boundaries
      vector<int> overset_boundary(6, 0); //xmin, xmax, ymin, ymax, zmin, zmax. 1 means yes.

      if(iod.mesh.bc_x0   == MeshData::OVERSET)  overset_boundary[0] = 1;
      if(iod.mesh.bc_xmax == MeshData::OVERSET)  overset_boundary[1] = 1;
      if(iod.mesh.bc_y0   == MeshData::OVERSET)  overset_boundary[2] = 1;
      if(iod.mesh.bc_ymax == MeshData::OVERSET)  overset_boundary[3] = 1;
      if(iod.mesh.bc_z0   == MeshData::OVERSET)  overset_boundary[4] = 1;
      if(iod.mesh.bc_zmax == MeshData::OVERSET)  overset_boundary[5] = 1;

      MPI_Send(overset_boundary.data(), overset_boundary.size(), MPI_INT, 0, 0, joint_comm);
      
    } 


    // -----------------------------------------------
    // Step 5. Receive nodes from follower procs
    // -----------------------------------------------
    recv_requests.clear();
    nNodes.assign(numFollowerProcs, 0);
    for(int proc=0; proc<numFollowerProcs; proc++) {
      recv_requests.push_back(MPI_Request());
      MPI_Irecv(&(nNodes[proc]), 1, MPI_INT, proc, proc, joint_comm,
                &recv_requests.back());
    }
    MPI_Waitall(recv_requests.size(), recv_requests.data(), MPI_STATUSES_IGNORE);
    recv_requests.clear();

    vector<Vec3D> export_points_all[numFollowerProcs];
    for(int proc=0; proc<numFollowerProcs; proc++) {
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
    // Step 6. Locate nodes and send a package to each follower proc
    // -----------------------------------------------
    export_points.resize(numFollowerProcs);
    for(int proc=0; proc<numFollowerProcs; proc++) {

      found[proc].clear();

      if(export_points_all[proc].empty())
        continue;

      Int3 ijk0;
      Vec3D xi;      
      for(int i=0; i<export_points_all[proc].size(); i++) {
        if(global_mesh->FindElementCoveringPoint(export_points_all[proc][i], ijk0, &xi, true)) {
          if(coordinates->IsHere(ijk0[0], ijk[1], ijk[2], true/*include_ghost*/) &&
             coordinates->IsHere(ijk0[0]+1,ijk[1]+1,ijk[2]+1, true/*include_ghost*/)) {
             found[proc].push_back(i);
             export_points[proc].push_back(GhostPoint(ijk0,xi));
          } 
        }
      }

      int count = export_points[proc].size();
      send_requests.push_back(MPI_Request());
      MPI_Isend(&count, 1, MPI_INT, proc, m2c_rank/*tag*/, joint_comm, &send_requests.back());
    }
    MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);
    send_requests.clear();

    for(int proc=0; proc<numLeaderProcs; proc++) {
      if(found[proc].size()>0) {
        send_requests.push_back(MPI_Request());
        MPI_Isend(found[proc].data(), found[proc].size(), proc, m2c_rank, joint_comm,
                  &send_requests.back());
      }
    }
    MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);
    send_requests.clear();


    // -----------------------------------------------
    // Step 7. Receive duplicates (if any) and remove them
    // -----------------------------------------------
    for(int proc=0; proc<numFollowerProcs; proc++)
      dups[proc].clear();
    vector<int> ndup(numFollowerProcs,0);
    recv_requests.clear();

    for(int proc=0; proc<numFollowerProcs; proc++) {
      if(export_points_all[proc].empty())
        continue;

      recv_requests.push_back(MPI_Request());
      MPI_Irecv(&ndup[proc], 1, MPI_INT, proc, proc, joint_comm, &recv_requests.back());
    }
    MPI_Waitall(recv_requests.size(), recv_requests.data(), MPI_STATUSES_IGNORE);
    recv_requests.clear();

    for(int proc=0; proc<numFollowerProcs; proc++) {
      if(ndup[proc]>0) {
        dups[proc].resize(ndup[proc]);

        recv_requests.push_back(MPI_Request());
        MPI_Irecv(dups[proc].data(), dups[proc].size(), MPI_INT, proc, proc,
                  joint_comm, &recv_requests.back());
      }
    }
    MPI_Waitall(recv_requests.size(), recv_requests.data(), MPI_STATUSES_IGNORE);
    recv_requests.clear();

    //find duplicates and remove them
    for(int proc=0; proc<numFollowerProcs; proc++) {
      for(auto&& dup : dups[proc]) {
        auto it = std::find(found[proc].begin(), found[proc].end(), dup);
        assert(it != found[proc].end());
        export_points[proc].erase(export_points[proc].begin() + (int)(it - found[proc].begin()));
      }
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
          if(coordinates->IsHere(ijk0[0],ijk0[1],ijk0[2], true/*include_ghost*/) &&
             coordinates->IsHere(ijk0[0]+1,ijk0[1]+1,ijk0[2]+1, true/*include_ghost*/)) {
            found[proc].push_back(i);
            export_points[proc].push_back(GhostPoint(ijk0, xi));
          }
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
      for(auto&& dup : dups[proc]) {
        auto it = std::find(found[proc].begin(), found[proc].end(), dup);
        assert(it != found[proc].end());
        export_points[proc].erase(export_points[proc].begin() + (int)(it - found[proc].begin()));
      }
    }
    
    //Tag nodes involved in interpolations ("red boxes" in KW's notes)
    double*** tag = TMP->GetDataPointer(); //value is 0 by default
    for(int proc = 0; proc < numLeaderProcs; proc++)
      for(auto&& gp : export_points[proc]) 
        for(int k=gp.ijk[2]; k<=gp.ijk[2]+1; k++)
          for(int j=gp.ijk[1]; j<=gp.ijk[1]+1; j++)
            for(int i=gp.ijk[0]; i<=gp.ijk[0]+1; i++)
              tag[k][j][i] = 1;
    TMP->RestoreDataPointerAndInsert();


    // -----------------------------------------------
    // Step 4. Get info about leader's mesh
    // -----------------------------------------------
    vector<int> overset_boundary(6, 0); //xmin, xmax, ymin, ymax, zmin, zmax, 1 means yes
  
    if(m2c_rank==0) { //receives data from leader proc #0, then broadcast
      int mysize;

      // x
      MPI_Recv(&mysize, 1, MPI_INT, 0, 0, joint_comm);
      assert(mysize>0);
      global_mesh_twin.x_glob.resize(mysize);
      MPI_Bcast(&mysize, 1, MPI_INT, 0, m2c_comm);

      MPI_Recv(global_mesh_twin.x_glob.data(), mysize, MPI_DOUBLE, 0, 0, joint_comm);
      MPI_Bcast(global_mesh_twin.x_glob.data(), mysize, MPI_DOUBLE, 0, m2c_comm);

      // dx
      MPI_Recv(&mysize, 1, MPI_INT, 0, 0, joint_comm);
      assert(mysize>0);
      global_mesh_twin.dx_glob.resize(mysize);
      MPI_Bcast(&mysize, 1, MPI_INT, 0, m2c_comm);

      MPI_Recv(global_mesh_twin.dx_glob.data(), mysize, MPI_DOUBLE, 0, 0, joint_comm);
      MPI_Bcast(global_mesh_twin.dx_glob.data(), mysize, MPI_DOUBLE, 0, m2c_comm);

      // y
      MPI_Recv(&mysize, 1, MPI_INT, 0, 0, joint_comm);
      assert(mysize>0);
      global_mesh_twin.y_glob.resize(mysize);
      MPI_Bcast(&mysize, 1, MPI_INT, 0, m2c_comm);

      MPI_Recv(global_mesh_twin.y_glob.data(), mysize, MPI_DOUBLE, 0, 0, joint_comm);
      MPI_Bcast(global_mesh_twin.y_glob.data(), mysize, MPI_DOUBLE, 0, m2c_comm);

      // dy
      MPI_Recv(&mysize, 1, MPI_INT, 0, 0, joint_comm);
      assert(mysize>0);
      global_mesh_twin.dy_glob.resize(mysize);
      MPI_Bcast(&mysize, 1, MPI_INT, 0, m2c_comm);

      MPI_Recv(global_mesh_twin.dy_glob.data(), mysize, MPI_DOUBLE, 0, 0, joint_comm);
      MPI_Bcast(global_mesh_twin.dy_glob.data(), mysize, MPI_DOUBLE, 0, m2c_comm);

      // z
      MPI_Recv(&mysize, 1, MPI_INT, 0, 0, joint_comm);
      assert(mysize>0);
      global_mesh_twin.z_glob.resize(mysize);
      MPI_Bcast(&mysize, 1, MPI_INT, 0, m2c_comm);

      MPI_Recv(global_mesh_twin.z_glob.data(), mysize, MPI_DOUBLE, 0, 0, joint_comm);
      MPI_Bcast(global_mesh_twin.z_glob.data(), mysize, MPI_DOUBLE, 0, m2c_comm);

      // dz
      MPI_Recv(&mysize, 1, MPI_INT, 0, 0, joint_comm);
      assert(mysize>0);
      global_mesh_twin.dz_glob.resize(mysize);
      MPI_Bcast(&mysize, 1, MPI_INT, 0, m2c_comm);

      MPI_Recv(global_mesh_twin.dz_glob.data(), mysize, MPI_DOUBLE, 0, 0, joint_comm);
      MPI_Bcast(global_mesh_twin.dz_glob.data(), mysize, MPI_DOUBLE, 0, m2c_comm);

      // receive the location(s) of the overset boundaries
      MPI_Recv(overset_boundary.data(), overset_boundary.size(), MPI_INT, 0, 0, joint_comm);
      MPI_Bcast(overset_boundary.data(), overset_boundary.size(), MPI_INT, 0, m2c_comm);

    } else {
      int mysize;

      // x
      MPI_Bcast(&mysize, 1, MPI_INT, 0, m2c_comm);
      assert(mysize>0);
      global_mesh_twin.x_glob.resize(mysize);
      MPI_Bcast(global_mesh_twin.x_glob.data(), mysize, MPI_DOUBLE, 0, m2c_comm);
       
      // dx
      MPI_Bcast(&mysize, 1, MPI_INT, 0, m2c_comm);
      assert(mysize>0);
      global_mesh_twin.dx_glob.resize(mysize);
      MPI_Bcast(global_mesh_twin.dx_glob.data(), mysize, MPI_DOUBLE, 0, m2c_comm);
       
      // y
      MPI_Bcast(&mysize, 1, MPI_INT, 0, m2c_comm);
      assert(mysize>0);
      global_mesh_twin.y_glob.resize(mysize);
      MPI_Bcast(global_mesh_twin.y_glob.data(), mysize, MPI_DOUBLE, 0, m2c_comm);
       
      // dy
      MPI_Bcast(&mysize, 1, MPI_INT, 0, m2c_comm);
      assert(mysize>0);
      global_mesh_twin.dy_glob.resize(mysize);
      MPI_Bcast(global_mesh_twin.dy_glob.data(), mysize, MPI_DOUBLE, 0, m2c_comm);
       
      // z
      MPI_Bcast(&mysize, 1, MPI_INT, 0, m2c_comm);
      assert(mysize>0);
      global_mesh_twin.z_glob.resize(mysize);
      MPI_Bcast(global_mesh_twin.z_glob.data(), mysize, MPI_DOUBLE, 0, m2c_comm);
       
      // dz
      MPI_Bcast(&mysize, 1, MPI_INT, 0, m2c_comm);
      assert(mysize>0);
      global_mesh_twin.dz_glob.resize(mysize);
      MPI_Bcast(global_mesh_twin.dz_glob.data(), mysize, MPI_DOUBLE, 0, m2c_comm);

      // receive the location(s) of the overset boundaries
      MPI_Bcast(overset_boundary.data(), overset_boundary.size(), MPI_INT, 0, m2c_comm);
    }

    // Now, global_mesh_twin has been created. overset_boundary has also been created.


    // -----------------------------------------------
    // Step 5. Identify import nodes ("green boxes" in KW notes) 
    //         and send all of them to all leader procs
    // -----------------------------------------------
    vector<Int3> import_all;
    vector<Vec3D> import_all_coords;

    double*** tag    = TMP->GetDataPointer();
    Vec3D***  coords = (Vec3D***)coordinates->GetDataPointer();
    Vec3D***  xx     = (Vec3D***)TMP3->GetDataPointer(); //value is 0 by default
    int i0, j0, k0, imax, jmax, kmax;
    coordinates->GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax); 

    // if the node is (1) inside the domain of leader's mesh AND (2) adjacent to either 
    // a red box or a real node outside leader's domain, then we tag it ("a green box")
    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++) {

          if(!global_mesh_twin.IsPointInDomain(coords[k][j][i], false))
            continue;

          if(tag[k][j][i] != 0) //already has a "status" (red or green box)
            continue;

          bool green = (tag[k][j][i-1]==1 || tag[k][j][i+1]==1 ||
                        tag[k][j-1][i]==1 || tag[k][j+1][i]==1 ||
                        tag[k-1][j][i]==1 || tag[k+1][j][i]==1); //adjacet to red box?
          if(!green) {
            if( (coordinates->IsHere(i-1,j,k) && 
                 !global_mesh_twin.IsPointInDomain(coords[k][j][i-1], false)) ||
                (coordinates->IsHere(i+1,j,k) && 
                 !global_mesh_twin.IsPointInDomain(coords[k][j][i+1], false)) ||
                (coordinates->IsHere(i,j-1,k) && 
                 !global_mesh_twin.IsPointInDomain(coords[k][j-1][i], false)) ||
                (coordinates->IsHere(i,j+1,k) && 
                 !global_mesh_twin.IsPointInDomain(coords[k][j+1][i], false)) ||
                (coordinates->IsHere(i,j,k-1) && 
                 !global_mesh_twin.IsPointInDomain(coords[k-1][j][i], false)) ||
                (coordinates->IsHere(i,j,k+1) && 
                 !global_mesh_twin.IsPointInDomain(coords[k+1][j][i], false)) )
              green = true; //adjacent to a real node outside leader's domain
          } 

          if(green) {
            tag[k][j][i] = 2; //register a green box

            import_all.push_back(Int3(i,j,k));
            import_all_coords.push_back(coords[k][j][i]);

            // block adjacent edges (6 of them) --- to verify that the green boxes 
            // form a closed surface
            xx[k][j][i] = Vec3D(1,1,1); 
            xx[k+1][j][i][2] = xx[k][j+1][i][1] = xx[k][j][i+1][1] = 1;
          }
     
        }

    TMP->RestoreDataPointerAndInsert();
    coordinates->RestoreDataPointerToLocalVector();
    TMP3->RestoreDataPointerAndInsert();

    // verify that the green boxes (tag = 2) form a closed interface in the outer
    // mesh that separates nodes that are "active" and "inactive"
    std::set<Int3> occluded_nodes; //not used
    int nReg = floodfiller->FillBasedOnEdgeObstructions(TMP3, 0/*non_obstruction_flag*/, 
                                                        occluded_nodes, Color);
    if(nReg != 2) {
      print_error("*** Error: Detected error in overset grids (M2C-M2C). Number of regions: %d."
                  " Should be 2.\n", nReg);
      exit(-1);
    }
 
    // Send green boxes to leader procs.
    int nNodes_all = import_all.size();
    send_requests.clear();
    for(int proc=0; proc<numLeaderProcs; proc++) {
      send_requests.push_back(MPI_Request());
      MPI_Isend(&nNodes_all, 1, MPI_INT, proc, m2c_rank/*tag*/, joint_comm,
                &send_requests.back());
    }    
    MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);
    send_requests.clear();

    if(nNodes_all>0) {
      for(int proc=0; proc<numLeaderProcs; proc++) {
        send_requests.push_back(MPI_Request());
        MPI_Isend((double*)import_all_coords.data(), 3*nNodes_all, MPI_DOUBLE,
                  proc, m2c_rank/*tag*/, joint_comm, &send_requests.back());
      }    
      MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);
      send_requests.clear();
    }

    // -----------------------------------------------
    // Step 6. Receive adoption info from leaders
    // -----------------------------------------------
    recv_requests.clear();
    nNodes.assign(numLeaderProcs,0);
    if(nNodes_all>0) {
      for(int proc=0; proc<numLeaderProcs; proc++) {
        recv_requests.push_back(MPI_Request());
        MPI_Irecv(&nNodes[proc], 1, MPI_INT, proc, proc, joint_comm, &recv_requests.back());
      }
    }
    MPI_Waitall(recv_requests.size(), recv_requests.data(), MPI_STATUSES_IGNORE);
    recv_requests.clear();

    for(int proc=0; proc<numLeaderProcs; proc++) {
      found[proc].clear();
      if(nNodes[proc]>0) {
        found[proc].resize(nNodes[proc]);
        recv_requests.push_back(MPI_Request()); 
        MPI_Irecv(found[proc].data(), found[proc].size(), proc, proc, joint_comm,
                  &recv_requests.back());
      }
    }
    MPI_Waitall(recv_requests.size(), recv_requests.data(), MPI_STATUSES_IGNORE);
    recv_requests.clear();


    // -----------------------------------------------
    // Step 7. Construct import_nodes; avoid duplication
    // -----------------------------------------------
    import_nodes.resize(numLeaderProcs);

    std::multimap<int,int> duplicates;
    vector<int> owner(import_all.size(),-1);
    for(int proc=0; proc<numLeaderProcs; proc++) {
      for(auto&& id : found[proc]) {
        if(owner[id]<0) {
          owner[id] = proc;
          import_nodes[proc].push_back(import_all[id]);
        } else //duplicates
          duplicates.insert(make_pair(proc, id));
      }
    }

    // make sure there are no orphans
    for(int i=0; i<owner.size(); i++) {
      if(owner[i] == -1) {
        fprintf(stderr,"\033[0;31m*** Error: [M2C-M2C] Node (%d,%d,%d) (%e,%e,%e) is not picked"
                       " up by any leader processor.\033[0m\n", import_all[i][0],
                       import_all[i][1], import_all[i][2], import_all_coords[i][0],
                       import_all_coords[i][1], import_all_coords[i][2]);
        exit(-1);
      }
    }

    // notify leaders about duplication, which should be erased
    send_requests.clear();
    for(int proc=0; proc<numLeaderProcs; proc++)
      dups[proc].clear();

    if(nNodes_all>0) {

      for(auto&& nod : duplicates)
        dups[nod.first].push_back(nod.second);

      for(int proc=0; proc<numLeaderProcs; proc++) {
        send_requests.push_back(MPI_Request());
        int ndup = dups[proc].size();
        MPI_Isend(&ndup, 1, MPI_INT, proc, m2c_rank, joint_comm, &send_requests.back()); 
      }
      MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);
      send_requests.clear();

      for(int proc=0; proc<numLeaderProcs; proc++) {
        if(dups[proc].size()>0) {
          send_requests.push_back(MPI_Request());
          MPI_Isend(dups[proc].data(), dups[proc].size(), proc, m2c_rank, joint_comm,
                    &send_requests.back());
        }
      }
      MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);
      send_requests.clear();
    }

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

