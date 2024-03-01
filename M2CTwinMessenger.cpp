/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<M2CTwinMessenger.h>
#include<trilinear_interpolation.h>
#include<cassert>
#include<climits>
#include<algorithm> //std::find


using std::vector;
using std::set;

extern int verbose;
extern int INACTIVE_MATERIAL_ID;

//---------------------------------------------------------------

M2CTwinMessenger::M2CTwinMessenger(IoData &iod_, MPI_Comm &m2c_comm_, MPI_Comm &joint_comm_,
                                   int status_)
              : iod(iod_), m2c_comm(m2c_comm_), joint_comm(joint_comm_),
                coordinates(NULL), ghost_nodes_inner(NULL), ghost_nodes_outer(NULL),
                global_mesh(NULL), TMP(NULL), TMP3(NULL), Color(NULL), floodfiller(NULL)
{

  MPI_Comm_rank(m2c_comm, &m2c_rank);
  MPI_Comm_size(m2c_comm, &m2c_size);

  assert(status_==1 || status_==2);
  if(status_==1)
    twinning_status = LEADER;
  else
    twinning_status = FOLLOWER;

  dt = -1.0;
  tmax = -1.0;

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
M2CTwinMessenger::CommunicateBeforeTimeStepping(SpaceVariable3D &coordinates_, DataManagers3D &dms_,
                                                vector<GhostPoint> &ghost_nodes_inner_,
                                                vector<GhostPoint> &ghost_nodes_outer_,
                                                GlobalMeshInfo &global_mesh_, SpaceVariable3D &V,
                                                SpaceVariable3D &ID, set<Int3> &spo_frozen_nodes)
{

  coordinates       = &coordinates_;
  ghost_nodes_inner = &ghost_nodes_inner_;
  ghost_nodes_outer = &ghost_nodes_outer_;
  global_mesh       = &global_mesh_;

  if(twinning_status == LEADER) {


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
    for(auto&& ghost : *ghost_nodes_outer)
      if(ghost.bcType == (int)MeshData::OVERSET && 
         ghost.type_projection == GhostPoint::FACE) {
        if(!coordinates->IsHereOrConnectedGhost(ghost.ijk[0], ghost.ijk[1], ghost.ijk[2]))
          continue; //Skip ghost nodes at corners/edges of the subdomain --- they are not needed
        import_all.push_back(ghost.ijk);
        import_all_coords.push_back(global_mesh->GetXYZ(ghost.ijk));
      }

    int nNodes_all = import_all.size();

    vector<MPI_Request> send_requests;
    for(int proc = 0; proc < numFollowerProcs; proc++) {
      send_requests.push_back(MPI_Request());
      //fprintf(stdout,"Leader proc %d sends nNodes_all = %d to Follower %d.\n", m2c_rank, nNodes_all, proc);
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
        MPI_Irecv(found[proc].data(), found[proc].size(), MPI_INT, proc, proc, joint_comm,
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
        assert(id>=0 && id<(int)import_all.size());
        if(owner[id]<0) {
          owner[id] = proc;
          import_nodes[proc].push_back(import_all[id]);  
        } else //duplicates
          duplicates.insert(std::make_pair(proc, id));
      }
    }

    // make sure there are no orphans
    for(int i=0; i<(int)owner.size(); i++) {
      if(owner[i] == -1) {
        fprintf(stdout,"\033[0;31m*** Error: [M2C-M2C] Node (%d,%d,%d) (%e,%e,%e) is not picked"
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
          MPI_Isend(dups[proc].data(), dups[proc].size(), MPI_INT, proc, m2c_rank, joint_comm,
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
      for(int i=0; i<(int)export_points_all[proc].size(); i++) {
        if(global_mesh->FindElementCoveringPoint(export_points_all[proc][i], ijk0, &xi, true)) {
          if(coordinates->IsHere(ijk0[0], ijk0[1], ijk0[2], true/*include_ghost*/) &&
             coordinates->IsHere(ijk0[0]+1,ijk0[1]+1,ijk0[2]+1, true/*include_ghost*/)) {
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

    for(int proc=0; proc<numFollowerProcs; proc++) {
      if(found[proc].size()>0) {
        send_requests.push_back(MPI_Request());
        MPI_Isend(found[proc].data(), found[proc].size(), MPI_INT, proc, m2c_rank, joint_comm,
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
      vector<int> local_dups;
      for(auto&& dup : dups[proc]) {
        auto it = std::find(found[proc].begin(), found[proc].end(), dup);
        assert(it != found[proc].end());
        local_dups.push_back(it - found[proc].begin()); //cannot directly erase!
      }
      std::sort(local_dups.begin(), local_dups.end());
      for(int n = local_dups.size()-1; n>=0; n--) //going backward to avoid wrong index 
        export_points[proc].erase(export_points[proc].begin() + local_dups[n]);
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
      for(int i=0; i<(int)export_points_all[proc].size(); i++) { 
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
      //fprintf(stdout,"Follower %d adopts %ld black points from Leader %d.\n", m2c_rank, found[proc].size(), proc); 
      if(found[proc].size()>0) {
        send_requests.push_back(MPI_Request());
        MPI_Isend(found[proc].data(), found[proc].size(), MPI_INT, proc, m2c_rank/*tag*/, joint_comm,
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
        dups[proc].resize(ndup[proc]);

        recv_requests.push_back(MPI_Request());
        MPI_Irecv(dups[proc].data(), dups[proc].size(), MPI_INT, proc, proc,
                  joint_comm, &recv_requests.back());
      }
    }
    MPI_Waitall(recv_requests.size(), recv_requests.data(), MPI_STATUSES_IGNORE);
    recv_requests.clear();
    
    //find duplicates and remove them
    for(int proc = 0; proc < numLeaderProcs; proc++) {
      vector<int> local_dups;
      for(auto&& dup : dups[proc]) {
        auto it = std::find(found[proc].begin(), found[proc].end(), dup);
        assert(it != found[proc].end());
        local_dups.push_back(it - found[proc].begin()); //can NOT directly erase (ordering!)
      }
      std::sort(local_dups.begin(), local_dups.end());
      for(int n = local_dups.size()-1; n>=0; n--) //going backward to avoid wrong index 
        export_points[proc].erase(export_points[proc].begin() + local_dups[n]);
    }
    
    //Tag nodes involved in interpolations ("red boxes" in KW's notes)
    double*** tag = TMP->GetDataPointer(); //value is 0 by default
    for(int proc = 0; proc < numLeaderProcs; proc++)
      for(auto&& gp : export_points[proc]) 
        for(int k=gp.ijk[2]; k<=gp.ijk[2]+1; k++)
          for(int j=gp.ijk[1]; j<=gp.ijk[1]+1; j++)
            for(int i=gp.ijk[0]; i<=gp.ijk[0]+1; i++)
              tag[k][j][i] = 1;
    TMP->RestoreDataPointerAndAdd(); //Note: "Add" allows us to keep track of ALL involved points,
                                     //      "Insert" would lead to some points being lost in the comm.
    tag = TMP->GetDataPointer();
    int ii0,jj0,kk0,iimax,jjmax,kkmax;
    coordinates->GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);
    for(int k=kk0; k<kkmax; k++)
      for(int j=jj0; j<jjmax; j++)
        for(int i=ii0; i<iimax; i++)
          if(tag[k][j][i]>0)
            tag[k][j][i] = 1; //tag should be 0 or 1
    TMP->RestoreDataPointerAndInsert();


    // -----------------------------------------------
    // Step 4. Get info about leader's mesh
    // -----------------------------------------------
    vector<int> overset_boundary(6, 0); //xmin, xmax, ymin, ymax, zmin, zmax, 1 means yes
  
    if(m2c_rank==0) { //receives data from leader proc #0, then broadcast
      int mysize;

      // x
      MPI_Recv(&mysize, 1, MPI_INT, 0, 0, joint_comm, MPI_STATUS_IGNORE);
      assert(mysize>0);
      global_mesh_twin.x_glob.resize(mysize);
      MPI_Bcast(&mysize, 1, MPI_INT, 0, m2c_comm);

      MPI_Recv(global_mesh_twin.x_glob.data(), mysize, MPI_DOUBLE, 0, 0, joint_comm, MPI_STATUS_IGNORE);
      MPI_Bcast(global_mesh_twin.x_glob.data(), mysize, MPI_DOUBLE, 0, m2c_comm);

      // dx
      MPI_Recv(&mysize, 1, MPI_INT, 0, 0, joint_comm, MPI_STATUS_IGNORE);
      assert(mysize>0);
      global_mesh_twin.dx_glob.resize(mysize);
      MPI_Bcast(&mysize, 1, MPI_INT, 0, m2c_comm);

      MPI_Recv(global_mesh_twin.dx_glob.data(), mysize, MPI_DOUBLE, 0, 0, joint_comm, MPI_STATUS_IGNORE);
      MPI_Bcast(global_mesh_twin.dx_glob.data(), mysize, MPI_DOUBLE, 0, m2c_comm);

      // y
      MPI_Recv(&mysize, 1, MPI_INT, 0, 0, joint_comm, MPI_STATUS_IGNORE);
      assert(mysize>0);
      global_mesh_twin.y_glob.resize(mysize);
      MPI_Bcast(&mysize, 1, MPI_INT, 0, m2c_comm);

      MPI_Recv(global_mesh_twin.y_glob.data(), mysize, MPI_DOUBLE, 0, 0, joint_comm, MPI_STATUS_IGNORE);
      MPI_Bcast(global_mesh_twin.y_glob.data(), mysize, MPI_DOUBLE, 0, m2c_comm);

      // dy
      MPI_Recv(&mysize, 1, MPI_INT, 0, 0, joint_comm, MPI_STATUS_IGNORE);
      assert(mysize>0);
      global_mesh_twin.dy_glob.resize(mysize);
      MPI_Bcast(&mysize, 1, MPI_INT, 0, m2c_comm);

      MPI_Recv(global_mesh_twin.dy_glob.data(), mysize, MPI_DOUBLE, 0, 0, joint_comm, MPI_STATUS_IGNORE);
      MPI_Bcast(global_mesh_twin.dy_glob.data(), mysize, MPI_DOUBLE, 0, m2c_comm);

      // z
      MPI_Recv(&mysize, 1, MPI_INT, 0, 0, joint_comm, MPI_STATUS_IGNORE);
      assert(mysize>0);
      global_mesh_twin.z_glob.resize(mysize);
      MPI_Bcast(&mysize, 1, MPI_INT, 0, m2c_comm);

      MPI_Recv(global_mesh_twin.z_glob.data(), mysize, MPI_DOUBLE, 0, 0, joint_comm, MPI_STATUS_IGNORE);
      MPI_Bcast(global_mesh_twin.z_glob.data(), mysize, MPI_DOUBLE, 0, m2c_comm);

      // dz
      MPI_Recv(&mysize, 1, MPI_INT, 0, 0, joint_comm, MPI_STATUS_IGNORE);
      assert(mysize>0);
      global_mesh_twin.dz_glob.resize(mysize);
      MPI_Bcast(&mysize, 1, MPI_INT, 0, m2c_comm);

      MPI_Recv(global_mesh_twin.dz_glob.data(), mysize, MPI_DOUBLE, 0, 0, joint_comm, MPI_STATUS_IGNORE);
      MPI_Bcast(global_mesh_twin.dz_glob.data(), mysize, MPI_DOUBLE, 0, m2c_comm);

      // receive the location(s) of the overset boundaries
      MPI_Recv(overset_boundary.data(), overset_boundary.size(), MPI_INT, 0, 0, joint_comm, MPI_STATUS_IGNORE);
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

    tag    = TMP->GetDataPointer();
    Vec3D***  coords = (Vec3D***)coordinates->GetDataPointer();
    Vec3D***  xx     = (Vec3D***)TMP3->GetDataPointer(); //value is 0 by default
    double*** id     = ID.GetDataPointer(); //set inactive id

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
 
            //add to "frozen nodes"
            spo_frozen_nodes.insert(Int3(i,j,k));

            // block adjacent edges (6 of them) --- to verify that the green boxes 
            // form a closed surface
            xx[k][j][i] = Vec3D(1,1,1); 
            xx[k+1][j][i][2] = xx[k][j+1][i][1] = xx[k][j][i+1][1] = 1;

          } else
            id[k][j][i] = INACTIVE_MATERIAL_ID;
     
        }

    TMP->RestoreDataPointerAndInsert();
    coordinates->RestoreDataPointerToLocalVector();
    TMP3->RestoreDataPointerAndInsert();
    ID.RestoreDataPointerAndInsert();

    // verify that the green boxes (tag = 2) form a closed interface in the outer
    // mesh that separates nodes that are "active" and "inactive"
    std::set<Int3> occluded_nodes; //to be filled w/ green boxes, including internal ghosts
    tag = TMP->GetDataPointer();
    for(int k=kk0; k<kkmax; k++)
      for(int j=jj0; j<jjmax; j++)
        for(int i=ii0; i<iimax; i++)
          if(tag[k][j][i]==2)
            occluded_nodes.insert(Int3(i,j,k)); // green box
    TMP->RestoreDataPointerToLocalVector();
    
    int nReg = floodfiller->FillBasedOnEdgeObstructions(*TMP3, 0/*non_obstruction_flag*/, 
                                                        occluded_nodes, *Color);
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
        MPI_Irecv(found[proc].data(), found[proc].size(), MPI_INT, proc, proc, joint_comm,
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
          duplicates.insert(std::make_pair(proc, id));
      }
    }

    // make sure there are no orphans
    for(int i=0; i<(int)owner.size(); i++) {
      if(owner[i] == -1) {
        fprintf(stdout,"\033[0;31m*** Error: [M2C-M2C] Node (%d,%d,%d) (%e,%e,%e) is not picked"
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
          MPI_Isend(dups[proc].data(), dups[proc].size(), MPI_INT, proc, m2c_rank, joint_comm,
                    &send_requests.back());
        }
      }
      MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);
      send_requests.clear();
    }

  }

  ExchangeData(V);

  if(twinning_status == FOLLOWER) {
    // in terms of time-stepping, the follower is one-step behind the leader --> O(dt) error
    dt   = 1.0e-40; //essentially, 0.
    tmax = DBL_MAX;
  }
}

//---------------------------------------------------------------

void
M2CTwinMessenger::FirstExchange(SpaceVariable3D &V, double dt0, double tmax0)
{
  Exchange(V,dt0,tmax0); 
}

//---------------------------------------------------------------

void
M2CTwinMessenger::Exchange(SpaceVariable3D &V, double dt0, double tmax0)
{
  ExchangeData(V);

  // leader sends dt and tmax to follower
  if(twinning_status == LEADER) {
    if(m2c_rank==0) {
      double buf[2] = {dt0, tmax0};
      MPI_Send(buf, 2, MPI_DOUBLE, 0, 0, joint_comm);
    }
  } else {
    double buf[2];
    if(m2c_rank==0)
      MPI_Recv(buf, 2, MPI_DOUBLE, 0, 0, joint_comm, MPI_STATUS_IGNORE);
    MPI_Bcast(buf, 2, MPI_DOUBLE, 0, m2c_comm);

    dt   = buf[0];
    tmax = buf[1];
  } 
 
}

//---------------------------------------------------------------

void
M2CTwinMessenger::FinalExchange(SpaceVariable3D &V)
{
  if(twinning_status == LEADER) {
    //when LEADER is here, the follower is still executing "Exchange(...)"
    
    ExchangeData(V);

    if(m2c_rank==0) {
      double dt0   = 1.0e-40;
      double tmax0 = 1.0e-40; //send a very small tmax to follower to trigger its termination
      double buf[2] = {dt0, tmax0};
      MPI_Send(buf, 2, MPI_DOUBLE, 0, 0, joint_comm);
    }
  }

  //FOLLOWER: When follower reaches here, the leader has finished.
}

//---------------------------------------------------------------

void
M2CTwinMessenger::ExchangeData(SpaceVariable3D &V)
{
  //First send data from follower to leader, then the opposite way (See KW's notes)
  InterpolateDataAndTransfer(V, FOLLOWER);
  InterpolateDataAndTransfer(V, LEADER);
}

//---------------------------------------------------------------

void
M2CTwinMessenger::InterpolateDataAndTransfer(SpaceVariable3D &V, TwinningStatus sender)
{

  int numTwinProcs(-1);
  MPI_Comm_remote_size(joint_comm, &numTwinProcs);

  // ----------------------------------------
  // Step 1. Check buffer size. Extend if needed.
  // ----------------------------------------
  int dim(-1);
  if((int)import_buffer.size() != numTwinProcs)
    import_buffer.resize(numTwinProcs);
  for(int proc=0; proc<numTwinProcs; proc++) {
    if(import_nodes[proc].size()>0) {
      dim = import_buffer[proc].size()/import_nodes[proc].size();
      break;
    }
  }
  if(dim<V.NumDOF()) {
    dim = V.NumDOF();
    for(int proc=0; proc<numTwinProcs; proc++)
      import_buffer[proc].resize(dim*import_nodes[proc].size());
  }

  dim = -1;
  if((int)export_buffer.size() != numTwinProcs)
    export_buffer.resize(numTwinProcs);
  for(int proc=0; proc<numTwinProcs; proc++) {
    if(export_points[proc].size()>0) {
      dim = export_buffer[proc].size()/export_points[proc].size();
      break;
    }
  }
  if(dim<V.NumDOF()) {
    dim = V.NumDOF();
    for(int proc=0; proc<numTwinProcs; proc++)
      export_buffer[proc].resize(dim*export_points[proc].size());
  }

  // ----------------------------------------
  // Step 2. Interpolate data and transfer
  // ----------------------------------------
  dim = V.NumDOF();

  if(twinning_status == sender) { 

    vector<MPI_Request> send_requests;
    double*** v = V.GetDataPointer();
    for(int proc=0; proc<numTwinProcs; proc++) {

      if(export_points[proc].size()==0)
        continue;

      for(int p=0; p<(int)export_points[proc].size(); p++) {
        int i = export_points[proc][p].ijk[0]; 
        int j = export_points[proc][p].ijk[1]; 
        int k = export_points[proc][p].ijk[2]; 
        for(int d=0; d<dim; d++)     
          export_buffer[proc][dim*p+d] 
              = MathTools::trilinear_interpolation(v[k][j][i*dim+d], v[k][j][(i+1)*dim+d], 
                    v[k][j+1][i*dim+d], v[k][j+1][(i+1)*dim+d], v[k+1][j][i*dim+d], 
                    v[k+1][j][(i+1)*dim+d], v[k+1][j+1][i*dim+d], v[k+1][j+1][(i+1)*dim+d],
                    (double*)export_points[proc][p].xi);
      }
      
      send_requests.push_back(MPI_Request());
      MPI_Isend(export_buffer[proc].data(), dim*export_points[proc].size(), MPI_DOUBLE, proc, 
                m2c_rank/*tag*/, joint_comm, &send_requests.back());
    }
    MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);
    send_requests.clear();

    V.RestoreDataPointerToLocalVector();
  
  } else { //receiver

    vector<MPI_Request> recv_requests;
    for(int proc=0; proc<numTwinProcs; proc++) {
      if(import_nodes[proc].size()>0) {
        recv_requests.push_back(MPI_Request());
        MPI_Irecv(import_buffer[proc].data(), dim*import_nodes[proc].size(), MPI_DOUBLE, 
                  proc, proc, joint_comm, &recv_requests.back()); 
      }
    }
    MPI_Waitall(recv_requests.size(), recv_requests.data(), MPI_STATUSES_IGNORE);
    recv_requests.clear();

    double*** v = V.GetDataPointer();
    for(int proc=0; proc<numTwinProcs; proc++) {
      for(int p=0; p<(int)import_nodes[proc].size(); p++) {
        Int3 &ijk(import_nodes[proc][p]);
        for(int d=0; d<dim; d++) 
          v[ijk[2]][ijk[1]][ijk[0]*dim+d] = import_buffer[proc][p*dim+d];
      }
    }
    V.RestoreDataPointerAndInsert();

  }
}

//---------------------------------------------------------------


//---------------------------------------------------------------

