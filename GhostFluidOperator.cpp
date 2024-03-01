/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<GhostFluidOperator.h>
#include<EmbeddedBoundaryFormula.h>
#include<EmbeddedBoundaryDataSet.h>
#include<Vector5D.h>
#include<Utils.h>
#include<cfloat> //DBL_MAX
#include<map>
#include<memory>

using std::vector;
using std::map;

extern int verbose;
extern int INACTIVE_MATERIAL_ID;

//-------------------------------------------------------------

GhostFluidOperator::GhostFluidOperator(MPI_Comm &comm_, DataManagers3D &dm_all_,
                                       GlobalMeshInfo &global_mesh_)
                  : comm(comm_), global_mesh(global_mesh_),
                    neicomm(comm, global_mesh_),
                    Tag(comm_, &(dm_all_.ghosted1_1dof))
{

  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);
  assert(mpi_rank>=0 && mpi_rank<mpi_size);

  Tag.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  Tag.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);
  Tag.GetGlobalSize(&NX,&NY,&NZ);

}

//-------------------------------------------------------------

GhostFluidOperator::~GhostFluidOperator()
{
}

//-------------------------------------------------------------

void
GhostFluidOperator::Destroy()
{
  Tag.Destroy(); 
}

//-------------------------------------------------------------

int
GhostFluidOperator::PopulateGhostNodesForViscosityOperator(SpaceVariable3D &V, 
                        SpaceVariable3D &ID,
                        vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS,
                        SpaceVariable3D &Velog)
{
  assert(EBDS && EBDS->size()>0);

  double*** id = ID.GetDataPointer();

  // ----------------------------------------------------------------------------
  // Step 1: Identify ghost nodes inside the subdomain that need to be populated
  // ----------------------------------------------------------------------------
  int number_of_ghosts = TagGhostNodes(id, 2); //2 layers; Tag and ghosts_ijk get computed 

  if(number_of_ghosts==0) {
    ID.RestoreDataPointerToLocalVector();
    return 0;
  }


  // ----------------------------------------------------------------------------
  // Step 2: Find image points and the needed geometric information
  // ----------------------------------------------------------------------------
  vector<Vec3D> ghosts_xyz;
  vector<Int3>  images_ijk;
  vector<Vec3D> images_xi;
  vector<Vec3D> images_xyz;
  vector<Vec3D> ghosts_xp; //position of projection points on the interface
  vector<Vec3D> ghosts_vp; //velocity at projection points on the interface

  vector<int>   ghosts_tag; //-1: undetermined, 0: const extrap should be performed
                            // 1: linear extrap should be performed, within this subdomain
                            // 2: linear extrap should be performed; requires neighbor info

  ProjectGhostsToInterface(EBDS, ghosts_xyz, ghosts_xp, ghosts_vp, 
                           images_ijk, images_xi, images_xyz, ghosts_tag);


  // ----------------------------------------------------------------------------
  // Step 3: Find a usable image for each ghost node, starting w/ the one found above
  //         Determine "ghosts_tag"
  // ----------------------------------------------------------------------------
  FindImagesForGhosts(id, ghosts_xyz, images_ijk, images_xi, images_xyz, ghosts_tag);


  // ----------------------------------------------------------------------------
  // Step 4: Populate ghost nodes
  // ----------------------------------------------------------------------------
  Vec5D*** v   = (Vec5D***)V.GetDataPointer();
  Vec3D*** vel = (Vec3D***)Velog.GetDataPointer();
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++) 
        for(int p=0; p<3; p++)
          vel[k][j][i][p] = v[k][j][i][1+p];
  V.RestoreDataPointerToLocalVector(); //V remains unchanged
  

  vector<Int3>           neighbor_requests;
  map<int, vector<int> > outside_images;
  vector<double>         neighbor_received; //will be fixed by neicomm

  // 4.1: Get necessary info from neighbors
  for(int n=0; n<(int)ghosts_ijk.size(); n++) {
    if(ghosts_tag[n]==2) {
      int i,j,k;
      for(int dk=0; dk<=1; dk++) {
        k = images_ijk[n][2]+dk;
        for(int dj=0; dj<=1; dj++) {
          j = images_ijk[n][1]+dj;
          for(int di=0; di<=1; di++) {
            i = images_ijk[n][0]+di;

            if(!Tag.IsHere(i,j,k,true)) { //get from neighbor
              neighbor_requests.push_back(Int3(i,j,k));
              if(outside_images.find(n)==outside_images.end())
                outside_images[n] = vector<int>();
              outside_images[n].push_back(neighbor_requests.size()-1);
            }
          }
        }
      }
    }
  }

  neicomm.Request((double***)vel, 3, neighbor_requests, neighbor_received, 0);
  assert(neighbor_received.size() == 3*neighbor_requests.size());

  // 4.2: Populate ghost nodes
  EmbeddedBoundaryFormula formula(EmbeddedBoundaryFormula::LINEAR_EXTRAPOLATION);

  for(int n=0; n<(int)ghosts_ijk.size(); n++) {

    if(ghosts_tag[n] == 0) {//constant extrapolation
      vel[ghosts_ijk[n][2]][ghosts_ijk[n][1]][ghosts_ijk[n][0]] = ghosts_vp[n];
    }
    else {//linear extrapolation

      double alpha = (ghosts_xyz[n]-ghosts_xp[n]).norm()
                   / (ghosts_xyz[n]-images_xyz[n]).norm();
      formula.BuildLinearExtrapolationFormula(ghosts_ijk[n], images_ijk[n], 
                                              images_xi[n], alpha); 

      if(ghosts_tag[n] == 1) {//within subdomain
        vel[ghosts_ijk[n][2]][ghosts_ijk[n][1]][ghosts_ijk[n][0]] 
            = formula.Evaluate3D(vel, ghosts_vp[n]);
      }
      else { //needs info outside subdomain

        assert(ghosts_tag[n]==2);

        vector<Int3>& nodes(formula.GetNodes());
        vector<Vec3D> vloc(nodes.size());

        for(int p=0; p<(int)nodes.size(); p++) {
          if(Tag.IsHere(nodes[p][0],nodes[p][1],nodes[p][2],true)) {
            vloc[p] = vel[nodes[p][2]][nodes[p][1]][nodes[p][0]];
          }
          else {
            auto it = outside_images.find(n);
            assert(it != outside_images.end());
            bool foundit = false;
            for(auto&& g : it->second) {
              if(neighbor_requests[g]==nodes[p]) {//found you!
                for(int dof=0; dof<3; dof++)
                  vloc[p][dof] = neighbor_received[3*g+dof];
                foundit = true;
                break;
              }
            } 
            assert(foundit);
          }
        }

        vel[ghosts_ijk[n][2]][ghosts_ijk[n][1]][ghosts_ijk[n][0]]
            = formula.Evaluate3D(vloc, ghosts_vp[n]); 
      }

    }
  }

  Velog.RestoreDataPointerAndInsert();

  ID.RestoreDataPointerToLocalVector();

  return number_of_ghosts;
}

//-------------------------------------------------------------

int
GhostFluidOperator::TagGhostNodes(double*** id, int nLayers)
{
  assert(nLayers>=1);

  ghosts_ijk.clear();

  double*** tag = Tag.GetDataPointer();

  // round 1
  auto current_layer = std::make_unique<std::set<Int3> >();
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        tag[k][j][i] = 0; //default
        if(id[k][j][i]==INACTIVE_MATERIAL_ID) {
          // check the neighbors
          for(int kk=k-1; kk<=k+1; kk++)
            for(int jj=j-1; jj<=j+1; jj++)
              for(int ii=i-1; ii<=i+1; ii++) {
                if(id[kk][jj][ii] != INACTIVE_MATERIAL_ID) {
                  ghosts_ijk.push_back(Int3(i,j,k)); //first-layer ghost
                  tag[k][j][i] = 1;
                  current_layer->insert(Int3(i,j,k));
                  goto found_layer_1_node;
                }
              }
        }
found_layer_1_node: { } //needs { } to avoid compilation error
      }

  Tag.RestoreDataPointerAndInsert();

  int total_ghosts = ghosts_ijk.size();
  MPI_Allreduce(MPI_IN_PLACE, &total_ghosts, 1, MPI_INT, MPI_SUM, comm);
  if(total_ghosts==0) {
    return total_ghosts;
  }

  // round 2, 3, ... 
  for(int layer=1; layer<nLayers; layer++) {

    auto previous_layer = std::move(current_layer);
    current_layer.reset(new std::set<Int3>());
    
    tag = Tag.GetDataPointer();
  
    int i,j,k;
    for(auto&& g : *previous_layer) {
      i = g[0];
      j = g[1];
      k = g[2];
      for(int kk=k-1; kk<=k+1; kk++)
        for(int jj=j-1; jj<=j+1; jj++)
          for(int ii=i-1; ii<=i+1; ii++) {
            if(!Tag.IsHere(ii,jj,kk,false))
              continue; //does not beyong to this subdomain
            if(id[kk][jj][ii]==INACTIVE_MATERIAL_ID && tag[kk][jj][ii]==0) {
              ghosts_ijk.push_back(Int3(ii,jj,kk)); //next-layer ghost
              tag[kk][jj][ii] = layer+1;
              current_layer->insert(Int3(ii,jj,kk));
            }
          }
    }

    Tag.RestoreDataPointerAndInsert();

    int this_layer_size = current_layer->size();
    MPI_Allreduce(MPI_IN_PLACE, &this_layer_size, 1, MPI_INT, MPI_SUM, comm);
    total_ghosts += this_layer_size;

    if(this_layer_size==0) {//no need to go to the next layer
      return total_ghosts;
    }
  }

  return total_ghosts;
}

//-------------------------------------------------------------

void
GhostFluidOperator::ProjectGhostsToInterface(vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS, 
                                             vector<Vec3D> &ghosts_xyz, vector<Vec3D> &ghosts_xp,
                                             vector<Vec3D> &ghosts_vp, vector<Int3> &images_ijk,
                                             vector<Vec3D> &images_xi, vector<Vec3D> &images_xyz,
                                             vector<int> &ghosts_tag)
{
  vector<double***> cpi;
  vector<vector<std::pair<Int3, ClosestPoint> > * > closest_points;

  for(auto&& data : *EBDS) {
    assert(data->ClosestPointIndex_ptr);
    cpi.push_back(data->ClosestPointIndex_ptr->GetDataPointer());
    assert(data->closest_points_ptr);
    closest_points.push_back(data->closest_points_ptr);
  }

  int nGhosts = ghosts_ijk.size();

  ghosts_xyz.resize(nGhosts);
  ghosts_xp.resize(nGhosts);
  ghosts_vp.resize(nGhosts);
  images_ijk.resize(nGhosts);
  images_xi.resize(nGhosts);
  images_xyz.resize(nGhosts);
  ghosts_tag.assign(nGhosts, -1);  //-1: undetermined, 0: const extrap should be performed
                                   // 1: linear extrap should be performed, within this subdomain
                                   // 2: linear extrap should be performed; requires neighbor info

  Int3 g;
  for(int n=0; n<nGhosts; n++) {

    g = ghosts_ijk[n];
    double min_dist = DBL_MAX;
    int mysurf = -1;
    for(int surf=0; surf<(int)cpi.size(); surf++) {
      int pid = cpi[surf][g[2]][g[1]][g[0]];
      if(pid<0)
        continue;
      ClosestPoint& cp((*closest_points[surf])[pid].second);
      if(cp.dist<min_dist) {
        min_dist = cp.dist; //this is unsigned distance
        mysurf = surf;
      }
    }

    if(mysurf<0) {
      fprintf(stdout,"\033[0;31m*** Error: Projection of node (%d,%d,%d) on embedded "
                     "surface not calculated.\033[0m\n", g[0], g[1], g[2]);
      exit(-1);
    }

    int pid = cpi[mysurf][g[2]][g[1]][g[0]];
    ClosestPoint& cp((*closest_points[mysurf])[pid].second);

    // all the nodes of the embedded surface closest to this ghost node
    vector<Vec3D> &X((*EBDS)[mysurf]->surface_ptr->X);
    vector<Vec3D> &Udot((*EBDS)[mysurf]->surface_ptr->Udot);

    Int3& tnodes((*EBDS)[mysurf]->surface_ptr->elems[cp.tid]);

    // get position at projeciton point
    ghosts_xp[n] = cp.xi[0]*X[tnodes[0]] + cp.xi[1]*X[tnodes[1]] 
                 + cp.xi[2]*X[tnodes[2]]; 

    // get velocity at projeciton point
    ghosts_vp[n] = cp.xi[0]*Udot[tnodes[0]] + cp.xi[1]*Udot[tnodes[1]]
                 + cp.xi[2]*Udot[tnodes[2]];

    ghosts_xyz[n] = global_mesh.GetXYZ(g);
    images_xyz[n] = 2.0*ghosts_xp[n] - ghosts_xyz[n]; //image point

    // get image point location w.r.t. mesh
    bool found = global_mesh.FindElementCoveringPoint(images_xyz[n], images_ijk[n], &images_xi[n],
                                                      true); //include ghost layer outside domain
    if(!found) {
      fprintf(stdout,"\033[0;35mWarning: Detected ghost node (%d,%d,%d) whose image is outside "
                     "physical domain.\033[0m\n", g[0], g[1], g[2]);
      ghosts_tag[n] = 0;
    }

  }

  for(auto&& data : *EBDS)
    data->ClosestPointIndex_ptr->RestoreDataPointerToLocalVector();

}

//-------------------------------------------------------------

void
GhostFluidOperator::FindImagesForGhosts(double*** id, vector<Vec3D> &ghosts_xyz,
                                        vector<Int3> &images_ijk, vector<Vec3D> &images_xi,
                                        vector<Vec3D> &images_xyz, vector<int> &ghosts_tag)
{

  vector<Int3>           neighbor_requests;
  map<int, vector<int> > outside_images;
  vector<double>         neighbor_received; //will be fixed by neicomm

  // Put the undertermined ghosts in a set. Clear this set gradually.
  std::set<int> undetermined;
  for(int n=0; n<(int)ghosts_ijk.size(); n++) {
    // if dist from ghost to interface is <15% of element size, apply const extrap
    if((images_xyz[n]-ghosts_xyz[n]).norm() < 0.3*global_mesh.GetMinDXYZ(ghosts_ijk[n])) 
      ghosts_tag[n] = 0;
    else if(ghosts_tag[n]<0)
      undetermined.insert(n);
  }
 
  // Run five iterations
  for(int iter=0; iter<5; iter++) {

    neighbor_requests.clear();
    outside_images.clear();
    neighbor_received.clear();
    
    std::set<int> found; //valid ghosts, to be removed from undetermined at the end

    //3.1 Figure out the need of information ("id") from neighbors
    for(auto&& n : undetermined) {
      int i,j,k;
      for(int dk=0; dk<=1; dk++) {
        k = images_ijk[n][2]+dk;
        for(int dj=0; dj<=1; dj++) {
          j = images_ijk[n][1]+dj;
          for(int di=0; di<=1; di++) {
            i = images_ijk[n][0]+di;

            if(!Tag.IsHere(i,j,k,true)) { //get from neighbor
              neighbor_requests.push_back(Int3(i,j,k));
              if(outside_images.find(n)==outside_images.end())
                outside_images[n] = vector<int>();
              outside_images[n].push_back(neighbor_requests.size()-1);
            }
          }
        }
      }
    }

    //3.2 Exchange info (id) w/ neighbors
    neicomm.Request(id, 1, neighbor_requests, neighbor_received, 0);
    assert(neighbor_received.size() == neighbor_requests.size());

    //3.3. Now, determine whether each ghost node is usable. If not, move further
    for(auto&& n : undetermined) {
      
      int i,j,k, counter=0;
      bool valid = true;
      for(int dk=0; dk<=1; dk++) {
        k = images_ijk[n][2]+dk;
        for(int dj=0; dj<=1; dj++) {
          j = images_ijk[n][1]+dj;
          for(int di=0; di<=1; di++) {
            i = images_ijk[n][0]+di;

            if(!Tag.IsHere(i,j,k,true)) { //get from neighbor
              int key = outside_images[n][counter++]; 
              if(neighbor_received[key]==INACTIVE_MATERIAL_ID &&
                 ghosts_ijk[n] != Int3(i,j,k)) 
                valid = false;
            }
            else {
              if(id[k][j][i]==INACTIVE_MATERIAL_ID &&
                 ghosts_ijk[n] != Int3(i,j,k))
                valid = false;
            }

            if(!valid)
              goto update_image; //no need to check others
          }
        }
      }

update_image:

      if(valid) { //found a valid image
        if(outside_images.find(n) != outside_images.end())
          ghosts_tag[n] = 2; //requires neighbor info
        else
          ghosts_tag[n] = 1; //do not need neighbor info
        found.insert(n);        
      }
      else { //image is invalid; move further
        Vec3D dir = images_xyz[n] - ghosts_xyz[n];
        images_xyz[n] = ghosts_xyz[n] + 1.2*dir;
        bool got = global_mesh.FindElementCoveringPoint(images_xyz[n], images_ijk[n], &images_xi[n],
                                                        true); //include ghost layer outside domain
        if(!got) {//constant extrapolation
          ghosts_tag[n] = 0; 
          found.insert(n);
        }
      }
    }     

    // remove "found" from "undetermined"
    for(auto&& n : found)
      assert(undetermined.erase(n)); //"erase" must return 1, i.e. "n" is found from undetermined

    if(undetermined.empty())
      break; //Hurray
  }

  int total_remaining = undetermined.size();
  MPI_Allreduce(MPI_IN_PLACE, &total_remaining, 1, MPI_INT, MPI_SUM, comm);
  
  if(total_remaining>0) {
    print_warning(comm, "Unable to find valid images for %d ghost nodes. "
                  "Apply const. extrapolation.\n", total_remaining);

    for(auto&& n : undetermined) 
      ghosts_tag[n] = 0;
  }

}

//-------------------------------------------------------------






//-------------------------------------------------------------

