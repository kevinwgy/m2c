/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<GhostFluidOperator.h>
#include<EmbeddedBoundaryFormula.h>
#include<EmbeddedBoundaryDataSet.h>
#include<cfloat> //DBL_MAX
#include<map>

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
                        SpaceVariable3D &Vgf)
{
  assert(EBDS && EBDS->size()>0);

  // ----------------------------------------------------------------------------
  // Step 1: Identify ghost nodes inside the subdomain that need to be populated
  // ----------------------------------------------------------------------------
  std::set<Int3> ghost_set;

  double*** id = ID.GetDataPointer();
  double*** tag = Tag.GetDataPointer();

  // round 1
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
                  ghost_set.insert(Int3(i,j,k)); //first-layer ghost
                  tag[k][j][i] = 1;
                  goto found_layer_1_node;
                }
              }
        }
found_layer_1_node: { } //needs { } to avoid compilation error
      }

  Tag.RestoreDataPointerAndInsert();

  tag = Tag.GetDataPointer();
  
  // round 2
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        if(tag[k][j][i]==1)
          continue;
        if(id[k][j][i]==INACTIVE_MATERIAL_ID) {
          for(int kk=k-1; kk<=k+1; kk++)
            for(int jj=j-1; jj<=j+1; jj++)
              for(int ii=i-1; ii<=i+1; ii++) {
                if(tag[kk][jj][ii]==1) {
                  ghost_set.insert(Int3(i,j,k)); //second-layer ghost
                  tag[k][j][i] = 2;
                  goto found_layer_2_node;
                }
              }
        }
found_layer_2_node: { } 
      }

  
  Tag.RestoreDataPointerAndInsert();




  // ----------------------------------------------------------------------------
  // Step 2: Find image points and the needed geometric information
  // ----------------------------------------------------------------------------
  vector<double***> cpi;
  vector<vector<std::pair<Int3, ClosestPoint> > * > closest_points;

  for(auto&& data : *EBDS) {
    assert(data->ClosestPointIndex_ptr);
    cpi.push_back(data->ClosestPointIndex_ptr->GetDataPointer());
    assert(data->closest_points_ptr);
    closest_points.push_back(data->closest_points_ptr);
  }

  vector<Int3> ghosts_ijk(ghost_set.begin(), ghost_set.end());
  vector<Int3> images_ijk;
  vector<Vec3D> images_xi;
  vector<Vec3D> images_xyz;
  vector<bool> const_extrap(ghosts_ijk.size(), false); //apply constant extrapolation
  vector<Vec3D> vp; //velocity at projection points on the interface
  for(auto&& g : ghosts_ijk) {
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

    assert(mysurf>=0);
    int pid = cpi[mysurf][g[2]][g[1]][g[0]];
    ClosestPoint& cp((*closest_points[mysurf])[pid].second);

    // all the nodes of the embedded surface closest to this ghost node
    vector<Vec3D> &X((*EBDS)[mysurf]->surface_ptr->X);
    vector<Vec3D> &Udot((*EBDS)[mysurf]->surface_ptr->Udot);

    Int3& tnodes((*EBDS)[mysurf]->surface_ptr->elems[cp.tid]);
    Vec3D xp = cp.xi[0]*X[tnodes[0]] + cp.xi[1]*X[tnodes[1]] 
             + cp.xi[2]*X[tnodes[2]]; 

    // add velocity at projeciton point
    vp.push_back(cp.xi[0]*Udot[tnodes[0]] + cp.xi[1]*Udot[tnodes[1]]
                                          + cp.xi[2]*Udot[tnodes[2]]);

    Vec3D xim = 2.0*xp - global_mesh.GetXYZ(g); //image point
    images_xyz.push_back(xim);

    // get image point location w.r.t. mesh
    images_ijk.push_back(Int3());
    images_xi.push_back(Vec3D());    
    bool found = global_mesh.FindElementCoveringPoint(xim, images_ijk.back(), &images_xi.back(),
                                                      true); //include ghost layer outside domain
    if(!found) {
      fprintf(stderr,"\033[0;35mWarning: Detected a ghost node (%d,%d,%d) whose image is outside "
                     "physical domain.\033[0m\n", g[0], g[1], g[2]);
      const_extrap[images_xyz.size()-1] = true;
    }

  }

  for(auto&& data : *EBDS)
    data->ClosestPointIndex_ptr->RestoreDataPointerToLocalVector();



//  vector<Int3> ghosts_ijk(ghost_set.begin(), ghost_set.end());
//  vector<Int3> images_ijk;
//  vector<Vec3D> images_xi;
  // ----------------------------------------------------------------------------
  // Step 3: Find a usable image for each ghost node, starting w/ the one found above
  // ----------------------------------------------------------------------------
  vector<Int3> neighbor_requests;
  map<int, vector<int> > outside_images;
  for(int iter=0; iter<5; iter++) {

    neighbor_requests.clear();
    outside_images.clear();

    // Figure out the need of information ("id") from neighbors
    for(int n=0; n<(int)ghosts_ijk.size(); n++) {

      if(const_extrap[n])
        continue; //this ghost node will be populated by constant extrapolation

      int i,j,k;
      for(int dk=0; dk<=1; dk++) {
        k = images_ijk[n][2]+dk;
        for(int dj=0; dj<=1; dj++) {
          j = images_ijk[n][1]+dj;
          for(int di=0; di<=1; di++) {
            i = images_ijk[n][0]+di;
            if(!Tag.IsHere(i,j,k,true)) {
              neighbor_requests.push_back(Int3(i,j,k));
              if(outside_images.find(n)==outside_images.end())
                outside_images[n] = vector<int>();
              outside_images[n].push_back(neighbor_requests.size()-1);
            }
          }
        }
      }

    }

    vector<double> neighbor_received; //will be fixed by neicomm
    neicomm.Request(id, 1, neighbor_requests, neighbor_received, 0);

  }
 
/*
  for(int n=0; n<(int)ghosts_ijk.size(); n++) {

    bool valid_image = false;

    if(const_extrap[n])
      continue;


      goto apply_const_extrapolation;

    Vec3D dir = images_xyz[n] - global_mesh.GetXYZ(ghost_ijk[n]);
    double dist = dir.norm();
    dir /= dist;

    // if dist from ghost to interface is less than 15% of element size, just apply
    // constant extrapolation
    if(dist<0.3*global_mesh.GetMinDXYZ(ghost_ijk[n]))
      goto apply_const_extrapolation;


    for(int iter=0; iter<5; iter++) {
      valid_image = CheckFeasibilityOfInterpolation(ghosts_ijk[n], images_ijk[n], images_xi[n], id);
      if(valid_image)
        break;
      // move a bit further into the domain 
      images_xyz[n] = global_mesh.GetXYZ(ghost_ijk[n]) + 1.2*dist*dir;
      bool found = global_mesh.FindElementCoveringPoint(images_xyz[n], images_ijk[n], &images_xi[n],
                                                        true); //include ghost layer outside domain
      if(!found)
        break; //apply const. extrapolation
    }





  for(int n=0; n<(int)ghosts_ijk.size(); n++) {

    bool valid_image = false;

    // if image is outside physical domain, apply constant extrapolation
    if(images_out[n])
      goto apply_const_extrapolation;

    Vec3D dir = images_xyz[n] - global_mesh.GetXYZ(ghost_ijk[n]);
    double dist = dir.norm();
    dir /= dist;

    // if dist from ghost to interface is less than 15% of element size, just apply
    // constant extrapolation
    if(dist<0.3*global_mesh.GetMinDXYZ(ghost_ijk[n]))
      goto apply_const_extrapolation;


    for(int iter=0; iter<5; iter++) {
      valid_image = CheckFeasibilityOfInterpolation(ghosts_ijk[n], images_ijk[n], images_xi[n], id);
      if(valid_image)
        break;
      // move a bit further into the domain 
      images_xyz[n] = global_mesh.GetXYZ(ghost_ijk[n]) + 1.2*dist*dir;
      bool found = global_mesh.FindElementCoveringPoint(images_xyz[n], images_ijk[n], &images_xi[n],
                                                        true); //include ghost layer outside domain
      if(!found)
        break; //apply const. extrapolation
    }

     
apply_const_extrapolation: 

    if(!valid_image) {

      for(int p=0; p<3; p++)
        v[ghost_ijk[n][2]][ghost_ijk[n][1]][ghost_ijk[n][0]][1+p] = vp[n][p];

      continue; //done with this one!
    }

    // I AM HERE
  }
*/


  ID.RestoreDataPointerToLocalVector();

  return 0;
}

//-------------------------------------------------------------
/*
bool
GhostFluidOperator::CheckFeasibilityOfInterpolation(Int3& ghost, Int3& image, Vec3D& image_xi, double*** id, *** REQUEST ***)
{
  for(int k=0; k<=1; k++)
    for(int j=0; j<=1; j++)
      for(int i=0; i<=1; i++) {

        if(global_mesh.Find
        image[2]+k

      }

}
*/
//-------------------------------------------------------------





//-------------------------------------------------------------

