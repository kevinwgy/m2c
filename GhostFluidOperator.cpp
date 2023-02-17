#include<GhostFluidOperator.h>
#include<EmbeddedBoundaryFormula.h>
#include<EmbeddedBoundaryDataSet.h>
#include <cfloat> //DBL_MAX

using std::vector;

extern int verbose;
extern int INACTIVE_MATERIAL_ID;

//-------------------------------------------------------------

GhostFluidOperator::GhostFluidOperator(MPI_Comm &comm_, DataManagers3D &dm_all_,
                                       GlobalMeshInfo &global_mesh_)
                  : comm(comm_), global_mesh(global_mesh_),
                    neicomm_ptr(NULL),
                    Tag(comm_, &(dm_all_.ghosted1_1dof))
{

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  assert(rank>=0 && rank<size);

  neicomm_ptr = new NeighborCommunicator(comm, global_mesh.GetAllNeighborsOfSub(rank),
                                         global_mesh.GetFaceEdgeNeighborsOfSub(rank),
                                         global_mesh.GetFaceNeighborsOfSub(rank));

  Tag.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  Tag.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);
  Tag.GetGlobalSize(&NX,&NY,&NZ);

}

//-------------------------------------------------------------

GhostFluidOperator::~GhostFluidOperator()
{
  if(neicomm_ptr)
    delete neicomm_ptr;
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

  
  ID.RestoreDataPointerToLocalVector();
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

    // get image point location w.r.t. mesh
    images_ijk.push_back(Int3());
    images_xi.push_back(Vec3D());    
    bool found = global_mesh.FindElementCoveringPoint(xim, images_ijk.back(), &images_xi.back(),
                                                      true); //include ghost layer outside domain
    if(!found) {
      fprintf(stderr,"\033[0;31m***Error: Detected a ghost node (%d,%d,%d) whose image is outside "
                     "physical domain.\033[0m\n", g[0], g[1], g[2]);
      exit(-1);
    }

  }

  for(auto&& data : *EBDS)
    data->ClosestPointIndex_ptr->RestoreDataPointerToLocalVector();



  // ----------------------------------------------------------------------------
  // Step 3: Apply linear extrapolation to populate ghost nodes
  // ----------------------------------------------------------------------------


  return 0;
}

//-------------------------------------------------------------


//-------------------------------------------------------------

