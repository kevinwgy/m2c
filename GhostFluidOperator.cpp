#include<GhostFluidOperator.h>
#include<EmbeddedBoundaryFormula.h>
#include<GhostPoint.h>

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
{ }

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
  std::set<Int3> ghost_ijk;

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
                  ghost_ijk.insert(Int3(i,j,k)); //first-layer ghost
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
                  ghost_ijk.insert(Int3(i,j,k)); //second-layer ghost
                  tag[k][j][i] = 2;
                  goto found_layer_2_node;
                }
              }
        }
found_layer_2_node: { } 
      }

  
  ID.RestoreDataPointerToLocalVector();
  Tag.RestoreDataPointerAndInsert();

/*


  // Step 2: Find image point
  vector<GhostPoint> ghosts;
  for(auto&& g : ghost_ijk) {


  }
*/


  return 0;
}

//-------------------------------------------------------------


//-------------------------------------------------------------

