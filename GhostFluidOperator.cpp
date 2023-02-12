#include<GhostFluidOperator.h>
#include<SpaceVariable.h>

using std::vector;
extern int verbose;
extern int INACTIVE_MATERIAL_ID;

//-------------------------------------------------------------

GhostFluidOperator::GhostFluidOperator(MPI_Comm &comm_, GlobalMeshInfo &global_mesh_)
                  : comm(comm_), global_mesh(global_mesh_),
                    neicomm_ptr(NULL)
{

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  assert(rank>=0 && rank<size);

  neicomm_ptr = new NeighborCommunicator(comm, global_mesh.GetAllNeighborsOfSub(rank),
                                         global_mesh.GetFaceEdgeNeighborsOfSub(rank),
                                         global_mesh.GetFaceNeighborsOfSub(rank));

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

  // Step 1: Tag ghost nodes that need to be populated
  vector<Int3> ghosts;

/*
  double*** id = (double***) ID.GetDataPointer();

  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        if(
  
      }
  
  ID.RestoreDataPointerToLocalVector();
*/
  return 0;
}

//-------------------------------------------------------------


//-------------------------------------------------------------

