#include<HyperelasticityOperator.h>

//------------------------------------------------------------

HyperelasticityOperator::HyperelasticityOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, 
                             IoData &iod_, SpaceVariable3D &coordinates_,
                             SpaceVariable3D &delta_xyz_, GlobalMeshInfo &global_mesh_,
                             std::vector<GhostPoint> &ghost_nodes_inner_,
                             std::vector<GhostPoint> &ghost_nodes_outer_)
                       : comm(comm_), iod(iod_), global_mesh(global_mesh_),
                         refmap(comm_, dm_all_, iod_, coordinates_, delta_xyz_,
                                global_mesh_, ghost_nodes_inner_, ghost_nodes_outer_)
{
  coordinates_.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates_.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);
}

//------------------------------------------------------------

HyperelasticityOperator::~HyperelasticityOperator()
{ }

//------------------------------------------------------------

void
HyperelasticityOperator::Destroy()
{
  refmap.Destroy();
}

//------------------------------------------------------------

void
HyperelasticityOperator::InitializeReferenceMap(SpaceVariable3D &Xi)
{
  refmap.SetInitialCondition(Xi);
}

//------------------------------------------------------------

void
HyperelasticityOperator::ApplyBoundaryConditionsToReferenceMap(SpaceVariable3D &Xi)
{
  refmap.ApplyBoundaryConditions(Xi);
}

//------------------------------------------------------------

void
HyperelasticityOperator::ComputeReferenceMapResidual(SpaceVariable3D &V, 
                                                     SpaceVariable3D &Xi, SpaceVariable3D &R)
{
  refmap.ComputeResidual(V,Xi,R);
}

//------------------------------------------------------------















//------------------------------------------------------------

