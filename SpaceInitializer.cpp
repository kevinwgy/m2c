/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <SpaceInitializer.h>
#include <SpaceOperator.h>
#include <LevelSetOperator.h>
#include <EmbeddedBoundaryDataSet.h>

using std::vector;

//---------------------------------------------------------------------

SpaceInitializer::SpaceInitializer(MPI_Comm &comm_, IoData &iod_, GlobalMeshInfo &global_mesh_)
                : comm(comm_), iod(iod_), global_mesh(global_mesh_)
{ }

//---------------------------------------------------------------------

SpaceInitializer::~SpaceInitialier()
{ }

//---------------------------------------------------------------------

void
SpaceInitializer::Destroy()
{ }

//---------------------------------------------------------------------

std::multimap<int, std::pair<int,int> >
SpaceInitializer::SetInitialCondition(SpaceVariable3D &V, SpaceVariable3D &ID, vector<SpaceVariable3D*> &Phi,
                                      SpaceOperator &spo, vector<SpaceVariable3D*> &lso,
                                      std::unique_ptr<vector<std::unique_ptr<EmbeddedBoundaryDataSet> > > EBDS)
{


}

//---------------------------------------------------------------------





















//---------------------------------------------------------------------

