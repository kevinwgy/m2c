#include<MultiPhaseOperator.h>
#include<SpaceOperator.h>

//-----------------------------------------------------

MultiPhaseOperator::MultiPhaseOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_,
                                       SpaceOperator &spo, vector<LevelSetOperator*> &lso)
                  : comm(comm_), iod(iod_),
                    coordinates(spo.GetMeshCoordinates()),
                    delta_xyz(spo.GetMeshDeltaXYZ()),
                    Tag(comm_, &(dm_all_.ghosted1_1dof))
{
  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);

  for(int i=0; i<lso.size(); i++)
    ls2matid[i] = lso[i]->GetMaterialID();
}

//-----------------------------------------------------

MultiPhaseOperator::~MultiPhaseOperator()
{ }

//-----------------------------------------------------

MultiPhaseOperator::Destroy()
{
  Tag.Destroy();
}

//-----------------------------------------------------

void 
MultiPhaseOperator::UpdateMaterialID(vector<SpaceVariable3D*> &Phi, SpaceVariable3D &ID)
{
  // reset tag to 0
  Tag.SetConstantValue(0.0, true/*workOnGhost*/);
  ID.SetConstantValue(0.0, true/*workOnGhost*/);
  int overlap = 0;

  double*** tag = (double***)Tag.GetDataPointer();
  double*** id  = (double***)ID.GetDataPointer();

  for(int ls = 0; ls<Phi.size(); ls++) {//loop through all the level set functions
  
    double*** phi = (double***)Phi[ls]->GetDataPointer();
    int matid = ls2matid[ls];

    for(int k=kk0; k<kmax; k++)
      for(int j=jj0; j<jjmax; j++)
        for(int i=ii0; i<iimax; i++)
          if(phi[k][j][i]<0) {
            if(id[k][j][i] != 0) {
              overlap++;
              tag[k][j][i] = 1;
            } 
            id[k][j][i] = matid; 
          }

    Phi[ls]->RestoreDataPointerToLocalVector(); //no changes made
  }

  MPI_Allreduce(MPI_IN_PLACE, &overlap, 1, MPI_INT, MPI_SUM, comm);


  if(overlap) {
    print_error("Error: Found overlapping material interfaces. Number of overlapped cells: %d.\n", overlap);
    exit_mpi();
  } 


  if(overlap) 
    Tag.RestoreDataPointerAndInsert();
  else
    Tag.RestoreDataPointerToLocalVector();
}

//-----------------------------------------------------

void
MultiPhaseOperator::UpdateStateVariablesAfterInterfaceMotion(SpaceVariable3D &IDn, 
                        SpaceVariable3D &ID, SpaceVariable3D &V)
{






}

//-----------------------------------------------------

