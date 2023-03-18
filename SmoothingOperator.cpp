/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<SmoothingOperator.h>
#include<cmath> //exp

//--------------------------------------------------------------------------

SmoothingOperator::SmoothingOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, SmoothingData &iod_smooth_,
                                     SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_,
                                     SpaceVariable3D &volume_)
                  : comm(comm_), iod_smooth(iod_smooth_), coordinates(coordinates_), volume(volume_),
                    delta_xyz(delta_xyz_),
                    V0(comm_, &(dm_all_.ghosted1_5dof))
{
  if(iod_smooth.type == SmoothingData::NONE || iod_smooth.iteration<1) {
    print_error("*** Error: Attempting to create SmoothingOperator with incorrect inputs. Type = %d, "
                "NumberOfIteration = %d.\n", iod_smooth.type, iod_smooth.iteration);
    exit_mpi();
  }
}

//--------------------------------------------------------------------------

SmoothingOperator::~SmoothingOperator()
{ }

//--------------------------------------------------------------------------

void SmoothingOperator::Destroy()
{
  V0.Destroy();
}

//--------------------------------------------------------------------------

void SmoothingOperator::ApplySmoothingFilter(SpaceVariable3D &V, SpaceVariable3D *ID)
{
  if(iod_smooth.type == SmoothingData::BOX)
    ApplyBoxFilter(V, ID);
  else if(iod_smooth.type == SmoothingData::GAUSSIAN)
    ApplyGausianFilter(V, ID);
}

//--------------------------------------------------------------------------

void SmoothingOperator::ApplyBoxFilter(SpaceVariable3D &V, SpaceVariable3D *ID)
{

  int Vdim = V.NumDOF(), V0dim = V0.NumDOF();
  if(Vdim>V0dim) {
    print_error("*** Error: Size of the internal variable in SmoothingOperator needs"
                " to be increased (%d vs. %d).\n", V0dim, Vdim);
    exit_mpi();
  } 

  int myid;
  double sum_weight = 0;

  //copy data from V to V0
  vector<int> ind;
  for(int i=0; i<Vdim; i++) 
    ind.push_back(i);
  V0.AXPlusBY(0.0, 1.0, V, ind, ind, true);

  //get mesh info
  int i0, j0, k0, imax, jmax, kmax;
  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);

  double*** v  = (double***)V.GetDataPointer();
  double*** v0 = (double***)V0.GetDataPointer();
  double*** id = ID ? (double***)ID->GetDataPointer() : NULL;
 
  //apply smoothing in the interior of the subdomain
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
            
        myid = ID ? id[k][j][i] : 0;

        // clear data
        sum_weight = 0.0;
        setValue(&v[k][j][i*Vdim], 0.0, Vdim);

        for(int neighk=k-1; neighk<=k+1; neighk++)
          for(int neighj=j-1; neighj<=j+1; neighj++)
            for(int neighi=i-1; neighi<=i+1; neighi++) {

              if(ID && myid != id[neighk][neighj][neighi])
                continue;

              //if(coordinates.OutsidePhysicalDomainAndUnpopulated(neighi,neighj,neighk))
              if(coordinates.OutsidePhysicalDomain(neighi,neighj,neighk))
                continue; //no valid data here.

              //fprintf(stdout,"(%d %d %d) %d : neigh (%d %d %d), weight = %e, rho = %e\n", i,j,k, (int)myid, neighi, neighj, neighk, 1.0, v0[neighk][neighj][neighi*V0dim+0]);

              sum_weight += 1.0;
              addarray(&v0[neighk][neighj][neighi*V0dim], &v[k][j][i*Vdim], Vdim);
            //  for(int p=0; p<Vdim; p++)
             //   v[k][j][i*Vdim+p] += v0[neighk][neighj][neighi*V0dim+p];
            }

        for(int p=0; p<Vdim; p++)
          v[k][j][i*Vdim+p] /= sum_weight;

      }

  V0.RestoreDataPointerToLocalVector();
  if(ID) ID->RestoreDataPointerToLocalVector();

  V.RestoreDataPointerAndInsert();

  if(iod_smooth.conservation == SmoothingData::ON)
    EnforceLocalConservation(V0, V, ID);

}

//--------------------------------------------------------------------------

void SmoothingOperator::ApplyGausianFilter(SpaceVariable3D &V, SpaceVariable3D *ID)
{

  int Vdim = V.NumDOF(), V0dim = V0.NumDOF();
  if(Vdim>V0dim) {
    print_error("*** Error: Size of the internal variable in SmoothingOperator needs"
                " to be increased (%d vs. %d).\n", V0dim, Vdim);
    exit_mpi();
  } 

  //get mesh info
  int i0, j0, k0, imax, jmax, kmax;
  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
  Vec3D*** dxyz = (Vec3D***)delta_xyz.GetDataPointer();
 
  int myid;
  double weight = 0, sum_weight = 0;
  double variance = 0.0; 

  //copy data from V to V0
  vector<int> ind;
  for(int i=0; i<Vdim; i++) 
    ind.push_back(i);
  V0.AXPlusBY(0.0, 1.0, V, ind, ind, true);

  double*** v  = (double***)V.GetDataPointer();
  double*** v0 = (double***)V0.GetDataPointer();
  double*** id = ID ? (double***)ID->GetDataPointer() : NULL;

  //apply smoothing in the interior of the subdomain
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
            
        myid = ID ? id[k][j][i] : 0;
        Vec3D &x(coords[k][j][i]);

        // clear data
        sum_weight = 0.0;
        setValue(&v[k][j][i*Vdim], 0.0, Vdim);

        // set sigma^2 (the variance of Gaussian distribution)
        variance = iod_smooth.sigma_factor*std::min(std::min(dxyz[k][j][i][0],dxyz[k][j][i][1]),dxyz[k][j][i][2]);
        variance = variance*variance;

        for(int neighk=k-1; neighk<=k+1; neighk++)
          for(int neighj=j-1; neighj<=j+1; neighj++)
            for(int neighi=i-1; neighi<=i+1; neighi++) {

              if(ID && myid != id[neighk][neighj][neighi])
                continue;

              //if(coordinates.OutsidePhysicalDomainAndUnpopulated(neighi,neighj,neighk))
              if(coordinates.OutsidePhysicalDomain(neighi,neighj,neighk))
                continue; //no valid data here.

              Vec3D &x1(coords[neighk][neighj][neighi]);

              weight = exp(-pow((x1-x).norm(),2)/(2.0*variance));

//              fprintf(stdout,"(%d %d %d) %d : neigh (%d %d %d), weight = %e, rho = %e\n", i,j,k, (int)myid, neighi, neighj, neighk, weight, v0[neighk][neighj][neighi*V0dim+0]);

              sum_weight += weight;
              for(int p=0; p<Vdim; p++)
                v[k][j][i*Vdim+p] += weight*v0[neighk][neighj][neighi*V0dim+p];
            }

        for(int p=0; p<Vdim; p++)
          v[k][j][i*Vdim+p] /= sum_weight;

      }

  coordinates.RestoreDataPointerToLocalVector();
  delta_xyz.RestoreDataPointerToLocalVector();
  V0.RestoreDataPointerToLocalVector();
  if(ID) ID->RestoreDataPointerToLocalVector();

  V.RestoreDataPointerAndInsert();

  if(iod_smooth.conservation == SmoothingData::ON)
    EnforceLocalConservation(V0, V, ID);
}

//--------------------------------------------------------------------------

void SmoothingOperator::EnforceLocalConservation(SpaceVariable3D &U0, SpaceVariable3D &U, SpaceVariable3D *ID)
{

  int Udim = U.NumDOF(), U0dim = U0.NumDOF();

  //get mesh info
  int i0, j0, k0, imax, jmax, kmax;
  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  double*** vol = (double***)volume.GetDataPointer();
 
  int myid;

  double*** u  = (double***)U.GetDataPointer();
  double*** u0 = (double***)U0.GetDataPointer();
  double*** id = ID ? (double***)ID->GetDataPointer() : NULL;

  double uold[Udim], unew[Udim];

  
  for(int iter=0;  iter<5;  iter++) {

    double maxdiff = 0.0;

    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++) {
              
          myid = ID ? id[k][j][i] : 0;

          // clear data
          for(int p=0; p<Udim; p++) {
            uold[p] = 0.0;
            unew[p] = 0.0;
          }

          for(int neighk=k-1; neighk<=k+1; neighk++)
            for(int neighj=j-1; neighj<=j+1; neighj++)
              for(int neighi=i-1; neighi<=i+1; neighi++) {

                if(ID && myid != id[neighk][neighj][neighi])
                  continue;

                //if(coordinates.OutsidePhysicalDomainAndUnpopulated(neighi,neighj,neighk))
                if(coordinates.OutsidePhysicalDomain(neighi,neighj,neighk))
                  continue; //no valid data here.

                for(int p=0; p<Udim; p++) {
                  uold[p] += vol[neighk][neighj][neighi]*u0[neighk][neighj][neighi*U0dim+p];
                  unew[p] += vol[neighk][neighj][neighi]*u[neighk][neighj][neighi*Udim+p];
                }

              }

          for(int p=0; p<Udim; p++) {
            if(unew[p]!=0.0) {
              double factor = uold[p]/unew[p];
              //double maxdifftmp = std::max(maxdiff, std::fabs((unew[p]-uold[p])/unew[p]));
              //fprintf(stdout,"(%d,%d,%d) p = %d, old = %e, new = %e, factor = %e, maxdiff = %e.\n", i,j,k, p, uold[p], unew[p], factor, maxdifftmp);
              factor = std::min(1.05, factor);
              factor = std::max(0.95, factor);
              u[k][j][Udim*i+p] *= factor; 
              maxdiff = std::max(maxdiff, std::fabs((unew[p]-uold[p])/unew[p]));
            }
          }
        }

    MPI_Allreduce(MPI_IN_PLACE, &maxdiff, 1, MPI_DOUBLE, MPI_MAX, comm);
    print("  o Enforcing conversation (It. %d): maxdiff = %e.\n", iter, maxdiff);
    if(maxdiff<iod_smooth.conservation_tol)
      break;
  }

  volume.RestoreDataPointerToLocalVector();
  U0.RestoreDataPointerToLocalVector();
  if(ID) ID->RestoreDataPointerToLocalVector();

  U.RestoreDataPointerAndInsert();

}

//--------------------------------------------------------------------------

