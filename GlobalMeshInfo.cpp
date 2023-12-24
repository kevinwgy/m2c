/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<GlobalMeshInfo.h>
#include<SpaceVariable.h>
#include<algorithm> //std::upper_bound

//------------------------------------------------------------------

GlobalMeshInfo::GlobalMeshInfo(std::vector<double> &x_glob_, std::vector<double> &y_glob_,
                               std::vector<double> &z_glob_, std::vector<double> &dx_glob_,
                               std::vector<double> &dy_glob_, std::vector<double> &dz_glob_)
{
  x_glob = x_glob_;   y_glob = y_glob_;   z_glob = z_glob_;
  dx_glob = dx_glob_; dy_glob = dy_glob_; dz_glob = dz_glob_;
  assert(x_glob.size() == dx_glob.size());
  assert(y_glob.size() == dy_glob.size());
  assert(z_glob.size() == dz_glob.size());

  xyz_min[0] = x_glob.front() - 0.5*dx_glob.front();
  xyz_max[0] = x_glob.back()  + 0.5*dx_glob.back();
  xyz_min[1] = y_glob.front() - 0.5*dy_glob.front();
  xyz_max[1] = y_glob.back()  + 0.5*dy_glob.back();
  xyz_min[2] = z_glob.front() - 0.5*dz_glob.front();
  xyz_max[2] = z_glob.back()  + 0.5*dz_glob.back();

  one_dimensional_mesh = (y_glob.size()==1) && (z_glob.size()==1);

  two_dimensional_mesh = (!one_dimensional_mesh) && (z_glob.size()==1);

  NX = x_glob.size();
  NY = y_glob.size();
  NZ = z_glob.size();

  domain_volume = ((x_glob[NX-1] + 0.5*dx_glob[NX-1]) - (x_glob[0] - 0.5*dx_glob[0]))
                * ((y_glob[NY-1] + 0.5*dy_glob[NY-1]) - (y_glob[0] - 0.5*dy_glob[0]))
                * ((z_glob[NZ-1] + 0.5*dz_glob[NZ-1]) - (z_glob[0] - 0.5*dz_glob[0]));
}

//------------------------------------------------------------------

GlobalMeshInfo::~GlobalMeshInfo()
{ }

//------------------------------------------------------------------

void
GlobalMeshInfo::FindSubdomainInfo(MPI_Comm& comm, DataManagers3D& dms)
{
  // Create a space variable only for use in this function. Destroyed at the end.
  SpaceVariable3D S(comm, &(dms.ghosted1_1dof));

  // Step 1: Construct subD_ijk_min/max
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  assert(rank>=0 && rank<size);

  int i0, j0, k0, imax, jmax, kmax;
  S.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);

  subD_ijk_min.resize(size);
  subD_ijk_min[rank] = Int3(i0,j0,k0);
  subD_ijk_max.resize(size);
  subD_ijk_max[rank] = Int3(imax,jmax,kmax);

  MPI_Allgather(MPI_IN_PLACE, 3, MPI_INT, (int*)subD_ijk_min.data(), 3, MPI_INT, comm); 
  MPI_Allgather(MPI_IN_PLACE, 3, MPI_INT, (int*)subD_ijk_max.data(), 3, MPI_INT, comm); 
  
  // Step 2: Construct subD_xyz_min/max
  subD_xyz_min.resize(size);
  subD_xyz_min[rank] = Vec3D(x_glob[i0]-0.5*dx_glob[i0], y_glob[j0]-0.5*dy_glob[j0],
                             z_glob[k0]-0.5*dz_glob[k0]);
  subD_xyz_max.resize(size);
  subD_xyz_max[rank] = Vec3D(x_glob[imax-1]+0.5*dx_glob[imax-1],
                             y_glob[jmax-1]+0.5*dy_glob[jmax-1],
                             z_glob[kmax-1]+0.5*dz_glob[kmax-1]);

  MPI_Allgather(MPI_IN_PLACE, 3, MPI_DOUBLE, (double*)subD_xyz_min.data(), 3, MPI_DOUBLE, comm); 
  MPI_Allgather(MPI_IN_PLACE, 3, MPI_DOUBLE, (double*)subD_xyz_max.data(), 3, MPI_DOUBLE, comm); 

  // Step 3: Find neighbours
  int none = -1;
  S.SetConstantValue(none, true); //ghosts outside physical domain will get 'none'

  double*** s = S.GetDataPointer();
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++)
        s[k][j][i] = rank;

  S.RestoreDataPointerAndInsert();

  s = S.GetDataPointer();
  std::vector<int> neighbors_all(size*27, none);

  int myi, myj, myk;
  for(int k=0; k<3; k++)
    for(int j=0; j<3; j++)
      for(int i=0; i<3; i++) {

        if(i==0)      myi = i0-1;
        else if(i==2) myi = imax;
        else          myi = i0; //doesn't have to be i0

        if(j==0)      myj = j0-1;
        else if(j==2) myj = jmax;
        else          myj = j0;

        if(k==0)      myk = k0-1;
        else if(k==2) myk = kmax;
        else          myk = k0;

        neighbors_all[27*rank + 9*k + 3*j + i] = (int)s[myk][myj][myi];
      }

  MPI_Allgather(MPI_IN_PLACE, 27, MPI_INT, neighbors_all.data(), 27, MPI_INT, comm); 
  S.RestoreDataPointerToLocalVector();

  // Step 4: Setup all the neighbour vectors
  subD_neighbors_27.resize(size);
  subD_neighbors_all.resize(size);
  subD_neighbors_19.resize(size);
  subD_neighbors_face_edge.resize(size);
  subD_neighbors_7.resize(size);
  subD_neighbors_face.resize(size);

  int neigh, neigh_type;
  for(int proc=0; proc<size; proc++) {
    for(int k=0; k<3; k++)
      for(int j=0; j<3; j++)
        for(int i=0; i<3; i++) {

          neigh = neighbors_all[27*proc+9*k+3*j+i];

          neigh_type = 0;
          if(k==0 || k==2) neigh_type++;
          if(j==0 || j==2) neigh_type++;
          if(i==0 || i==2) neigh_type++;

          subD_neighbors_27[proc].push_back(neigh);
          if(neigh_type!=0 && neigh!=none)
            subD_neighbors_all[proc].push_back(neigh);

          if(neigh_type!=3) {
            subD_neighbors_19[proc].push_back(neigh);
            if(neigh_type!=0 && neigh!=none)
              subD_neighbors_face_edge[proc].push_back(neigh);
          }

          if(neigh_type!=3 && neigh_type!=2) {
            subD_neighbors_7[proc].push_back(neigh);
            if(neigh_type!=0 && neigh!=none)
              subD_neighbors_face[proc].push_back(neigh);
          } 
        }

    assert(subD_neighbors_27[proc].size()==27);
    assert(subD_neighbors_all[proc].size()<=26);
    assert(subD_neighbors_19[proc].size()==19);
    assert(subD_neighbors_face_edge[proc].size()<=18);
    assert(subD_neighbors_7[proc].size()==7);
    assert(subD_neighbors_face[proc].size()<=6);
  }


  S.Destroy();

}

//------------------------------------------------------------------

double
GlobalMeshInfo::GetX(int i) {
  if(i<0)
    return x_glob[0] + i*dx_glob[0];
  else if(i>=(int)x_glob.size())
    return x_glob.back() + (i-x_glob.size()+1)*dx_glob.back();
  else
    return x_glob[i];
}

//------------------------------------------------------------------

double
GlobalMeshInfo::GetY(int j) {
  if(j<0)
    return y_glob[0] + j*dy_glob[0];
  else if(j>=(int)y_glob.size())
    return y_glob.back() + (j-y_glob.size()+1)*dy_glob.back();
  else
    return y_glob[j];
}

//------------------------------------------------------------------

double
GlobalMeshInfo::GetZ(int k) {
  if(k<0)
    return z_glob[0] + k*dz_glob[0];
  else if(k>=(int)z_glob.size())
    return z_glob.back() + (k-z_glob.size()+1)*dz_glob.back();
  else
    return z_glob[k];
}

//------------------------------------------------------------------

double
GlobalMeshInfo::GetDx(int i) {
  if(i<0)
    return dx_glob[0];
  else if(i>=(int)dx_glob.size())
    return dx_glob.back();
  else
    return dx_glob[i];
}

//------------------------------------------------------------------

double
GlobalMeshInfo::GetDy(int j) {
  if(j<0)
    return dy_glob[0];
  else if(j>=(int)dy_glob.size())
    return dy_glob.back();
  else
    return dy_glob[j];
}

//------------------------------------------------------------------

double
GlobalMeshInfo::GetDz(int k) {
  if(k<0)
    return dz_glob[0];
  else if(k>=(int)dz_glob.size())
    return dz_glob.back();
  else
    return dz_glob[k];
}

//------------------------------------------------------------------

double
GlobalMeshInfo::GetX(Int3 ijk)
{
  return GetX(ijk[0]);
}

//------------------------------------------------------------------

double
GlobalMeshInfo::GetY(Int3 ijk)
{
  return GetY(ijk[1]);
}

//------------------------------------------------------------------

double
GlobalMeshInfo::GetZ(Int3 ijk)
{
  return GetZ(ijk[2]);
}

//------------------------------------------------------------------

double
GlobalMeshInfo::GetDx(Int3 ijk)
{
  return GetDx(ijk[0]);
}

//------------------------------------------------------------------

double
GlobalMeshInfo::GetDy(Int3 ijk)
{
  return GetDy(ijk[1]);
}

//------------------------------------------------------------------

double
GlobalMeshInfo::GetDz(Int3 ijk)
{
  return GetDz(ijk[2]);
}

//------------------------------------------------------------------

Vec3D
GlobalMeshInfo::GetXYZ(Int3 ijk)
{
  return Vec3D(GetX(ijk[0]), GetY(ijk[1]), GetZ(ijk[2]));
}

//------------------------------------------------------------------

Vec3D
GlobalMeshInfo::GetXYZ(int i, int j, int k)
{
  return Vec3D(GetX(i), GetY(j), GetZ(k));
}

//------------------------------------------------------------------

Vec3D
GlobalMeshInfo::GetDXYZ(Int3 ijk)
{
  return Vec3D(GetDx(ijk[0]), GetDy(ijk[1]), GetDz(ijk[2]));
}

//------------------------------------------------------------------

Vec3D
GlobalMeshInfo::GetDXYZ(int i, int j, int k)
{
  return Vec3D(GetDx(i), GetDy(j), GetDz(k));
}

//------------------------------------------------------------------

double
GlobalMeshInfo::GetMinDXYZ(Int3 ijk)
{
  return two_dimensional_mesh ? std::min(GetDx(ijk[0]), GetDy(ijk[1]))
                              : std::min(std::min(GetDx(ijk[0]), GetDy(ijk[1])), 
                                         GetDz(ijk[2]));
}

//------------------------------------------------------------------

double
GlobalMeshInfo::GetMaxDXYZ(Int3 ijk)
{
  return two_dimensional_mesh ? std::max(GetDx(ijk[0]), GetDy(ijk[1]))
                              : std::max(std::max(GetDx(ijk[0]), GetDy(ijk[1])),
                                         GetDz(ijk[2]) );
}

//------------------------------------------------------------------

bool
GlobalMeshInfo::IsPointInDomain(Vec3D &p, bool include_ghost_layer)
{

  int Nx = x_glob.size();
  int Ny = y_glob.size();
  int Nz = z_glob.size();

  double shift = include_ghost_layer ? 1.5 : 0.5;

  if(p[0] < x_glob[0]    - shift*dx_glob[0])     return false;
  if(p[0] > x_glob[Nx-1] + shift*dx_glob[Nx-1])  return false;
  if(p[1] < y_glob[0]    - shift*dy_glob[0])     return false;
  if(p[1] > y_glob[Ny-1] + shift*dy_glob[Ny-1])  return false;
  if(p[2] < z_glob[0]    - shift*dz_glob[0])     return false;
  if(p[2] > z_glob[Nz-1] + shift*dz_glob[Nz-1])  return false;

  return true;

}

//------------------------------------------------------------------

bool
GlobalMeshInfo::IsPointInNodalMesh(Vec3D &p, bool include_ghost_layer)
{

  int Nx = x_glob.size();
  int Ny = y_glob.size();
  int Nz = z_glob.size();

  double shift = include_ghost_layer ? 1.0 : 0.0;

  if(p[0] < x_glob[0]    - shift*dx_glob[0])     return false;
  if(p[0] > x_glob[Nx-1] + shift*dx_glob[Nx-1])  return false;
  if(p[1] < y_glob[0]    - shift*dy_glob[0])     return false;
  if(p[1] > y_glob[Ny-1] + shift*dy_glob[Ny-1])  return false;
  if(p[2] < z_glob[0]    - shift*dz_glob[0])     return false;
  if(p[2] > z_glob[Nz-1] + shift*dz_glob[Nz-1])  return false;

  return true;

}

//------------------------------------------------------------------

Int3
GlobalMeshInfo::FindClosestNodeToPoint(Vec3D &p, bool include_ghost_layer)
{
  Int3 ijk;
  double d1, d2;

  int i = std::upper_bound(x_glob.begin(), x_glob.end(), p[0]) - x_glob.begin();
  d1 = (i==(int)x_glob.size()) ? fabs(x_glob[i-1] + dx_glob[i-1] - p[0]) : fabs(x_glob[i] - p[0]);
  d2 = (i==0)                  ? fabs(x_glob[i] - dx_glob[i] - p[0])     : fabs(x_glob[i-1] - p[0]);
  ijk[0] = d1<d2 ? i : i-1;

  if(!include_ghost_layer) {
    if(ijk[0]==-1)                      ijk[0]++;
    else if(ijk[0]==(int)x_glob.size()) ijk[0]--;
  }

  int j = std::upper_bound(y_glob.begin(), y_glob.end(), p[1]) - y_glob.begin();
  d1 = (j==(int)y_glob.size()) ? fabs(y_glob[j-1] + dy_glob[j-1] - p[1]) : fabs(y_glob[j] - p[1]);
  d2 = (j==0)                  ? fabs(y_glob[j] - dy_glob[j] - p[1])     : fabs(y_glob[j-1] - p[1]);
  ijk[1] = d1<d2 ? j : j-1;

  if(!include_ghost_layer) {
    if(ijk[1]==-1)                      ijk[1]++;
    else if(ijk[1]==(int)y_glob.size()) ijk[1]--;
  }

  int k = std::upper_bound(z_glob.begin(), z_glob.end(), p[2]) - z_glob.begin();
  d1 = (k==(int)z_glob.size()) ? fabs(z_glob[k-1] + dz_glob[k-1] - p[2]) : fabs(z_glob[k] - p[2]);
  d2 = (k==0)                  ? fabs(z_glob[k] - dz_glob[k] - p[2])     : fabs(z_glob[k-1] - p[2]);
  ijk[2] = d1<d2 ? k : k-1;

  if(!include_ghost_layer) {
    if(ijk[1]==-1)                      ijk[1]++;
    else if(ijk[1]==(int)y_glob.size()) ijk[1]--;
  }

  return ijk;
}

//------------------------------------------------------------------

bool
GlobalMeshInfo::FindCellCoveringPoint(Vec3D &p, Int3 &ijk, bool include_ghost_layer)
{
  if(!IsPointInDomain(p, include_ghost_layer))
    return false;

  double shift = 0.5;

  if(include_ghost_layer && p[0]<x_glob[0]-shift*dx_glob[0]) {
    ijk[0] = -1;  
  } else {
    bool found = false;
    for(int i=0; i<(int)x_glob.size(); i++) {
      if(p[0]<x_glob[i]+shift*dx_glob[i]) {
        ijk[0] = i;
        found = true;
        break;
      }
    }
    if(!found)
      ijk[0] = include_ghost_layer ? x_glob.size() : x_glob.size() - 1;
  }


  if(include_ghost_layer && p[1]<y_glob[0]-shift*dy_glob[0]) {
    ijk[1] = -1;  
  } else {
    bool found = false;
    for(int i=0; i<(int)y_glob.size(); i++) {
      if(p[1]<y_glob[i]+shift*dy_glob[i]) {
        ijk[1] = i;
        found = true;
        break;
      }
    }
    if(!found)
      ijk[1] = include_ghost_layer ? y_glob.size() : y_glob.size() - 1;
  }


  if(include_ghost_layer && p[2]<z_glob[0]-shift*dz_glob[0]) {
    ijk[2] = -1;  
  } else {
    bool found = false;
    for(int i=0; i<(int)z_glob.size(); i++) {
      if(p[2]<z_glob[i]+shift*dz_glob[i]) {
        ijk[2] = i;
        found = true;
        break;
      }
    }
    if(!found)
      ijk[2] = include_ghost_layer ? z_glob.size() : z_glob.size() - 1;
  }

  return true;
}

//------------------------------------------------------------------

bool
GlobalMeshInfo::FindElementCoveringPoint(Vec3D &p, Int3 &ijk0, Vec3D *xi, bool include_ghost_layer)
{
  if(!IsPointInNodalMesh(p, include_ghost_layer))
    return false;

  ijk0[0] = int(std::upper_bound(x_glob.begin(), x_glob.end(), p[0]) - x_glob.begin()) - 1;
  ijk0[1] = int(std::upper_bound(y_glob.begin(), y_glob.end(), p[1]) - y_glob.begin()) - 1;
  ijk0[2] = int(std::upper_bound(z_glob.begin(), z_glob.end(), p[2]) - z_glob.begin()) - 1;

  if(xi) {
    double s0 = GetX(ijk0[0]);
    (*xi)[0] = (p[0] - s0)/(GetX(ijk0[0]+1) - s0);
    s0 = GetY(ijk0[1]);
    (*xi)[1] = (p[1] - s0)/(GetY(ijk0[1]+1) - s0);
    s0 = GetZ(ijk0[2]);
    (*xi)[2] = (p[2] - s0)/(GetZ(ijk0[2]+1) - s0);
  }

  return true;
}

//------------------------------------------------------------------

int
GlobalMeshInfo::GetOwnerOfCell(int i, int j, int k, bool include_ghost_layer)
{

  assert(!subD_ijk_min.empty()); //otherwise, this vector has not been created!

  if(include_ghost_layer) { //"pull" the node into the domain interior
    if(i==-1)                  i = 0;
    if(i==(int)x_glob.size())  i = x_glob.size()-1;
    if(j==-1)                  j = 0;
    if(j==(int)y_glob.size())  j = y_glob.size()-1;
    if(k==-1)                  k = 0;
    if(k==(int)z_glob.size())  k = z_glob.size()-1;
  }

  for(int proc=0; proc<(int)subD_ijk_min.size(); proc++) {
    if(k >= subD_ijk_min[proc][2] && k < subD_ijk_max[proc][2] &&
       j >= subD_ijk_min[proc][1] && j < subD_ijk_max[proc][1] &&
       i >= subD_ijk_min[proc][0] && i < subD_ijk_max[proc][0])
      return proc; //got you!
  }

  return -1; //not found!

}

//------------------------------------------------------------------

int
GlobalMeshInfo::GetOwnerOfPoint(Vec3D &p, bool include_ghost_layer)
{
  Int3 ijk;
  bool found = FindCellCoveringPoint(p, ijk, include_ghost_layer);
  if(!found)
    return -1; //not found!

  return GetOwnerOfCell(ijk[0], ijk[1], ijk[2], include_ghost_layer);
}

//------------------------------------------------------------------

bool
GlobalMeshInfo::IsCellInSubdomain(int i, int j, int k, int sub, 
                                  bool include_ext_ghost_layer)
{
  assert(!subD_ijk_min.empty()); //otherwise, this vector has not been created!

  if(sub<0 || sub>=(int)subD_ijk_min.size()){
    fprintf(stdout,"\033[0;31m*** Error: Calling IsCellInSubdomain with incorrect "
                   "subdomain id (%d).\033[0m\n", sub);
    exit(-1);
  }
    
  if(include_ext_ghost_layer) {
    if(i==-1)                  i = 0;
    if(i==(int)x_glob.size())  i = x_glob.size()-1;
    if(j==-1)                  j = 0;
    if(j==(int)y_glob.size())  j = y_glob.size()-1;
    if(k==-1)                  k = 0;
    if(k==(int)z_glob.size())  k = z_glob.size()-1;
  }

  if(k >= subD_ijk_min[sub][2] && k < subD_ijk_max[sub][2] &&
     j >= subD_ijk_min[sub][1] && j < subD_ijk_max[sub][1] &&
     i >= subD_ijk_min[sub][0] && i < subD_ijk_max[sub][0])
    return true;

  return false;
}

//------------------------------------------------------------------

bool
GlobalMeshInfo::IsPointInSubdomain(Vec3D &p, int sub, bool include_ext_ghost_layer)
{
  Int3 ijk;
  bool found = FindCellCoveringPoint(p, ijk, include_ext_ghost_layer);
  if(!found)
    return false;

  return IsCellInSubdomain(ijk[0], ijk[1], ijk[2], sub, include_ext_ghost_layer);
}

//------------------------------------------------------------------

std::vector<int> &
GlobalMeshInfo::GetAllNeighborsOfSub(int sub) {
  assert(!subD_neighbors_all.empty()); //otherwise, it has not been created!
  assert(sub>=0 && sub<(int)subD_neighbors_all.size());
  return subD_neighbors_all[sub];
}

//------------------------------------------------------------------

std::vector<int> &
GlobalMeshInfo::GetFaceEdgeNeighborsOfSub(int sub) {
  assert(!subD_neighbors_face_edge.empty()); //otherwise, it has not been created!
  assert(sub>=0 && sub<(int)subD_neighbors_face_edge.size());
  return subD_neighbors_face_edge[sub];
}

//------------------------------------------------------------------

std::vector<int> &
GlobalMeshInfo::GetFaceNeighborsOfSub(int sub) {
  assert(!subD_neighbors_face.empty()); //otherwise, it has not been created!
  assert(sub>=0 && sub<(int)subD_neighbors_face.size());
  return subD_neighbors_face[sub];
}

//------------------------------------------------------------------

std::vector<int> &
GlobalMeshInfo::Get27NeighborhoodOfSub(int sub) {
  assert(!subD_neighbors_27.empty()); //otherwise, it has not been created!
  assert(sub>=0 && sub<(int)subD_neighbors_27.size());
  return subD_neighbors_27[sub];
}

//------------------------------------------------------------------

std::vector<int> &
GlobalMeshInfo::Get19NeighborhoodOfSub(int sub) {
  assert(!subD_neighbors_19.empty()); //otherwise, it has not been created!
  assert(sub>=0 && sub<(int)subD_neighbors_19.size());
  return subD_neighbors_19[sub];
}

//------------------------------------------------------------------

std::vector<int> &
GlobalMeshInfo::Get7NeighborhoodOfSub(int sub) {
  assert(!subD_neighbors_7.empty()); //otherwise, it has not been created!
  assert(sub>=0 && sub<(int)subD_neighbors_7.size());
  return subD_neighbors_7[sub];
}

//------------------------------------------------------------------


//------------------------------------------------------------------





