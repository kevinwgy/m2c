#include<GlobalMeshInfo.h>
#include<Vector3D.h>
#include<algorithm> //std::upper_bound

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
  d1 = (i==x_glob.size()) ? fabs(x_glob[i-1] + dx_glob[i-1] - p[0]) : fabs(x_glob[i] - p[0]);
  d2 = (i==0)             ? fabs(x_glob[i] - dx_glob[i] - p[0])     : fabs(x_glob[i-1] - p[0]);
  ijk[0] = d1<d2 ? i : i-1;

  if(!include_ghost_layer) {
    if(ijk[0]==-1)                 ijk[0]++;
    else if(ijk[0]==x_glob.size()) ijk[0]--;
  }

  int j = std::upper_bound(y_glob.begin(), y_glob.end(), p[1]) - y_glob.begin();
  d1 = (j==y_glob.size()) ? fabs(y_glob[j-1] + dy_glob[j-1] - p[1]) : fabs(y_glob[j] - p[1]);
  d2 = (j==0)             ? fabs(y_glob[j] - dy_glob[j] - p[1])     : fabs(y_glob[j-1] - p[1]);
  ijk[1] = d1<d2 ? j : j-1;

  if(!include_ghost_layer) {
    if(ijk[1]==-1)                 ijk[1]++;
    else if(ijk[1]==y_glob.size()) ijk[1]--;
  }

  int k = std::upper_bound(z_glob.begin(), z_glob.end(), p[2]) - z_glob.begin();
  d1 = (k==z_glob.size()) ? fabs(z_glob[k-1] + dz_glob[k-1] - p[2]) : fabs(z_glob[k] - p[2]);
  d2 = (k==0)             ? fabs(z_glob[k] - dz_glob[k] - p[2])     : fabs(z_glob[k-1] - p[2]);
  ijk[2] = d1<d2 ? k : k-1;

  if(!include_ghost_layer) {
    if(ijk[1]==-1)                 ijk[1]++;
    else if(ijk[1]==y_glob.size()) ijk[1]--;
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
    for(int i=0; i<x_glob.size(); i++) {
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
    for(int i=0; i<y_glob.size(); i++) {
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
    for(int i=0; i<z_glob.size(); i++) {
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
GlobalMeshInfo::FindElementCoveringPoint(Vec3D &p, Int3 &ijk0, bool include_ghost_layer)
{
  if(!IsPointInNodalMesh(p, include_ghost_layer))
    return false;

  ijk0[0] = int(std::upper_bound(x_glob.begin(), x_glob.end(), p[0]) - x_glob.begin()) - 1;
  ijk0[1] = int(std::upper_bound(y_glob.begin(), y_glob.end(), p[1]) - y_glob.begin()) - 1;
  ijk0[2] = int(std::upper_bound(z_glob.begin(), z_glob.end(), p[2]) - z_glob.begin()) - 1;

  return true;
}

//------------------------------------------------------------------


//------------------------------------------------------------------






