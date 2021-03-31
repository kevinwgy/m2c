#include <MeshGenerator.h>

//-----------------------------------------------------

void MeshGenerator::ComputeMeshCoordinatesAndDeltas(MeshData &iod_mesh,
                                                    vector<double> &x, vector<double> &y,
                                                    vector<double> &z, vector<double> &dx,
                                                    vector<double> &dy, vector<double> &dz)
{

  // read input info
  double x0 = iod_mesh.x0, xmax = iod_mesh.xmax;
  double y0 = iod_mesh.y0, ymax = iod_mesh.ymax;
  double z0 = iod_mesh.z0, zmax = iod_mesh.zmax;
  if(xmax<x0 || ymax<y0 || zmax<z0) {
    print_error("*** Error: Dimensions of the computational domain are incorrect. [%e,%e,%e]->[%e,%e,%e]\n",
                x0,y0,z0,xmax,ymax,zmax);
    exit_mpi();
  }

  int Nx = iod_mesh.Nx, Ny = iod_mesh.Ny, Nz = iod_mesh.Nz;

  vector<pair<double,double> > xpoints;
  for(auto it = iod_mesh.xpoints_map.dataMap.begin(); it != iod_mesh.xpoints_map.dataMap.end(); it++)
    xpoints.push_back(make_pair(it->second->coord, it->second->h));
  vector<pair<double,double> > ypoints;
  for(auto it = iod_mesh.ypoints_map.dataMap.begin(); it != iod_mesh.ypoints_map.dataMap.end(); it++)
    ypoints.push_back(make_pair(it->second->coord, it->second->h));
  vector<pair<double,double> > zpoints;
  for(auto it = iod_mesh.zpoints_map.dataMap.begin(); it != iod_mesh.zpoints_map.dataMap.end(); it++)
    zpoints.push_back(make_pair(it->second->coord, it->second->h));
  
  // check for error
  if ( (Nx<=0 && xpoints.size()==0) || (Ny<=0 && ypoints.size()==0) || (Nz<=0 && zpoints.size()==0) ) {
    print_error("*** Error: Unable to create mesh due to insufficient input information.\n");
    exit_mpi();
  }
  if ( (Nx>0 && xpoints.size()>0) || (Ny>0 && ypoints.size()>0) || (Nz>0 && zpoints.size()>0) ) {
    print_error("*** Error: Unable to create mesh due to conflicting (or redundant) input information.\n");
    exit_mpi();
  }

  // sort by coord
  sort(xpoints.begin(), xpoints.end());
  sort(ypoints.begin(), ypoints.end());
  sort(zpoints.begin(), zpoints.end());

  // insert min and max
  if(xpoints.size() != 0) {
    if(xpoints.front.first > x0)
      xpoints.insert(xpoints.front, make_pair(x0, xpoints.front.second));
    else if(xpoints.front.first < x0) {
      print_error("*** Error: Detected mesh control point (x-dir) outside physical domain.\n");
      exit_mpi();
    }
      
    if(xpoints.back.first < xmax)
      xpoints.push_back(make_pair(xmax, xpoints.back.second));
    else if (xpoints.back.first > xmax) {
      print_error("*** Error: Detected mesh control point (x-dir) outside physical domain.\n");
      exit_mpi();
    }
  }

  if(ypoints.size() != 0) {
    if(ypoints.front.first > y0)
      ypoints.insert(ypoints.front, make_pair(y0, ypoints.front.second));
    else if(ypoints.front.first < y0) {
      print_error("*** Error: Detected mesh control point (y-dir) outside physical domain.\n");
      exit_mpi();
    }

    if(ypoints.back.first < ymax)
      ypoints.push_back(make_pair(ymax, ypoints.back.second));
    else if(ypoints.back.first > ymax) { 
      print_error("*** Error: Detected mesh control point (y-dir) outside physical domain.\n");
      exit_mpi();
    }
  }

  if(zpoints.size() != 0) {
    if(zpoints.front.first > z0)
      zpoints.insert(zpoints.front, make_pair(z0, zpoints.front.second));
    else if(zpoints.front.first < z0) {
      print_error("*** Error: Detected mesh control point (z-dir) outside physical domain.\n");
      exit_mpi();
    }

    if(zpoints.back.first < zmax)
      zpoints.push_back(make_pair(zmax, zpoints.back.second));
    else if(zpoints.back.first > zmax) {
      print_error("*** Error: Detected mesh control point (z-dir) outside physical domain.\n");
      exit_mpi();
    }
  }

  // generate mesh in x-dir
  if(xpoints.size()==0)
    ComputeMesh1DUniform(x0, xmax, Nx, x, dx);
  else 
    ComputeMesh1DNonUniform(x0, xmax, xpoints, x, dx);

  // generate mesh in y-dir
  if(ypoints.size()==0)
    ComputeMesh1DUniform(y0, ymax, Ny, y, dy); 
  else
    ComputeMesh1DNonUniform(y0, ymax, ypoints, y, dy);

  // generate mesh in z-dir
  if(zpoints.size()==0)
    ComputeMesh1DUniform(z0, zmax, Nz, z, dz);
  else
    ComputeMesh1DNonUniform(z0, zmax, zpoints, z, dz);

}

//-----------------------------------------------------

void MeshGenerator::ComputeMesh1DUniform(double x0, double xmax, double Nx, 
                                         vector<dobule> &x, vector<double> &dx)
{
  x.clear();
  dx.clear();

  double h = (xmax-xmin)/Nx;
  dx.resize(Nx, h);

  x.resize(Nx);
  for(int i=0; i<Nx; i++)
    x[i] = x0 + 0.5*h + i*h;
}

//-----------------------------------------------------

// Ref. Alfio Quarteroni, Numerical Models for Differential Problems, Chapter 6
// Applying linear interpolation to get local cell size. (TODO: can be extended if needed)
void MeshGenerator::ComputeMesh1DNonUniform(double x0, double xmax, 
                                            vector<pair<double,double> > &xpoints,
                                            vector<double> &x, vector<double> &dx)
{
  //The cell size at x0 and xmax must be defined in xpoints. So the coords of the first & last
  //nodes (cell centers) are specified
  int Np = xpoints.size();
  x0   = x0 + 0.5*xpoints[0].second();
  xmax = xmax - 0.5*xpoints[Np-1].second();

  //Calculate the number nodes / cells (in the real domain)
  double h = 
  double Nreal = 0.0; //integration of 1/h from x0 to xmax
  for(int i=0; i<Np; i++)
    
    

}

//-----------------------------------------------------




