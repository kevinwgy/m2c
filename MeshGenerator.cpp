/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <MeshGenerator.h>
using std::pair;
using std::vector;

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
    xpoints.push_back(std::make_pair(it->second->coord, it->second->h));
  vector<pair<double,double> > ypoints;
  for(auto it = iod_mesh.ypoints_map.dataMap.begin(); it != iod_mesh.ypoints_map.dataMap.end(); it++)
    ypoints.push_back(std::make_pair(it->second->coord, it->second->h));
  vector<pair<double,double> > zpoints;
  for(auto it = iod_mesh.zpoints_map.dataMap.begin(); it != iod_mesh.zpoints_map.dataMap.end(); it++)
    zpoints.push_back(std::make_pair(it->second->coord, it->second->h));
  
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
    if(xpoints.front().first > x0)
      xpoints.insert(xpoints.begin(), std::make_pair(x0, xpoints.front().second));
    else if(xpoints.front().first < x0) {
      print_error("*** Error: Detected mesh control point (x-dir) outside physical domain.\n");
      exit_mpi();
    }
      
    if(xpoints.back().first < xmax)
      xpoints.push_back(std::make_pair(xmax, xpoints.back().second));
    else if (xpoints.back().first > xmax) {
      print_error("*** Error: Detected mesh control point (x-dir) outside physical domain.\n");
      exit_mpi();
    }
  }

  if(ypoints.size() != 0) {
    if(ypoints.front().first > y0)
      ypoints.insert(ypoints.begin(), std::make_pair(y0, ypoints.front().second));
    else if(ypoints.front().first < y0) {
      print_error("*** Error: Detected mesh control point (y-dir) outside physical domain.\n");
      exit_mpi();
    }

    if(ypoints.back().first < ymax)
      ypoints.push_back(std::make_pair(ymax, ypoints.back().second));
    else if(ypoints.back().first > ymax) { 
      print_error("*** Error: Detected mesh control point (y-dir) outside physical domain.\n");
      exit_mpi();
    }
  }

  if(zpoints.size() != 0) {
    if(zpoints.front().first > z0)
      zpoints.insert(zpoints.begin(), std::make_pair(z0, zpoints.front().second));
    else if(zpoints.front().first < z0) {
      print_error("*** Error: Detected mesh control point (z-dir) outside physical domain.\n");
      exit_mpi();
    }

    if(zpoints.back().first < zmax)
      zpoints.push_back(std::make_pair(zmax, zpoints.back().second));
    else if(zpoints.back().first > zmax) {
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


  // Print mesh statistics to the screen
  print("\n- Mesh Statistics:\n");
  print("  X-Direction: [%e, %e], %d nodes/cells, dx_min = %e, dx_max = %e.\n", x0, xmax, x.size(), 
         *std::min_element(dx.begin(), dx.end()), *std::max_element(dx.begin(), dx.end()));
  print("  Y-Direction: [%e, %e], %d nodes/cells, dy_min = %e, dy_max = %e.\n", y0, ymax, y.size(), 
         *std::min_element(dy.begin(), dy.end()), *std::max_element(dy.begin(), dy.end()));
  print("  Z-Direction: [%e, %e], %d nodes/cells, dz_min = %e, dz_max = %e.\n", z0, zmax, z.size(), 
         *std::min_element(dz.begin(), dz.end()), *std::max_element(dz.begin(), dz.end()));
  print("  Total number of nodes/cells: %d.\n", x.size()*y.size()*z.size());

  if(iod_mesh.type == MeshData::CYLINDRICAL) 
    print("  Imposing cylindrical symmetry: x ~ axial coordinate, y ~ radial coordinate.\n");
  else if(iod_mesh.type == MeshData::SPHERICAL)
    print("  Imposing spherical symmetry: x ~ radial coordinate.\n");

/*
  for(int i=0; i<x.size(); i++)
    print("x[%d] = %e.\n", i, x[i]);
*/
}

//-----------------------------------------------------

void MeshGenerator::ComputeMesh1DUniform(double x0, double xmax, double Nx, 
                                         vector<double> &x, vector<double> &dx)
{
  x.clear();
  dx.clear();

  double h = (xmax-x0)/Nx;
  dx.resize(Nx, h);

  x.resize(Nx);
  for(int i=0; i<Nx; i++)
    x[i] = x0 + 0.5*h + i*h;
}

//-----------------------------------------------------

// Ref. (1) KW's notes. (2) Alfio Quarteroni, Numerical Models for Differential Problems, Chapter 6
// Applying linear interpolation to get local cell size. (TODO: can be extended if needed)
void MeshGenerator::ComputeMesh1DNonUniform(double x0, double xmax, 
                                            vector<pair<double,double> > &xpoints,
                                            vector<double> &x, vector<double> &dx)
{
  //The cell size at x0 and xmax must be defined in xpoints. So the coords of the first & last
  //nodes (cell centers) are specified
  int Np = xpoints.size();

  //Calculate the number nodes / cells (in the real domain)
  double Nreal = 0.0;
  double a, b, xp1, xp2, h1, h2;
  for(int i=0; i<Np-1; i++) {
    xp1 = xpoints[i].first;
    h1  = xpoints[i].second;
    xp2 = xpoints[i+1].first;
    h2  = xpoints[i+1].second;
    a   = (h2 - h1)/(xp2 - xp1);
    b   = (h1*xp2 - h2*xp1)/(xp2 - xp1);
    Nreal += (a==0.0) ? (xp2-xp1)/b : 1.0/a*log((a*xp2+b)/(a*xp1+b));
  }  
  int N = std::max(1, (int)(Nreal+1e-12)); //(int)Nreal --> integer floor of Nreal (since Nreal>0)


  //Calculating the coordinates of the cell boundaries
  double kappa = N/Nreal;
  double xi[N+1]; //these are the cell boundaries

  xi[0] = x0;

  int i = 0;
  xp1 = xpoints[i].first;
  h1  = xpoints[i].second;
  xp2 = xpoints[i+1].first;
  h2  = xpoints[i+1].second;
  a   = (h2 - h1)/(xp2 - xp1);
  b   = (h1*xp2 - h2*xp1)/(xp2 - xp1);

  for(int k=0; k<N; k++) {

    xi[k+1] = (a==0.0) ? xi[k] + b/kappa : (xi[k] + b/a)*exp(a/kappa) - b/a;

    if(xi[k+1] >= xp2) {
      if(k+1==N && i+1 == (int)xpoints.size()-1) {//done
        xi[k+1] = xp2;
        break;
      } else if(k+1==N || i+1 == (int)xpoints.size()-1) {
        print_error("*** Error: Cannot generate mesh (possibly a software bug)\n");
        exit_mpi();
      }

      double C = (a==0.0) ? (xp2-xi[k])/b : 1.0/a*log((a*xp2+b)/(a*xi[k]+b));

      i++;
      xp1 = xp2; 
      h1  = h2; 
      xp2 = xpoints[i+1].first;
      h2  = xpoints[i+1].second;
      a   = (h2 - h1)/(xp2 - xp1);
      b   = (h1*xp2 - h2*xp1)/(xp2 - xp1);

      xi[k+1] = (a==0.0) ? xp1 + b*(1/kappa - C) : (xp1 + b/a)*exp(a*(1/kappa - C)) - b/a;

      if(xi[k+1]>xp2) {
        print_error("*** Error: Detected conflicts in mesh control points and cell widths. Cannot generate mesh.\n");
        exit_mpi();
      }
    }

    //fprintf(stdout,"xi[%d] = %e\n", k+1, xi[k+1]);

  }
  
  if(xi[N]<xmax)
    xi[N] = xmax;


  // Now calculate x and dx
  x.resize(N);
  dx.resize(N);

  for(int i=0; i<N; i++)
    dx[i] = xi[i+1]-xi[i];

  x[0] = x0 + 0.5*dx[0];
  for(int i=0; i<N-1; i++)
    x[i+1] = x[i] + 0.5*(dx[i]+dx[i+1]);
}

//-----------------------------------------------------




