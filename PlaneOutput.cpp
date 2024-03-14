/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<PlaneOutput.h>
#include<Intersections.h>
#include<trilinear_interpolation.h>
#include<numeric> //std::accumulate

using std::vector;
using std::array;
using std::pair;
using std::string;

//-------------------------------------------------------------------------

int WORKER_PROC = 0;

//-------------------------------------------------------------------------

PlaneOutput::PlaneOutput(MPI_Comm &comm_, OutputData &iod_output_, PlanePlot &iod_pplot_, 
                         std::vector<VarFcnBase*> &vf_, GlobalMeshInfo &global_mesh_,
                         IonizationOperator *ion_)
           : comm(comm_), iod_pplot(iod_pplot_), vf(vf_), global_mesh(global_mesh_),
             ion(ion_)
{

  iFrame = 0;

  frequency = iod_pplot.frequency;
  frequency_dt = iod_pplot.frequency_dt;

  last_snapshot_time = -1.0;

  Vec3D p(iod_pplot.x0, iod_pplot.y0, iod_pplot.z0);
  Vec3D n(iod_pplot.normal_x, iod_pplot.normal_y, iod_pplot.normal_z);
  if(n.norm()==0) {
    print_error("*** Error: Detected invalid normal for user-specified cut-plane (norm = 0).\n");
    exit_mpi();
  }
  n /= n.norm();
  print("\n- [Cut-Plane]: Coords = (%e, %e, %e), Normal = (%e, %e, %e).\n", p[0], p[1], p[2],
        n[0], n[1], n[2]);


  if(iod_pplot.mesh[0] == 0) {
    print_error("*** Error: Cut-plane output: Mesh file name is not provided.\n");
    exit(-1);
  } 

  // ----------------------------------------------
  // Get proc number
  // ----------------------------------------------
  int mpi_rank(-1);
  MPI_Comm_rank(comm, &mpi_rank);


  // ----------------------------------------------
  // set up files (Similar to ProbeOutput constructor)
  // ----------------------------------------------
  filename_mesh = string(iod_output_.prefix) + string(iod_pplot.mesh);

  if (iod_pplot.density[0] != 0) {
    filename[PlanePlot::DENSITY] = string(iod_output_.prefix) + string(iod_pplot.density);
    if(mpi_rank == WORKER_PROC) {
      FILE *file = fopen(filename[PlanePlot::DENSITY].c_str(), "w");
      assert(file); //if file is not opened, the pointer would be NULL
      fprintf(file, "Scalar Density under load for SurfaceNodes\n");
      fclose(file);
    }
  }

  if (iod_pplot.velocity[0] != 0) {
    filename[PlanePlot::VELOCITY] = string(iod_output_.prefix) + string(iod_pplot.velocity);
    if(mpi_rank == WORKER_PROC) {
      FILE *file = fopen(filename[PlanePlot::VELOCITY].c_str(), "w");
      assert(file); //if file is not opened, the pointer would be NULL
      fprintf(file, "Vector Velocity under load for SurfaceNodes\n");
      fclose(file);
    }
  }

  if (iod_pplot.pressure[0] != 0) {
    filename[PlanePlot::PRESSURE] = string(iod_output_.prefix) + string(iod_pplot.pressure);
    if(mpi_rank == WORKER_PROC) {
      FILE *file = fopen(filename[PlanePlot::PRESSURE].c_str(), "w");
      assert(file); //if file is not opened, the pointer would be NULL
      fprintf(file, "Scalar Pressure under load for SurfaceNodes\n");
      fclose(file);
    }
  }

  if (iod_pplot.temperature[0] != 0) {
    filename[PlanePlot::TEMPERATURE] = string(iod_output_.prefix) + string(iod_pplot.temperature);
    if(mpi_rank == WORKER_PROC) {
      FILE *file = fopen(filename[PlanePlot::TEMPERATURE].c_str(), "w");
      assert(file); //if file is not opened, the pointer would be NULL
      fprintf(file, "Scalar Temperature under load for SurfaceNodes\n");
      fclose(file);
    }
  }

  if (iod_pplot.delta_temperature[0] != 0) {
    filename[PlanePlot::DELTA_TEMPERATURE] = string(iod_output_.prefix) + string(iod_pplot.delta_temperature);
    if(mpi_rank == WORKER_PROC) {
      FILE *file = fopen(filename[PlanePlot::DELTA_TEMPERATURE].c_str(), "w");
      assert(file); //if file is not opened, the pointer would be NULL
      fprintf(file, "Scalar DeltaTemperature under load for SurfaceNodes\n");
      fclose(file);
    }
  }

  if (iod_pplot.materialid[0] != 0) {
    filename[PlanePlot::MATERIALID] = string(iod_output_.prefix) + string(iod_pplot.materialid);
    if(mpi_rank == WORKER_PROC) {
      FILE *file = fopen(filename[PlanePlot::MATERIALID].c_str(), "w");
      assert(file); //if file is not opened, the pointer would be NULL
      fprintf(file, "Scalar MaterialID under load for SurfaceNodes\n");
      fclose(file);
    }
  }

  if (iod_pplot.laser_radiance[0] != 0) {
    filename[PlanePlot::LASERRADIANCE] = string(iod_output_.prefix) + string(iod_pplot.laser_radiance);
    if(mpi_rank == WORKER_PROC) {
      FILE *file = fopen(filename[PlanePlot::LASERRADIANCE].c_str(), "w");
      assert(file); //if file is not opened, the pointer would be NULL
      fprintf(file, "Scalar LaserIrradiance under load for SurfaceNodes\n");
      fclose(file);
    }
  }

  if (iod_pplot.levelset0[0] != 0) {
    filename[PlanePlot::LEVELSET0] = string(iod_output_.prefix) + string(iod_pplot.levelset0);
    if(mpi_rank == WORKER_PROC) {
      FILE *file = fopen(filename[PlanePlot::LEVELSET0].c_str(), "w");
      assert(file); //if file is not opened, the pointer would be NULL
      fprintf(file, "Scalar LevelSet0 under load for SurfaceNodes\n");
      fclose(file);
    }
  }

  if (iod_pplot.levelset1[0] != 0) {
    filename[PlanePlot::LEVELSET1] = string(iod_output_.prefix) + string(iod_pplot.levelset1);
    if(mpi_rank == WORKER_PROC) {
      FILE *file = fopen(filename[PlanePlot::LEVELSET1].c_str(), "w");
      assert(file); //if file is not opened, the pointer would be NULL
      fprintf(file, "Scalar LevelSet1 under load for SurfaceNodes\n");
      fclose(file);
    }
  }

  if (iod_pplot.levelset2[0] != 0) {
    filename[PlanePlot::LEVELSET2] = string(iod_output_.prefix) + string(iod_pplot.levelset2);
    if(mpi_rank == WORKER_PROC) {
      FILE *file = fopen(filename[PlanePlot::LEVELSET2].c_str(), "w");
      assert(file); //if file is not opened, the pointer would be NULL
      fprintf(file, "Scalar LevelSet2 under load for SurfaceNodes\n");
      fclose(file);
    }
  }

  if (iod_pplot.levelset3[0] != 0) {
    filename[PlanePlot::LEVELSET3] = string(iod_output_.prefix) + string(iod_pplot.levelset3);
    if(mpi_rank == WORKER_PROC) {
      FILE *file = fopen(filename[PlanePlot::LEVELSET3].c_str(), "w");
      assert(file); //if file is not opened, the pointer would be NULL
      fprintf(file, "Scalar LevelSet3 under load for SurfaceNodes\n");
      fclose(file);
    }
  }

  if (iod_pplot.levelset4[0] != 0) {
    filename[PlanePlot::LEVELSET4] = string(iod_output_.prefix) + string(iod_pplot.levelset4);
    if(mpi_rank == WORKER_PROC) {
      FILE *file = fopen(filename[PlanePlot::LEVELSET4].c_str(), "w");
      assert(file); //if file is not opened, the pointer would be NULL
      fprintf(file, "Scalar LevelSet4 under load for SurfaceNodes\n");
      fclose(file);
    }
  }

  if (iod_pplot.ionization_result[0] != 0) {
    if(!ion) {
      print_error("*** Error: Requested ionization result, "
                  "without specifying an ionization model.\n");
      exit_mpi();
    }
    filename[PlanePlot::IONIZATION] = string(iod_output_.prefix) + string(iod_pplot.ionization_result);
    if(mpi_rank == WORKER_PROC) {
      FILE *file = fopen(filename[PlanePlot::IONIZATION].c_str(), "w");
      assert(file); //if file is not opened, the pointer would be NULL
      fprintf(file, "Scalar IonizationResult under load for SurfaceNodes\n");
      fclose(file);
    }
  }

}

//-------------------------------------------------------------------------

PlaneOutput::~PlaneOutput()
{
}

//-------------------------------------------------------------------------

void
PlaneOutput::InitializeOutput(SpaceVariable3D &coordinates)
{

  points.clear();
  ijk.clear();
  ijk_valid.clear();
  trilinear_coords.clear();
  sol1_global.clear();
  sol3_global.clear();

  Vec3D pop(iod_pplot.x0, iod_pplot.y0, iod_pplot.z0);
  Vec3D normal(iod_pplot.normal_x, iod_pplot.normal_y, iod_pplot.normal_z);

  int i0, j0, k0, imax, jmax, kmax;
  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);

  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();

  // loop through cells; cut each cell by the plane
  Vec3D Vmin, Vmax, dxyz_over_two;
  vector<Vec3D> xlocal(6);
  vector<array<int,4> > elems;
  vector<Vec3D> points_dup;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        dxyz_over_two = 0.5*global_mesh.GetDXYZ(Int3(i,j,k));
        Vmin = coords[k][j][i] - dxyz_over_two;
        Vmax = coords[k][j][i] + dxyz_over_two;
        int np = GeoTools::PlaneCuttingAxisAlignedBox(pop,normal,Vmin,Vmax,&xlocal,0);

        if(np<3)
          continue;

        //got at least a triangle

        int n0 = points_dup.size();
        points_dup.insert(points_dup.end(), xlocal.begin(), xlocal.begin()+np);

        assert(np<=6);
        switch(np) {
          case 3:
            elems.push_back(array<int,4>{n0,n0+1,n0+2,n0+2});
            break;
          case 4: 
            elems.push_back(array<int,4>{n0,n0+1,n0+2,n0+3});
            break;
          case 5: 
            elems.push_back(array<int,4>{n0,n0+1,n0+2,n0+3});
            elems.push_back(array<int,4>{n0,n0+3,n0+4,n0+4});
            break;
          case 6:
            elems.push_back(array<int,4>{n0,n0+1,n0+2,n0+3});
            elems.push_back(array<int,4>{n0,n0+3,n0+4,n0+5});
            break;
        }  
     }


  // eliminate duplicated nodes (within each subdomain)
  vector<int> nodemap;
  if(!points_dup.empty()) {
    nodemap.resize(points_dup.size());
    for(int n=0; n<(int)points_dup.size(); n++) {
      bool found = false;
      for(int i=0; i<(int)points.size(); i++) {
        if((points_dup[n]-points[i]).norm()<1e-12) {
          nodemap[n] = i;
          found = true;
          break;
        }
      }
      if(!found) {
        points.push_back(points_dup[n]);
        nodemap[n] = points.size()-1;
      }
    }
  }

  for(auto&& elem : elems)
    for(int n=0; n<4; n++)
      elem[n] = nodemap[elem[n]];


  // get interpolation info
  Int3 ijk0;
  Vec3D xi0;
  for(int n=0; n<(int)points.size(); n++) {

    bool found = global_mesh.FindElementCoveringPoint(points[n], ijk0, &xi0, true);
    assert(found);
    assert(coordinates.IsHere(ijk0[0],ijk0[1],ijk0[2],true) &&
           coordinates.IsHere(ijk0[0]+1,ijk0[1]+1,ijk0[2]+1,true));

    ijk.push_back(ijk0);
    trilinear_coords.push_back(xi0);

    ijk_valid.push_back(std::make_pair(8, array<bool,8>{true,true,true,true,
                                                        true,true,true,true}));
    int ii(ijk0[0]), jj(ijk0[1]), kk(ijk0[2]);
    if(coordinates.OutsidePhysicalDomainAndUnpopulated(ii,jj,kk)) {//c000
      ijk_valid.back().first--;
      ijk_valid.back().second[0] = false;
    }
    if(coordinates.OutsidePhysicalDomainAndUnpopulated(ii+1,jj,kk)) {//c100
      ijk_valid.back().first--;
      ijk_valid.back().second[1] = false;
    }
    if(coordinates.OutsidePhysicalDomainAndUnpopulated(ii,jj+1,kk)) {//c010
      ijk_valid.back().first--;
      ijk_valid.back().second[2] = false;
    }
    if(coordinates.OutsidePhysicalDomainAndUnpopulated(ii+1,jj+1,kk)) {//c110
      ijk_valid.back().first--;
      ijk_valid.back().second[3] = false;
    }
    if(coordinates.OutsidePhysicalDomainAndUnpopulated(ii,jj,kk+1)) {//c001
      ijk_valid.back().first--;
      ijk_valid.back().second[4] = false;
    }
    if(coordinates.OutsidePhysicalDomainAndUnpopulated(ii+1,jj,kk+1)) {//c101
      ijk_valid.back().first--;
      ijk_valid.back().second[5] = false;
    }
    if(coordinates.OutsidePhysicalDomainAndUnpopulated(ii,jj+1,kk+1)) {//c011
      ijk_valid.back().first--;
      ijk_valid.back().second[6] = false;
    }
    if(coordinates.OutsidePhysicalDomainAndUnpopulated(ii+1,jj+1,kk+1)) {//c111
      ijk_valid.back().first--;
      ijk_valid.back().second[7] = false;
    } 
    assert(ijk_valid.back().first>0);

  }


  coordinates.RestoreDataPointerToLocalVector();

  int nNodes = points.size();
  int nElems = elems.size();

  sol1.resize(nNodes);
  sol3.resize(nNodes);

  // ---------------
  // MPI comms
  // ---------------
  int WORKER_PROC = 0;
  int mpi_rank(-1), mpi_size(-1);
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  vector<int> nNodes_all;
  vector<int> nElems_all;
  if(mpi_rank==WORKER_PROC) { //lead proc

    nNodes_all.resize(mpi_size);
    nElems_all.resize(mpi_size);
    MPI_Gather(&nNodes, 1, MPI_INT, nNodes_all.data(), 1, MPI_INT, WORKER_PROC, comm);
    MPI_Gather(&nElems, 1, MPI_INT, nElems_all.data(), 1, MPI_INT, WORKER_PROC, comm);

    int totalNodes = std::accumulate(nNodes_all.begin(), nNodes_all.end(), 0);
    sol1_global.resize(totalNodes);
    sol3_global.resize(totalNodes);

    
    package_disp1.resize(mpi_size);
    package_size1.resize(mpi_size);
    package_disp3.resize(mpi_size);
    package_size3.resize(mpi_size);
    package_disp1[0] = 0;
    package_size1[0] = nNodes_all[0];
    package_disp3[0] = 0;
    package_size3[0] = 3*nNodes_all[0];
    for(int proc=1; proc<mpi_size; proc++) {
      package_disp1[proc] = package_disp1[proc-1] + package_size1[proc-1];
      package_size1[proc] = nNodes_all[proc];
      package_disp3[proc] = package_disp3[proc-1] + package_size3[proc-1];
      package_size3[proc] = 3*nNodes_all[proc];
    }

  } else {
    MPI_Gather(&nNodes, 1, MPI_INT, NULL, 0, MPI_INT, WORKER_PROC, comm);
    MPI_Gather(&nElems, 1, MPI_INT, NULL, 0, MPI_INT, WORKER_PROC, comm);
  }

  // WORKER_PROC gathers all the nodes
  MPI_Gatherv((double*)points.data(), 3*points.size(), MPI_DOUBLE,
              (double*)sol3_global.data(), package_size3.data(),
              package_disp3.data(), MPI_DOUBLE, WORKER_PROC, comm);


  // WORKER_PROC gathers all the elements
  vector<int> elems_size(mpi_size);
  vector<int> elems_disp(mpi_size);
  vector<array<int,4> > elems_global;
  if(mpi_rank==WORKER_PROC) {
    elems_disp[0] = 0;
    elems_size[0] = 4*nElems_all[0];
    for(int proc=1; proc<mpi_size; proc++) {
      elems_disp[proc] = elems_disp[proc-1] + elems_size[proc-1];
      elems_size[proc] = 4*nElems_all[proc]; 
    }
    int totalElems = std::accumulate(nElems_all.begin(), nElems_all.end(), 0);
    elems_global.resize(totalElems);
  }
  MPI_Gatherv((int*)elems.data(), 4*elems.size(), MPI_INT,
              (int*)elems_global.data(), elems_size.data(),
              elems_disp.data(), MPI_INT, WORKER_PROC, comm);

  // WORKER_PROC updates node ids in elements
  if(mpi_rank==WORKER_PROC) {
    for(int proc=1; proc<mpi_size; proc++) {
      int n0 = package_disp1[proc];
      int first_elem = elems_disp[proc]/4;
      for(int e=0; e<nElems_all[proc]; e++)
        for(int n=0; n<4; n++)
          elems_global[first_elem+e][n] += n0;
    }
  }

  // WORKER_PROC prints mesh to file 
  if(mpi_rank==WORKER_PROC) {

    FILE *file = fopen(filename_mesh.c_str(), "w");
    assert(file); //if file is not opened, the pointer would be NULL
    fprintf(file, "Nodes CutPlaneNodes\n");
    for(int i=0; i<(int)sol3_global.size(); i++)
      fprintf(file, "%10d    %16.8e    %16.8e    %16.8e\n",
              i+1, sol3_global[i][0], sol3_global[i][1], sol3_global[i][2]);
    fprintf(file, "Elements CutPlane using CutPlaneNodes\n");
    for(int i=0; i<(int)elems_global.size(); i++) {
      if(elems_global[i][2]==elems_global[i][3]) //triangle
        fprintf(file, "%10d    4    %10d    %10d    %10d\n", i+1,
                elems_global[i][0]+1, elems_global[i][1]+1, elems_global[i][2]+1);
      else //quad
        fprintf(file, "%10d    2    %10d    %10d    %10d    %10d\n", i+1,
                elems_global[i][0]+1, elems_global[i][1]+1, elems_global[i][2]+1,
                elems_global[i][3]+1);
    }
    fclose(file);

  }

  // WORKER_PROC writes number of nodes to solution file(s)
  if(mpi_rank==WORKER_PROC) {

    for(int i=0; i<PlanePlot::SIZE; i++)
      if(!filename[i].empty()) {
        FILE *file = fopen(filename[i].c_str(), "a");
        assert(file); //if file is not opened, the pointer would be NULL
        fprintf(file,"%ld\n", sol1_global.size());
        fclose(file);
      }
  }
  
}

//-------------------------------------------------------------------------

void
PlaneOutput::WriteSolutionOnPlane(double time, double dt, int time_step, SpaceVariable3D &V,
                                  SpaceVariable3D &ID, std::vector<SpaceVariable3D*> &Phi,
                                  SpaceVariable3D* L, bool force_write)
{

  if(!isTimeToWrite(time,dt,time_step,frequency_dt,frequency,last_snapshot_time,force_write))
    return;

  double*** v  = V.GetDataPointer(); //note: dim = 5

  if(!filename[PlanePlot::DENSITY].empty()) 
    GetScalarSolutionAndWrite(time, v, 5, 0, filename[PlanePlot::DENSITY]);

  if(!filename[PlanePlot::VELOCITY].empty()) 
    GetVec3SolutionAndWrite(time, v, 5, 123, filename[PlanePlot::VELOCITY]);

  if(!filename[PlanePlot::PRESSURE].empty()) 
    GetScalarSolutionAndWrite(time, v, 5, 4, filename[PlanePlot::PRESSURE]);

  if(!filename[PlanePlot::TEMPERATURE].empty()) {
    double*** id = ID.GetDataPointer();
    GetTemperatureAndWrite(time, v, id, filename[PlanePlot::TEMPERATURE]);
    ID.RestoreDataPointerToLocalVector();
  }

  if(!filename[PlanePlot::DELTA_TEMPERATURE].empty()) {
    double*** id = ID.GetDataPointer();
    GetDeltaTemperatureAndWrite(time, v, id, filename[PlanePlot::DELTA_TEMPERATURE]);
    ID.RestoreDataPointerToLocalVector();
  }

  if(!filename[PlanePlot::MATERIALID].empty()) {
    double*** id = ID.GetDataPointer();
    GetScalarSolutionAndWrite(time, id, 1, 0, filename[PlanePlot::MATERIALID]);
    ID.RestoreDataPointerToLocalVector();
  }

  if(!filename[PlanePlot::LASERRADIANCE].empty()) {
    if(L == NULL) {
      print_error("*** Error: Requested laser radiance output, but laser source is not specified.\n");
      exit_mpi();
    }
    double*** l  = (double***)L->GetDataPointer();
    GetScalarSolutionAndWrite(time, l, 1, 0, filename[PlanePlot::LASERRADIANCE]);
    L->RestoreDataPointerToLocalVector();
  }

  if(!filename[PlanePlot::LEVELSET0].empty() && Phi.size()>=1) {
    double*** phi = (double***)Phi[0]->GetDataPointer();
    GetScalarSolutionAndWrite(time, phi, 1, 0, filename[PlanePlot::LEVELSET0]);
    Phi[0]->RestoreDataPointerToLocalVector();
  }

  if(!filename[PlanePlot::LEVELSET1].empty() && Phi.size()>=2) {
    double*** phi = (double***)Phi[1]->GetDataPointer();
    GetScalarSolutionAndWrite(time, phi, 1, 0, filename[PlanePlot::LEVELSET1]);
    Phi[1]->RestoreDataPointerToLocalVector();
  }

  if(!filename[PlanePlot::LEVELSET2].empty() && Phi.size()>=3) {
    double*** phi = (double***)Phi[2]->GetDataPointer();
    GetScalarSolutionAndWrite(time, phi, 1, 0, filename[PlanePlot::LEVELSET2]);
    Phi[2]->RestoreDataPointerToLocalVector();
  }

  if(!filename[PlanePlot::LEVELSET3].empty() && Phi.size()>=4) {
    double*** phi = (double***)Phi[3]->GetDataPointer();
    GetScalarSolutionAndWrite(time, phi, 1, 0, filename[PlanePlot::LEVELSET3]);
    Phi[3]->RestoreDataPointerToLocalVector();
  }

  if(!filename[PlanePlot::LEVELSET4].empty() && Phi.size()>=5) {
    double*** phi = (double***)Phi[4]->GetDataPointer();
    GetScalarSolutionAndWrite(time, phi, 1, 0, filename[PlanePlot::LEVELSET4]);
    Phi[4]->RestoreDataPointerToLocalVector();
  }

 if(!filename[PlanePlot::IONIZATION].empty()) {
    if(ion == NULL) {
      print_error("*** Error: Output of ionization result is requested, but it is not computed.\n");
      exit_mpi();
    }
    double*** id = ID.GetDataPointer();
    GetIonizationAndWrite(time, v, id, filename[PlanePlot::IONIZATION]);
    ID.RestoreDataPointerToLocalVector();
  }

  V.RestoreDataPointerToLocalVector();

}

//-------------------------------------------------------------------------

void
PlaneOutput::GetScalarSolutionAndWrite(double time, double*** v, int dim, int p, string& fname)
{

  // interpolate solution
  int i,j,k;
  double c[8]; //c000,c100,c010,c110,c001,c101,c011,c111;
  for(int node=0; node<(int)points.size(); node++) {
    i = ijk[node][0];
    j = ijk[node][1];
    k = ijk[node][2];
    c[0] = ijk_valid[node].second[0] ? v[k][j][i*dim+p]         : 0.0;
    c[1] = ijk_valid[node].second[1] ? v[k][j][(i+1)*dim+p]     : 0.0;
    c[2] = ijk_valid[node].second[2] ? v[k][j+1][i*dim+p]       : 0.0;
    c[3] = ijk_valid[node].second[3] ? v[k][j+1][(i+1)*dim+p]   : 0.0;
    c[4] = ijk_valid[node].second[4] ? v[k+1][j][i*dim+p]       : 0.0;
    c[5] = ijk_valid[node].second[5] ? v[k+1][j][(i+1)*dim+p]   : 0.0;
    c[6] = ijk_valid[node].second[6] ? v[k+1][j+1][i*dim+p]     : 0.0;
    c[7] = ijk_valid[node].second[7] ? v[k+1][j+1][(i+1)*dim+p] : 0.0;

    if(ijk_valid[node].first<8) {//fill invalid slots with average value
      double c_avg = std::accumulate(c,c+8,0)/ijk_valid[node].first;
      for(int s=0; s<8; s++)
        if(!ijk_valid[node].second[s])  c[s] = c_avg;
    }

    sol1[node] = MathTools::trilinear_interpolation(trilinear_coords[node],c[0],c[1],c[2],
                                                    c[3],c[4],c[5],c[6],c[7]);
  }
  
  //gathers solution
  MPI_Gatherv(sol1.data(), sol1.size(), MPI_DOUBLE, sol1_global.data(), package_size1.data(),
              package_disp1.data(), MPI_DOUBLE, WORKER_PROC, comm);
  
  //writes to file
  WriteScalarSolutionToFile(time, sol1_global, fname);

}

//-------------------------------------------------------------------------

void
PlaneOutput::GetVec3SolutionAndWrite(double time, double*** v, int dim, int p, string& fname)
{

  assert(p==123); //other options can be added later (easily)

  // interpolate solution
  int i,j,k;
  Vec3D c[8]; //c000,c100,c010,c110,c001,c101,c011,c111;
  for(int node=0; node<(int)points.size(); node++) {
    i = ijk[node][0];
    j = ijk[node][1];
    k = ijk[node][2];
    for(int q=1; q<4; q++) {
      c[0][q-1] = ijk_valid[node].second[0] ? v[k][j][i*dim+q]         : 0.0;
      c[1][q-1] = ijk_valid[node].second[1] ? v[k][j][(i+1)*dim+q]     : 0.0;
      c[2][q-1] = ijk_valid[node].second[2] ? v[k][j+1][i*dim+q]       : 0.0;
      c[3][q-1] = ijk_valid[node].second[3] ? v[k][j+1][(i+1)*dim+q]   : 0.0;
      c[4][q-1] = ijk_valid[node].second[4] ? v[k+1][j][i*dim+q]       : 0.0;
      c[5][q-1] = ijk_valid[node].second[5] ? v[k+1][j][(i+1)*dim+q]   : 0.0;
      c[6][q-1] = ijk_valid[node].second[6] ? v[k+1][j+1][i*dim+q]     : 0.0;
      c[7][q-1] = ijk_valid[node].second[7] ? v[k+1][j+1][(i+1)*dim+q] : 0.0;
    }
    if(ijk_valid[node].first<8) {//fill invalid slots with average value
      Vec3D c_avg = std::accumulate(c,c+8,Vec3D(0.0))/(double)ijk_valid[node].first;
      for(int s=0; s<8; s++)
        if(!ijk_valid[node].second[s])  c[s] = c_avg;
    }

    sol3[node] = MathTools::trilinear_interpolation(trilinear_coords[node],c[0],c[1],c[2],
                                                    c[3],c[4],c[5],c[6],c[7]);
  }
  
  int mpi_rank(-1);
  MPI_Comm_rank(comm, &mpi_rank);

  //gathers solution
  MPI_Gatherv((double*)sol3.data(), 3*sol3.size(), MPI_DOUBLE, (double*)sol3_global.data(),
              package_size3.data(), package_disp3.data(), MPI_DOUBLE, WORKER_PROC, comm);
  
  //writes to file
  WriteVec3SolutionToFile(time, sol3_global, fname);

}

//-------------------------------------------------------------------------

void
PlaneOutput::GetTemperatureAndWrite(double time, double*** v, double*** id, string& fname)
{

  // interpolate solution
  int i,j,k;
  double c[8]; //c000,c100,c010,c110,c001,c101,c011,c111;
  double rho, p, e;
  int dim = 5, myid;
  for(int node=0; node<(int)points.size(); node++) {

    i = ijk[node][0];
    j = ijk[node][1];
    k = ijk[node][2];
 
    for(int s=0; s<8; s++)
      c[s] = 0.0;

    // c000
    if(ijk_valid[node].second[0]) {
      myid = id[k][j][i];
      rho  =  v[k][j][i*dim];
      p    =  v[k][j][i*dim+4];
      e    = vf[myid]->GetInternalEnergyPerUnitMass(rho,p);
      c[0] = vf[myid]->GetTemperature(rho,e);
    }

    // c100
    if(ijk_valid[node].second[1]) {
      myid = id[k][j][i+1];
      rho  =  v[k][j][(i+1)*dim];
      p    =  v[k][j][(i+1)*dim+4];
      e    = vf[myid]->GetInternalEnergyPerUnitMass(rho,p);
      c[1] = vf[myid]->GetTemperature(rho,e);
    }

    // c010
    if(ijk_valid[node].second[2]) {
      myid = id[k][j+1][i];
      rho  =  v[k][j+1][i*dim];
      p    =  v[k][j+1][i*dim+4];
      e    = vf[myid]->GetInternalEnergyPerUnitMass(rho,p);
      c[2] = vf[myid]->GetTemperature(rho,e);
    }

    // c110
    if(ijk_valid[node].second[3]) {
      myid = id[k][j+1][i+1];
      rho  =  v[k][j+1][(i+1)*dim];
      p    =  v[k][j+1][(i+1)*dim+4];
      e    = vf[myid]->GetInternalEnergyPerUnitMass(rho,p);
      c[3] = vf[myid]->GetTemperature(rho,e);
    }

    // c001
    if(ijk_valid[node].second[4]) {
      myid = id[k+1][j][i];
      rho  =  v[k+1][j][i*dim];
      p    =  v[k+1][j][i*dim+4];
      e    = vf[myid]->GetInternalEnergyPerUnitMass(rho,p);
      c[4] = vf[myid]->GetTemperature(rho,e);
    }

    // c101
    if(ijk_valid[node].second[5]) {
      myid = id[k+1][j][i+1];
      rho  =  v[k+1][j][(i+1)*dim];
      p    =  v[k+1][j][(i+1)*dim+4];
      e    = vf[myid]->GetInternalEnergyPerUnitMass(rho,p);
      c[5] = vf[myid]->GetTemperature(rho,e);
    }

    // c011
    if(ijk_valid[node].second[6]) {
      myid = id[k+1][j+1][i];
      rho  =  v[k+1][j+1][i*dim];
      p    =  v[k+1][j+1][i*dim+4];
      e    = vf[myid]->GetInternalEnergyPerUnitMass(rho,p);
      c[6] = vf[myid]->GetTemperature(rho,e);
    }

    // c111
    if(ijk_valid[node].second[7]) {
      myid = id[k+1][j+1][i+1];
      rho  =  v[k+1][j+1][(i+1)*dim];
      p    =  v[k+1][j+1][(i+1)*dim+4];
      e    = vf[myid]->GetInternalEnergyPerUnitMass(rho,p);
      c[7] = vf[myid]->GetTemperature(rho,e);
    }


    if(ijk_valid[node].first<8) {//fill invalid slots with average value
      double c_avg = std::accumulate(c,c+8,0)/ijk_valid[node].first;
      for(int s=0; s<8; s++)
        if(!ijk_valid[node].second[s])  c[s] = c_avg;
    }

    sol1[node] = MathTools::trilinear_interpolation(trilinear_coords[node],c[0],c[1],c[2],
                                                    c[3],c[4],c[5],c[6],c[7]);
  }
  
  int mpi_rank(-1);
  MPI_Comm_rank(comm, &mpi_rank);

  //gathers solution
  MPI_Gatherv(sol1.data(), sol1.size(), MPI_DOUBLE, sol1_global.data(), package_size1.data(),
              package_disp1.data(), MPI_DOUBLE, WORKER_PROC, comm);
  
  //writes to file
  WriteScalarSolutionToFile(time, sol1_global, fname);

}

//-------------------------------------------------------------------------

void
PlaneOutput::GetDeltaTemperatureAndWrite(double time, double*** v, double*** id, string& fname)
{

  // interpolate solution
  int i,j,k;
  double c[8]; //c000,c100,c010,c110,c001,c101,c011,c111;
  double rho, p, e;
  int dim = 5, myid;
  for(int node=0; node<(int)points.size(); node++) {

    i = ijk[node][0];
    j = ijk[node][1];
    k = ijk[node][2];
 
    for(int s=0; s<8; s++)
      c[s] = 0.0;

    // c000
    if(ijk_valid[node].second[0]) {
      myid = id[k][j][i];
      rho  =  v[k][j][i*dim];
      p    =  v[k][j][i*dim+4];
      e    = vf[myid]->GetInternalEnergyPerUnitMass(rho,p);
      c[0] = vf[myid]->GetTemperature(rho,e) - vf[myid]->GetReferenceTemperature();
    }

    // c100
    if(ijk_valid[node].second[1]) {
      myid = id[k][j][i+1];
      rho  =  v[k][j][(i+1)*dim];
      p    =  v[k][j][(i+1)*dim+4];
      e    = vf[myid]->GetInternalEnergyPerUnitMass(rho,p);
      c[1] = vf[myid]->GetTemperature(rho,e) - vf[myid]->GetReferenceTemperature();
    }

    // c010
    if(ijk_valid[node].second[2]) {
      myid = id[k][j+1][i];
      rho  =  v[k][j+1][i*dim];
      p    =  v[k][j+1][i*dim+4];
      e    = vf[myid]->GetInternalEnergyPerUnitMass(rho,p);
      c[2] = vf[myid]->GetTemperature(rho,e) - vf[myid]->GetReferenceTemperature();
    }

    // c110
    if(ijk_valid[node].second[3]) {
      myid = id[k][j+1][i+1];
      rho  =  v[k][j+1][(i+1)*dim];
      p    =  v[k][j+1][(i+1)*dim+4];
      e    = vf[myid]->GetInternalEnergyPerUnitMass(rho,p);
      c[3] = vf[myid]->GetTemperature(rho,e) - vf[myid]->GetReferenceTemperature();
    }

    // c001
    if(ijk_valid[node].second[4]) {
      myid = id[k+1][j][i];
      rho  =  v[k+1][j][i*dim];
      p    =  v[k+1][j][i*dim+4];
      e    = vf[myid]->GetInternalEnergyPerUnitMass(rho,p);
      c[4] = vf[myid]->GetTemperature(rho,e) - vf[myid]->GetReferenceTemperature();
    }

    // c101
    if(ijk_valid[node].second[5]) {
      myid = id[k+1][j][i+1];
      rho  =  v[k+1][j][(i+1)*dim];
      p    =  v[k+1][j][(i+1)*dim+4];
      e    = vf[myid]->GetInternalEnergyPerUnitMass(rho,p);
      c[5] = vf[myid]->GetTemperature(rho,e) - vf[myid]->GetReferenceTemperature();
    }

    // c011
    if(ijk_valid[node].second[6]) {
      myid = id[k+1][j+1][i];
      rho  =  v[k+1][j+1][i*dim];
      p    =  v[k+1][j+1][i*dim+4];
      e    = vf[myid]->GetInternalEnergyPerUnitMass(rho,p);
      c[6] = vf[myid]->GetTemperature(rho,e) - vf[myid]->GetReferenceTemperature();
    }

    // c111
    if(ijk_valid[node].second[7]) {
      myid = id[k+1][j+1][i+1];
      rho  =  v[k+1][j+1][(i+1)*dim];
      p    =  v[k+1][j+1][(i+1)*dim+4];
      e    = vf[myid]->GetInternalEnergyPerUnitMass(rho,p);
      c[7] = vf[myid]->GetTemperature(rho,e) - vf[myid]->GetReferenceTemperature();
    }


    if(ijk_valid[node].first<8) {//fill invalid slots with average value
      double c_avg = std::accumulate(c,c+8,0)/ijk_valid[node].first;
      for(int s=0; s<8; s++)
        if(!ijk_valid[node].second[s])  c[s] = c_avg;
    }

    sol1[node] = MathTools::trilinear_interpolation(trilinear_coords[node],c[0],c[1],c[2],
                                                    c[3],c[4],c[5],c[6],c[7]);
  }
  
  int mpi_rank(-1);
  MPI_Comm_rank(comm, &mpi_rank);

  //gathers solution
  MPI_Gatherv(sol1.data(), sol1.size(), MPI_DOUBLE, sol1_global.data(), package_size1.data(),
              package_disp1.data(), MPI_DOUBLE, WORKER_PROC, comm);
  
  //writes to file
  WriteScalarSolutionToFile(time, sol1_global, fname);

}

//-------------------------------------------------------------------------

void
PlaneOutput::GetIonizationAndWrite(double time, double*** v, double*** id, string& fname)
{

  assert(ion);

  // interpolate solution
  int i,j,k;
  Vec3D c[8]; //c000,c100,c010,c110,c001,c101,c011,c111;
  double rho, p;
  int dim = 5, myid;
  for(int node=0; node<(int)points.size(); node++) {

    i = ijk[node][0];
    j = ijk[node][1];
    k = ijk[node][2];
 
    for(int s=0; s<8; s++)
      c[s] = 0.0;

    // c000
    if(ijk_valid[node].second[0]) {
      myid = id[k][j][i];
      rho  =  v[k][j][i*dim];
      p    =  v[k][j][i*dim+4];
      c[0] = ion->ComputeIonizationAtOnePoint(myid, rho, p);
    }

    // c100
    if(ijk_valid[node].second[1]) {
      myid = id[k][j][i+1];
      rho  =  v[k][j][(i+1)*dim];
      p    =  v[k][j][(i+1)*dim+4];
      c[1] = ion->ComputeIonizationAtOnePoint(myid, rho, p);
    }

    // c010
    if(ijk_valid[node].second[2]) {
      myid = id[k][j+1][i];
      rho  =  v[k][j+1][i*dim];
      p    =  v[k][j+1][i*dim+4];
      c[2] = ion->ComputeIonizationAtOnePoint(myid, rho, p);
    }

    // c110
    if(ijk_valid[node].second[3]) {
      myid = id[k][j+1][i+1];
      rho  =  v[k][j+1][(i+1)*dim];
      p    =  v[k][j+1][(i+1)*dim+4];
      c[3] = ion->ComputeIonizationAtOnePoint(myid, rho, p);
    }

    // c001
    if(ijk_valid[node].second[4]) {
      myid = id[k+1][j][i];
      rho  =  v[k+1][j][i*dim];
      p    =  v[k+1][j][i*dim+4];
      c[4] = ion->ComputeIonizationAtOnePoint(myid, rho, p);
    }

    // c101
    if(ijk_valid[node].second[5]) {
      myid = id[k+1][j][i+1];
      rho  =  v[k+1][j][(i+1)*dim];
      p    =  v[k+1][j][(i+1)*dim+4];
      c[5] = ion->ComputeIonizationAtOnePoint(myid, rho, p);
    }

    // c011
    if(ijk_valid[node].second[6]) {
      myid = id[k+1][j+1][i];
      rho  =  v[k+1][j+1][i*dim];
      p    =  v[k+1][j+1][i*dim+4];
      c[6] = ion->ComputeIonizationAtOnePoint(myid, rho, p);
    }

    // c111
    if(ijk_valid[node].second[7]) {
      myid = id[k+1][j+1][i+1];
      rho  =  v[k+1][j+1][(i+1)*dim];
      p    =  v[k+1][j+1][(i+1)*dim+4];
      c[7] = ion->ComputeIonizationAtOnePoint(myid, rho, p);
    }


    if(ijk_valid[node].first<8) {//fill invalid slots with average value
      Vec3D c_avg = std::accumulate(c,c+8,Vec3D(0.0))/(double)ijk_valid[node].first;
      for(int s=0; s<8; s++)
        if(!ijk_valid[node].second[s])  c[s] = c_avg;
    }

    sol3[node] = MathTools::trilinear_interpolation(trilinear_coords[node],c[0],c[1],c[2],
                                                    c[3],c[4],c[5],c[6],c[7]);
  }
  
  int mpi_rank(-1);
  MPI_Comm_rank(comm, &mpi_rank);

  //gathers solution
  MPI_Gatherv((double*)sol3.data(), sol3.size(), MPI_DOUBLE, (double*)sol3_global.data(),
              package_size3.data(), package_disp3.data(), MPI_DOUBLE, WORKER_PROC, comm);
  
  //writes to file
  WriteVec3SolutionToFile(time, sol3_global, fname);

}

//-------------------------------------------------------------------------

void
PlaneOutput::WriteScalarSolutionToFile(double time, vector<double>& S, string& fname)
{
  int mpi_rank(-1);
  MPI_Comm_rank(comm, &mpi_rank);

  //writes to file
  if(mpi_rank==WORKER_PROC) {
    FILE *file = fopen(fname.c_str(), "a");
    assert(file); //if file is not opened, the pointer would be NULL
    fprintf(file,"%e\n", time);
    for(auto&& s : S)
      fprintf(file,"%16.8e\n", s);
    fclose(file);
  }
}

//-------------------------------------------------------------------------

void
PlaneOutput::WriteVec3SolutionToFile(double time, vector<Vec3D>& S, string& fname)
{
  int mpi_rank(-1);
  MPI_Comm_rank(comm, &mpi_rank);

  //writes to file
  if(mpi_rank==WORKER_PROC) {
    FILE *file = fopen(fname.c_str(), "a");
    assert(file); //if file is not opened, the pointer would be NULL
    fprintf(file,"%e\n", time);
    for(auto&& s : S)
      fprintf(file,"%16.8e    %16.8e    %16.8e\n", s[0], s[1], s[2]);
    fclose(file);
  }
}

//-------------------------------------------------------------------------







//-------------------------------------------------------------------------

