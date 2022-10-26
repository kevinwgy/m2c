#include<PlaneOutput.h>
#include<Intersections.h>
#include<numeric> //std::accumulate

using std::vector;
using std::array;
using std::pair;
using std::string;

//-------------------------------------------------------------------------

PlaneOutput::PlaneOutput(MPI_Comm &comm_, OutputData &iod_output_, PlanePlot &iod_pplot_, 
                         std::vector<VarFcnBase*> &vf_, GlobalMeshInfo &global_mesh_,
                         IonizationOperator *ion_)
           : comm(comm_), iod_pplot(iod_pplot_), vf(vf_), global_mesh(global_mesh_),
             probe_util(comm_, iod_output_, vf_, ion_)
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
  print("- [Cut-Plane] Coords = (%e, %e, %e), Normal = (%e, %e, %e).\n", p[0], p[1], p[2],
        n[0], n[1], n[2]);


  if(iod_pplot.mesh[0] == 0) {
    print_error("*** Error: [Cut-Plane] Mesh file name is not provided.\n");
    exit(-1);
  } 

  // ----------------------------------------------
  // RETURN IF NOT PROC 0 (ONLY PROC 0 CONTINUES)
  // ----------------------------------------------
  int worker = 0;
  int mpi_rank(-1);
  MPI_Comm_rank(comm, &mpi_rank);
  if(mpi_rank != worker)
    return;


  // ----------------------------------------------
  // set up files (Similar to ProbeOutput constructor)
  // ----------------------------------------------
  filename_mesh = string(iod_output_.prefix) + string(iod_pplot.mesh);

  if (iod_pplot.density[0] != 0) {
    filename[PlanePlot::DENSITY] = string(iod_output_.prefix) + string(iod_pplot.density);
    FILE *file = fopen(filename[PlanePlot::DENSITY].c_str(), "w");
    fprintf(file, "Scalar Density under load for SurfaceNodes\n");
    fclose(file);
  }

  if (iod_pplot.velocity[0] != 0) {
    filename[PlanePlot::VELOCITY] = string(iod_output_.prefix) + string(iod_pplot.velocity);
    FILE *file = fopen(filename[PlanePlot::VELOCITY].c_str(), "w");
    fprintf(file, "Vector Velocity under load for SurfaceNodes\n");
    fclose(file);
  }

  if (iod_pplot.pressure[0] != 0) {
    filename[PlanePlot::PRESSURE] = string(iod_output_.prefix) + string(iod_pplot.pressure);
    FILE *file = fopen(filename[PlanePlot::PRESSURE].c_str(), "w");
    fprintf(file, "Scalar Pressure under load for SurfaceNodes\n");
    fclose(file);
  }

  if (iod_pplot.temperature[0] != 0) {
    filename[PlanePlot::TEMPERATURE] = string(iod_output_.prefix) + string(iod_pplot.temperature);
    FILE *file = fopen(filename[PlanePlot::TEMPERATURE].c_str(), "w");
    fprintf(file, "Scalar Temperature under load for SurfaceNodes\n");
    fclose(file);
  }

  if (iod_pplot.delta_temperature[0] != 0) {
    filename[PlanePlot::DELTA_TEMPERATURE] = string(iod_output_.prefix) + string(iod_pplot.delta_temperature);
    FILE *file = fopen(filename[PlanePlot::DELTA_TEMPERATURE].c_str(), "w");
    fprintf(file, "Scalar DeltaTemperature under load for SurfaceNodes\n");
    fclose(file);
  }

  if (iod_pplot.materialid[0] != 0) {
    filename[PlanePlot::MATERIALID] = string(iod_output_.prefix) + string(iod_pplot.materialid);
    FILE *file = fopen(filename[PlanePlot::MATERIALID].c_str(), "w");
    fprintf(file, "Scalar MaterialID under load for SurfaceNodes\n");
    fclose(file);
  }

  if (iod_pplot.laser_radiance[0] != 0) {
    filename[PlanePlot::LASERRADIANCE] = string(iod_output_.prefix) + string(iod_pplot.laser_radiance);
    FILE *file = fopen(filename[PlanePlot::LASERRADIANCE].c_str(), "w");
    fprintf(file, "Scalar LaserRadiance under load for SurfaceNodes\n");
    fclose(file);
  }


  if (iod_pplot.levelset0[0] != 0) {
    filename[PlanePlot::LEVELSET0] = string(iod_output_.prefix) + string(iod_pplot.levelset0);
    FILE *file = fopen(filename[PlanePlot::LEVELSET0].c_str(), "w");
    fprintf(file, "Scalar LevelSet0 under load for SurfaceNodes\n");
    fclose(file);
  }

  if (iod_pplot.levelset1[0] != 0) {
    filename[PlanePlot::LEVELSET1] = string(iod_output_.prefix) + string(iod_pplot.levelset1);
    FILE *file = fopen(filename[PlanePlot::LEVELSET1].c_str(), "w");
    fprintf(file, "Scalar LevelSet1 under load for SurfaceNodes\n");
    fclose(file);
  }

  if (iod_pplot.levelset2[0] != 0) {
    filename[PlanePlot::LEVELSET2] = string(iod_output_.prefix) + string(iod_pplot.levelset2);
    FILE *file = fopen(filename[PlanePlot::LEVELSET2].c_str(), "w");
    fprintf(file, "Scalar LevelSet2 under load for SurfaceNodes\n");
    fclose(file);
  }

  if (iod_pplot.levelset3[0] != 0) {
    filename[PlanePlot::LEVELSET3] = string(iod_output_.prefix) + string(iod_pplot.levelset3);
    FILE *file = fopen(filename[PlanePlot::LEVELSET3].c_str(), "w");
    fprintf(file, "Scalar LevelSet3 under load for SurfaceNodes\n");
    fclose(file);
  }

  if (iod_pplot.levelset4[0] != 0) {
    filename[PlanePlot::LEVELSET4] = string(iod_output_.prefix) + string(iod_pplot.levelset4);
    FILE *file = fopen(filename[PlanePlot::LEVELSET4].c_str(), "w");
    fprintf(file, "Scalar LevelSet4 under load for SurfaceNodes\n");
    fclose(file);
  }

  if (iod_pplot.ionization_result[0] != 0) {
    if(!ion_) {
      fprintf(stderr,"\033[0;31m*** Error: Requested ionization result, "
              "without specifying an ionization model.\033[0m\n");
      exit(-1);
    }
    filename[PlanePlot::IONIZATION] = string(iod_output_.prefix) + string(iod_pplot.ionization_result);
    FILE *file = fopen(filename[PlanePlot::IONIZATION].c_str(), "w");
    fprintf(file, "Scalar IonizationResult under load for SurfaceNodes\n");
    fclose(file);
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

        int n0 = points.size();
        points.insert(points.end(), xlocal.begin(), xlocal.begin()+np);


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
            elems.push_back(array<int,4>{n0,n0+3,n0+4});
            break;
          case 6:
            elems.push_back(array<int,4>{n0,n0+1,n0+2,n0+3});
            elems.push_back(array<int,4>{n0,n0+3,n0+4,n0+5});
            break;
        }  

        // interpolation info
        Int3 ijk0;
        Vec3D xi0;
        for(int n=0; n<np; n++) {

          bool found = global_mesh.FindElementCoveringPoint(xlocal[n], ijk0, &xi0, true);
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

      }

  coordinates.RestoreDataPointerToLocalVector();

  int nNodes = points.size();
  int nElems = elems.size();

  // ---------------
  // MPI comms
  // ---------------
  int worker = 0;
  int mpi_rank(-1), mpi_size(-1);
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  vector<int> nNodes_all;
  vector<int> nElems_all;
  if(mpi_rank==worker) { //lead proc

    nNodes_all.resize(mpi_size);
    nElems_all.resize(mpi_size);
    MPI_Gather(&nNodes, 1, MPI_INT, nNodes_all.data(), 1, MPI_INT, worker, comm);
    MPI_Gather(&nElems, 1, MPI_INT, nElems_all.data(), 1, MPI_INT, worker, comm);

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
    MPI_Gather(&nNodes, 1, MPI_INT, NULL, 0, MPI_INT, worker, comm);
    MPI_Gather(&nElems, 1, MPI_INT, NULL, 0, MPI_INT, worker, comm);
  }

  // worker gathers all the nodes
  MPI_Gatherv((double*)points.data(), 3*points.size(), MPI_DOUBLE,
              (double*)sol3_global.data(), package_size3.data(),
              package_disp3.data(), MPI_DOUBLE, worker, comm);


  // worker gathers all the elements
  vector<int> elems_size(mpi_size);
  vector<int> elems_disp(mpi_size);
  vector<array<int,4> > elems_global;
  if(mpi_rank==worker) {
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
              elems_disp.data(), MPI_INT, worker, comm);

  // worker updates node ids in elements
  if(mpi_rank==worker) {
    for(int proc=1; proc<mpi_size; proc++) {
      int n0 = package_disp1[proc];
      for(int e=elems_disp[proc]/4; e<nElems_all[proc]; e++)
        for(int n=0; n<4; n++)
          elems_global[e][n] += n0;
    }
  }

  // worker prints mesh to file 
  if(mpi_rank==worker) {

    FILE *file = fopen(filename_mesh.c_str(), "a");
    fprintf(file, "Nodes CutPlaneNodes\n");
    for(int i=0; i<sol3_global.size(); i++)
      fprintf(file, "%10d    %16.8e    %16.8e    %16.8e\n",
              i+1, sol3_global[i][0], sol3_global[i][1], sol3_global[i][2]);
    fprintf(file, "Elements CutPlane using CutPlaneNodes\n");
    for(int i=0; i<elems_global.size(); i++) {
      if(elems_global[i][2]==elems_global[i][3]) //triangle
        fprintf(file, "%10d    4    %10d    %10d    %10d\n", i+1,
                elems_global[i][0], elems_global[i][1], elems_global[i][2]);
      else //quad
        fprintf(file, "%10d    2    %10d    %10d    %10d    %10d\n", i+1,
                elems_global[i][0], elems_global[i][1], elems_global[i][2],
                elems_global[i][3]);
    }
    fclose(file);

  }


}

//-------------------------------------------------------------------------

void
PlaneOutput::WriteSolutionOnPlane(double time, double dt, int time_step, SpaceVariable3D &V,
                                  SpaceVariable3D &ID, std::vector<SpaceVariable3D*> &Phi,
                                  SpaceVariable3D* L, bool force_write)
{




}

//-------------------------------------------------------------------------


