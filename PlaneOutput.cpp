#include<PlaneOutput.h>

using std::vector;
using std::array;
using std::pair;

//-------------------------------------------------------------------------

PlaneOutput::PlaneOutput(MPI_Comm &comm_, OutputData &iod_output_, PlanePlot &iod_pplot_, 
                         std::vector<VarFcnBase*> &vf_, GlobalMeshInfo &global_mesh_,
                         IonizationOperator *ion_)
           : comm(comm_), iod_pplot(iod_pplot_), vf(vf_), 
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


  for(int i=0; i<PlanePlot::SIZE; i++)
    file[i] = NULL;


  // ----------------------------------------------
  // Set up mesh and interpolation formulas (and print to a "top" file)
  // ----------------------------------------------
  int total_nodes = SetupMeshAndInterpolation(p, n); //only proc #0 will print


  // ----------------------------------------------
  // RETURN IF NOT PROC 0 (ONLY PROC 0 CONTINUES)
  // ----------------------------------------------
  int mpi_rank(-1);
  MPI_Comm_rank(comm, &mpi_rank);
  if(mpi_rank)
    return;


  // ----------------------------------------------
  // set up files (Similar to ProbeOutput constructor)
  // ----------------------------------------------
  int spn = strlen(iod_output.prefix) + 1;

  if (iod_pplot.density[0] != 0) {
    char *filename = new char[spn + strlen(iod_pplot.density)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_pplot.density);
    file[PlanePlot::DENSITY] = fopen(filename, "w");
    fprintf(file[PlanePlot::DENSITY], "Scalar Density under load for SurfaceNodes\n%d\n",
            total_nodes);
    fclose(file[PlanePlot::DENSITY]);
    delete [] filename;
  }

  if (iod_pplot.velocity[0] != 0) {
    char *filename = new char[spn + strlen(iod_pplot.velocity)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_pplot.velocity);
    file[PlanePlot::VELOCITY] = fopen(filename, "w");
    fprintf(file[PlanePlot::VELOCITY], "Vector Velocity under load for SurfaceNodes\n%d\n",
            total_nodes);
    fclose(file[PlanePlot::VELOCITY]);
    delete [] filename;
  }

  if (iod_pplot.pressure[0] != 0) {
    char *filename = new char[spn + strlen(iod_pplot.pressure)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_pplot.pressure);
    file[PlanePlot::PRESSURE] = fopen(filename, "w");
    fprintf(file[PlanePlot::PRESSURE], "Scalar Pressure under load for SurfaceNodes\n%d\n",
            total_nodes);
    fclose(file[PlanePlot::PRESSURE]);
    delete [] filename;
  }

  if (iod_pplot.temperature[0] != 0) {
    char *filename = new char[spn + strlen(iod_pplot.temperature)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_pplot.temperature);
    file[PlanePlot::TEMPERATURE] = fopen(filename, "w");
    fprintf(file[PlanePlot::TEMPERATURE], "Scalar Temperature under load for SurfaceNodes\n%d\n",
            total_nodes);
    fclose(file[PlanePlot::TEMPERATURE]);
    delete [] filename;
  }

  if (iod_pplot.delta_temperature[0] != 0) {
    char *filename = new char[spn + strlen(iod_pplot.delta_temperature)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_pplot.delta_temperature);
    file[PlanePlot::DELTA_TEMPERATURE] = fopen(filename, "w");
    fprintf(file[PlanePlot::DELTA_TEMPERATURE], "Scalar DeltaTemperature under load for SurfaceNodes\n%d\n",
            total_nodes);
    fclose(file[PlanePlot::DELTA_TEMPERATURE]);
    delete [] filename;
  }

  if (iod_pplot.materialid[0] != 0) {
    char *filename = new char[spn + strlen(iod_pplot.materialid)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_pplot.materialid);
    file[PlanePlot::MATERIALID] = fopen(filename, "w");
    fprintf(file[PlanePlot::MATERIALID], "Scalar MaterialID under load for SurfaceNodes\n%d\n",
            total_nodes);
    fclose(file[PlanePlot::MATERIALID]);
    delete [] filename;
  }

  if (iod_pplot.laser_radiance[0] != 0) {
    char *filename = new char[spn + strlen(iod_pplot.materialid)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_pplot.laser_radiance);
    file[PlanePlot::LASERRADIANCE] = fopen(filename, "w");
    fprintf(file[PlanePlot::LASERRADIANCE], "Scalar LaserRadiance under load for SurfaceNodes\n%d\n",
            total_nodes);
    fclose(file[PlanePlot::LASERRADIANCE]);
    delete [] filename;
  }


  if (iod_pplot.levelset0[0] != 0) {
    char *filename = new char[spn + strlen(iod_pplot.levelset0)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_pplot.levelset0);
    file[PlanePlot::LEVELSET0] = fopen(filename, "w");
    fprintf(file[PlanePlot::LEVELSET0], "Scalar LevelSet0 under load for SurfaceNodes\n%d\n",
            total_nodes);
    fclose(file[PlanePlot::LEVELSET0]);
    delete [] filename;
  }

  if (iod_pplot.levelset1[0] != 0) {
    char *filename = new char[spn + strlen(iod_pplot.levelset1)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_pplot.levelset1);
    file[PlanePlot::LEVELSET1] = fopen(filename, "w");
    fprintf(file[PlanePlot::LEVELSET1], "Scalar LevelSet1 under load for SurfaceNodes\n%d\n",
            total_nodes);
    fclose(file[PlanePlot::LEVELSET1]);
    delete [] filename;
  }

  if (iod_pplot.levelset2[0] != 0) {
    char *filename = new char[spn + strlen(iod_pplot.levelset2)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_pplot.levelset2);
    file[PlanePlot::LEVELSET2] = fopen(filename, "w");
    fprintf(file[PlanePlot::LEVELSET2], "Scalar LevelSet2 under load for SurfaceNodes\n%d\n",
            total_nodes);
    fclose(file[PlanePlot::LEVELSET2]);
    delete [] filename;
  }

  if (iod_pplot.levelset3[0] != 0) {
    char *filename = new char[spn + strlen(iod_pplot.levelset3)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_pplot.levelset3);
    file[PlanePlot::LEVELSET3] = fopen(filename, "w");
    fprintf(file[PlanePlot::LEVELSET3], "Scalar LevelSet3 under load for SurfaceNodes\n%d\n",
            total_nodes);
    fclose(file[PlanePlot::LEVELSET3]);
    delete [] filename;
  }

  if (iod_pplot.levelset4[0] != 0) {
    char *filename = new char[spn + strlen(iod_pplot.levelset4)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_pplot.levelset4);
    file[PlanePlot::LEVELSET4] = fopen(filename, "w");
    fprintf(file[PlanePlot::LEVELSET4], "Scalar LevelSet4 under load for SurfaceNodes\n%d\n",
            total_nodes);
    fclose(file[PlanePlot::LEVELSET4]);
    delete [] filename;
  }


  if (iod_pplot.ionization_result[0] != 0) {
    if(!ion) {
      fprintf(stderr,"*** Error: Requested ionization result, without specifying an ionization model.\n");
      exit(-1);
    }
    char *filename = new char[spn + strlen(iod_pplot.ionization_result)];
    sprintf(filename, "%s%s", iod_output.prefix, iod_pplot.ionization_result);
    file[PlanePlot::IONIZATION] = fopen(filename, "w");
    fprintf(file[PlanePlot::IONIZATION], "Vector IonizationResult under load for SurfaceNodes\n%d\n",
            total_nodes);
    fclose(file[PlanePlot::IONIZATION]);
    delete [] filename;
  }

}

//-------------------------------------------------------------------------

PlaneOutput::~PlaneOutput()
{
  for(int i=0; i<PlanePlot::SIZE; i++)
    if(file[i]) fclose(file[i]);
}

//-------------------------------------------------------------------------

int
PlaneOutput::SetupMeshAndInterpolation(Vec3D& point, Vec3D& normal)
{
  for(



}

//-------------------------------------------------------------------------

void
PlaneOutput::WriteSolutionOnPlane(double time, double dt, int time_step, SpaceVariable3D &V,
                                  SpaceVariable3D &ID, std::vector<SpaceVariable3D*> &Phi,
                                  SpaceVariable3D* L, bool force_write)
{




}

//-------------------------------------------------------------------------


