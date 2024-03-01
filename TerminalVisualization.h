/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _TERMINAL_VISUALIZATION_H_
#define _TERMINAL_VISUALIZATION_H_

#include<GlobalMeshInfo.h>
#include<IonizationOperator.h>
#include<fstream>
#include<time.h>

/*****************************************************************
 * class TerminalVisualization is reponsible for printing a user
 * specified solution field to either a file or the screen using
 * ANSI color codes (256 colors). This is a feature for quick
 * diagnostic of simulations.
 ****************************************************************/

class TerminalVisualization
{

  MPI_Comm &comm;
  int mpi_rank, mpi_size;

  int number_colors;

  TerminalVisualizationData &iod_terminal;
  std::vector<VarFcnBase*> &vf;

  //! global mesh
  GlobalMeshInfo &global_mesh;

  //! post-processor
  IonizationOperator* ion;

  //! output stream
  std::ofstream outfile;
  std::ostream *myout;

  int iFrame;
  double last_snapshot_time;
  bool use_clocktime;
  double last_snapshot_clocktime;

  //! colormap
  std::vector<Vec3D> turbo_rgb;
  std::vector<int>   turbo_ANSI;
  std::vector<int>   gray_ANSI;
  int under_color_ANSI, over_color_ANSI; //!< for values below/over the min/max.

  //! visualization grid
  int nrows, ncols;
  Vec3D xyzmin, xyzmax;
  double dx, dh, dv;
  std::vector<Int3> ijk;
  std::vector<double> sol; 
  

public:

  TerminalVisualization(MPI_Comm &comm_, TerminalVisualizationData &iod_terminal_, GlobalMeshInfo &global_mesh_,
                        vector<VarFcnBase*> &vf_, IonizationOperator* ion_ = NULL);

  ~TerminalVisualization();

  void PrintSolutionSnapshot(double time, double dt, int time_step, SpaceVariable3D &V, SpaceVariable3D &ID,
                             std::vector<SpaceVariable3D*> &Phi, SpaceVariable3D *L, bool force_write);

private:

  void SetupGrayColorMap();
  void SetupTurboColorMap();
  void SetupTurboRGB();

};

#endif
