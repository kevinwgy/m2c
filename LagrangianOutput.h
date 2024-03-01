/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _LAGRANGIAN_OUTPUT_H_
#define _LAGRANGIAN_OUTPUT_H_

#include<IoData.h>
#include<Utils.h>

/****************************************
 * class LagrangianOutput is responsible
 * for printing results on a Lagrangian
 * mesh to files in the ASCII "xpost" 
 * format. This class does not store 
 * any result. But it is responsible for
 * checking IoData as to what and when 
 * to write.
 ***************************************/

class LagrangianOutput {

  LagrangianMeshOutputData& iod_lag;

  MPI_Comm& comm;

  int iFrame; //!< frame id
  double last_snapshot_time; //!< latest time when solution snapshot is written to file

  FILE* disp_file;
  FILE* sol_file;
  FILE* sol2_file;

public:

  LagrangianOutput(MPI_Comm &comm_, LagrangianMeshOutputData &iod_lag_);
  ~LagrangianOutput();

  void OutputTriangulatedMesh(std::vector<Vec3D>& X0, std::vector<Int3>& elems);
  void OutputResults(double t, double dt, int time_step, std::vector<Vec3D>& X0, std::vector<Vec3D>& X, 
                     std::vector<Vec3D>& F, std::vector<Vec3D>* F2_ptr, bool force_write); 
                
private:

  void AppendResultToFile(FILE* file, double time, int N, int dim, double* S);

};

#endif
