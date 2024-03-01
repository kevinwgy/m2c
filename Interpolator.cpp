/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<Interpolator.h>

//--------------------------------------------------------------------------

InterpolatorLinear::InterpolatorLinear(MPI_Comm &comm_, DataManagers3D &dm_all_,
                                       SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_)
                   : InterpolatorBase(coordinates_,delta_xyz_),
                     Cl(comm_, &(dm_all_.ghosted1_1dof)),
                     Cr(comm_, &(dm_all_.ghosted1_1dof)),
                     Ct(comm_, &(dm_all_.ghosted1_1dof)),
                     Cb(comm_, &(dm_all_.ghosted1_1dof)),
                     Ck(comm_, &(dm_all_.ghosted1_1dof)),
                     Cf(comm_, &(dm_all_.ghosted1_1dof))
{
  Cl.SetConstantValue(0.0,true);
  Cr.SetConstantValue(0.0,true);
  Ct.SetConstantValue(0.0,true);
  Cb.SetConstantValue(0.0,true);
  Ck.SetConstantValue(0.0,true);
  Cf.SetConstantValue(0.0,true);

  CalculateInterpolationCoefficients();
}

//--------------------------------------------------------------------------

void InterpolatorLinear::InterpolateAtCellInterfaces(int dir/*0~x,1~y,2~z*/,
                             SpaceVariable3D &Vin, std::vector<int> &input_dof,
                             SpaceVariable3D &Vout, std::vector<int> &output_dof)
{
  //check for common errors
  if(input_dof.size() != output_dof.size()) {
    print_error("*** Error: Interpolation failed due to inconsistent number of dof's (%d vs. %d)\n",
                input_dof.size(), output_dof.size());
    exit_mpi();
  }
  if(*std::min_element(input_dof.begin(), input_dof.end())<0 ||
     *std::max_element(input_dof.begin(), input_dof.end())>=Vin.NumDOF()) {
    print_error("*** Error: Interpolation failed due to incorrect dof number(s) (min:%d, max:%d).\n",
                *std::min_element(input_dof.begin(), input_dof.end()),
                *std::max_element(input_dof.begin(), input_dof.end()));
    exit_mpi();
  }
  if(*std::min_element(output_dof.begin(), output_dof.end())<0 ||
     *std::max_element(output_dof.begin(), output_dof.end())>=Vout.NumDOF()) {
    print_error("*** Error: Interpolation failed due to incorrect dof number(s) (min:%d, max:%d).\n",
                *std::min_element(output_dof.begin(), output_dof.end()),
                *std::max_element(output_dof.begin(), output_dof.end()));
    exit_mpi();
  }

  // Get Data
  double*** vin  = (double***)Vin.GetDataPointer();
  double*** vout = (double***)Vout.GetDataPointer();
  int DOFin = Vin.NumDOF();
  int DOFout = Vout.NumDOF();

  // Loop through the domain interior, and the right, top, and front ghost layers. For each cell, calculate
  // interpolation at the left, lower, and back cell boundaries/interfaces
  if(dir == 0) {
    // Interpolate Vin at i-1/2 and i+1/2 within the physical domain (not the ghost boundary outside the physical
    // domain). Specifically, Vout[k][j][i] is the interpolated value at (i-1/2, j, k)
    int pin, pout;
    double*** cl = (double***)Cl.GetDataPointer();
    double*** cr = (double***)Cr.GetDataPointer();
    for(int k=k0; k<kkmax; k++)
      for(int j=j0; j<jjmax; j++)
        for(int i=i0; i<iimax; i++) {

          if(k==kkmax-1 || j==jjmax-1) 
            continue; //either done in another subdomain, or not needed (outside physical domain)

          for(int id=0; id<(int)input_dof.size(); id++) {
            pin  = input_dof[id];
            pout = output_dof[id];
            vout[k][j][i*DOFout+pout] = cl[k][j][i-1]*vin[k][j][(i-1)*DOFin+pin] 
                                      + cr[k][j][i]*vin[k][j][i*DOFin+pin];
          }

        }
    Cl.RestoreDataPointerAndInsert();
    Cr.RestoreDataPointerAndInsert();
  }
  else if(dir == 1) {
    // Interpolate Vin at j-1/2 and j+1/2 within the physical domain (not the ghost boundary outside the physical
    // domain). Specifically, Vout[k][j][i] is the interpolated value at (i, j-1/2, k)
    int pin, pout;
    double*** cb = (double***)Cb.GetDataPointer();
    double*** ct = (double***)Ct.GetDataPointer();
    for(int k=k0; k<kkmax; k++)
      for(int j=j0; j<jjmax; j++)
        for(int i=i0; i<iimax; i++) {

          if(i==iimax-1 || k==kkmax-1)
            continue;

          for(int id=0; id<(int)input_dof.size(); id++) {
            pin  = input_dof[id];
            pout = output_dof[id];
            vout[k][j][i*DOFout+pout] = cb[k][j-1][i]*vin[k][j-1][i*DOFin+pin] 
                                      + ct[k][j][i]*vin[k][j][i*DOFin+pin];
          } 

        }
    Cb.RestoreDataPointerAndInsert();
    Ct.RestoreDataPointerAndInsert();
  }
  else if(dir == 2) {
    // Interpolate Vin at k-1/2 and k+1/2 within the physical domain (not the ghost boundary outside the physical
    // domain). Specifically, Vout[k][j][i] is the interpolated value at (i, j, k-1/2)
    int pin, pout;
    double*** ck = (double***)Ck.GetDataPointer();
    double*** cf = (double***)Cf.GetDataPointer();
    for(int k=k0; k<kkmax; k++)
      for(int j=j0; j<jjmax; j++)
        for(int i=i0; i<iimax; i++) {

          if(i==iimax-1 || j==jjmax-1)
            continue;

          for(int id=0; id<(int)input_dof.size(); id++) {
            pin  = input_dof[id];
            pout = output_dof[id];
            vout[k][j][i*DOFout+pout] = ck[k-1][j][i]*vin[k-1][j][i*DOFin+pin] 
                                      + cf[k][j][i]*vin[k][j][i*DOFin+pin];
          } 

        }
    Ck.RestoreDataPointerAndInsert();
    Cf.RestoreDataPointerAndInsert();
  }

  Vin.RestoreDataPointerAndInsert();
  Vout.RestoreDataPointerAndInsert();
}

//--------------------------------------------------------------------------

void InterpolatorLinear::CalculateInterpolationCoefficients()
{

  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
  Vec3D*** dxyz   = (Vec3D***)delta_xyz.GetDataPointer();

  double*** cl = (double***)Cl.GetDataPointer();
  double*** cr = (double***)Cr.GetDataPointer();
  double*** cb = (double***)Cb.GetDataPointer();
  double*** ct = (double***)Ct.GetDataPointer();
  double*** ck = (double***)Ck.GetDataPointer();
  double*** cf = (double***)Cf.GetDataPointer();

  // Loop through the real cells and the right, upper, and front ghost cells
  // Note: The notation "i/j/k_minus_half" just indicates the cell interface
  //       between i/j/k and i/j/k-1. It is NOT at the middle if the two cells have
  //       different widths.
  for(int k=k0; k<kkmax; k++)
    for(int j=j0; j<jjmax; j++)
      for(int i=i0; i<iimax; i++) {

        double x_i_minus_half = coords[k][j][i][0] - 0.5*dxyz[k][j][i][0];
        double dx = coords[k][j][i][0] - coords[k][j][i-1][0];
        cl[k][j][i-1] = (coords[k][j][i][0] - x_i_minus_half)/dx;
        cr[k][j][i]   = (x_i_minus_half - coords[k][j][i-1][0])/dx;
        
        double y_j_minus_half = coords[k][j][i][1] - 0.5*dxyz[k][j][i][1];
        double dy = coords[k][j][i][1] - coords[k][j-1][i][1];
        cb[k][j-1][i] = (coords[k][j][i][1] - y_j_minus_half)/dy;
        ct[k][j][i]   = (y_j_minus_half - coords[k][j-1][i][1])/dy;

        double z_k_minus_half = coords[k][j][i][2] - 0.5*dxyz[k][j][i][2];
        double dz = coords[k][j][i][2] - coords[k-1][j][i][2];
        ck[k-1][j][i] = (coords[k][j][i][2] - z_k_minus_half)/dz;
        cf[k][j][i]   = (z_k_minus_half - coords[k-1][j][i][2])/dz;

      }

  Cl.RestoreDataPointerAndInsert();
  Cr.RestoreDataPointerAndInsert();
  Cb.RestoreDataPointerAndInsert();
  Ct.RestoreDataPointerAndInsert();
  Ck.RestoreDataPointerAndInsert();
  Cf.RestoreDataPointerAndInsert();

  delta_xyz.RestoreDataPointerToLocalVector(); //no changes
  coordinates.RestoreDataPointerToLocalVector(); //no changes

}

//--------------------------------------------------------------------------





