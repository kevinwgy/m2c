/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <GradientCalculatorFD3.h>
#include <algorithm>
#include <cassert>
#include <Eigen/Dense> //linear solve

using std::vector;

//--------------------------------------------------------------------------

GradientCalculatorFD3::GradientCalculatorFD3(MPI_Comm &comm_, DataManagers3D &dm_all_,
                               SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_,
                               int shift_)
                          : GradientCalculatorBase(coordinates_, delta_xyz_),
                            Cx(comm_, &(dm_all_.ghosted1_4dof)),
                            Cy(comm_, &(dm_all_.ghosted1_4dof)),
                            Cz(comm_, &(dm_all_.ghosted1_4dof)),
                            Var(comm_, &(dm_all_.ghosted1_5dof)),
                            coordinatesG2(comm_, &(dm_all_.ghosted2_3dof)),
                            shift(shift_)
                   
{
  Cx.SetConstantValue(0.0,true);
  Cy.SetConstantValue(0.0,true);
  Cz.SetConstantValue(0.0,true);
  Var.SetConstantValue(0.0,true);

  //build coordinatesG2. Still, only one layer of nodes outside physical domain is populated
  Vec3D*** coordsG2 = (Vec3D***)coordinatesG2.GetDataPointer();
  Vec3D*** coords   = (Vec3D***)coordinates.GetDataPointer();
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++) 
        coordsG2[k][j][i] = coords[k][j][i];
  coordinates.RestoreDataPointerToLocalVector();
  coordinatesG2.RestoreDataPointerAndInsert();

  if(shift == -1) 
    CalculateCoefficientsLeftBiased();
  else if(shift == 1)
    CalculateCoefficientsRightBiased();
  else {
    print_error("*** Error: Cannot initialize GradientCalculatorFD3 for shift = %d.\n", shift);
    exit_mpi();
  }
}

//--------------------------------------------------------------------------

void GradientCalculatorFD3::CalculateFirstDerivativeAtNodes(int dir/*0~d/dx,1~d/dy,2~d/dz*/,
                                SpaceVariable3D &V, std::vector<int> &input_dof,
                                SpaceVariable3D &DV, std::vector<int> &output_dof)
{
  //check for common errors
  if(input_dof.size() != output_dof.size()) {
    print_error("*** Error: CalculateFirstDerivativeAtNodes failed due to inconsistent number of dof's (%d vs. %d)\n",
                input_dof.size(), output_dof.size());
    exit_mpi();
  }
  if(*std::min_element(input_dof.begin(), input_dof.end())<0 ||
     *std::max_element(input_dof.begin(), input_dof.end())>=V.NumDOF()) {
    print_error("*** Error: CalculateFirstDerivativeAtNodes failed due to incorrect dof number(s) (min:%d, max:%d).\n",
                *std::min_element(input_dof.begin(), input_dof.end()),
                *std::max_element(input_dof.begin(), input_dof.end()));
    exit_mpi();
  }
  if(*std::min_element(output_dof.begin(), output_dof.end())<0 ||
     *std::max_element(output_dof.begin(), output_dof.end())>=DV.NumDOF()) {
    print_error("*** Error: CalculateFirstDerivativeAtNodes failed due to incorrect dof number(s) (min:%d, max:%d).\n",
                *std::min_element(output_dof.begin(), output_dof.end()),
                *std::max_element(output_dof.begin(), output_dof.end()));
    exit_mpi();
  }
  if(V.NumGhostLayers()<2) {
    print_error("*** Error: CalculateFirstDerivativeAtNodes(FD3) failed. Input variable should have at least 2 ghost layers.\n");
    exit_mpi();
  }

  // Get data
  double*** v = (double***)V.GetDataPointer();
  double*** dv = (double***)DV.GetDataPointer(); 

  int DOFin = V.NumDOF(), DOFout = DV.NumDOF();

  int q0;
  assert(shift==-1 || shift==1);
  if(shift == -1)
    q0 = 2;
  else
    q0 = 1;
  
    
  //-------------------------------------------------------------
  // Loop through the domain interior 
  //-------------------------------------------------------------

  if(dir==0) { //calculating d/dx
    int pin, pout;
    double*** cx = (double***)Cx.GetDataPointer();

    // loop through domain interior
    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++) { 

          for(int id=0; id<(int)input_dof.size(); id++) {
            pin  = input_dof[id];
            pout = output_dof[id];
            dv[k][j][i*DOFout+pout] = cx[k][j][i*4]*v[k][j][(i-q0)*DOFin+pin]
                                    + cx[k][j][i*4+1]*v[k][j][(i-q0+1)*DOFin+pin]
                                    + cx[k][j][i*4+2]*v[k][j][(i-q0+2)*DOFin+pin]
                                    + cx[k][j][i*4+3]*v[k][j][(i-q0+3)*DOFin+pin];
          }
        }

    Cx.RestoreDataPointerToLocalVector();
  }  

  else if(dir==1) { //calculating d/dy
    int pin, pout;
    double*** cy = (double***)Cy.GetDataPointer();

    // loop through domain interior 
    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++) {

          for(int id=0; id<(int)input_dof.size(); id++) {
            pin  = input_dof[id];
            pout = output_dof[id];
            dv[k][j][i*DOFout+pout] = cy[k][j][i*4]*v[k][j-q0][i*DOFin+pin]
                                    + cy[k][j][i*4+1]*v[k][j-q0+1][i*DOFin+pin]
                                    + cy[k][j][i*4+2]*v[k][j-q0+2][i*DOFin+pin]
                                    + cy[k][j][i*4+3]*v[k][j-q0+3][i*DOFin+pin];
          }
        }

    Cy.RestoreDataPointerToLocalVector();
  }

  else if(dir==2) { //calculating d/dz
    int pin, pout;
    double*** cz = (double***)Cz.GetDataPointer();

    // loop through domain interior 
    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++) {

          for(int id=0; id<(int)input_dof.size(); id++) {
            pin  = input_dof[id];
            pout = output_dof[id];
            dv[k][j][i*DOFout+pout] = cz[k][j][i*4]*v[k-q0][j][i*DOFin+pin]
                                    + cz[k][j][i*4+1]*v[k-q0+1][j][i*DOFin+pin]
                                    + cz[k][j][i*4+2]*v[k-q0+2][j][i*DOFin+pin]
                                    + cz[k][j][i*4+3]*v[k-q0+3][j][i*DOFin+pin];
          }
        }

    Cz.RestoreDataPointerToLocalVector();
  }

  // merge data 
  DV.RestoreDataPointerAndInsert();

  V.RestoreDataPointerToLocalVector();
}

//--------------------------------------------------------------------------

void 
GradientCalculatorFD3::CalculateFirstDerivativeAtSelectedNodes(int dir/*0~d/dx,1~d/dy,2~d/dz*/, 
                           vector<Int3> &nodes, SpaceVariable3D &V, std::vector<int> &input_dof,
                           SpaceVariable3D &DV, std::vector<int> &output_dof)
{
  //check for common errors
  if(input_dof.size() != output_dof.size()) {
    print_error("*** Error: CalculateFirstDerivativeAtSelectedNodes failed due to inconsistent number of dof's (%d vs. %d)\n",
                input_dof.size(), output_dof.size());
    exit_mpi();
  }
  if(*std::min_element(input_dof.begin(), input_dof.end())<0 ||
     *std::max_element(input_dof.begin(), input_dof.end())>=V.NumDOF()) {
    print_error("*** Error: CalculateFirstDerivativeAtSelectedNodes failed due to incorrect dof number(s) (min:%d, max:%d).\n",
                *std::min_element(input_dof.begin(), input_dof.end()),
                *std::max_element(input_dof.begin(), input_dof.end()));
    exit_mpi();
  }
  if(*std::min_element(output_dof.begin(), output_dof.end())<0 ||
     *std::max_element(output_dof.begin(), output_dof.end())>=DV.NumDOF()) {
    print_error("*** Error: CalculateFirstDerivativeAtSelectedNodes failed due to incorrect dof number(s) (min:%d, max:%d).\n",
                *std::min_element(output_dof.begin(), output_dof.end()),
                *std::max_element(output_dof.begin(), output_dof.end()));
    exit_mpi();
  }
  if(V.NumGhostLayers()<2) {
    print_error("*** Error: CalculateFirstDerivativeAtSelectedNodes(FD3) failed. Input variable should have at least 2 ghost layers.\n");
    exit_mpi();
  }

  // Get data
  double*** v = (double***)V.GetDataPointer();
  double*** dv = (double***)DV.GetDataPointer(); 

  int DOFin = V.NumDOF(), DOFout = DV.NumDOF();

  int q0;
  assert(shift==-1 || shift==1);
  if(shift == -1)
    q0 = 2;
  else
    q0 = 1;
  
    
  //-------------------------------------------------------------
  // Loop through the domain interior 
  //-------------------------------------------------------------

  if(dir==0) { //calculating d/dx
    int pin, pout;
    double*** cx = (double***)Cx.GetDataPointer();

    for(auto it = nodes.begin(); it != nodes.end(); it++) {

      int i((*it)[0]), j((*it)[1]), k((*it)[2]); 
      if(!coordinates.IsHere(i,j,k,false)) // only work inside the subdomain
        continue;

      for(int id=0; id<(int)input_dof.size(); id++) {
        pin  = input_dof[id];
        pout = output_dof[id];
        dv[k][j][i*DOFout+pout] = cx[k][j][i*4]*v[k][j][(i-q0)*DOFin+pin]
                                + cx[k][j][i*4+1]*v[k][j][(i-q0+1)*DOFin+pin]
                                + cx[k][j][i*4+2]*v[k][j][(i-q0+2)*DOFin+pin]
                                + cx[k][j][i*4+3]*v[k][j][(i-q0+3)*DOFin+pin];
      }
    }

    Cx.RestoreDataPointerToLocalVector();
  }  

  else if(dir==1) { //calculating d/dy
    int pin, pout;
    double*** cy = (double***)Cy.GetDataPointer();

    for(auto it = nodes.begin(); it != nodes.end(); it++) {

      int i((*it)[0]), j((*it)[1]), k((*it)[2]); 
      if(!coordinates.IsHere(i,j,k,false)) // only work inside the subdomain
        continue;

      for(int id=0; id<(int)input_dof.size(); id++) {
        pin  = input_dof[id];
        pout = output_dof[id];
        dv[k][j][i*DOFout+pout] = cy[k][j][i*4]*v[k][j-q0][i*DOFin+pin]
                                + cy[k][j][i*4+1]*v[k][j-q0+1][i*DOFin+pin]
                                + cy[k][j][i*4+2]*v[k][j-q0+2][i*DOFin+pin]
                                + cy[k][j][i*4+3]*v[k][j-q0+3][i*DOFin+pin];
      }
    }

    Cy.RestoreDataPointerToLocalVector();
  }

  else if(dir==2) { //calculating d/dz
    int pin, pout;
    double*** cz = (double***)Cz.GetDataPointer();

    for(auto it = nodes.begin(); it != nodes.end(); it++) {

      int i((*it)[0]), j((*it)[1]), k((*it)[2]); 
      if(!coordinates.IsHere(i,j,k,false)) // only work inside the subdomain
        continue;

      for(int id=0; id<(int)input_dof.size(); id++) {
        pin  = input_dof[id];
        pout = output_dof[id];
        dv[k][j][i*DOFout+pout] = cz[k][j][i*4]*v[k-q0][j][i*DOFin+pin]
                                + cz[k][j][i*4+1]*v[k-q0+1][j][i*DOFin+pin]
                                + cz[k][j][i*4+2]*v[k-q0+2][j][i*DOFin+pin]
                                + cz[k][j][i*4+3]*v[k-q0+3][j][i*DOFin+pin];
      }
    }

    Cz.RestoreDataPointerToLocalVector();
  }

  // merge data 
  DV.RestoreDataPointerAndInsert();

  V.RestoreDataPointerToLocalVector();
}

//--------------------------------------------------------------------------

void GradientCalculatorFD3::CalculateCoefficientsLeftBiased()
{
  // Calculate the coefficients Cx0, Cx1, Cx2, Cx3, Cy0, ..., which are for
  // calculating the derivative at cell centers using finite differencing,
  // For example, the partial derivative du/dx at (i,j,k) should be computed as
  // du/dx(i,j,k) = Cx[k][j][i][0]*u[k][j][i-2] + Cx[k][j][i][1]*u[k][j][i-1]
  //              + Cx[k][j][i][2]*u[k][j][i] + Cx[k][j][i][3]*u[k][j][i+1]
  
  Vec3D*** coordsG2 = (Vec3D***)coordinatesG2.GetDataPointer();

  double*** cx = (double***)Cx.GetDataPointer();
  double*** cy = (double***)Cy.GetDataPointer();
  double*** cz = (double***)Cz.GetDataPointer();

  //Loop through the real cells
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        if(i>0) { //third-order
          double dx0 = coordsG2[k][j][i-2][0] - coordsG2[k][j][i][0];
          double dx1 = coordsG2[k][j][i-1][0] - coordsG2[k][j][i][0];
          double dx2 = 0.0; //coordsG2[k][j][i][0] - coordsG2[k][j][i][0]
          double dx3 = coordsG2[k][j][i+1][0] - coordsG2[k][j][i][0];
          CalculateCoefficientsByTaylerSeries(dx0, dx1, dx2, dx3, &cx[k][j][i*4]); 
        } 
        else { //first-order
          double dx = coordsG2[k][j][i][0] - coordsG2[k][j][i-1][0];
          cx[k][j][i*4]   = 0.0;
          cx[k][j][i*4+1] = -1.0/dx;
          cx[k][j][i*4+2] = 1.0/dx;
          cx[k][j][i*4+3] = 0.0;
        } 

        if(j>0) {
          double dy0 = coordsG2[k][j-2][i][1] - coordsG2[k][j][i][1];
          double dy1 = coordsG2[k][j-1][i][1] - coordsG2[k][j][i][1];
          double dy2 = 0.0; 
          double dy3 = coordsG2[k][j+1][i][1] - coordsG2[k][j][i][1];
          CalculateCoefficientsByTaylerSeries(dy0, dy1, dy2, dy3, &cy[k][j][i*4]);
        } 
        else { 
          double dy = coordsG2[k][j][i][1] - coordsG2[k][j-1][i][1];
          cy[k][j][i*4]   = 0.0;
          cy[k][j][i*4+1] = -1.0/dy;
          cy[k][j][i*4+2] = 1.0/dy;
          cy[k][j][i*4+3] = 0.0;
        } 

        if(k>0) {
          double dz0 = coordsG2[k-2][j][i][2] - coordsG2[k][j][i][2];
          double dz1 = coordsG2[k-1][j][i][2] - coordsG2[k][j][i][2];
          double dz2 = 0.0; 
          double dz3 = coordsG2[k+1][j][i][2] - coordsG2[k][j][i][2];
          CalculateCoefficientsByTaylerSeries(dz0, dz1, dz2, dz3, &cz[k][j][i*4]);
        } 
        else { 
          double dz = coordsG2[k][j][i][2] - coordsG2[k-1][j][i][2];
          cz[k][j][i*4]   = 0.0;
          cz[k][j][i*4+1] = -1.0/dz;
          cz[k][j][i*4+2] = 1.0/dz;
          cz[k][j][i*4+3] = 0.0;
        } 

      }

  Cx.RestoreDataPointerAndInsert();
  Cy.RestoreDataPointerAndInsert();
  Cz.RestoreDataPointerAndInsert();

  coordinatesG2.RestoreDataPointerToLocalVector();

}

//--------------------------------------------------------------------------

void GradientCalculatorFD3::CalculateCoefficientsRightBiased()
{
  // Calculate the coefficients Cx0, Cx1, Cx2, Cx3, Cy0, ..., which are for
  // calculating the derivative at cell centers using finite differencing,
  // For example, the partial derivative du/dx at (i,j,k) should be computed as
  // du/dx(i,j,k) = Cx[k][j][i][0]*u[k][j][i-1] + Cx[k][j][i][1]*u[k][j][i]
  //              + Cx[k][j][i][2]*u[k][j][i+1] + Cx[k][j][i][3]*u[k][j][i+2]
  
  int NX, NY, NZ;
  coordinates.GetGlobalSize(&NX, &NY, &NZ);

  Vec3D*** coordsG2 = (Vec3D***)coordinatesG2.GetDataPointer();

  double*** cx = (double***)Cx.GetDataPointer();
  double*** cy = (double***)Cy.GetDataPointer();
  double*** cz = (double***)Cz.GetDataPointer();

  //Loop through the real cells
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        if(i<NX-1) { //third-order
          double dx0 = coordsG2[k][j][i-1][0] - coordsG2[k][j][i][0];
          double dx1 = 0.0; //coordsG2[k][j][i][0] - coordsG2[k][j][i][0]
          double dx2 = coordsG2[k][j][i+1][0] - coordsG2[k][j][i][0];
          double dx3 = coordsG2[k][j][i+2][0] - coordsG2[k][j][i][0];
          CalculateCoefficientsByTaylerSeries(dx0, dx1, dx2, dx3, &cx[k][j][i*4]);
        } 
        else { //first-order
          double dx = coordsG2[k][j][i+1][0] - coordsG2[k][j][i][0];
          cx[k][j][i*4]   = 0.0;
          cx[k][j][i*4+1] = -1.0/dx;
          cx[k][j][i*4+2] = 1.0/dx;
          cx[k][j][i*4+3] = 0.0;
        } 

        if(i<NY-1) {
          double dy0 = coordsG2[k][j-1][i][1] - coordsG2[k][j][i][1];
          double dy1 = 0.0; 
          double dy2 = coordsG2[k][j+1][i][1] - coordsG2[k][j][i][1];
          double dy3 = coordsG2[k][j+2][i][1] - coordsG2[k][j][i][1];
          CalculateCoefficientsByTaylerSeries(dy0, dy1, dy2, dy3, &cy[k][j][i*4]);
        } 
        else { 
          double dy = coordsG2[k][j+1][i][1] - coordsG2[k][j][i][1];
          cy[k][j][i*4]   = 0.0;
          cy[k][j][i*4+1] = -1.0/dy;
          cy[k][j][i*4+2] = 1.0/dy;
          cy[k][j][i*4+3] = 0.0;
        } 

        if(k<NZ-1) {
          double dz0 = coordsG2[k-1][j][i][2] - coordsG2[k][j][i][2];
          double dz1 = 0.0; 
          double dz2 = coordsG2[k+1][j][i][2] - coordsG2[k][j][i][2];
          double dz3 = coordsG2[k+2][j][i][2] - coordsG2[k][j][i][2];
          CalculateCoefficientsByTaylerSeries(dz0, dz1, dz2, dz3, &cz[k][j][i*4]);
        } 
        else { 
          double dz = coordsG2[k+1][j][i][2] - coordsG2[k][j][i][2];
          cz[k][j][i*4]   = 0.0;
          cz[k][j][i*4+1] = -1.0/dz;
          cz[k][j][i*4+2] = 1.0/dz;
          cz[k][j][i*4+3] = 0.0;
        } 
      }

  Cx.RestoreDataPointerAndInsert();
  Cy.RestoreDataPointerAndInsert();
  Cz.RestoreDataPointerAndInsert();

  coordinatesG2.RestoreDataPointerToLocalVector();

}


//--------------------------------------------------------------------------

void GradientCalculatorFD3::CalculateCoefficientsByTaylerSeries(double dx0, double dx1, double dx2, double dx3, double *coeffs)
{
  Eigen::Matrix4d A;
  Eigen::Vector4d b;
  A << 1.0,         1.0,         1.0,         1.0,
       dx0,         dx1,         dx2,         dx3,
       dx0*dx0,     dx1*dx1,     dx2*dx2,     dx3*dx3,
       dx0*dx0*dx0, dx1*dx1*dx1, dx2*dx2*dx2, dx3*dx3*dx3;
  b << 0.0, 1.0, 0.0, 0.0;

  Eigen::Vector4d x = A.colPivHouseholderQr().solve(b); //solve Ax = b

  for(int i=0; i<4; i++)
    coeffs[i] = x(i);
}

//--------------------------------------------------------------------------


