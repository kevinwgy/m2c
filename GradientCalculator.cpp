#include <GradientCalculator.h>

//--------------------------------------------------------------------------

GradientCalculatorCentral::GradientCalculatorCentral(MPI_Comm &comm_, DataManagers3D &dm_all_,
                               SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_,
                               Interpolator &interpolator_)
                          : GradientCalculatorBase(coordinates_, delta_xyz_),
                            Cx(comm_, &(dm_all_.ghosted1_3dof)),
                            Cy(comm_, &(dm_all_.ghosted1_3dof)),
                            Cz(comm_, &(dm_all_.ghosted1_3dof)),
                            Var(comm_, &(dm_all_.ghosted1_5dof)),
                            interpolator(interpolator_)
{
  Cx.SetConstantValue(0.0,true);
  Cy.SetConstantValue(0.0,true);
  Cz.SetConstantValue(0.0,true);
  Var.SetConstantValue(0.0,true);

  CalculateCoefficients();
}

//--------------------------------------------------------------------------

void GradientCalculatorCentral::CalculateFirstDerivativeAtNodes(int dir/*0~d/dx,1~d/dy,2~d/dz*/,
                                    SpaceVariable3D &V, std::vector<int> &input_dof,
                                    SpaceVariable3D &DV, std::vector<int> &output_dof);
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

  // Get data
  double*** v = (double***)V.GetDataPointer();
  double*** dv = (double***)DV.GetDataPointer(); 
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();

  int DOFin = V.NumDOF(), DOFout = DV.NumDOF();

  //-------------------------------------------------------------
  // Loop through the domain interior and relevant ghost layer
  // The solution is valid in the ghost layer outside the physical domain
  // except edges and corners
  //-------------------------------------------------------------
  int NX, NY, NZ;
  V.GetGlobalSize(&NX, &NY, &NZ);

  if(dir==0) { //calculating d/dx
    double pin, pout;
    double*** cx = (double***)Cx.GetDataPointer();

    // loop through domain interior and the relevant ghost region
    for(int k=kk0; k<kkmax; k++)
      for(int j=jj0; j<jjmax; j++)
        for(int i=i0; i<imax; i++) { //here, we do not visit ii0 and iimax-1

          //We don't interpolate at edges and corners
          if(V.BoundaryType(i,j,k)>=2)
            continue;

          //For example, v[k][-1][-1] is generally invalid (not populated).
          //In this type of cases, switch to one-sided difference (1st order)
          if(i==0 && (j==-1 || j==NY || k==-1 || k==NZ)) {
            for(int id=0; id<input_dof.size(); id++) {
              pin  = input_dof[id];
              pout = output_dof[id];
              dv[k][j][i*DOFout+pout] = (v[k][j][(i+1)*DOFin+pin] - v[k][j][i*DOFin+pin])
                                       /(coords[k][j][i+1][0] - coords[k][j][i][0]);
            }
          }
          else if(i==NX-1 && (j==-1 || j==NY || k==-1 || k==NZ)) {
            for(int id=0; id<input_dof.size(); id++) {
              pin  = input_dof[id];
              pout = output_dof[id];
              dv[k][j][i*DOFout+pout] = (v[k][j][i*DOFin+pin] - v[k][j][(i-1)*DOFin+pin])
                                       /(coords[k][j][i][0] - coords[k][j][i-1][0]);
            }
          }
          else { //the normal central differencing scheme
            for(int id=0; id<input_dof.size(); id++) {
              pin  = input_dof[id];
              pout = output_dof[id];
              dv[k][j][i*DOFout+pout] = cx[k][j][i][0]*v[k][j][(i-1)*DOFin+pin]
                                      + cx[k][j][i][1]*v[k][j][i*DOFin+pin]
                                      + cx[k][j][i][2]*v[k][j][(i+1)*DOFin+pin];
            }
          }

        }
    Cx.RestoreDataPointerToLocalVector();
  }  

  else if(dir==1) { //calculating d/dy
    double pin, pout;
    double*** cy = (double***)Cy.GetDataPointer();

    // loop through domain interior and the relevant ghost region
    for(int k=kk0; k<kkmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=ii0; i<iimax; i++) {

          //We don't interpolate at edges and corners
          if(V.BoundaryType(i,j,k)>=2)
            continue;

          if(j==0 && (i==-1 || i==NX || k==-1 || k==NZ)) {
            for(int id=0; id<input_dof.size(); id++) {
              pin  = input_dof[id];
              pout = output_dof[id];
              dv[k][j][i*DOFout+pout] = (v[k][j+1][i*DOFin+pin] - v[k][j][i*DOFin+pin])
                                       /(coords[k][j+1][i][1] - coords[k][j][i][1]);
            }
          }
          else if(j==NY-1 && (i==-1 || i==NX || k==-1 || k==NZ)) {
            for(int id=0; id<input_dof.size(); id++) {
              pin  = input_dof[id];
              pout = output_dof[id];
              dv[k][j][i*DOFout+pout] = (v[k][j][i*DOFin+pin] - v[k][j-1][i*DOFin+pin])
                                       /(coords[k][j][i][1] - coords[k][j-1][i][1]);
            }
          }
          else { //the normal central differencing scheme
            for(int id=0; id<input_dof.size(); id++) {
              pin  = input_dof[id];
              pout = output_dof[id];
              dv[k][j][i*DOFout+pout] = cy[k][j][i][0]*v[k][j-1][i*DOFin+pin]
                                      + cy[k][j][i][1]*v[k][j][i*DOFin+pin]
                                      + cy[k][j][i][2]*v[k][j+1][i*DOFin+pin];
            }
          }

        }
    Cy.RestoreDataPointerToLocalVector();
  }

  else if(dir==2) { //calculating d/dz
    double pin, pout;
    double*** cz = (double***)Cz.GetDataPointer();

    // loop through domain interior and the relevant ghost region
    for(int k=k0; k<kmax; k++)
      for(int j=jj0; j<jjmax; j++)
        for(int i=ii0; i<iimax; i++) {

          //We don't interpolate at edges and corners
          if(V.BoundaryType(i,j,k)>=2)
            continue;

          if(k==0 && (i==-1 || i==NX || j==-1 || j==NY)) {
            for(int id=0; id<input_dof.size(); id++) {
              pin  = input_dof[id];
              pout = output_dof[id];
              dv[k][j][i*DOFout+pout] = (v[k+1][j][i*DOFin+pin] - v[k][j][i*DOFin+pin])
                                       /(coords[k+1][j][i][2] - coords[k][j][i][2]);
            }
          }
          else if(k==NZ-1 && (i==-1 || i==NX || j==-1 || j==NY)) {
            for(int id=0; id<input_dof.size(); id++) {
              pin  = input_dof[id];
              pout = output_dof[id];
              dv[k][j][i*DOFout+pout] = (v[k][j][i*DOFin+pin] - v[k-1][j][i*DOFin+pin])
                                       /(coords[k][j][i][2] - coords[k-1][j][i][2]);
            }
          }
          else { //the normal central differencing scheme
            for(int id=0; id<input_dof.size(); id++) {
              pin  = input_dof[id];
              pout = output_dof[id];
              dv[k][j][i*DOFout+pout] = cz[k][j][i][0]*v[k-1][j][i*DOFin+pin]
                                      + cz[k][j][i][1]*v[k][j][i*DOFin+pin]
                                      + cz[k][j][i][2]*v[k+1][j][i*DOFin+pin];
            }
          }

        }
    Cz.RestoreDataPointerToLocalVector();
  }

  // merge data 
  DV.RestoreDataPointerAndInsert();

  V.RestoreDataPointerToLocalVector();
  coordinates.RestoreDataPointerToLocalVector();
}

//--------------------------------------------------------------------------
//! calculates x-,y-,or z-derivative at i +/- 1/2, j +/- 1/2, or k +/- 1/2
void GradientCalculatorCentral::CalculateFirstDerivativeAtCellInterfaces(int dir/*0~d/dx,1~d/dy,2~d/dz*/,
                                    int inter/*0~i +/- 1/2, 1~j +/- 1/2, 2~k +/- 1/2*/,
                                    SpaceVariable3D &V, std::vector<int> &input_dof,
                                    SpaceVariable3D &DV, std::vector<int> &output_dof)
{
  //check for common errors
  if(input_dof.size() != output_dof.size()) {
    print_error("*** Error: CalculateFirstDerivativeAtIMinusHalf failed due to inconsistent "
                "number of dof's (%d vs. %d)\n",
                input_dof.size(), output_dof.size());
    exit_mpi();
  }
  if(*std::min_element(input_dof.begin(), input_dof.end())<0 ||
     *std::max_element(input_dof.begin(), input_dof.end())>=V.NumDOF()) {
    print_error("*** Error: CalculateFirstDerivativeAtIMinusHalf failed due to incorrect "
                "dof number(s) (min:%d, max:%d).\n",
                *std::min_element(input_dof.begin(), input_dof.end()),
                *std::max_element(input_dof.begin(), input_dof.end()));
    exit_mpi();
  }
  if(*std::min_element(output_dof.begin(), output_dof.end())<0 ||
     *std::max_element(output_dof.begin(), output_dof.end())>=DV.NumDOF()) {
    print_error("*** Error: CalculateFirstDerivativeAtIMinusHalf failed due to incorrect "
                "dof number(s) (min:%d, max:%d).\n",
                *std::min_element(output_dof.begin(), output_dof.end()),
                *std::max_element(output_dof.begin(), output_dof.end()));
    exit_mpi();
  }

  //call internal functions
  if(dir==inter)
    CentralDifferencingAtCellInterfaces(dir, V, input_dof, DV, output_dof);
  else {
    CalculateFirstDerivativeAtNodes(dir, V, input_dof, Var, input_dof);
    interpolator.InterpolateAtCellInterfaces(inter, Var, input_dof, DV, output_dof);
  }

}

//--------------------------------------------------------------------------

void GradientCalculatorCentral::CalculateCoefficients()
{
  // Calculate the coefficients Cx0, Cx1, Cx2, Cy0, ..., which are for
  // calculating the derivative at cell centers using central differencing,
  // i.e. centered quadratic interpolation + differentiation
  // For example, the partial derivative du/dx at (i,j,k) should be computed as
  // du/dx(i,j,k) = Cx0[k][j][i-1]*u[k][j][i-1] + Cx1[k][j][i]*u[k][j][i]
  //              + Cx2[k][j][i+1]*u[k][j][i+1]
  // In other words, Cx0[k][j][i] stores the coefficient when cell [k][j][i] is
  // used as the left cell.
  
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
  Vec3D*** dxyz   = (Vec3D***)delta_xyz.GetDataPointer();

  Vec3D*** cx = (double***)Cx.GetDataPointer();
  Vec3D*** cy = (double***)Cy.GetDataPointer();
  Vec3D*** cz = (double***)Cz.GetDataPointer();

  //Loop through the real cells
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        //du/dx[k][j][i] = cx[k][j][i][0]*u[k][j][i-1] + cx[k][j][i][1]*u[k][j][i]
        //               + cx[k][j][i][2]*u[k][j][i+1]
        double dx0 = coords[k][j][i][0] - coords[k][j][i-1][0];
        double dx1 = coords[k][j][i+1][0] - coords[k][j][i][0];
        double dx2 = dx0 + dx1;
        cx[k][j][i][0] = -dx1/(dx0*dx2);
        cx[k][j][i][1] = 1.0/dx0 - 1.0/dx1;
        cx[k][j][i][2] = dx0/(dx1*dx2); 

        //du/dy[k][j][i] = cy[k][j][i][0]*u[k][j-1][i] + cy[k][j][i][1]*u[k][j][i]
        //               + cy[k][j][i][2]*u[k][j+1][i]
        double dy0 = coords[k][j][i][1] - coords[k][j-1][i][1];
        double dy1 = coords[k][j+1][i][1] - coords[k][j][i][1];
        double dy2 = dy0 + dy1;
        cy[k][j][i][0] = -dy1/(dy0*dy2);
        cy[k][j][i][1] = 1.0/dy0 - 1.0/dy1;
        cy[k][j][i][2] = dy0/(dy1*dy2); 

        //du/dz[k][j][i] = cz[k][j][i][0]*u[k-1][j][i] + cz[k][j][i][1]*u[k][j][i]
        //               + cz[k][j][i][2]*u[k+1][j][i]
        double dz0 = coords[k][j][i][2] - coords[k-1][j][i][2];
        double dz1 = coords[k+1][j][i][2] - coords[k][j][i][2];
        double dz2 = dz0 + dz1;
        cz[k][j][i][0] = -dz1/(dz0*dz2);
        cz[k][j][i][1] = 1.0/dz0 - 1.0/dz1;
        cz[k][j][i][2] = dz0/(dz1*dz2); 

      }

  Cx.RestoreDataPointerAndInsert();
  Cy.RestoreDataPointerAndInsert();
  Cz.RestoreDataPointerAndInsert();

  coordinates.RestoreDataPointerToLocalVector();
  delta_xyz.RestoreDataPointerToLocalVector();

}

//--------------------------------------------------------------------------
// calculate dV/dx at i +/- 1/2, dV/dy at j +/- 1/2, or dV/dz at k +/- 1/2
void GradientCalculatorCentral::CentralDifferencingAtCellInterfaces(int dir/*0~d/dx,1~d/dy,2~d/dz*/,
                                    SpaceVariable3D &V, std::vector<int> &input_dof,
                                    SpaceVariable3D &DV, std::vector<int> &output_dof)
{
  // Get data
  double*** v = (double***)V.GetDataPointer();
  double*** dv = (double***)DV.GetDataPointer(); 
  int DOFin = V.NumDOF(), DOFout = DV.NumDOF();

  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();

  double pin, pout;

  if(dir==0) {
    for(int k=k0; k<kkmax; k++)
      for(int j=j0; j<jjmax; j++)
        for(int i=i0; i<iimax; i++) {

          if(k==kkmax-1 || j==jjmax-1)
            continue;

          double dx = coords[k][j][i][0] - coords[k][j][i-1][0];
          for(int id=0; id<input_dof.size(); id++) {
            pin  = input_dof[id];
            pout = output_dof[id];
            dv[k][j][i*DOFout+pout] = (v[k][j][i*DOFin+pin] - v[k][j][(i-1)*DOFin+pin])/dx;
          }
        }
  }
  else if(dir==1) {
    for(int k=k0; k<kkmax; k++)
      for(int j=j0; j<jjmax; j++)
        for(int i=i0; i<iimax; i++) {

          if(i==iimax-1 || k==kkmax-1)
            continue;

          double dy = coords[k][j][i][1] - coords[k][j-1][i][1];
          for(int id=0; id<input_dof.size(); id++) {
            pin  = input_dof[id];
            pout = output_dof[id];
            dv[k][j][i*DOFout+pout] = (v[k][j][i*DOFin+pin] - v[k][j-1][i*DOFin+pin])/dy;
          }
        }
  }
  else if(dir==2) {
    for(int k=k0; k<kkmax; k++)
      for(int j=j0; j<jjmax; j++)
        for(int i=i0; i<iimax; i++) {

          if(i==iimax-1 || j==jjmax-1)
            continue;

          double dz = coords[k][j][i][2] - coords[k-1][j][i][2];
          for(int id=0; id<input_dof.size(); id++) {
            pin  = input_dof[id];
            pout = output_dof[id];
            dv[k][j][i*DOFout+pout] = (v[k][j][i*DOFin+pin] - v[k-1][j][i*DOFin+pin])/dz;
          }
        }
  }

  DV.RestoreDataPointerAndInsert();
  V.RestoreDataPointerToLocalVector();
  coordinates.RestoreDataPointerToLocalVector();

}

//--------------------------------------------------------------------------




