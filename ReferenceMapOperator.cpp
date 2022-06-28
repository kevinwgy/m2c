#include<ReferenceMapOperator.h>

//-------------------------------------------------------------------

ReferenceMapOperator::ReferenceMapOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_,
                                           SpaceVariable3D &coordinates_,
                                           SpaceVariable3D &delta_xyz_, 
                                           std::vector<GhostPoint> &ghost_nodes_inner_,
                                           std::vector<GhostPoint> &ghost_nodes_outer_)
                    : comm(comm_), iod(iod_), coordinates(coordinates_),
                      delta_xyz(delta_xyz_), ghost_nodes_inner(ghost_nodes_inner_),
                      ghost_nodes_outer(ghost_nodes_outer_)
{

  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);

  //set up gradient calculator
  if(iod.refmap.fd == ReferenceMapData::UPWIND_CENTRAL_3) { //currently the only option
    grad_minus = new GradientCalculatorFD3(comm, dm_all_, coordinates, delta_xyz, -1);
    grad_plus  = new GradientCalculatorFD3(comm, dm_all_, coordinates, delta_xyz,  1);
  }

}

//-------------------------------------------------------------------

ReferenceMapOperator::~ReferenceMapOperator()
{
  if(grad_minus) delete grad_minus;
  if(grad_plus)  delete grad_plus;
}

//-------------------------------------------------------------------
// Initialize Xi to coords
void
ReferenceMapOperator::SetInitialCondition(SpaceVariable3D &Xi)
{
  Xi.AXPlusBY(0.0, 1.0, coordinates, true);
}

//-------------------------------------------------------------------

// Apply boundary conditions by populating ghost cells of Xi
void
ReferenceMapOperator::ApplyBoundaryConditions(SpaceVariable3D &Xi, double time)
{

  Vec3D*** xi = (Vec3D***)Xi.GetDataPointer();

}

//-------------------------------------------------------------------

void
ReferenceMapOperator::ComputeResidual(SpaceVariable3D &V, SpaceVariable3D &Xi, SpaceVariable3D &R)
{
  Vec5D*** v   = (Vec5D***) V.GetDataPointer();
  Vec3D*** xi  = (Vec3D***) Xi.GetDataPointer();
  Vec3D*** res = (Vec3D***) R.GetDataPointer();

  //***************************************************************
  // Step 1: Calculate partial derivatives of xi 
  //***************************************************************
  vector<int> ind0{0,1,2};

  Vec3D*** s = scalarG2.GetDataPointer();
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++)
        s[k][j][i] = xi[k][j][i];
  scalarG2.RestoreDataPointerAndInsert(); //need to exchange

  grad_minus->CalculateFirstDerivativeAtNodes(0/*x*/, scalarG2, ind0, Xil, ind0);
  grad_plus->CalculateFirstDerivativeAtNodes(0/*x*/, scalarG2, ind0, Xir, ind0);
  grad_minus->CalculateFirstDerivativeAtNodes(1/*y*/, scalarG2, ind0, Xib, ind0);
  grad_plus->CalculateFirstDerivativeAtNodes(1/*y*/, scalarG2, ind0, Xit, ind0);
  grad_minus->CalculateFirstDerivativeAtNodes(2/*z*/, scalarG2, ind0, Xik, ind0);
  grad_plus->CalculateFirstDerivativeAtNodes(2/*z*/, scalarG2, ind0, Xif, ind0); 


  //***************************************************************
  // Step 2: Loop through active nodes and compute residual
  //***************************************************************
  int NX, NY, NZ;
  coordinates.GetGlobalSize(&NX, &NY, &NZ);

  Vec3D*** xil = Xil.GetDataPointer(); //d(Xi)/dx, left-biased
  Vec3D*** xir = Xir.GetDataPointer(); //d(Xi)/dx, right-biased
  Vec3D*** xib = Xib.GetDataPointer();
  Vec3D*** xit = Xit.GetDataPointer();
  Vec3D*** xik = Xik.GetDataPointer();
  Vec3D*** xif = Xif.GetDataPointer();
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();

  double a, b, c, d, e, f;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        a = (i-2>=-1) ? xil[k][j][i] : (xi[k][j][i]-xi[k][j][i-1])/(coords[k][j][i][0]-coords[k][j][i-1][0]);
        b = (i+2<=NX) ? xir[k][j][i] : (xi[k][j][i+1]-xi[k][j][i])/(coords[k][j][i+1][0]-coords[k][j][i][0]);
        c = (j-2>=-1) ? xib[k][j][i] : (xi[k][j][i]-xi[k][j-1][i])/(coords[k][j][i][1]-coords[k][j-1][i][1]);
        d = (j+2<=NY) ? xit[k][j][i] : (xi[k][j+1][i]-xi[k][j][i])/(coords[k][j+1][i][1]-coords[k][j][i][1]);
        e = (k-2>=-1) ? xik[k][j][i] : (xi[k][j][i]-xi[k-1][j][i])/(coords[k][j][i][2]-coords[k-1][j][i][2]);
        f = (k+2<=NZ) ? xif[k][j][i] : (xi[k+1][j][i]-xi[k][j][i])/(coords[k+1][j][i][2]-coords[k][j][i][2]);

        res[k][j][i] = (v[k][j][i][1]>=0 ? a*v[k][j][i][1] : b*v[k][j][i][1])
                     + (v[k][j][i][2]>=0 ? c*v[k][j][i][2] : d*v[k][j][i][2])
                     + (v[k][j][i][3]>=0 ? e*v[k][j][i][3] : f*v[k][j][i][3]);

        res[k][j][i] *= -1.0; //moves the residual to the RHS

      }


  Xil.RestoreDataPointerToLocalVector();
  Xir.RestoreDataPointerToLocalVector();
  Xib.RestoreDataPointerToLocalVector();
  Xit.RestoreDataPointerToLocalVector();
  Xik.RestoreDataPointerToLocalVector();
  Xif.RestoreDataPointerToLocalVector();
  coordinates.RestoreDataPointerToLocalVector();

  V.RestoreDataPointerToLocalVector();
  Xi.RestoreDataPointerToLocalVector();

  R.RestoreDataPointerAndInsert();

}

//-------------------------------------------------------------------


//-------------------------------------------------------------------




