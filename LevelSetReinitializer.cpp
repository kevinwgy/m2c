#include <LevelSetReinitializer.h>

extern int verbose;

//--------------------------------------------------------------------------

LevelSetReinitializer::LevelSetReinitializer(MPI_Comm &comm_, DataManagers3D &dm_all_, 
                           LevelSetSchemeData &iod_ls_, SpaceVariable3D &coordinates_, 
                           SpaceVariable3D &delta_xyz_)
                     : comm(comm_), iod_ls(iod_ls_), coordinates(coordinates_),
                       delta_xyz(delta_xyz_), interp(NULL), grad(NULL),
                       Tag(comm_, &(dm_all_.ghosted1_1dof)),
                       R(comm_, &(dm_all_.ghosted1_1dof)),
                       Phi1(comm_, &(dm_all_.ghosted1_1dof)),
                       Sign(comm_, &(dm_all_.ghosted1_1dof)),
                       PhiG2(comm_, &(dm_all_.ghosted2_1dof))
{
  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);

  if(true) //can add other interpolation methods later
    interp = new InterpolatorLinear(comm_, dm_all_, coordinates_, delta_xyz_);
  if(true) //can add other differentiation methods later
    grad = new GradientCalculatorCentral(comm_, dm_all_, coordinates_, delta_xyz_, *interp);
}


//--------------------------------------------------------------------------

LevelSetReinitializer::~LevelSetReinitializer()
{
  if(interp) delete interp;
  if(grad)   delete grad;
}

//--------------------------------------------------------------------------

void
LevelSetReinitializer::Destroy()
{
  if(interp)
    interp->Destroy();
  if(grad)
    grad->Destroy();

  Tag.Destroy();
  R.Destroy();
  Phi1.Destroy();
  Sign.Destroy();
  PhiG2.Destroy();
}

//--------------------------------------------------------------------------

void
LevelSetReinitializer::Reinitialize(SpaceVariable3D &Phi)
{

  // Step 1: Prep: Tag first layer nodes & store the sign function
  TagFirstLayerNodes(Phi); 
  EvaluateSignFunction(Phi);

  // Step 2: Reinitialize first layer nodes (no iterations needed)
  if(iod_ls.reinit.updateFirstLayer == LevelSetReinitializationData::ON) {
    Phi1.AXPlusBY(0.0, 1.0, Phi, true); //Phi1 = Phi
    ReinitializeFirstLayerNodes(Phi1, Phi); //Phi is updated
    ApplyBoundaryConditions(Phi);
  }

  // Step 3: Main loop -- 3rd-order Runge-Kutta w/ spatially varying dt
  double cfl = 0.9;
  double residual = 0.0, residual0 = 0.0;
  int iter;
  for(iter = 0; iter < iod_ls.reinit.maxIts; iter++) {

    //************** Step 1 of RK3 *****************
    residual = ComputeResidual(Phi, R, cfl);  //R = R(Phi)
    if(iter==0)
      residual0 = residual;
    if(verbose>=1)
      print("  o Iter. %d: Residual = %e, Rel. Residual = %e, Tol = %e.\n", iter, residual, residual/residual0, 
            iod_ls.reinit.convergence_tolerance);
    if(residual/residual0 < iod_ls.reinit.convergence_tolerance ||
       residual < iod_ls.reinit.convergence_tolerance) //residual itself is nondimensional
      break;

    Phi1.AXPlusBY(0.0, 1.0, Phi); //Phi1 = Phi
    Phi1.AXPlusBY(1.0, 1.0, R);   //Phi1 = Phi + R(Phi)
    ApplyBoundaryConditions(Phi1);
    //*********************************************


    //************** Step 2 of RK3 *****************
    ComputeResidual(Phi1, R, cfl);
    Phi1.AXPlusBY(0.25, 0.75, Phi);
    Phi1.AXPlusBY(1.0, 0.25, R);
    ApplyBoundaryConditions(Phi1);
    //*********************************************

    //************** Step 3 of RK3 *****************
    ComputeResidual(Phi1, R, cfl);
    Phi.AXPlusBY(1.0/3.0, 2.0/3.0, Phi1);
    Phi.AXPlusBY(1.0, 2.0/3.0, R);
    ApplyBoundaryConditions(Phi);
    //*********************************************
    
  }

  if(iter==iod_ls.reinit.maxIts)
    print("  o Warning: Failed to converge. Residual = %e, Rel. Residual = %e, Tol = %e.\n", 
          residual, residual/residual0, iod_ls.reinit.convergence_tolerance);
}

//--------------------------------------------------------------------------

void
LevelSetReinitializer::TagFirstLayerNodes(SpaceVariable3D &Phi)
{

  firstLayer.resize(0);

  double*** phi   = (double***)Phi.GetDataPointer();
  double*** tag   = (double***)Tag.GetDataPointer();

  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++) {

        tag[k][j][i] = 0; //default

        if(Phi.OutsidePhysicalDomainAndUnpopulated(i,j,k))
          continue;

        if(i-1>=ii0 && phi[k][j][i]*phi[k][j][i-1]<=0) {
          tag[k][j][i] = 1;
          if(Phi.IsHere(i,j,k,false)) //only push nodes in the subdomain interior into firstLayer
            firstLayer.push_back(Int3(i,j,k));
          continue;
        }

        if(i+1<iimax && phi[k][j][i]*phi[k][j][i+1]<=0) {
          tag[k][j][i] = 1;
          if(Phi.IsHere(i,j,k,false)) //only push nodes in the subdomain interior into firstLayer
            firstLayer.push_back(Int3(i,j,k));
          continue;
        }

        if(j-1>=jj0 && phi[k][j][i]*phi[k][j-1][i]<=0) {
          tag[k][j][i] = 1;
          if(Phi.IsHere(i,j,k,false)) //only push nodes in the subdomain interior into firstLayer
            firstLayer.push_back(Int3(i,j,k));
          continue;
        }

        if(j+1<jjmax && phi[k][j][i]*phi[k][j+1][i]<=0) {
          tag[k][j][i] = 1;
          if(Phi.IsHere(i,j,k,false)) //only push nodes in the subdomain interior into firstLayer
            firstLayer.push_back(Int3(i,j,k));
          continue;
        }

        if(k-1>=kk0 && phi[k][j][i]*phi[k-1][j][i]<=0) {
          tag[k][j][i] = 1;
          if(Phi.IsHere(i,j,k,false)) //only push nodes in the subdomain interior into firstLayer
            firstLayer.push_back(Int3(i,j,k));
          continue;
        }
        
        if(k+1<kkmax && phi[k][j][i]*phi[k+1][j][i]<=0) {
          tag[k][j][i] = 1;
          if(Phi.IsHere(i,j,k,false)) //only push nodes in the subdomain interior into firstLayer
            firstLayer.push_back(Int3(i,j,k));
          continue;
        }
      }

  Tag.RestoreDataPointerAndInsert();

  Phi.RestoreDataPointerToLocalVector(); //no changes made
}

//--------------------------------------------------------------------------

void
LevelSetReinitializer::EvaluateSignFunction(SpaceVariable3D &Phi)
{
  double*** phi   = (double***)Phi.GetDataPointer();
  double*** sign  = (double***)Sign.GetDataPointer();

  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++) {
        if(phi[k][j][i]>0)
          sign[k][j][i] = 1.0;
        else if(phi[k][j][i]<0)
          sign[k][j][i] = -1.0;
        else
          sign[k][j][i] = 0.0;
      }

  Phi.RestoreDataPointerToLocalVector();
  Sign.RestoreDataPointerAndInsert(); 
}

//--------------------------------------------------------------------------

void
LevelSetReinitializer::ReinitializeFirstLayerNodes(SpaceVariable3D &Phi0, SpaceVariable3D &Phi)
{

  PopulatePhiG2(Phi0);

  int NX, NY, NZ;
  Phi.GetGlobalSize(&NX, &NY, &NZ);

  double*** phi   = (double***)Phi.GetDataPointer();
  double*** phi0  = (double***)Phi0.GetDataPointer();
  double*** phig  = (double***)PhiG2.GetDataPointer();
  double*** tag   = (double***)Tag.GetDataPointer();
  Vec3D*** dxyz   = (Vec3D***)delta_xyz.GetDataPointer();
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();

  int i,j,k;
  Vec3D gradphi;
  double gradphi_norm;
  double epsx, epsy, epsz;
  for(int n=0; n<firstLayer.size(); n++) {

    // This must be a node in the interior of the subdomain (see how firstLayer is populated)
    i = firstLayer[n][0]; 
    j = firstLayer[n][1]; 
    k = firstLayer[n][2]; 

    // set phi = phi / |grad(phi)|, using tagged nodes to calculate the derivatives 
    // ref: Hartmann et al. 2008, Eqs.(20)(21) (which cites another paper)
    
    //Eq.(21a) of Hartmann et al., 2008
    epsx = 1.0e-3*dxyz[k][j][i][0];
    epsy = 1.0e-3*dxyz[k][j][i][1];
    epsz = 1.0e-3*dxyz[k][j][i][2];
    gradphi[0] = DifferentiateInFirstLayer(coords[k][j][i-1][0], coords[k][j][i][0], coords[k][j][i+1][0],
                                              tag[k][j][i-1],       tag[k][j][i],       tag[k][j][i+1],
                                             phig[k][j][i-1],      phig[k][j][i],      phig[k][j][i+1], 
                                           (i-2>-1) ? phig[k][j][i-2] : phig[k][j][i-1],
                                           (i+2<=NX) ? phig[k][j][i+2] : phig[k][j][i+1], epsx);
    gradphi[1] = DifferentiateInFirstLayer(coords[k][j-1][i][1], coords[k][j][i][1], coords[k][j+1][i][1],
                                              tag[k][j-1][i],       tag[k][j][i],       tag[k][j+1][i],
                                             phig[k][j-1][i],      phig[k][j][i],      phig[k][j+1][i],
                                           (j-2>-1) ? phig[k][j-2][i] : phig[k][j-1][i],
                                           (j+2<=NY) ? phig[k][j+2][i] : phig[k][j+2][i], epsy);
    gradphi[2] = DifferentiateInFirstLayer(coords[k-1][j][i][2], coords[k][j][i][2], coords[k+1][j][i][2],
                                              tag[k-1][j][i],       tag[k][j][i],       tag[k+1][j][i],
                                             phig[k-1][j][i],      phig[k][j][i],      phig[k+1][j][i],
                                           (k-2>-1) ? phig[k-2][j][i] : phig[k-1][j][i],
                                           (k+2<=NZ) ? phig[k+2][j][i] : phig[k+1][j][i], epsz);
    gradphi_norm = gradphi.norm();

    if(gradphi_norm == 0) {
      fprintf(stderr,"Warning: (%d,%d,%d)(%e,%e,%e): Updating first layer node led to zero gradient.\n",
              i,j,k,coords[k][j][i][0],coords[k][j][i][1],coords[k][j][i][2]);
      phi[k][j][i] = phig[k][j][i];
    } else
      phi[k][j][i] = phig[k][j][i]/gradphi_norm;
  }


  Phi.RestoreDataPointerAndInsert();

  Phi0.RestoreDataPointerToLocalVector();
  PhiG2.RestoreDataPointerToLocalVector();
  Tag.RestoreDataPointerToLocalVector();
  coordinates.RestoreDataPointerToLocalVector();
  delta_xyz.RestoreDataPointerToLocalVector();
  
}

//--------------------------------------------------------------------------

void
LevelSetReinitializer::PopulatePhiG2(SpaceVariable3D &Phi0)
{
  PhiG2.SetConstantValue(0.0, true);

  double*** phig2 = (double***)PhiG2.GetDataPointer();
  double*** phi0  = (double***)Phi0.GetDataPointer();

  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++) {
        phig2[k][j][i] = phi0[k][j][i]; 
      }

  Phi0.RestoreDataPointerToLocalVector();
  PhiG2.RestoreDataPointerAndInsert();
}

//--------------------------------------------------------------------------

//Eq.(21a) of Hartmann et al., 2008, simplified
double
LevelSetReinitializer::DifferentiateInFirstLayer(double x0, double x1, double x2, 
                                                 double tag0, double tag1, double tag2,
                                                 double phi0, double phi1, double phi2, 
                                                 double phi00, double phi3, double eps)
{
  bool phi0_useful = true, phi2_useful = true;


  //identify cases in which phi0 is not useful
  if(tag0==0) phi0_useful = false;
  if(tag2==0) phi2_useful = false;
  if(!phi0_useful && !phi2_useful)
    return 0.0;


  double dphi0 = phi1 - phi0;
  double dphi1 = phi2 - phi1;

  bool condA0 = (phi0*phi2<0) && fabs(dphi0 + eps)<fabs(dphi1);
  bool condA2 = (phi0*phi2<0) && fabs(dphi1 + eps)<fabs(dphi0);
  bool condB = dphi0*dphi1<0 || phi0*phi00<0 || phi2*phi3<0;
  
  if(phi0_useful) {
    if(condA0 && condB)
      phi0_useful = false;
  }
  if(phi2_useful) {
    if(condA2 && condB)
      phi2_useful = false;
  }
  if(!phi0_useful && !phi2_useful)
    return 0.0;


  if(phi0_useful) {
    if(phi2_useful) { //central differencing
      double c0 = -(x2-x1)/((x1-x0)*(x2-x0));
      double c1 = 1.0/(x1-x0) - 1.0/(x2-x1);
      double c2 = (x1-x0)/((x2-x0)*(x2-x1));
      return c0*phi0 + c1*phi1 + c2*phi2;
    }
    else  //backward difference
      return (phi1-phi0)/(x1-x0);
  }
  else if(phi2_useful) { //forward difference
    return (phi2-phi1)/(x2-x1);
  }
  else //neigher phi0 nor phi2 useful
    return 0.0;    

  return 0.0;
}

//--------------------------------------------------------------------------

double
LevelSetReinitializer::ComputeResidual(SpaceVariable3D &Phi, SpaceVariable3D &R, double cfl)
{

  int NX, NY, NZ;
  Phi.GetGlobalSize(&NX, &NY, &NZ);

  // get data
  double*** tag    = (double***)Tag.GetDataPointer();
  double*** sign   = (double***)Sign.GetDataPointer();
  double*** phi    = (double***)Phi.GetDataPointer();
  Vec3D***  dxyz   = (Vec3D***)delta_xyz.GetDataPointer();
  Vec3D***  coords = (Vec3D***)coordinates.GetDataPointer();
  double*** res    = (double***)R.GetDataPointer();

  // loop through the interior of the subdomain
  double max_residual = 0.0, dt, dx, local_res;
  double a,b,c,d,e,f, ap,am, bp,bm, cp,cm, dp,dm, ep,em, fp,fm;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        
        if(tag[k][j][i]!=0) {//fixed node (first layer)
          res[k][j][i] = 0.0;
          continue;
        }

        dx = std::min(dxyz[k][j][i][0], dxyz[k][j][i][1]);
        //dx = std::min(dxyz[k][j][i][0], std::min(dxyz[k][j][i][1], dxyz[k][j][i][2]));
        dt = cfl*dx;
        
        a = (phi[k][j][i]-phi[k][j][i-1])/(coords[k][j][i][0]-coords[k][j][i-1][0]);
        b = (phi[k][j][i+1]-phi[k][j][i])/(coords[k][j][i+1][0]-coords[k][j][i][0]);
        c = (phi[k][j][i]-phi[k][j-1][i])/(coords[k][j][i][1]-coords[k][j-1][i][1]);
        d = (phi[k][j+1][i]-phi[k][j][i])/(coords[k][j+1][i][1]-coords[k][j][i][1]);
        e = (phi[k][j][i]-phi[k-1][j][i])/(coords[k][j][i][2]-coords[k-1][j][i][2]);
        f = (phi[k+1][j][i]-phi[k][j][i])/(coords[k+1][j][i][2]-coords[k][j][i][2]);
        ap = std::max(a,0.0), am = std::min(a,0.0);
        bp = std::max(b,0.0), bm = std::min(b,0.0);
        cp = std::max(c,0.0), cm = std::min(c,0.0);
        dp = std::max(d,0.0), dm = std::min(d,0.0);
        ep = std::max(e,0.0), em = std::min(e,0.0);
        fp = std::max(f,0.0), fm = std::min(f,0.0);

        if(phi[k][j][i]>=0) {
          local_res = sqrt(std::max(ap*ap, bm*bm) + std::max(cp*cp, dm*dm) + std::max(ep*ep, fm*fm)) - 1.0;
        } else {
          local_res = sqrt(std::max(am*am, bp*bp) + std::max(cm*cm, dp*dp) + std::max(em*em, fp*fp)) - 1.0;
        }

        res[k][j][i] = -dt*sign[k][j][i]*local_res;

        max_residual = std::max(max_residual, fabs(local_res));
      }

  // get max residual over all the processors
  MPI_Allreduce(MPI_IN_PLACE, &max_residual, 1, MPI_DOUBLE, MPI_MAX, comm);

  // restore data
  Tag.RestoreDataPointerToLocalVector();
  Sign.RestoreDataPointerToLocalVector();
  Phi.RestoreDataPointerToLocalVector();
  delta_xyz.RestoreDataPointerToLocalVector();
  coordinates.RestoreDataPointerToLocalVector();

  R.RestoreDataPointerAndInsert();

  return max_residual;
}

//--------------------------------------------------------------------------

// Apply boundary conditions by populating ghost cells of Phi
// TODO:KW(07/2021) IDENTICAL TO THE FUNCTION IN LEVELSETOPERATOR!! (Should do it in a better way)
void
LevelSetReinitializer::ApplyBoundaryConditions(SpaceVariable3D &Phi)
{
  double*** phi = (double***) Phi.GetDataPointer();
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();

  int NX, NY, NZ;
  Phi.GetGlobalSize(&NX, &NY, &NZ);

  double r, r1, r2, f1, f2;
  //! Left boundary
  if(ii0==-1) { 
    switch (iod_ls.bc_x0) {

      case LevelSetSchemeData::ZERO_NEUMANN :
        for(int k=k0; k<kmax; k++)
          for(int j=j0; j<jmax; j++)
            phi[k][j][ii0] = phi[k][j][ii0+1];
        break;

      case LevelSetSchemeData::LINEAR_EXTRAPOLATION :
        for(int k=k0; k<kmax; k++)
          for(int j=j0; j<jmax; j++) {
            if(ii0+2<NX) { //ii0+2 is within physical domain
              r  = coords[k][j][ii0][0];
              r1 = coords[k][j][ii0+1][0];  f1 = phi[k][j][ii0+1];
              r2 = coords[k][j][ii0+2][0];  f2 = phi[k][j][ii0+2];
              phi[k][j][ii0] = f1 + (f2-f1)/(r2-r1)*(r-r1);
            } else
              phi[k][j][ii0] = phi[k][j][ii0+1];
          }
        break;

      case LevelSetSchemeData::NONE :
        break;

      default :
        print_error("*** Error: Level set boundary condition at x=x0 cannot be specified!\n");
        exit_mpi();
    }
  }

  //! Right boundary
  if(iimax==NX+1) { 
    switch (iod_ls.bc_xmax) {

      case LevelSetSchemeData::ZERO_NEUMANN :
        for(int k=k0; k<kmax; k++)
          for(int j=j0; j<jmax; j++)
            phi[k][j][iimax-1] = phi[k][j][iimax-2];
        break;

      case LevelSetSchemeData::LINEAR_EXTRAPOLATION :
        for(int k=k0; k<kmax; k++)
          for(int j=j0; j<jmax; j++) {
            if(iimax-3>=0) { //iimax-3 is within physical domain
              r  = coords[k][j][iimax-1][0];
              r1 = coords[k][j][iimax-2][0];  f1 = phi[k][j][iimax-2];
              r2 = coords[k][j][iimax-3][0];  f2 = phi[k][j][iimax-3];
              phi[k][j][iimax-1] = f1 + (f2-f1)/(r2-r1)*(r-r1);
            } else
              phi[k][j][iimax-1] = phi[k][j][iimax-2];
          }
        break;
 
      case LevelSetSchemeData::NONE :
        break;

      default :
        print_error("*** Error: Level set boundary condition at x=xmax cannot be specified!\n");
        exit_mpi();
    }
  }

  //! Bottom boundary
  if(jj0==-1) { 
    switch (iod_ls.bc_y0) {

      case LevelSetSchemeData::ZERO_NEUMANN :
        for(int k=k0; k<kmax; k++)
          for(int i=i0; i<imax; i++)
            phi[k][jj0][i] = phi[k][jj0+1][i];
        break;

      case LevelSetSchemeData::LINEAR_EXTRAPOLATION :
        for(int k=k0; k<kmax; k++)
          for(int i=i0; i<imax; i++) {
            if(jj0+2<NY) { //jj0+2 is within physical domain
              r  = coords[k][jj0][i][1];
              r1 = coords[k][jj0+1][i][1];  f1 = phi[k][jj0+1][i];
              r2 = coords[k][jj0+2][i][1];  f2 = phi[k][jj0+2][i];
              phi[k][jj0][i] = f1 + (f2-f1)/(r2-r1)*(r-r1);
            } else
              phi[k][jj0][i] = phi[k][jj0+1][i];
          }
        break;

      case LevelSetSchemeData::NONE :
        break;

      default :
        print_error("*** Error: Level set boundary condition at y=y0 cannot be specified!\n");
        exit_mpi();
    }
  }

  //! Top boundary
  if(jjmax==NY+1) { 
    switch (iod_ls.bc_ymax) {

      case LevelSetSchemeData::ZERO_NEUMANN :
        for(int k=k0; k<kmax; k++)
          for(int i=i0; i<imax; i++)
            phi[k][jjmax-1][i] = phi[k][jjmax-2][i]; 
        break;

      case LevelSetSchemeData::LINEAR_EXTRAPOLATION :
        for(int k=k0; k<kmax; k++)
          for(int i=i0; i<imax; i++) {
            if(jjmax-3>=0) { //jjmax-3 is within physical domain
              r  = coords[k][jjmax-1][i][1];
              r1 = coords[k][jjmax-2][i][1];  f1 = phi[k][jjmax-2][i];
              r2 = coords[k][jjmax-3][i][1];  f2 = phi[k][jjmax-3][i];
              phi[k][jjmax-1][i] = f1 + (f2-f1)/(r2-r1)*(r-r1);
            } else
              phi[k][jjmax-1][i] = phi[k][jjmax-2][i];
          }
        break;
 
      case LevelSetSchemeData::NONE :
        break;

      default :
        print_error("*** Error: Level set boundary condition at y=ymax cannot be specified!\n");
        exit_mpi();
    }
  }

  //! Back boundary (z min)
  if(kk0==-1) { 
    switch (iod_ls.bc_z0) {

      case LevelSetSchemeData::ZERO_NEUMANN :
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++)
            phi[kk0][j][i] = phi[kk0+1][j][i]; 
        break;

      case LevelSetSchemeData::LINEAR_EXTRAPOLATION :
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {
            if(kk0+2<NZ) { //kk0+2 is within physical domain
              r  = coords[kk0][j][i][2];
              r1 = coords[kk0+1][j][i][2];  f1 = phi[kk0+1][j][i];
              r2 = coords[kk0+2][j][i][2];  f2 = phi[kk0+2][j][i];
              phi[kk0][j][i] = f1 + (f2-f1)/(r2-r1)*(r-r1);
            } else
              phi[kk0][j][i] = phi[kk0+1][j][i]; 
          }
        break;

      case LevelSetSchemeData::NONE :
        break;

      default :
        print_error("*** Error: Level set boundary condition at z=z0 cannot be specified!\n");
        exit_mpi();
    }
  }

  //! Front boundary (z max)
  if(kkmax==NZ+1) { 
    switch (iod_ls.bc_zmax) {

      case LevelSetSchemeData::ZERO_NEUMANN :
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++)
            phi[kkmax-1][j][i] = phi[kkmax-2][j][i];
        break;

      case LevelSetSchemeData::LINEAR_EXTRAPOLATION :
        for(int j=j0; j<jmax; j++)
          for(int i=i0; i<imax; i++) {
            if(kkmax-3>=0) { //kkmax-3 is within physical domain
              r  = coords[kkmax-1][j][i][2];
              r1 = coords[kkmax-2][j][i][2];  f1 = phi[kkmax-2][j][i];
              r2 = coords[kkmax-3][j][i][2];  f2 = phi[kkmax-3][j][i];
              phi[kkmax-1][j][i] = f1 + (f2-f1)/(r2-r1)*(r-r1);
            } else
              phi[kkmax-1][j][i] = phi[kkmax-2][j][i];
          }
        break;

      case LevelSetSchemeData::NONE :
        break;

      default :
        print_error("*** Error: Boundary condition at z=zmax cannot be specified!\n");
        exit_mpi();
    }
  }

  Phi.RestoreDataPointerAndInsert();

  coordinates.RestoreDataPointerToLocalVector();
}

//--------------------------------------------------------------------------


//--------------------------------------------------------------------------












