/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<TimeIntegratorSemiImp.h>

//----------------------------------------------------------------------------
// SIMPLE 
//----------------------------------------------------------------------------

TimeIntegratorSIMPLE::TimeIntegratorSIMPLE(MPI_Comm &comm_, IoData& iod_, DataManagers3D& dms_,
                                           SpaceOperator& spo_, IncompressibleOperator &inco_,
                                           vector<LevelSetOperator*>& lso_, MultiPhaseOperator &mpo_,
                                           LaserAbsorptionSolver* laser_, EmbeddedBoundaryOperator* embed_,
                                           HyperelasticityOperator* heo_, PrescribedMotionOperator* pmo_)
                    : TimeIntegratorBase(comm_, iod_, dms_, spo_, lso_, mpo_, laser_, embed_, heo_, pmo_),
                      inco(inco_), Homo(comm_, &(dms_.ghosted1_1dof)), VXstar(comm_, &(dms_.ghosted1_1dof)),
                      VYstar(comm_, &(dms_.ghosted1_1dof)), VZstar(comm_, &(dms_.ghosted1_1dof)),
                      Pprime(comm_, &(dms_.ghosted1_1dof)), B(comm_, &(dms_.ghosted1_1dof)),
                      DX(comm_, &(dms_.ghosted1_1dof)), DY(comm_, &(dms_.ghosted1_1dof)),
                      DZ(comm_, &(dms_.ghosted1_1dof)),
                      vlin_solver(comm_, dms_.ghosted1_1dof, iod.ts.semi_impl.velocity_linear_solver,
                      "velocity"),
                      plin_solver(comm_, dms_.ghosted1_1dof, iod.ts.semi_impl.pressure_linear_solver,
                      "pressure")
{
  type = SIMPLE;

  Homo.SetConstantValue(1, true); //default

  if(iod.ts.semi_impl.E<=0.0) {
    print_error("*** Error: In the SIMPLE family of methods, E must be set to a positive value.\n");
    exit_mpi();
  }
  Efactor = iod.ts.semi_impl.E;
  
  if(iod.ts.semi_impl.alphaP<=0.0) {
    print_error("*** Error: In the SIMPLE family of methods, alphaP must be set to a positive value "
                "(usually less than 1).\n");
    exit_mpi();
  }
  alphaP = iod.ts.semi_impl.alphaP;

  ijk_zero_p = FindCornerFixedPressure();


  // screen outputs
  print("- Setting up a semi-implicit time integrator.\n");
  string ksp_type, pc_type;
  vlin_solver.GetSolverType(&ksp_type, &pc_type);
  print("  o Linear solver for velocity: %s, Preconditioner: %s.\n",
        ksp_type.c_str(), pc_type.c_str());
  plin_solver.GetSolverType(&ksp_type, &pc_type);
  print("  o Linear solver for pressure: %s, Preconditioner: %s.\n",
        ksp_type.c_str(), pc_type.c_str());
}

//----------------------------------------------------------------------------

TimeIntegratorSIMPLE::~TimeIntegratorSIMPLE()
{ }

//----------------------------------------------------------------------------

void
TimeIntegratorSIMPLE::Destroy()
{
  VXstar.Destroy();
  VYstar.Destroy();
  VZstar.Destroy();
  Pprime.Destroy();
  B.Destroy();
  Homo.Destroy();
  DX.Destroy();
  DY.Destroy();
  DZ.Destroy();

  vlin_solver.Destroy();
  plin_solver.Destroy();

  TimeIntegratorBase::Destroy();
}

//----------------------------------------------------------------------------

void
TimeIntegratorSIMPLE::AdvanceOneTimeStep(SpaceVariable3D &V, SpaceVariable3D &ID,
                                         vector<SpaceVariable3D*>& Phi, vector<SpaceVariable3D*> &NPhi,
                                         vector<SpaceVariable3D*> &KappaPhi,
                                         SpaceVariable3D *L, SpaceVariable3D *Xi, SpaceVariable3D *LocalDt,
                                         [[maybe_unused]] double time, double dt, 
                                         int time_step, int subcycle, double dts)
{

  if(mpo.NumberOfMaterials()>1) {
    print_error("*** Error: Need to update homogeneity. Currently, the incompressible flow solver does not allow"
                " more than one material.\n");
    exit_mpi();
  }
  if(Phi.size()>0 || NPhi.size()>0 || KappaPhi.size()>0 || L || Xi || subcycle>0 || dts != dt) {
    print_error("*** Error: Problem setup is not supported by TimeIntegratorSIMPLE(or SIMPLEC).\n");
    exit_mpi();
  }


  GlobalMeshInfo &global_mesh(spo.GetGlobalMeshInfo());

  double*** id = ID.GetDataPointer();
  double*** homo = Homo.GetDataPointer();

  int iter, maxIter = time_step == 1 ? 10*iod.ts.semi_impl.maxIts : iod.ts.semi_impl.maxIts;
  vector<double> lin_rnorm; 
  bool lin_success, converged(false);
  double rel_err(10000.0);

  if(type == SIMPLEC)
    print("  o Running the iterative SIMPLEC procedure (E = %e).\n", Efactor);
  else
    print("  o Running the iterative SIMPLE procedure (E = %e, alphaP = %e).\n", Efactor, alphaP);

  for(iter = 0; iter < maxIter; iter++) {

    Vec5D*** v = (Vec5D***)V.GetDataPointer();

    print("Here.\n");
    ExtractVariableComponents(v, &VXstar, &VYstar, &VZstar, NULL);

    //-----------------------------------------------------
    // Step 1: Solve the momentum equations for u*, v*, w*
    //-----------------------------------------------------

    print("Here 2.\n");
    // Solve the x-momentum equation
    inco.BuildVelocityEquationSIMPLE(0, v, id, homo, vlin_rows, B, DX, type==SIMPLEC, Efactor, dt, LocalDt);
    vlin_solver.SetLinearOperator(vlin_rows);
    lin_success = vlin_solver.Solve(B, VXstar, NULL, NULL, &lin_rnorm);
    if(!lin_success) {
      print_warning("  x Warning: Linear solver for the x-momentum equation failed to converge.\n");
      for(int i=0; i<(int)lin_rnorm.size(); i++)
        print_warning("    > It. %d: residual = %e.\n", i+1, lin_rnorm[i]);
    }

    print("Here 3.\n");
    // Solve the y-momentum equation
    if(!global_mesh.IsMesh1D()) {
      inco.BuildVelocityEquationSIMPLE(1, v, id, homo, vlin_rows, B, DY, type==SIMPLEC, Efactor, dt, LocalDt);
      vlin_solver.SetLinearOperator(vlin_rows);
      lin_success = vlin_solver.Solve(B, VYstar, NULL, NULL, &lin_rnorm);
      if(!lin_success) {
        print_warning("  x Warning: Linear solver for the y-momentum equation failed to converge.\n");
        for(int i=0; i<(int)lin_rnorm.size(); i++)
          print_warning("      > It. %d: residual = %e.\n", i+1, lin_rnorm[i]);
      }
    }

    print("Here 4.\n");
    // Solve the z-momentum equation
    if(!global_mesh.IsMesh1D() && !global_mesh.IsMesh2D()) {
      inco.BuildVelocityEquationSIMPLE(2, v, id, homo, vlin_rows, B, DZ, type==SIMPLEC, Efactor, dt, LocalDt);
      vlin_solver.SetLinearOperator(vlin_rows);
      lin_success = vlin_solver.Solve(B, VZstar, NULL, NULL, &lin_rnorm);
      if(!lin_success) {
        print_warning("  x Warning: Linear solver for the z-momentum equation failed to converge.\n");
        for(int i=0; i<(int)lin_rnorm.size(); i++)
          print_warning("    > It. %d: residual = %e.\n", i+1, lin_rnorm[i]);
      }
    }

    print("Good for now.\n");
    exit_mpi();

    
    //-----------------------------------------------------
    // Step 2: Solve the p' equation
    //-----------------------------------------------------
    inco.BuildPressureEquationSIMPLE(v, homo, VXstar, VYstar, VZstar, DX, DY, DZ, plin_rows, B, &ijk_zero_p);
    plin_solver.SetLinearOperator(plin_rows);
    Pprime.SetConstantValue(0.0, true); //!< This is p *correction*. Set init guess to 0 (Patankar 6.7-4)
    lin_success = plin_solver.Solve(B, Pprime, NULL, NULL, &lin_rnorm);
    if(!lin_success) {
      print_warning("  x Warning: Linear solver for the pressure correction equation failed to converge.\n");
      for(int i=0; i<(int)lin_rnorm.size(); i++)
        print_warning("    > It. %d: residual = %e.\n", i+1, lin_rnorm[i]);
    }


    //-----------------------------------------------------
    // Step 3: Update p, u, v, w, and compute relative error in velocity
    //-----------------------------------------------------
    rel_err = UpdateStates(v, Pprime, DX, DY, DZ, VXstar, VYstar, VZstar, alphaP); 

    V.RestoreDataPointerAndInsert();

    if(rel_err<iod.ts.semi_impl.convergence_tolerance) {
      converged = true;
      break; 
    }

    print("  o It. %d: Relative error in velocity (2-norm): %e.\n", iter+1, rel_err);

  }

  if(converged)
    print("  o Converged after %d iterations. Relative error in velocity (2-norm): %e.\n", iter+1, rel_err);
  else
    print_warning("  o Failed to converge. Relative error in velocity (2-norm): %e.\n", rel_err);
    

  ID.RestoreDataPointerToLocalVector();
  Homo.RestoreDataPointerToLocalVector();
}

//----------------------------------------------------------------------------

void
TimeIntegratorSIMPLE::ExtractVariableComponents(Vec5D*** v, SpaceVariable3D *VX_ptr, SpaceVariable3D *VY_ptr,
                                                SpaceVariable3D *VZ_ptr, SpaceVariable3D *P_ptr)
{
  double*** vx = VX_ptr ? VX_ptr->GetDataPointer() : NULL;
  double*** vy = VY_ptr ? VY_ptr->GetDataPointer() : NULL;
  double*** vz = VZ_ptr ? VZ_ptr->GetDataPointer() : NULL;
  double*** pp = P_ptr  ? P_ptr->GetDataPointer()  : NULL;

  int ii0, jj0, kk0, iimax, jjmax, kkmax;
  VXstar.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);

  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++) {
        if(vx) vx[k][j][i] = v[k][j][i][1];
        if(vy) vy[k][j][i] = v[k][j][i][2];
        if(vz) vz[k][j][i] = v[k][j][i][3];
        if(pp) pp[k][j][i] = v[k][j][i][4];
      }

  if(VX_ptr) VX_ptr->RestoreDataPointerToLocalVector(); //no need to exchange, as we have covered ghost layers
  if(VY_ptr) VY_ptr->RestoreDataPointerToLocalVector();
  if(VZ_ptr) VZ_ptr->RestoreDataPointerToLocalVector();
  if(P_ptr)   P_ptr->RestoreDataPointerToLocalVector();
}

//----------------------------------------------------------------------------

Int3
TimeIntegratorSIMPLE::FindCornerFixedPressure()
{
  GlobalMeshInfo &global_mesh(spo.GetGlobalMeshInfo());
  return Int3(global_mesh.NX-1, global_mesh.NY-1, global_mesh.NZ-1);
}

//----------------------------------------------------------------------------

double
TimeIntegratorSIMPLE::UpdateStates(Vec5D*** v, SpaceVariable3D &Pprime, SpaceVariable3D &DX,
                                   SpaceVariable3D &DY, SpaceVariable3D &DZ,
                                   SpaceVariable3D &VX, SpaceVariable3D &VY,
                                   SpaceVariable3D &VZ, double prelax)
{
  GlobalMeshInfo &global_mesh(spo.GetGlobalMeshInfo());

  double*** diagx = DX.GetDataPointer();
  double*** diagy = DY.GetDataPointer();
  double*** diagz = DZ.GetDataPointer();
  double*** ustar = VX.GetDataPointer();
  double*** vstar = VY.GetDataPointer();
  double*** wstar = VZ.GetDataPointer();
  double*** pp    = Pprime.GetDataPointer();

  double uerr  = 0.0;
  double unorm = 0.0;
  double ucorr, vcorr, wcorr, unew, vnew, wnew;
  double dz, dydz, dxdydz;

  int i0, j0, k0, imax, jmax, kmax;
  DX.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);

  for(int k=k0; k<kmax; k++) {
    dz = global_mesh.GetDz(k);
    for(int j=j0; j<jmax; j++) {
      dydz = dz*global_mesh.GetDy(j);
      for(int i=i0; i<imax; i++) {
        dxdydz = dydz*global_mesh.GetDx(i);

        ucorr = i>0 ? diagx[k][j][i]*(pp[k][j][i-1] - pp[k][j][i]) : 0.0;
        vcorr = j>0 ? diagy[k][j][i]*(pp[k][j-1][i] - pp[k][j][i]) : 0.0;
        wcorr = k>0 ? diagz[k][j][i]*(pp[k-1][j][i] - pp[k][j][i]) : 0.0;
     
        if(i>0) v[k][j][i][1] = ustar[k][j][i] + ucorr;
        if(j>0) v[k][j][i][2] = vstar[k][j][i] + vcorr;
        if(k>0) v[k][j][i][3] = wstar[k][j][i] + wcorr;
        v[k][j][i][4] += prelax*pp[k][j][i];

        unew = v[k][j][i][1];
        vnew = v[k][j][i][2];
        wnew = v[k][j][i][3];

        unorm += (unew*unew + vnew*vnew + wnew*wnew)*dxdydz;
        uerr  += (ucorr*ucorr + vcorr*vcorr + wcorr*wcorr)*dxdydz;
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &unorm, 1, MPI_DOUBLE, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, &uerr, 1, MPI_DOUBLE, MPI_SUM, comm);

  DX.RestoreDataPointerToLocalVector();
  DY.RestoreDataPointerToLocalVector();
  DZ.RestoreDataPointerToLocalVector();
  VX.RestoreDataPointerToLocalVector();
  VY.RestoreDataPointerToLocalVector();
  VZ.RestoreDataPointerToLocalVector();
  Pprime.RestoreDataPointerToLocalVector();
   
  return sqrt(uerr/unorm);

}


//----------------------------------------------------------------------------



//----------------------------------------------------------------------------
// SIMPLER
//----------------------------------------------------------------------------

TimeIntegratorSIMPLER::TimeIntegratorSIMPLER(MPI_Comm &comm_, IoData& iod_, DataManagers3D& dms_,
                                             SpaceOperator& spo_, IncompressibleOperator &inco_,
                                             vector<LevelSetOperator*>& lso_, MultiPhaseOperator &mpo_,
                                             LaserAbsorptionSolver* laser_, EmbeddedBoundaryOperator* embed_,
                                             HyperelasticityOperator* heo_, PrescribedMotionOperator* pmo_)
                     : TimeIntegratorSIMPLE(comm_, iod_, dms_, spo_, inco_, lso_, mpo_, laser_, embed_, 
                                            heo_, pmo_),
                       Bu(comm_, &(dms_.ghosted1_1dof)), Bv(comm_, &(dms_.ghosted1_1dof)),
                       Bw(comm_, &(dms_.ghosted1_1dof)), P(comm_, &(dms_.ghosted1_1dof)),
                       ulin_solver(comm_, dms_.ghosted1_1dof, iod.ts.semi_impl.velocity_linear_solver),
                       wlin_solver(comm_, dms_.ghosted1_1dof, iod.ts.semi_impl.pressure_linear_solver)
{
  type = SIMPLER;
  alphaP = 0.0; //not used
}

//----------------------------------------------------------------------------

TimeIntegratorSIMPLER::~TimeIntegratorSIMPLER()
{ }

//----------------------------------------------------------------------------

void
TimeIntegratorSIMPLER::Destroy()
{
  Bu.Destroy();
  Bv.Destroy();
  Bw.Destroy();
  P.Destroy();
  TimeIntegratorSIMPLE::Destroy();
}

//----------------------------------------------------------------------------

void
TimeIntegratorSIMPLER::AdvanceOneTimeStep(SpaceVariable3D &V, SpaceVariable3D &ID,
                                          vector<SpaceVariable3D*>& Phi, vector<SpaceVariable3D*> &NPhi,
                                          vector<SpaceVariable3D*> &KappaPhi,
                                          SpaceVariable3D *L, SpaceVariable3D *Xi, SpaceVariable3D *LocalDt,
                                          [[maybe_unused]] double time, double dt, int time_step, 
                                          int subcycle, double dts)
{

  if(mpo.NumberOfMaterials()>1) {
    print_error("*** Error: Need to update homogeneity. Currently, the incompressible flow solver does not allow"
                " more than one material.\n");
    exit_mpi();
  }
  if(Phi.size()>0 || NPhi.size()>0 || KappaPhi.size()>0 || L || Xi || subcycle>0 || dts != dt) {
    print_error("*** Error: Problem setup is not supported by TimeIntegratorSIMPLEC.\n");
    exit_mpi();
  }


  GlobalMeshInfo &global_mesh(spo.GetGlobalMeshInfo());

  double*** id = ID.GetDataPointer();
  double*** homo = Homo.GetDataPointer();

  int iter, maxIter = time_step == 1 ? 10*iod.ts.semi_impl.maxIts : iod.ts.semi_impl.maxIts;
  vector<double> lin_rnorm; 
  bool lin_success, converged(false);
  double rel_err(10000.0);

  print("  o Running the iterative SIMPLER procedure (E = %e).\n", Efactor);

  for(iter = 0; iter < maxIter; iter++) {

    Vec5D*** v = (Vec5D***)V.GetDataPointer();

    ExtractVariableComponents(v, NULL, NULL, NULL, &P);
    
    //-----------------------------------------------------
    // Step 1: Solve the pressure (not p') equation
    //-----------------------------------------------------
    B.SetConstantValue(0.0);
    Bu.SetConstantValue(0.0);
    Bv.SetConstantValue(0.0);
    Bw.SetConstantValue(0.0);
    // ulin_rows, vlin_rows, and wlin_rows will be used in Step 2
    inco.CalculateCoefficientsSIMPLER(0, v, id, homo, ulin_rows, Bu, VXstar, DX, Efactor, dt, LocalDt); //"Uhat"
    inco.CalculateCoefficientsSIMPLER(1, v, id, homo, vlin_rows, Bv, VYstar, DY, Efactor, dt, LocalDt); //"Vhat"
    inco.CalculateCoefficientsSIMPLER(2, v, id, homo, wlin_rows, Bw, VZstar, DZ, Efactor, dt, LocalDt); //"What"
    inco.BuildPressureEquationSIMPLE(v, homo, VXstar, VYstar, VZstar, DX, DY, DZ, plin_rows, B, &ijk_zero_p);
    plin_solver.SetLinearOperator(plin_rows);
    lin_success = plin_solver.Solve(B, P, NULL, NULL, &lin_rnorm);
    if(!lin_success) {
      print_warning("  x Warning: Linear solver for the pressure equation failed to converge.\n");
      for(int i=0; i<(int)lin_rnorm.size(); i++)
        print_warning("    > It. %d: residual = %e.\n", i+1, lin_rnorm[i]);
    }


    //-----------------------------------------------------
    // Step 2: Solve the momentum equations for u*, v*, w*
    //-----------------------------------------------------
    ExtractVariableComponents(v, &VXstar, &VYstar, &VZstar, NULL); //initial guesses

    // Solve the x-momentum equation
    inco.UpdateVelocityEquationRHS_SIMPLER(0, P, Bu);
    vlin_solver.SetLinearOperator(ulin_rows);
    lin_success = vlin_solver.Solve(Bu, VXstar, NULL, NULL, &lin_rnorm);
    if(!lin_success) {
      print_warning("  x Warning: Linear solver for the x-momentum equation failed to converge.\n");
      for(int i=0; i<(int)lin_rnorm.size(); i++)
        print_warning("    > It. %d: residual = %e.\n", i+1, lin_rnorm[i]);
    }

    // Solve the y-momentum equation
    if(!global_mesh.IsMesh1D()) {
      inco.UpdateVelocityEquationRHS_SIMPLER(1, P, Bv);
      vlin_solver.SetLinearOperator(vlin_rows);
      lin_success = vlin_solver.Solve(Bv, VYstar, NULL, NULL, &lin_rnorm);
      if(!lin_success) {
        print_warning("  x Warning: Linear solver for the y-momentum equation failed to converge.\n");
        for(int i=0; i<(int)lin_rnorm.size(); i++)
          print_warning("      > It. %d: residual = %e.\n", i+1, lin_rnorm[i]);
      }
    }

    // Solve the z-momentum equation
    if(!global_mesh.IsMesh1D() && !global_mesh.IsMesh2D()) {
      inco.UpdateVelocityEquationRHS_SIMPLER(2, P, Bw);
      vlin_solver.SetLinearOperator(wlin_rows);
      lin_success = vlin_solver.Solve(Bw, VZstar, NULL, NULL, &lin_rnorm);
      if(!lin_success) {
        print_warning("  x Warning: Linear solver for the z-momentum equation failed to converge.\n");
        for(int i=0; i<(int)lin_rnorm.size(); i++)
          print_warning("    > It. %d: residual = %e.\n", i+1, lin_rnorm[i]);
      }
    }

    
    //-----------------------------------------------------
    // Step 3: Solve the p' equation
    //-----------------------------------------------------
    inco.BuildPressureEquationRHS_SIMPLER(v, homo, VXstar, VYstar, VZstar, B, &ijk_zero_p);
    plin_solver.UsePreviousPreconditioner(true); //The matrix A is still the same
    Pprime.SetConstantValue(0.0, true); //!< This is p *correction*. Set init guess to 0 (Patankar 6.7-4)
    lin_success = plin_solver.Solve(B, Pprime, NULL, NULL, &lin_rnorm);
    if(!lin_success) {
      print_warning("  x Warning: Linear solver for the pressure correction equation failed to converge.\n");
      for(int i=0; i<(int)lin_rnorm.size(); i++)
        print_warning("    > It. %d: residual = %e.\n", i+1, lin_rnorm[i]);
    }


    //-----------------------------------------------------
    // Step 4: Update p, u, v, w, and compute relative error in velocity
    //-----------------------------------------------------
    rel_err = UpdateStates(v, P, Pprime, DX, DY, DZ, VXstar, VYstar, VZstar); 

    V.RestoreDataPointerAndInsert();

    if(rel_err<iod.ts.semi_impl.convergence_tolerance) {
      converged = true;
      break; 
    }

    print("  o It. %d: Relative error in velocity (2-norm): %e.\n", iter+1, rel_err);

  }

  if(converged)
    print("  o Converged after %d iterations. Relative error in velocity (2-norm): %e.\n", iter+1, rel_err);
  else
    print_warning("  o Failed to converge. Relative error in velocity (2-norm): %e.\n", rel_err);
    

  ID.RestoreDataPointerToLocalVector();
  Homo.RestoreDataPointerToLocalVector();
}

//----------------------------------------------------------------------------

double
TimeIntegratorSIMPLER::UpdateStates(Vec5D*** v, SpaceVariable3D &P, SpaceVariable3D &Pprime,
                                    SpaceVariable3D &DX, SpaceVariable3D &DY, SpaceVariable3D &DZ,
                                    SpaceVariable3D &VX, SpaceVariable3D &VY, SpaceVariable3D &VZ)
{
  GlobalMeshInfo &global_mesh(spo.GetGlobalMeshInfo());

  double*** diagx = DX.GetDataPointer();
  double*** diagy = DY.GetDataPointer();
  double*** diagz = DZ.GetDataPointer();
  double*** ustar = VX.GetDataPointer();
  double*** vstar = VY.GetDataPointer();
  double*** wstar = VZ.GetDataPointer();
  double*** p     = P.GetDataPointer();
  double*** pp    = Pprime.GetDataPointer();

  double uerr  = 0.0;
  double unorm = 0.0;
  double ucorr, vcorr, wcorr, unew, vnew, wnew;
  double dz, dydz, dxdydz;

  int i0, j0, k0, imax, jmax, kmax;
  DX.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);

  for(int k=k0; k<kmax; k++) {
    dz = global_mesh.GetDz(k);
    for(int j=j0; j<jmax; j++) {
      dydz = dz*global_mesh.GetDy(j);
      for(int i=i0; i<imax; i++) {
        dxdydz = dydz*global_mesh.GetDx(i);

        ucorr = i>0 ? diagx[k][j][i]*(pp[k][j][i-1] - pp[k][j][i]) : 0.0;
        vcorr = j>0 ? diagy[k][j][i]*(pp[k][j-1][i] - pp[k][j][i]) : 0.0;
        wcorr = k>0 ? diagz[k][j][i]*(pp[k-1][j][i] - pp[k][j][i]) : 0.0;
     
        if(i>0) v[k][j][i][1] = ustar[k][j][i] + ucorr;
        if(j>0) v[k][j][i][2] = vstar[k][j][i] + vcorr;
        if(k>0) v[k][j][i][3] = wstar[k][j][i] + wcorr;
        v[k][j][i][4] = p[k][j][i];

        unew = v[k][j][i][1];
        vnew = v[k][j][i][2];
        wnew = v[k][j][i][3];

        unorm += (unew*unew + vnew*vnew + wnew*wnew)*dxdydz;
        uerr  += (ucorr*ucorr + vcorr*vcorr + wcorr*wcorr)*dxdydz;
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &unorm, 1, MPI_DOUBLE, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, &uerr, 1, MPI_DOUBLE, MPI_SUM, comm);

  DX.RestoreDataPointerToLocalVector();
  DY.RestoreDataPointerToLocalVector();
  DZ.RestoreDataPointerToLocalVector();
  VX.RestoreDataPointerToLocalVector();
  VY.RestoreDataPointerToLocalVector();
  VZ.RestoreDataPointerToLocalVector();
  P.RestoreDataPointerToLocalVector();
  Pprime.RestoreDataPointerToLocalVector();
   
  return sqrt(uerr/unorm);

}




//----------------------------------------------------------------------------
// SIMPLEC
//----------------------------------------------------------------------------

TimeIntegratorSIMPLEC::TimeIntegratorSIMPLEC(MPI_Comm &comm_, IoData& iod_, DataManagers3D& dms_,
                                             SpaceOperator& spo_, IncompressibleOperator &inco_,
                                             vector<LevelSetOperator*>& lso_, MultiPhaseOperator &mpo_,
                                             LaserAbsorptionSolver* laser_, EmbeddedBoundaryOperator* embed_,
                                             HyperelasticityOperator* heo_, PrescribedMotionOperator* pmo_)
                     : TimeIntegratorSIMPLE(comm_, iod_, dms_, spo_, inco_, lso_, mpo_, laser_, embed_, 
                                            heo_, pmo_)
{
  type = SIMPLEC;
  alphaP = 1.0; //!< fixed to 1.0 in the SIMPLEC algorithm
}

//----------------------------------------------------------------------------

TimeIntegratorSIMPLEC::~TimeIntegratorSIMPLEC()
{ }

//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
// PISO 
//----------------------------------------------------------------------------

TimeIntegratorPISO::TimeIntegratorPISO(MPI_Comm &comm_, IoData& iod_, DataManagers3D& dms_,
                                       SpaceOperator& spo_, IncompressibleOperator &inco_,
                                       vector<LevelSetOperator*>& lso_, MultiPhaseOperator &mpo_,
                                       LaserAbsorptionSolver* laser_, EmbeddedBoundaryOperator* embed_,
                                       HyperelasticityOperator* heo_, PrescribedMotionOperator* pmo_)
                     : TimeIntegratorSIMPLE(comm_, iod_, dms_, spo_, inco_, lso_, mpo_, laser_, embed_, 
                                            heo_, pmo_),
                       VXprime(comm_, &(dms_.ghosted1_1dof)), VYprime(comm_, &(dms_.ghosted1_1dof)),
                       VZprime(comm_, &(dms_.ghosted1_1dof)), Pstar(comm_, &(dms_.ghosted1_1dof)),
                       VXtildep(comm_, &(dms_.ghosted1_1dof)), VYtildep(comm_, &(dms_.ghosted1_1dof)),
                       VZtildep(comm_, &(dms_.ghosted1_1dof))
{
  type = PISO;
  Efactor = 1.0e8; // essentially, no relaxation
  alphaP  = 1.0; //no relaxation
}
//----------------------------------------------------------------------------

TimeIntegratorPISO::~TimeIntegratorPISO()
{ }

//----------------------------------------------------------------------------

void
TimeIntegratorPISO::Destroy()
{ 
  VXprime.Destroy();
  VYprime.Destroy();
  VZprime.Destroy();
  VXtildep.Destroy();
  VYtildep.Destroy();
  VZtildep.Destroy();
  Pstar.Destroy();

  TimeIntegratorSIMPLE::Destroy();
}

//----------------------------------------------------------------------------

void
TimeIntegratorPISO::AdvanceOneTimeStep(SpaceVariable3D &V, SpaceVariable3D &ID,
                                       vector<SpaceVariable3D*>& Phi, vector<SpaceVariable3D*> &NPhi,
                                       vector<SpaceVariable3D*> &KappaPhi,
                                       SpaceVariable3D *L, SpaceVariable3D *Xi, SpaceVariable3D *LocalDt,
                                       [[maybe_unused]] double time, double dt, int time_step, int subcycle,
                                       double dts)
{

  if(mpo.NumberOfMaterials()>1) {
    print_error("*** Error: Need to update homogeneity. Currently, the incompressible flow solver does not allow"
                " more than one material.\n");
    exit_mpi();
  }
  if(Phi.size()>0 || NPhi.size()>0 || KappaPhi.size()>0 || L || Xi || subcycle>0 || dts != dt) {
    print_error("*** Error: Problem setup is not supported by TimeIntegratorSIMPLE(or SIMPLEC).\n");
    exit_mpi();
  }


  GlobalMeshInfo &global_mesh(spo.GetGlobalMeshInfo());

  double*** id = ID.GetDataPointer();
  double*** homo = Homo.GetDataPointer();

  int iter(0), maxIter = time_step == 1 ? 10*iod.ts.semi_impl.maxIts : iod.ts.semi_impl.maxIts;
  vector<double> lin_rnorm; 
  bool lin_success, converged(false);
  double rel_err;

  print("  o Running the PISO procedure.\n");

  Vec5D*** v = (Vec5D***)V.GetDataPointer();

  ExtractVariableComponents(v, &VXstar, &VYstar, &VZstar, &Pstar);

  //-----------------------------------------------------
  // Step 1: Solve the momentum equations for u*, v*, w*
  //-----------------------------------------------------

  // Solve the x-momentum equation
  inco.BuildVelocityEquationSIMPLE(0, v, id, homo, vlin_rows, B, DX, false, Efactor, dt, LocalDt);
  vlin_solver.SetLinearOperator(vlin_rows);
  lin_success = vlin_solver.Solve(B, VXstar, NULL, NULL, &lin_rnorm);
  if(!lin_success) {
    print_warning("  x Warning: Linear solver for the x-momentum equation failed to converge.\n");
    for(int i=0; i<(int)lin_rnorm.size(); i++)
      print_warning("    > It. %d: residual = %e.\n", i+1, lin_rnorm[i]);
  }

  // Solve the y-momentum equation
  if(!global_mesh.IsMesh1D()) {
    inco.BuildVelocityEquationSIMPLE(1, v, id, homo, vlin_rows, B, DY, false, Efactor, dt, LocalDt);
    vlin_solver.SetLinearOperator(vlin_rows);
    lin_success = vlin_solver.Solve(B, VYstar, NULL, NULL, &lin_rnorm);
    if(!lin_success) {
      print_warning("  x Warning: Linear solver for the y-momentum equation failed to converge.\n");
      for(int i=0; i<(int)lin_rnorm.size(); i++)
        print_warning("      > It. %d: residual = %e.\n", i+1, lin_rnorm[i]);
    }
  }

  // Solve the z-momentum equation
  if(!global_mesh.IsMesh1D() && !global_mesh.IsMesh2D()) {
    inco.BuildVelocityEquationSIMPLE(2, v, id, homo, vlin_rows, B, DZ, false, Efactor, dt, LocalDt);
    vlin_solver.SetLinearOperator(vlin_rows);
    lin_success = vlin_solver.Solve(B, VZstar, NULL, NULL, &lin_rnorm);
    if(!lin_success) {
      print_warning("  x Warning: Linear solver for the z-momentum equation failed to converge.\n");
      for(int i=0; i<(int)lin_rnorm.size(); i++)
        print_warning("    > It. %d: residual = %e.\n", i+1, lin_rnorm[i]);
    }
  }

 
   
  //-----------------------------------------------------
  // Step 2: Solve the p' equation (First Corrector Step)
  //-----------------------------------------------------
  inco.BuildPressureEquationSIMPLE(v, homo, VXstar, VYstar, VZstar, DX, DY, DZ, plin_rows, B, &ijk_zero_p);
  plin_solver.SetLinearOperator(plin_rows);
  Pprime.SetConstantValue(0.0, true); //!< This is p *correction*. Set init guess to 0 (Patankar 6.7-4)
  lin_success = plin_solver.Solve(B, Pprime, NULL, NULL, &lin_rnorm);
  if(!lin_success) {
    print_warning("  x Warning: Linear solver for the pressure correction equation failed to converge.\n");
    for(int i=0; i<(int)lin_rnorm.size(); i++)
      print_warning("    > It. %d: residual = %e.\n", i+1, lin_rnorm[i]);
  }

  
  //-----------------------------------------------------
  // Step 3: Update p, u, v, w, and compute relative error in velocity (First Corrector Step)
  //-----------------------------------------------------
  rel_err = UpdateStatesCorrector(Pstar, Pprime, DX, DY, DZ, VXstar, VYstar, VZstar, 
                                  VXprime, VYprime, VZprime); 
  print("  o Corrector Step %d: Relative error in velocity (2-norm): %e.\n", iter+1, rel_err);
   
  if(rel_err<iod.ts.semi_impl.convergence_tolerance) {
    converged = true;
    goto END_CORRECTORS;
  }

  //-----------------------------------------------------
  // Step 4: Update p, u, v, w, and compute relative error in velocity (Subsequent Corrector Steps)
  //-----------------------------------------------------
  for(iter = 1; iter < maxIter; iter++) {

    inco.CalculateVelocityTildePISO(0, v, id, homo, VXprime, VXtildep, Efactor, dt, LocalDt);
    inco.CalculateVelocityTildePISO(1, v, id, homo, VYprime, VYtildep, Efactor, dt, LocalDt);
    inco.CalculateVelocityTildePISO(2, v, id, homo, VZprime, VZtildep, Efactor, dt, LocalDt);

    inco.BuildPressureEquationRHS_SIMPLER(v, homo, VXtildep, VYtildep, VZtildep, B, &ijk_zero_p);
    plin_solver.UsePreviousPreconditioner(true); //The matrix A is still the same
    Pprime.SetConstantValue(0.0, true); //!< This is p *correction*. Set init guess to 0 (Patankar 6.7-4)
    lin_success = plin_solver.Solve(B, Pprime, NULL, NULL, &lin_rnorm);
    if(!lin_success) {
      print_warning("  x Warning: Linear solver for the pressure correction equation failed to converge.\n");
      for(int i=0; i<(int)lin_rnorm.size(); i++)
        print_warning("    > It. %d: residual = %e.\n", i+1, lin_rnorm[i]);
    }

    rel_err = UpdateStatesCorrector(Pstar, Pprime, DX, DY, DZ, VXstar, VYstar, VZstar, 
                                    VXprime, VYprime, VZprime, &VXtildep, &VYtildep, &VZtildep); 

    if(rel_err<iod.ts.semi_impl.convergence_tolerance) {
      converged = true;
      break; 
    }

    print("  o It. %d: Relative error in velocity (2-norm): %e.\n", iter+1, rel_err);

  }

END_CORRECTORS:

  UpdateStatesFinal(v, Pstar, VXstar, VYstar, VZstar); // update v
  V.RestoreDataPointerAndInsert();

  if(converged)
    print("  o Converged after %d iterations. Relative error in velocity (2-norm): %e.\n", iter+1, rel_err);
  else
    print_warning("  o Failed to converge. Relative error in velocity (2-norm): %e.\n", rel_err);


  ID.RestoreDataPointerToLocalVector();
  Homo.RestoreDataPointerToLocalVector();

}

//----------------------------------------------------------------------------

double
TimeIntegratorPISO::UpdateStatesCorrector(SpaceVariable3D &Pstar, SpaceVariable3D &Pprime,
                        SpaceVariable3D &DX, SpaceVariable3D &DY, SpaceVariable3D &DZ,
                        SpaceVariable3D &VXstar, SpaceVariable3D &VYstar, SpaceVariable3D &VZstar,
                        SpaceVariable3D &VXprime, SpaceVariable3D &VYprime, SpaceVariable3D &VZprime,
                        SpaceVariable3D *VXtildep, SpaceVariable3D *VYtildep,
                        SpaceVariable3D *VZtildep)
{

  GlobalMeshInfo &global_mesh(spo.GetGlobalMeshInfo());

  //inputs: DX, DY, DZ, VXstar, VYstar, VZstar, Pstar, Pprime
  //        (VXtildep, VYtildep, VZtildep for 2nd, 3rd, etc. correctors)
  //outputs: VXstar, VYstar, VZstar, Pstar, VXprime, VYprime, VZprime
  double*** diagx = DX.GetDataPointer();
  double*** diagy = DY.GetDataPointer();
  double*** diagz = DZ.GetDataPointer();
  double*** ustar = VXstar.GetDataPointer();
  double*** vstar = VYstar.GetDataPointer();
  double*** wstar = VZstar.GetDataPointer();

  double*** utildep = VXtildep ? VXtildep->GetDataPointer() : NULL;
  double*** vtildep = VYtildep ? VYtildep->GetDataPointer() : NULL;
  double*** wtildep = VZtildep ? VZtildep->GetDataPointer() : NULL;
  
  double*** up    = VXprime.GetDataPointer();
  double*** vp    = VYprime.GetDataPointer();
  double*** wp    = VZprime.GetDataPointer();
  double*** pstar = Pstar.GetDataPointer();
  double*** pp    = Pprime.GetDataPointer();

  double uerr  = 0.0;
  double unorm = 0.0;
  double ucorr, vcorr, wcorr, unew, vnew, wnew;
  double dz, dydz, dxdydz;

  int i0, j0, k0, imax, jmax, kmax;
  DX.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);

  for(int k=k0; k<kmax; k++) {
    dz = global_mesh.GetDz(k);
    for(int j=j0; j<jmax; j++) {
      dydz = dz*global_mesh.GetDy(j);
      for(int i=i0; i<imax; i++) {
        dxdydz = dydz*global_mesh.GetDx(i);

        ucorr = i>0 ? (utildep ? utildep[k][j][i] : 0.0) +
                      diagx[k][j][i]*(pp[k][j][i-1] - pp[k][j][i])
                    : 0.0;
        vcorr = j>0 ? (vtildep ? vtildep[k][j][i] : 0.0) +
                      diagy[k][j][i]*(pp[k][j-1][i] - pp[k][j][i])
                    : 0.0;
        wcorr = k>0 ? (wtildep ? wtildep[k][j][i] : 0.0) +
                      diagz[k][j][i]*(pp[k-1][j][i] - pp[k][j][i])
                    : 0.0;
     
        // Updates
        ustar[k][j][i] += ucorr;
        vstar[k][j][i] += vcorr;
        wstar[k][j][i] += wcorr;
        up[k][j][i]     = ucorr;
        vp[k][j][i]     = vcorr;
        wp[k][j][i]     = wcorr;
        pstar[k][j][i] += pp[k][j][i];

        unew = ustar[k][j][i];
        vnew = vstar[k][j][i];
        wnew = wstar[k][j][i];
        unorm += (unew*unew + vnew*vnew + wnew*wnew)*dxdydz;
        uerr  += (ucorr*ucorr + vcorr*vcorr + wcorr*wcorr)*dxdydz;
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &unorm, 1, MPI_DOUBLE, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, &uerr, 1, MPI_DOUBLE, MPI_SUM, comm);

  DX.RestoreDataPointerToLocalVector();
  DY.RestoreDataPointerToLocalVector();
  DZ.RestoreDataPointerToLocalVector();
  Pprime.RestoreDataPointerToLocalVector();
  VXstar.RestoreDataPointerAndInsert();
  VYstar.RestoreDataPointerAndInsert();
  VZstar.RestoreDataPointerAndInsert();
  VXprime.RestoreDataPointerAndInsert();
  VYprime.RestoreDataPointerAndInsert();
  VZprime.RestoreDataPointerAndInsert();
  Pstar.RestoreDataPointerAndInsert();
   
  if(VXtildep) VXtildep->RestoreDataPointerToLocalVector();
  if(VYtildep) VYtildep->RestoreDataPointerToLocalVector();
  if(VZtildep) VZtildep->RestoreDataPointerToLocalVector();

  return sqrt(uerr/unorm);

}

//----------------------------------------------------------------------------

void
TimeIntegratorPISO::UpdateStatesFinal(Vec5D*** v, SpaceVariable3D &P, SpaceVariable3D &VX,
                                      SpaceVariable3D &VY, SpaceVariable3D &VZ)
{

  double*** pp = P.GetDataPointer();
  double*** uu = VX.GetDataPointer();
  double*** vv = VY.GetDataPointer();
  double*** ww = VZ.GetDataPointer();

  int i0, j0, k0, imax, jmax, kmax;
  P.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);

  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        if(i>0) v[k][j][i][1] = uu[k][j][i];
        if(j>0) v[k][j][i][2] = vv[k][j][i];
        if(k>0) v[k][j][i][3] = ww[k][j][i];
        v[k][j][i][4] = pp[k][j][i];
      }

  VX.RestoreDataPointerToLocalVector();
  VY.RestoreDataPointerToLocalVector();
  VZ.RestoreDataPointerToLocalVector();
  P.RestoreDataPointerToLocalVector();
   
}

//----------------------------------------------------------------------------







