/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<TimeIntegratorSemiImp.h>
extern int verbose;

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
                      DZ(comm_, &(dms_.ghosted1_1dof)), V0(comm_, &(dms_.ghosted1_5dof)),
                      Tmp1(comm_, &(dms_.ghosted1_1dof)), Tmp5(comm_, &(dms_.ghosted1_5dof)),
                      vlin_solver(comm_, dms_.ghosted1_1dof, iod.ts.semi_impl.velocity_linear_solver,
                      "velocity"),
                      plin_solver(comm_, dms_.ghosted1_1dof, iod.ts.semi_impl.pressure_linear_solver,
                      "pressure"),
                      vturb_lin_solver(comm_, dms_.ghosted1_1dof, iod.ts.semi_impl.turbulence_linear_solver,
                      "turbulence"),
                      Vturb0_ptr(NULL), R3_ptr(NULL)
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

  if(iod.rans.sa_cw1_reduction_factor != 1.0) {
    if(iod.rans.model != RANSTurbulenceModelData::SPALART_ALLMARAS) {
      print_error("*** Error: CW1ReductionFactor is used with Spalart-Allmaras only.\n");
      exit_mpi();
    }
    if(iod.rans.sa_cw1_reduction_factor <= 0.0) {
      print_error("*** Error: CW1ReductionFactor must be positive.\n");
      exit_mpi();
    }
    if(iod.rans.sa_cw1_reduction_t1<0 || iod.rans.sa_cw1_reduction_t2<iod.rans.sa_cw1_reduction_t1) {
      print_error("*** Error: Detected error in CW1ReductionTime1 and/or CW1ReductionTime2.\n");
      exit_mpi();
    }
  }


  fix_pressure_at_one_corner = iod.ts.semi_impl.fix_pressure_at_one_corner == SemiImplicitTsData::YES;
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
  if(iod.rans.model != RANSTurbulenceModelData::NONE) {
    vturb_lin_solver.GetSolverType(&ksp_type, &pc_type);
    print("  o Linear solver for turbulence closure equations: %s, Preconditioner: %s.\n",
          ksp_type.c_str(), pc_type.c_str());
  }

  if(sso)
    R3_ptr = new SpaceVariable3D(comm_, &(dms_.ghosted1_3dof));

  if(iod.rans.model == RANSTurbulenceModelData::SPALART_ALLMARAS) {
    if(iod.ts.semi_impl.type != SemiImplicitTsData::SIMPLE &&
       iod.ts.semi_impl.type != SemiImplicitTsData::SIMPLEC) {
      print_error("*** Error: Currently, only SIMPLE and SIMPLEC support RANS Spalart-Allmaras.\n");
      exit_mpi();
    }
    Vturb0_ptr = new SpaceVariable3D(comm_, &(dms_.ghosted1_1dof));
  }

}

//----------------------------------------------------------------------------

TimeIntegratorSIMPLE::~TimeIntegratorSIMPLE()
{
  if(R3_ptr)
    delete R3_ptr;

  if(Vturb0_ptr)
    delete Vturb0_ptr;
}

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
  V0.Destroy();
  Tmp1.Destroy();
  Tmp5.Destroy();

  vlin_solver.Destroy();
  plin_solver.Destroy();
  vturb_lin_solver.Destroy();

  if(R3_ptr)
    R3_ptr->Destroy();

  if(Vturb0_ptr)
    Vturb0_ptr->Destroy();

  TimeIntegratorBase::Destroy();
}

//----------------------------------------------------------------------------

void
TimeIntegratorSIMPLE::AdvanceOneTimeStep(SpaceVariable3D &V, SpaceVariable3D &ID,
                                         vector<SpaceVariable3D*>& Phi, vector<SpaceVariable3D*> &NPhi,
                                         vector<SpaceVariable3D*> &KappaPhi,
                                         SpaceVariable3D *L, SpaceVariable3D *Xi, SpaceVariable3D *Vturb,
                                         SpaceVariable3D *LocalDt,
                                         double time, double dt, 
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
  vector<int>    lin_rnorm_its;
  bool lin_success, converged(false);
  int nLinIts(0);
  double rel_err_prev(100000), rel_err(10000.0);

  if(type == SIMPLEC)
    print("  o Running the iterative SIMPLEC procedure (E = %e).\n", Efactor);
  else
    print("  o Running the iterative SIMPLE procedure (E = %e, alphaP = %e).\n", Efactor, alphaP);

  // Store V0
  V0.AXPlusBY(0.0, 1.0, V); //V0 = V (only copy nodes inside physical domain)
  Vec5D*** v0 = (Vec5D***)V0.GetDataPointer();

  // Store Vturb0 (if provided)
  double*** vturb0 = NULL;
  if(Vturb) {//with turbulence
    assert(Vturb0_ptr);
    Vturb0_ptr->AXPlusBY(0.0, 1.0, *Vturb); //Vturb0 = Vturb
    vturb0 = Vturb0_ptr->GetDataPointer();
  }


  bool repetition = false;
  double plin_rtol, plin_abstol, plin_div;
  int plin_maxIts;
  plin_solver.GetTolerances(&plin_rtol, &plin_abstol, &plin_div, &plin_maxIts);

  for(iter = 0; iter < maxIter; iter++) {

    // Backup
    Tmp5.AXPlusBY(0.0, 1.0, V, true); //Tmp = V (include ghost nodes)
    if(Vturb)
      Tmp1.AXPlusBY(0.0, 1.0, *Vturb, true); //Tmp = vturb (include ghost nodes)

    Vec5D***      v = (Vec5D***)V.GetDataPointer();
    double*** vturb = Vturb ? Vturb->GetDataPointer() : NULL;

    ExtractVariableComponents(v, &VXstar, &VYstar, &VZstar, NULL);

    //-----------------------------------------------------
    // Step 1: Solve the momentum equations for u*, v*, w*
    //-----------------------------------------------------

    // Solve the x-momentum equation
    if(global_mesh.x_glob.size()>1) {
      inco.BuildVelocityEquationSIMPLE(0, v0, v, id, vturb, homo, vlin_rows, B, DX, type==SIMPLEC, Efactor,
                                       dt, LocalDt);
      vlin_solver.SetLinearOperator(vlin_rows);

      //print("X-Momentum condition number is: %f\n",vlin_solver.EstimateConditionNumber()); //prints out the condition number

      lin_success = vlin_solver.Solve(B, VXstar, NULL, &nLinIts, &lin_rnorm, &lin_rnorm_its);
      if(!lin_success) {
        print_warning("    x Warning: Linear solver for the x-momentum equation failed to converge."
                      " Residual: %e -> %e.\n", lin_rnorm.front(), lin_rnorm.back());
        if(verbose>=2)
          for(int i=0; i<(int)lin_rnorm.size(); i++)
            print_warning("      > It. %d: residual = %e.\n", lin_rnorm_its[i]+1, lin_rnorm[i]);
      } else {
        if(verbose>=1)
          print("    * Solver of the x-momentum equation converged in %d iterations. Residual: %e -> %e.\n",
                nLinIts, lin_rnorm.front(), lin_rnorm.back());
      }
    }


    // Solve the y-momentum equation
    if(global_mesh.y_glob.size()>1) {
      inco.BuildVelocityEquationSIMPLE(1, v0, v, id, vturb, homo, vlin_rows, B, DY, type==SIMPLEC, Efactor,
                                       dt, LocalDt);
      vlin_solver.SetLinearOperator(vlin_rows);

      //print("Y-Momentum condition number is: %f\n",vlin_solver.EstimateConditionNumber()); //prints out the condition number

      lin_success = vlin_solver.Solve(B, VYstar, NULL, &nLinIts, &lin_rnorm, &lin_rnorm_its);
      if(!lin_success) {
        print_warning("    x Warning: Linear solver for the y-momentum equation failed to converge."
                      " Residual: %e -> %e.\n", lin_rnorm.front(), lin_rnorm.back());
        if(verbose>=2)
          for(int i=0; i<(int)lin_rnorm.size(); i++)
            print_warning("        > It. %d: residual = %e.\n", lin_rnorm_its[i]+1, lin_rnorm[i]);
      } else {
        if(verbose>=1)
          print("    * Solver of the y-momentum equation converged in %d iterations. Residual: %e -> %e.\n",
                nLinIts, lin_rnorm.front(), lin_rnorm.back());
      }
    }


    // Solve the z-momentum equation
    if(global_mesh.z_glob.size()>1) {
      inco.BuildVelocityEquationSIMPLE(2, v0, v, id, vturb, homo, vlin_rows, B, DZ, type==SIMPLEC, Efactor,
                                       dt, LocalDt);
      vlin_solver.SetLinearOperator(vlin_rows);
      lin_success = vlin_solver.Solve(B, VZstar, NULL, &nLinIts, &lin_rnorm, &lin_rnorm_its);
      if(!lin_success) {
        print_warning("    x Warning: Linear solver for the z-momentum equation failed to converge."
                      " Residual: %e -> %e.\n", lin_rnorm.front(), lin_rnorm.back());
        if(verbose>=2)
          for(int i=0; i<(int)lin_rnorm.size(); i++)
            print_warning("      > It. %d: residual = %e.\n", lin_rnorm_its[i]+1, lin_rnorm[i]);
      } else {
        if(verbose>=1)
          print("    * Solver of the z-momentum equation converged in %d iterations. Residual: %e -> %e.\n",
                nLinIts, lin_rnorm.front(), lin_rnorm.back());
      }
    }



    //-----------------------------------------------------
    // Step 2: Solve the p' equation
    //-----------------------------------------------------
    inco.BuildPressureEquationSIMPLE(v, homo, VXstar, VYstar, VZstar, DX, DY, DZ, plin_rows, B,
                                     fix_pressure_at_one_corner ? &ijk_zero_p : NULL);
    plin_solver.SetLinearOperator(plin_rows);

    //print("Pressure condition number is: %f\n",plin_solver.EstimateConditionNumber()); //prints out the condition number

    Pprime.SetConstantValue(0.0, true); //!< This is p *correction*. Set init guess to 0 (Patankar 6.7-4)

    lin_success = plin_solver.Solve(B, Pprime, NULL, &nLinIts, &lin_rnorm, &lin_rnorm_its);
    if(!lin_success) {
      print_warning("    x Warning: Linear solver for the pressure correction equation failed to converge."
                    " Residual: %e -> %e.\n", lin_rnorm.front(), lin_rnorm.back());
/*
      plin_solver.WriteToMatlabFile("Mat.m", "A");
      B.WriteToMatlabFile("B.m", "b");
      Pprime.WriteToMatlabFile("P.m","x");
      double sigma_min, sigma_max;
      plin_solver.CalculateExtremeSingularValues(sigma_max, sigma_min);
      print("sigma_max = %e, sigma_min = %e.\n", sigma_max, sigma_min);
*/
      if(verbose>=2)
        for(int i=0; i<(int)lin_rnorm.size(); i++)
          print_warning("      > It. %d: residual = %e.\n", lin_rnorm_its[i]+1, lin_rnorm[i]);
    } else {
      if(verbose>=1)
        print("    * Solver of the pressure correction equation converged in %d iterations."
              " Residual: %e -> %e.\n", nLinIts, lin_rnorm.front(), lin_rnorm.back());
    }


    //plin_solver.WriteToMatlabFile("PMatrix.m", "PMat");
    //B.WriteToMatlabFile("PB.m", "PB");
    //Pprime.WriteToVTRFile("Pprime.vtr","Pprime");


    //-----------------------------------------------------
    // Step 3: Update p, u, v, w, and compute relative error in velocity
    //-----------------------------------------------------
    rel_err_prev = rel_err;
    rel_err = UpdateStates(v, Pprime, DX, DY, DZ, VXstar, VYstar, VZstar, alphaP); 

    V.RestoreDataPointerAndInsert();
    inco.ApplyBoundaryConditions(V);


    //-----------------------------------------------------
    // Step 4: Solve turbulence closure equations
    //-----------------------------------------------------
    if(Vturb) {
      if(iod.rans.model != RANSTurbulenceModelData::SPALART_ALLMARAS) {
        print_error("*** Error: Found unknown turbulence closure model.\n");
        exit_mpi();
      }
      v = (Vec5D***)V.GetDataPointer();
      double cw1_reduction = 1.0; //Ramping function imposed here to reduce cw1 
      if(time<iod.rans.sa_cw1_reduction_t1) {
        cw1_reduction = iod.rans.sa_cw1_reduction_factor;
      } else if(time<iod.rans.sa_cw1_reduction_t2) {
        double t1 = iod.rans.sa_cw1_reduction_t1;
        double t2 = iod.rans.sa_cw1_reduction_t2;
        double factor = iod.rans.sa_cw1_reduction_factor;
        cw1_reduction = factor - (factor - 1.0)/(t2 - t1)*(time-t1);
      }

      inco.BuildSATurbulenceEquationSIMPLE(v, id, vturb0, vturb, vturb_lin_rows, B, Efactor, cw1_reduction, dt, LocalDt);
      V.RestoreDataPointerToLocalVector();
      Vturb->RestoreDataPointerAndInsert();
      vturb_lin_solver.SetLinearOperator(vturb_lin_rows);

      //print("SA condition number is: %f\n",vturb_lin_solver.EstimateConditionNumber()); //prints out the condition number

      lin_success = vturb_lin_solver.Solve(B, *Vturb, NULL, &nLinIts, &lin_rnorm, &lin_rnorm_its);
      if(!lin_success) {
        print_warning("    x Warning: Linear solver for the turbulence closure equation failed to converge."
                      " Residual: %e -> %e.\n", lin_rnorm.front(), lin_rnorm.back());
        if(verbose>=2)
          for(int i=0; i<(int)lin_rnorm.size(); i++)
            print_warning("      > It. %d: residual = %e.\n", lin_rnorm_its[i]+1, lin_rnorm[i]);
      } else {
        if(verbose>=1)
          print("    * Solver of the turbulence closure equation converged in %d iterations."
                " Residual: %e -> %e.\n", nLinIts, lin_rnorm.front(), lin_rnorm.back());
      }

      inco.ApplyBoundaryConditionsTurbulenceVariables(*Vturb);
    }



    print("  o It. %d: Relative error in velocity (2-norm): %e.\n", iter+1, rel_err);
    if(rel_err/rel_err_prev>3.0) { 
      print_warning("  x Warning: Detected drastic increase of error (%eX). Repeating.\n",
                    rel_err/rel_err_prev);
      repetition = true;
      rel_err = rel_err_prev;
      iter--;

      double current_rtol;
      plin_solver.GetTolerances(&current_rtol, NULL, NULL, NULL);
      plin_solver.SetTolerances(10.0*current_rtol, plin_abstol, plin_div, plin_maxIts);

      V.AXPlusBY(0.0, 1.0, Tmp5, true); 
      if(Vturb)
        Vturb->AXPlusBY(0.0, 1.0, Tmp1, true); //Tmp = vturb (include ghost nodes)
    } else {
      if(repetition == true) {
        plin_solver.SetTolerances(plin_rtol, plin_abstol, plin_div, plin_maxIts);
        repetition = false;
      }
    }


    //TODO: For the moment, only check convergence of the N-S equations, not turbulence closure.
    if(rel_err<iod.ts.semi_impl.convergence_tolerance) {
      converged = true;
      break; 
    }
  }

  if(converged)
    print("  o Converged after %d iterations. Relative error in velocity (2-norm): %e.\n", iter+1, rel_err);
  else
    print_warning("  o Failed to converge. Relative error in velocity (2-norm): %e.\n", rel_err);
    

  if(sso) {
    assert(R3_ptr);
    inco.CalculateMomentumChanges(v0, V, id, *R3_ptr);
  }

  V0.RestoreDataPointerToLocalVector();

  if(Vturb0_ptr)
    Vturb0_ptr->RestoreDataPointerToLocalVector();

  ID.RestoreDataPointerToLocalVector();
  Homo.RestoreDataPointerToLocalVector();

  if(sso)
    sso->MonitorConvergence(*R3_ptr, ID);
    
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

  if(VX_ptr) VX_ptr->RestoreDataPointerAndInsert();
  if(VY_ptr) VY_ptr->RestoreDataPointerAndInsert();
  if(VZ_ptr) VZ_ptr->RestoreDataPointerAndInsert();
  if(P_ptr)   P_ptr->RestoreDataPointerAndInsert();
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
   
  if(fabs(unorm)==0.0) {//new velocity is zero everywhere (return abs. error in this case)
    unorm = 1.0;
  }

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
                                          SpaceVariable3D *L, SpaceVariable3D *Xi, SpaceVariable3D *Vturb,
                                          SpaceVariable3D *LocalDt,
                                          [[maybe_unused]] double time, double dt, int time_step, 
                                          int subcycle, double dts)
{

  if(mpo.NumberOfMaterials()>1) {
    print_error("*** Error: Need to update homogeneity. Currently, the incompressible flow solver does not allow"
                " more than one material.\n");
    exit_mpi();
  }
  if(Phi.size()>0 || NPhi.size()>0 || KappaPhi.size()>0 || L || Xi || subcycle>0 || dts != dt) {
    print_error("*** Error: Problem setup is not supported by TimeIntegratorSIMPLER.\n");
    exit_mpi();
  }

  if(Vturb) {
    print_error("*** Error: Turbulence models are not supported by TimeIntegratorSIMPLER.\n");
    exit_mpi();
  }

  GlobalMeshInfo &global_mesh(spo.GetGlobalMeshInfo());

  double*** id = ID.GetDataPointer();
  double*** homo = Homo.GetDataPointer();

  int iter, maxIter = time_step == 1 ? 10*iod.ts.semi_impl.maxIts : iod.ts.semi_impl.maxIts;
  vector<double> lin_rnorm; 
  vector<int>    lin_rnorm_its; 
  bool lin_success, converged(false);
  int nLinIts(0);
  double rel_err(10000.0);

  print("  o Running the iterative SIMPLER procedure (E = %e).\n", Efactor);

  // Store V0
  V0.AXPlusBY(0.0, 1.0, V); //V0 = V (only copy nodes inside physical domain)
  Vec5D*** v0 = (Vec5D***)V0.GetDataPointer();

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
    if(global_mesh.x_glob.size()>1)
      inco.CalculateCoefficientsSIMPLER(0, v0, v, id, homo, ulin_rows, Bu, VXstar, DX, Efactor, dt, LocalDt); //"Uhat"
    if(global_mesh.y_glob.size()>1)
      inco.CalculateCoefficientsSIMPLER(1, v0, v, id, homo, vlin_rows, Bv, VYstar, DY, Efactor, dt, LocalDt); //"Vhat"
    if(global_mesh.z_glob.size()>1)
      inco.CalculateCoefficientsSIMPLER(2, v0, v, id, homo, wlin_rows, Bw, VZstar, DZ, Efactor, dt, LocalDt); //"What"
    inco.BuildPressureEquationSIMPLE(v, homo, VXstar, VYstar, VZstar, DX, DY, DZ, plin_rows, B, 
                                     fix_pressure_at_one_corner ? &ijk_zero_p : NULL);
    plin_solver.SetLinearOperator(plin_rows);
    lin_success = plin_solver.Solve(B, P, NULL, &nLinIts, &lin_rnorm, &lin_rnorm_its);
    if(!lin_success) {
      print_warning("    x Warning: Linear solver for the pressure equation failed to converge."
                    " Residual: %e - > %e.\n", lin_rnorm.front(), lin_rnorm.back());
      if(verbose>=2)
        for(int i=0; i<(int)lin_rnorm.size(); i++)
          print_warning("      > It. %d: residual = %e.\n", lin_rnorm_its[i]+1, lin_rnorm[i]);
    } else {
      if(verbose>=1)
        print("    * Solver of the pressure equation converged in %d iterations. Residual: %e -> %e.\n",
              nLinIts, lin_rnorm.front(), lin_rnorm.back());
    }


    //-----------------------------------------------------
    // Step 2: Solve the momentum equations for u*, v*, w*
    //-----------------------------------------------------
    ExtractVariableComponents(v, &VXstar, &VYstar, &VZstar, NULL); //initial guesses

    // Solve the x-momentum equation
    if(global_mesh.x_glob.size()>1) {
      inco.UpdateVelocityEquationRHS_SIMPLER(0, P, Bu);
      vlin_solver.SetLinearOperator(ulin_rows);
      lin_success = vlin_solver.Solve(Bu, VXstar, NULL, &nLinIts, &lin_rnorm, &lin_rnorm_its);
      if(!lin_success) {
        print_warning("    x Warning: Linear solver for the x-momentum equation failed to converge."
                      " Residual: %e -> %e.\n", lin_rnorm.front(), lin_rnorm.back());
        if(verbose>=2)
          for(int i=0; i<(int)lin_rnorm.size(); i++)
            print_warning("      > It. %d: residual = %e.\n", lin_rnorm_its[i]+1, lin_rnorm[i]);
      } else {
        if(verbose>=1)
          print("    * Solver of the x-momentum equation converged in %d iterations. Residual: %e -> %e.\n",
                nLinIts, lin_rnorm.front(), lin_rnorm.back());
      }
    }

    // Solve the y-momentum equation
    if(global_mesh.y_glob.size()>1) {
      inco.UpdateVelocityEquationRHS_SIMPLER(1, P, Bv);
      vlin_solver.SetLinearOperator(vlin_rows);
      lin_success = vlin_solver.Solve(Bv, VYstar, NULL, &nLinIts, &lin_rnorm, &lin_rnorm_its);
      if(!lin_success) {
        print_warning("    x Warning: Linear solver for the y-momentum equation failed to converge."
                      " Residual: %e -> %e.\n", lin_rnorm.front(), lin_rnorm.back());
        if(verbose>=2)
          for(int i=0; i<(int)lin_rnorm.size(); i++)
            print_warning("      > It. %d: residual = %e.\n", lin_rnorm_its[i]+1, lin_rnorm[i]);
      } else {
        if(verbose>=1)
          print("    * Solver of the y-momentum equation converged in %d iterations. Residual: %e -> %e.\n",
                nLinIts, lin_rnorm.front(), lin_rnorm.back());
      }
    }

    // Solve the z-momentum equation
    if(global_mesh.z_glob.size()>1) {
      inco.UpdateVelocityEquationRHS_SIMPLER(2, P, Bw);
      vlin_solver.SetLinearOperator(wlin_rows);
      lin_success = vlin_solver.Solve(Bw, VZstar, NULL, &nLinIts, &lin_rnorm, &lin_rnorm_its);
      if(!lin_success) {
        print_warning("    x Warning: Linear solver for the z-momentum equation failed to converge."
                      " Residual: %e -> %e.\n", lin_rnorm.front(), lin_rnorm.back());
        if(verbose>=2)
          for(int i=0; i<(int)lin_rnorm.size(); i++)
            print_warning("      > It. %d: residual = %e.\n", lin_rnorm_its[i]+1, lin_rnorm[i]);
      } else {
        if(verbose>=1)
          print("    * Solver of the z-momentum equation converged in %d iterations. Residual: %e -> %e.\n",
                nLinIts, lin_rnorm.front(), lin_rnorm.back());
      }
    }

    
    //-----------------------------------------------------
    // Step 3: Solve the p' equation
    //-----------------------------------------------------
    inco.BuildPressureEquationRHS_SIMPLER(v, homo, VXstar, VYstar, VZstar, B,
                                          fix_pressure_at_one_corner ? &ijk_zero_p : NULL);
    plin_solver.UsePreviousPreconditioner(true); //The matrix A is still the same
    Pprime.SetConstantValue(0.0, true); //!< This is p *correction*. Set init guess to 0 (Patankar 6.7-4)
    lin_success = plin_solver.Solve(B, Pprime, NULL, &nLinIts, &lin_rnorm, &lin_rnorm_its);
    if(!lin_success) {
      print_warning("    x Warning: Linear solver for the pressure correction equation failed to converge."
                    " Residual: %e -> %e.\n", lin_rnorm.front(), lin_rnorm.back());
      if(verbose>=2)
        for(int i=0; i<(int)lin_rnorm.size(); i++)
          print_warning("      > It. %d: residual = %e.\n", lin_rnorm_its[i]+1, lin_rnorm[i]);
    } else {
      if(verbose>=1)
        print("    * Solver of the pressure correction equation converged in %d iterations. Residual: %e -> %e.\n",
              nLinIts, lin_rnorm.front(), lin_rnorm.back());
    }


    //-----------------------------------------------------
    // Step 4: Update p, u, v, w, and compute relative error in velocity
    //-----------------------------------------------------
    rel_err = UpdateStates(v, P, Pprime, DX, DY, DZ, VXstar, VYstar, VZstar); 

    V.RestoreDataPointerAndInsert();
    inco.ApplyBoundaryConditions(V);


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
 

  if(sso) {
    assert(R3_ptr);
    inco.CalculateMomentumChanges(v0, V, id, *R3_ptr);
  }
   
  V0.RestoreDataPointerToLocalVector();

  ID.RestoreDataPointerToLocalVector();
  Homo.RestoreDataPointerToLocalVector();

  if(sso)
    sso->MonitorConvergence(*R3_ptr, ID);
    
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
   
  if(fabs(unorm)==0.0) {//new velocity is zero everywhere (return abs. error in this case)
    unorm = 1.0;
  }

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
                                       SpaceVariable3D *L, SpaceVariable3D *Xi, SpaceVariable3D *Vturb,
                                       SpaceVariable3D *LocalDt,
                                       [[maybe_unused]] double time, double dt, int time_step, int subcycle,
                                       double dts)
{

  if(mpo.NumberOfMaterials()>1) {
    print_error("*** Error: Need to update homogeneity. Currently, the incompressible flow solver does not allow"
                " more than one material.\n");
    exit_mpi();
  }
  if(Phi.size()>0 || NPhi.size()>0 || KappaPhi.size()>0 || L || Xi || subcycle>0 || dts != dt) {
    print_error("*** Error: Problem setup is not supported by TimeIntegratorPISO.\n");
    exit_mpi();
  }

  if(Vturb) {
    print_error("*** Error: Turbulence models are not supported by TimeIntegratorPISO.\n");
    exit_mpi();
  }
  double*** vturb = NULL; //not supporting turbulence right now

  GlobalMeshInfo &global_mesh(spo.GetGlobalMeshInfo());

  double*** id = ID.GetDataPointer();
  double*** homo = Homo.GetDataPointer();

  int iter(0), maxIter = time_step == 1 ? 10*iod.ts.semi_impl.maxIts : iod.ts.semi_impl.maxIts;
  vector<double> lin_rnorm; 
  vector<int>    lin_rnorm_its; 
  bool lin_success, converged(false);
  int nLinIts(0);
  double rel_err;

  print("  o Running the PISO procedure.\n");

  // Store V0
  V0.AXPlusBY(0.0, 1.0, V); //V0 = V (only copy nodes inside physical domain)
  Vec5D*** v0 = (Vec5D***)V0.GetDataPointer();

  Vec5D*** v = (Vec5D***)V.GetDataPointer();

  ExtractVariableComponents(v, &VXstar, &VYstar, &VZstar, &Pstar);

  //-----------------------------------------------------
  // Step 1: Solve the momentum equations for u*, v*, w*
  //-----------------------------------------------------

  // Solve the x-momentum equation
  if(global_mesh.x_glob.size()>1) {
    inco.BuildVelocityEquationSIMPLE(0, v0, v, id, vturb, homo, vlin_rows, B, DX, false, Efactor, dt, LocalDt);
    vlin_solver.SetLinearOperator(vlin_rows);
    lin_success = vlin_solver.Solve(B, VXstar, NULL, &nLinIts, &lin_rnorm, &lin_rnorm_its);
    if(!lin_success) {
      print_warning("    x Warning: Linear solver for the x-momentum equation failed to converge."
                    " Residual: %e -> %e.\n", lin_rnorm.front(), lin_rnorm.back());
      if(verbose>=2)
        for(int i=0; i<(int)lin_rnorm.size(); i++)
          print_warning("      > It. %d: residual = %e.\n", lin_rnorm_its[i]+1, lin_rnorm[i]);
    } else {
      if(verbose>=1)
        print("    * Solver of the x-momentum equation converged in %d iterations. Residual: %e -> %e.\n",
              nLinIts, lin_rnorm.front(), lin_rnorm.back());
    }
  }

  // Solve the y-momentum equation
  if(global_mesh.y_glob.size()>1) {
    inco.BuildVelocityEquationSIMPLE(1, v0, v, id, vturb, homo, vlin_rows, B, DY, false, Efactor, dt, LocalDt);
    vlin_solver.SetLinearOperator(vlin_rows);
    lin_success = vlin_solver.Solve(B, VYstar, NULL, &nLinIts, &lin_rnorm, &lin_rnorm_its);
    if(!lin_success) {
      print_warning("    x Warning: Linear solver for the y-momentum equation failed to converge."
                    " Residual: %e -> %e.\n", lin_rnorm.front(), lin_rnorm.back());
      if(verbose>=2)
        for(int i=0; i<(int)lin_rnorm.size(); i++)
          print_warning("      > It. %d: residual = %e.\n", lin_rnorm_its[i]+1, lin_rnorm[i]);
    } else {
      if(verbose>=1)
        print("    * Solver of the y-momentum equation converged in %d iterations. Residual: %e -> %e.\n",
              nLinIts, lin_rnorm.front(), lin_rnorm.back());
    }
  }

  // Solve the z-momentum equation
  if(global_mesh.z_glob.size()>1) {
    inco.BuildVelocityEquationSIMPLE(2, v0, v, id, vturb, homo, vlin_rows, B, DZ, false, Efactor, dt, LocalDt);
    vlin_solver.SetLinearOperator(vlin_rows);
    lin_success = vlin_solver.Solve(B, VZstar, NULL, &nLinIts, &lin_rnorm, &lin_rnorm_its);
    if(!lin_success) {
      print_warning("    x Warning: Linear solver for the z-momentum equation failed to converge."
                    " Residual: %e -> %e.\n", lin_rnorm.front(), lin_rnorm.back());
      if(verbose>=2)
        for(int i=0; i<(int)lin_rnorm.size(); i++)
          print_warning("      > It. %d: residual = %e.\n", lin_rnorm_its[i]+1, lin_rnorm[i]);
    } else {
      if(verbose>=1)
        print("    * Solver of the z-momentum equation converged in %d iterations. Residual: %e -> %e.\n",
              nLinIts, lin_rnorm.front(), lin_rnorm.back());
    }
  }

 
   
  //-----------------------------------------------------
  // Step 2: Solve the p' equation (First Corrector Step)
  //-----------------------------------------------------
  inco.BuildPressureEquationSIMPLE(v, homo, VXstar, VYstar, VZstar, DX, DY, DZ, plin_rows, B,
                                   fix_pressure_at_one_corner ? &ijk_zero_p : NULL);
  plin_solver.SetLinearOperator(plin_rows);
  Pprime.SetConstantValue(0.0, true); //!< This is p *correction*. Set init guess to 0 (Patankar 6.7-4)
  lin_success = plin_solver.Solve(B, Pprime, NULL, &nLinIts, &lin_rnorm, &lin_rnorm_its);
  if(!lin_success) {
    print_warning("    x Warning: Linear solver for the pressure correction equation failed to converge."
                  " Residual: %e -> %e.\n", lin_rnorm.front(), lin_rnorm.back());
    if(verbose>=2)
      for(int i=0; i<(int)lin_rnorm.size(); i++)
        print_warning("      > It. %d: residual = %e.\n", lin_rnorm_its[i]+1, lin_rnorm[i]);
  } else {
    if(verbose>=1)
      print("    * Solver of the pressure correction equation converged in %d iterations. Residual: %e -> %e.\n",
            nLinIts, lin_rnorm.front(), lin_rnorm.back());
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

    inco.CalculateVelocityTildePISO(0, v0, v, id, homo, VXprime, VXtildep, Efactor, dt, LocalDt);
    inco.CalculateVelocityTildePISO(1, v0, v, id, homo, VYprime, VYtildep, Efactor, dt, LocalDt);
    inco.CalculateVelocityTildePISO(2, v0, v, id, homo, VZprime, VZtildep, Efactor, dt, LocalDt);

    inco.BuildPressureEquationRHS_SIMPLER(v, homo, VXtildep, VYtildep, VZtildep, B,
                                          fix_pressure_at_one_corner ? &ijk_zero_p : NULL);
    plin_solver.UsePreviousPreconditioner(true); //The matrix A is still the same
    Pprime.SetConstantValue(0.0, true); //!< This is p *correction*. Set init guess to 0 (Patankar 6.7-4)
    lin_success = plin_solver.Solve(B, Pprime, NULL, &nLinIts, &lin_rnorm, &lin_rnorm_its);
    if(!lin_success) {
      print_warning("    x Warning: Linear solver for the pressure correction equation failed to converge."
                    " Residual: %e -> %e.\n", lin_rnorm.front(), lin_rnorm.back());
      if(verbose>=2)
        for(int i=0; i<(int)lin_rnorm.size(); i++)
          print_warning("      > It. %d: residual = %e.\n", lin_rnorm_its[i]+1, lin_rnorm[i]);
    } else {
      if(verbose>=1)
        print("    * Solver of the pressure correction equation converged in %d iterations."
              " Residual: %e -> %e.\n", nLinIts, lin_rnorm.front(), lin_rnorm.back());
    }

    rel_err = UpdateStatesCorrector(Pstar, Pprime, DX, DY, DZ, VXstar, VYstar, VZstar, 
                                    VXprime, VYprime, VZprime, &VXtildep, &VYtildep, &VZtildep); 

    if(rel_err<iod.ts.semi_impl.convergence_tolerance) {
      converged = true;
      break; 
    }

    print("  o It. %d: Relative error in velocity (2-norm): %e.\n", iter+1, rel_err);

    //NOTE: inco.ApplyBoundaryConditions not called. I think it is needed only for turbulent flows
  }

END_CORRECTORS:

  UpdateStatesFinal(v, Pstar, VXstar, VYstar, VZstar); // update v
  V.RestoreDataPointerAndInsert();
  inco.ApplyBoundaryConditions(V); 

  if(converged)
    print("  o Converged after %d iterations. Relative error in velocity (2-norm): %e.\n", iter+1, rel_err);
  else
    print_warning("  o Failed to converge. Relative error in velocity (2-norm): %e.\n", rel_err);


  if(sso) {
    assert(R3_ptr);
    inco.CalculateMomentumChanges(v0, V, id, *R3_ptr);
  }
   
  V0.RestoreDataPointerToLocalVector();

  ID.RestoreDataPointerToLocalVector();
  Homo.RestoreDataPointerToLocalVector();

  if(sso)
    sso->MonitorConvergence(*R3_ptr, ID);
    
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

  if(fabs(unorm)==0.0) {//new velocity is zero everywhere (return abs. error in this case)
    unorm = 1.0;
  }

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







