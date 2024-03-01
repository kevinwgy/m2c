/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<ExactRiemannSolverBase.h>
#include<array>
#include<utility> //std::pair
#include<bits/stdc++.h> //std::swap
#include<iostream>

//#include <chrono> // for timing

#ifndef WITHOUT_BOOST
#include<boost/math/tools/roots.hpp>
using namespace boost::math::tools;
#endif

using std::pair;
using std::cout;
using std::endl;

// for timing
//using std::chrono::high_resolution_clock;
//using std::chrono::duration_cast;
//using std::chrono::duration;
//using std::chrono::milliseconds;

extern int verbose;

#define INVALID_MATERIAL_ID -1
//-----------------------------------------------------

ExactRiemannSolverBase::ExactRiemannSolverBase(std::vector<VarFcnBase*> &vf_, 
    ExactRiemannSolverData &iod_riemann_) : vf(vf_), iod_riemann(iod_riemann_)
{
  maxIts_main          = iod_riemann.maxIts_main;
  maxIts_bracket       = iod_riemann.maxIts_bracket;
  maxIts_shock         = iod_riemann.maxIts_shock;
  numSteps_rarefaction = iod_riemann.numSteps_rarefaction;
  tol_main             = iod_riemann.tol_main; //applied to both pressure and velocity!
  tol_shock            = iod_riemann.tol_shock; //a density tolerance
  tol_rarefaction      = iod_riemann.tol_rarefaction; //a non-dimensional tolerance
  min_pressure         = iod_riemann.min_pressure;
  failure_threshold    = iod_riemann.failure_threshold;
  pressure_at_failure  = iod_riemann.pressure_at_failure;
  surface_tension      = iod_riemann.surface_tension == ExactRiemannSolverData::YES;
  integrationPath1.reserve(500);
  integrationPath3.reserve(500);
}

//-----------------------------------------------------
/** Solves the one-dimensional Riemann problem. Extension of Kamm 2015 
 * to Two Materials. See KW's notes for details
 * Returns an integer error code
 * 0: no errors
 * 1: riemann solver failed to find a bracketing interval
 */
int
ExactRiemannSolverBase::ComputeRiemannSolution(double *dir, 
    double *Vm, int idl /*"left" state*/, 
    double *Vp, int idr /*"right" state*/, 
    double *Vs, int &id /*solution at xi = 0 (i.e. x=0) */,
    double *Vsm /*left 'star' solution*/,
    double *Vsp /*right 'star' solution*/,
    double curvature)
{
  assert(curvature == 0.0); //the base class does not handle curvature!

  //std::cout << "ExactRiemannSolverBase::ComputeRiemannSolution: this is the base version!" << std::endl;
  // Convert to a 1D problem (i.e. One-Dimensional Riemann)
  double rhol  = Vm[0];
  double ul    = Vm[1]*dir[0] + Vm[2]*dir[1] + Vm[3]*dir[2];
  double pl    = Vm[4];
  double rhor  = Vp[0];
  double ur    = Vp[1]*dir[0] + Vp[2]*dir[1] + Vp[3]*dir[2];
  double pr    = Vp[4];
  //fprintf(stdout,"1DRiemann: left = %e %e %e (%d) : right = %e %e %e (%d)\n", rhol, ul, pl, idl, rhor, ur, pr, idr);

  integrationPath1.clear();
  integrationPath3.clear();
  std::vector<double> vectL{pl, rhol, ul};
  std::vector<double> vectR{pr, rhor, ur};
  integrationPath1.push_back(vectL);
  integrationPath3.push_back(vectR); 

#if PRINT_RIEMANN_SOLUTION == 1
  std::cout << "Left State (rho, u, p): " << rhol << ", " << ul << ", " << pl << "." << std::endl;
  std::cout << "Right State (rho, u, p): " << rhor << ", " << ur << ", " << pr << "." << std::endl;
#endif

  double el = vf[idl]->GetInternalEnergyPerUnitMass(rhol, pl);
  double cl = vf[idl]->ComputeSoundSpeedSquare(rhol, el);

  if(rhol<=0 || cl<0) {
    fprintf(stdout,"*** Error: Negative density or c^2 (square of sound speed) in ComputeRiemannSolution(l)." 
	" rho = %e, u = %e, p = %e, e = %e, c^2 = %e, ID = %d.\n",
	rhol, ul, pl, el, cl, idl);
    exit(-1);
  } else
    cl = sqrt(cl);

  double er = vf[idr]->GetInternalEnergyPerUnitMass(rhor, pr);
  double cr = vf[idr]->ComputeSoundSpeedSquare(rhor, er);

  if(rhor<=0 || cr<0) {
    fprintf(stdout,"*** Error: Negative density or c^2 (square of sound speed) in ComputeRiemannSolution(r)."
	" rho = %e, u = %e, p = %e, e = %e, c^2 = %e, ID = %d.\n",
	rhor, ur, pr, er, cr, idr);
    exit(-1);
  } else
    cr = sqrt(cr);


  // Declare variables in the "star region"
  double p0(DBL_MIN), ul0(0.0), ur0(0.0), rhol0(DBL_MIN), rhor0(DBL_MIN);
  double p1(DBL_MIN), ul1(0.0), ur1(0.0), rhol1(DBL_MIN), rhor1(DBL_MIN); //Secant Method ("k-1","k" in Kamm, (19))
  double p2(DBL_MIN), ul2(0.0), ur2(0.0), rhol2(DBL_MIN), rhor2(DBL_MIN); // "k+1"
  double f0(0.0), f1(0.0), f2(0.0); //difference between ul and ur


  // -------------------------------
  // Now, Solve The Riemann Problem
  // -------------------------------
  bool success = true;

  // monitor if the solution involves a transonic rarefaction. This is special as the solution 
  // at xi = x = 0 is within the rarefaction fan.
  bool trans_rare = false;
  double Vrare_x0[3]; //rho, u, and p at x = 0, in the case of a transonic rarefaction


  // A Trivial Case
  if(ul == ur && pl == pr) {
    FinalizeSolution(dir, Vm, Vp, rhol, ul, pl, idl, rhor, ur, pr, idr, rhol, rhor, ul, pl, 
	trans_rare, Vrare_x0, //inputs
	Vs, id, Vsm, Vsp/*outputs*/);
    return 0;
  }

  // -------------------------------
  // Step 1: Initialization
  //         (find initial interval [p0, p1])
  // -------------------------------
  success = FindInitialInterval(rhol, ul, pl, el, cl, idl, rhor, ur, pr, er, cr, idr, /*inputs*/
      p0, rhol0, rhor0, ul0, ur0,
      p1, rhol1, rhor1, ul1, ur1/*outputs*/);
  /* our convention is that p0 < p1 */

  if(!success) { //failed to find a bracketing interval. Output the state corresponding smallest "f"

    // get sol1d, trans_rare and Vrare_x0
#if PRINT_RIEMANN_SOLUTION == 1
    sol1d.clear();
#endif
    success = ComputeRhoUStar(1, integrationPath1, rhol, ul, pl, p1, idl, rhol0, rhol0*1.1, rhol2, ul2,
	&trans_rare, Vrare_x0/*filled only if found a trans. rarefaction*/);
    success = success && ComputeRhoUStar(3, integrationPath3, rhor, ur, pr, p1, idr, rhor0, rhol0*1.1, rhor2, ur2,
	&trans_rare, Vrare_x0/*filled only if found a trans. rarefaction*/);

    if(!success) {
      cout << "Warning: Riemann solver failed to find an initial bracketing interval (Returning initial states as solution)." << endl;
      for(int i=0; i<5; i++) {
	Vsm[i] = Vm[i];
	Vsp[i] = Vp[i];
      }
      double u_avg = 0.5*(ul+ur);
      if(u_avg >= 0) {
	id = idl;
	for(int i=0; i<5; i++)
	  Vs[i] = Vm[i];
      } else {
	id = idr;
	for(int i=0; i<5; i++)
	  Vs[i] = Vp[i];
      }

      ExactRiemannSolverNonAdaptive riemannNonAdaptive(vf, iod_riemann); 
      int retryRiemann = riemannNonAdaptive.ComputeRiemannSolution(dir, Vm, idl, Vp, idr, Vs, id, Vsm, Vsp);
      if(verbose>=1)
	cout << "Warning: Riemann solver failed to find an initial bracketing interval. Activated the non-adaptive version." << endl; 
      return retryRiemann;
    }

    FinalizeSolution(dir, Vm, Vp, rhol, ul, pl, idl, rhor, ur, pr, idr, rhol2, rhor2, 0.5*(ul2+ur2), p1,
	trans_rare, Vrare_x0, /*inputs*/
	Vs, id, Vsm, Vsp /*outputs*/);

    ExactRiemannSolverNonAdaptive riemannNonAdaptive(vf, iod_riemann); 
    int retryRiemann = riemannNonAdaptive.ComputeRiemannSolution(dir, Vm, idl, Vp, idr, Vs, id, Vsm, Vsp);
    if(verbose>=1)
      cout << "Warning: Riemann solver failed to find an initial bracketing interval. Activated the non-adaptive version." << endl;
    return retryRiemann;
  }

  f0 = ul0 - ur0;
  f1 = ul1 - ur1;

#if PRINT_RIEMANN_SOLUTION == 1
  cout << "Found initial interval: p0 = " << p0 << ", f0 = " << f0 << ", p1 = " << p1 << ", f1 = " << f1 << endl;
#endif

  // -------------------------------
  // Step 2: Main Loop (Secant Method, Safeguarded) 
  // -------------------------------
  double denom = 0; 

  int iter = 0;
  double err_p = 1.0, err_u = 1.0;

  //p2 (and f2) is always the latest one
  p2 = p1; 
  f2 = f1;

#if PRINT_RIEMANN_SOLUTION == 1
  sol1d.clear();
#endif

  for(iter=0; iter<maxIts_main; iter++) {

    // 2.1: Update p using the Brent method (safeguarded secant method)
    denom = f1 - f0;
    if(denom == 0) {
      cout << "Warning: Division-by-zero while using the secant method to solve the Riemann problem." << endl;
      cout << "         left state: " << rhol << ", " << ul << ", " << pl << ", " << idl << " | right: " 
	<< rhor << ", " << ur << ", " << pr << ", " << idr << endl;
      cout << "         dir = " << dir[0] << "," << dir[1] << "," << dir[2]
	<< ", f0 = " << f0 << ", f1 = " << f1 << endl;
      p2 = 0.5*(p0+p1);
    } else {
      p2 = p2 - f2*(p1-p0)/denom;  // update p2
      if(p2<=p0 || p2>=p1) //discard and switch to bisection
	p2 = 0.5*(p0+p1);
    }

    //fprintf(stdout,"iter = %d, p0 = %e, p1 = %e, p2 = %e, f0 = %e, f1 = %e.\n", iter, p0, p1, p2, f0, f1);


try_again:

    // 2.2: Calculate ul2, ur2 
    success = ComputeRhoUStar(1, integrationPath1, rhol, ul, pl, p2, idl/*inputs*/, 
	rhol0, rhol1/*initial guesses for Hugo. eq.*/,
	rhol2, ul2/*outputs*/, 
	&trans_rare, Vrare_x0/*filled only if found a trans. rarefaction*/);

    if(!success) {
      //fprintf(stdout,"*** Error: Exact Riemann solver failed. left: %e %e %e (%d) | right: %e %e %e (%d).\n", 
      //        rhol, ul, pl, idl, rhor, ur, pr, idr);
      p2 = 0.5*(p2+p0); //move closer to p0

      iter++;
      if(iter<maxIts_main)
	goto try_again;
      else
	break;
    }

    success = ComputeRhoUStar(3, integrationPath3, rhor, ur, pr,  p2, idr/*inputs*/, 
	rhor0, rhor1/*initial guesses for Hugo. erq.*/,
	rhor2, ur2/*outputs*/,
	&trans_rare, Vrare_x0/*filled only if found a trans. rarefaction*/);

    if(!success) {
      //fprintf(stdout,"*** Error: Exact Riemann solver failed (2). left: %e %e %e (%d) | right: %e %e %e (%d).\n", 
      //        rhol, ul, pl, idl, rhor, ur, pr, idr);
      p2 = 0.5*(p2+p1); //move closer to p1

      iter++;
      if(iter<maxIts_main)
	goto try_again;
      else
	break;
    }

    f2 = ul2 - ur2;


    // 2.3: Update for the next iteration
    if(f0*f2<0.0) {
      p1 = p2;
      f1 = f2;
      rhol1 = rhol2;
      rhor1 = rhor2;
    } else {
      p0 = p2;
      f0 = f2;
      rhol0 = rhol2;
      rhor0 = rhor2;
    }


    // 2.4: Check stopping criterion
    err_p = fabs(p1 - p0)/std::max(fabs(pl + 0.5*rhol*ul*ul), fabs(pr + 0.5*rhor*ur*ur));
    err_u = fabs(f2)/std::max(cl, cr);

#if PRINT_RIEMANN_SOLUTION == 1
    cout << "Iter " << iter << ": p-interval = [" << p0 << ", " << p1 << "], p2 = " << p2 << ", err_p = " << err_p 
      << ", err_u = " << err_u << "." << endl;
#endif

    if( (err_p < tol_main && err_u < tol_main) || (err_p < tol_main*1e-3) || (err_u < tol_main*1e-3) )
      break; // converged

    trans_rare = false; //reset

#if PRINT_RIEMANN_SOLUTION == 1
    sol1d.clear();
#endif

  }


  // -------------------------------
  // Step 3: Find state at xi = x = 0 (for output)
  // -------------------------------
  double u2 = 0.5*(ul2 + ur2);
  FinalizeSolution(dir, Vm, Vp, rhol, ul, pl, idl, rhor, ur, pr, idr, rhol2, rhor2, u2, p2, 
      trans_rare, Vrare_x0, //inputs
      Vs, id, Vsm, Vsp/*outputs*/);


  if(iter == maxIts_main) {
    if(verbose>=1) {
      cout << "Warning: Exact Riemann solver failed to converge. err_p = " << err_p
	<< ", err_u = " << err_u << "." << endl;
      cout << "    Vm = [" << Vm[0] << ", " << Vm[1] << ", " << Vm[2] << ", " << Vm[3] << ", " << Vm[4] << "] (" << idl
	<< "),  Vp = [" << Vp[0] << ", " << Vp[1] << ", " << Vp[2] << ", " << Vp[3] << ", " << Vp[4] << "] (" << idr
	<< ")" << endl;
    }

    ExactRiemannSolverNonAdaptive riemannNonAdaptive(vf, iod_riemann);
    int retryRiemann = riemannNonAdaptive.ComputeRiemannSolution(dir, Vm, idl, Vp, idr, Vs, id, Vsm, Vsp);
    if(verbose>=1)
      cout << "Warning: Exact Riemann solver (adaptive) failed to converge. Activated the non-adaptive version." << endl;
    return retryRiemann;
  }

#if PRINT_RIEMANN_SOLUTION == 1
  std::cout << "Star State: (rhols, rhors, us, ps): " << rhol2 << ", " << rhor2 << ", " << u2 << ", " << p2 << "." << std::endl;
#endif


  //success!
  return 0;

}

//-----------------------------------------------------

  void
ExactRiemannSolverBase::FinalizeSolution(double *dir, double *Vm, double *Vp,
    double rhol, double ul, double pl, int idl, 
    double rhor, double ur, double pr, int idr, 
    double rhol2, double rhor2, double u2, double p2,
    bool trans_rare, double Vrare_x0[3], /*inputs*/
    double *Vs, int &id, double *Vsm, double *Vsp /*outputs*/)
{
  // find tangential velocity from input
  double utanl[3] = {Vm[1]-ul*dir[0], Vm[2]-ul*dir[1], Vm[3]-ul*dir[2]};
  double utanr[3] = {Vp[1]-ur*dir[0], Vp[2]-ur*dir[1], Vp[3]-ur*dir[2]};

  // find material id at xi = x = 0
  if(u2>=0)
    id = idl;
  else
    id = idr;

#if PRINT_RIEMANN_SOLUTION == 1
  // the 2-wave
  sol1d.push_back(vector<double>{u2 - std::max(1e-6, 0.001*fabs(u2)), rhol2, u2, p2, (double)idl});
  sol1d.push_back(vector<double>{u2, rhor2, u2, p2, (double)idr});
  // 1- and 3- waves
  integrationPath1.clear();
  integrationPath3.clear();
  std::vector<double> vectL{pl, rhol, ul};
  std::vector<double> vectR{pr, rhor, ur};
  integrationPath1.push_back(vectL);
  integrationPath3.push_back(vectR); 
 
  bool success;
  double ul2_tmp, ur2_tmp, rhol2_tmp, rhor2_tmp;
  success = ComputeRhoUStar(1, integrationPath1, rhol, ul, pl, p2, idl/*inputs*/, 
      rhol2, 0.9*rhol2/*initial guesses for Hugo. eq.*/,
      rhol2_tmp, ul2_tmp/*outputs*/, 
      &trans_rare, Vrare_x0/*filled only if found a trans. rarefaction*/);
  if (!success) {
    std::cout <<  "*** Error: ComputeRhoUStar(1) failed when finalizng the solution." << std::endl;
    exit(-1);
  } 
  success = ComputeRhoUStar(3, integrationPath3, rhor, ur, pr,  p2, idr/*inputs*/, 
      rhor2, 0.9*rhor2/*initial guesses for Hugo. erq.*/,
      rhor2_tmp, ur2_tmp/*outputs*/,
      &trans_rare, Vrare_x0/*filled only if found a trans. rarefaction*/);
  if (!success) {
    std::cout <<  "*** Error: ComputeRhoUStar(3) failed when finalizng the solution." << std::endl;
    exit(-1);
  } 
#endif

  Vs[0] = Vs[1] = Vs[2] = Vs[3] = Vs[4] = 0.0;

  if(trans_rare) {
    Vs[0] = Vrare_x0[0];
    Vs[1] = Vrare_x0[1]*dir[0];
    Vs[2] = Vrare_x0[1]*dir[1];
    Vs[3] = Vrare_x0[1]*dir[2];
    Vs[4] = Vrare_x0[2];
  }
  else { 
    //find state variables at xi = x = 0

    if(u2>=0) { //either Vl or Vlstar --- check the 1-wave

      bool is_star_state = false;

      if(pl >= p2) {//1-wave is rarefaction
	double el2 = vf[idl]->GetInternalEnergyPerUnitMass(rhol2, p2);
	double cl2 = vf[idl]->ComputeSoundSpeedSquare(rhol2, el2);

	if(rhol2<=0 || cl2<0) {
	  fprintf(stdout,"*** Error: Negative density or c^2 (square of sound speed) in ComputeRiemannSolution(l2)."
	      " rho = %e, p = %e, e = %e, c^2 = %e, id = %d.\n",
	      rhol2, pl, el2, cl2, idl);
	  exit(-1);
	} else
	  cl2 = sqrt(cl2);

	if(u2 - cl2 <= 0) //rarefaction tail speed
	  is_star_state = true;
      } 
      else {//1-wave is shock
	double us = (rhol2*u2 - rhol*ul)/(rhol2 - rhol); //shock speed
	if(us <= 0)      
	  is_star_state = true;
      }

      if(is_star_state) {
	Vs[0] = rhol2;
	Vs[1] = u2*dir[0];
	Vs[2] = u2*dir[1];
	Vs[3] = u2*dir[2];
	Vs[4] = p2;
      } else {
	Vs[0] = rhol;
	Vs[1] = ul*dir[0];
	Vs[2] = ul*dir[1];
	Vs[3] = ul*dir[2];
	Vs[4] = pl;
      }

    } else { //either Vr or Vrstar --- check the 3-wave

      bool is_star_state = false;

      if(pr >= p2) {//3-wave is rarefaction
	double er2 = vf[idr]->GetInternalEnergyPerUnitMass(rhor2, p2);
	double cr2 = vf[idr]->ComputeSoundSpeedSquare(rhor2, er2);

	if(rhor2<=0 || cr2<0) {
	  fprintf(stdout,"*** Error: Negative density or c^2 (square of sound speed) in ComputeRiemannSolution(r2)." 
	      " rho = %e, p = %e, e = %e, c^2 = %e, id = %d.\n",
	      rhor2, p2, er2, cr2, idr);
	  exit(-1);
	} else
	  cr2 = sqrt(cr2);

	if(u2 - cr2 >= 0)
	  is_star_state = true;
      }
      else {//3-wave is shock
	double us = (rhor2*u2 - rhor*ur)/(rhor2 - rhor); //shock speed
	if(us >= 0)
	  is_star_state = true;
      }

      if(is_star_state) {
	Vs[0] = rhor2;
	Vs[1] = u2*dir[0];
	Vs[2] = u2*dir[1];
	Vs[3] = u2*dir[2];
	Vs[4] = p2;
      } else {
	Vs[0] = rhor;
	Vs[1] = ur*dir[0];
	Vs[2] = ur*dir[1];
	Vs[3] = ur*dir[2];
	Vs[4] = pr;
      }

    }
  }

  // determine the tangential components of velocity -- upwinding
  if(u2>0) {
    for(int i=1; i<=3; i++)
      Vs[i] += utanl[i-1];
  } 
  else if(u2<0) {
    for(int i=1; i<=3; i++)
      Vs[i] += utanr[i-1];
  } 
  else {//u2 == 0
    for(int i=1; i<=3; i++)
      Vs[i] += 0.5*(utanl[i-1]+utanr[i-1]);
  }


  // determine Vsm and Vsp, i.e. the star states on the minus and plus sides of the contact discontinuity
  Vsm[0] = rhol2;
  Vsm[1] = utanl[0] + u2*dir[0];
  Vsm[2] = utanl[1] + u2*dir[1];
  Vsm[3] = utanl[2] + u2*dir[2];
  Vsm[4] = p2;
  Vsp[0] = rhor2;
  Vsp[1] = utanr[0] + u2*dir[0];
  Vsp[2] = utanr[1] + u2*dir[1];
  Vsp[3] = utanr[2] + u2*dir[2];
  Vsp[4] = p2;



#if PRINT_RIEMANN_SOLUTION == 1
  std::sort(sol1d.begin(), sol1d.end(), 
      [](vector<double> v1, vector<double> v2){return v1[0]<v2[0];});
  int last = sol1d.size()-1;
  double xi_span = sol1d[last][0] - sol1d[0][0];
  sol1d.insert(sol1d.begin(), vector<double>{sol1d[0][0]-xi_span, sol1d[0][1], sol1d[0][2], sol1d[0][3], sol1d[0][4]});
  last++;
  sol1d.push_back(vector<double>{sol1d[last][0]+xi_span, sol1d[last][1], sol1d[last][2], sol1d[last][3], sol1d[last][4]});

  FILE* solFile = fopen("RiemannSolution.txt", "w");
  print(solFile, "## One-Dimensional Riemann Problem.\n");
  print(solFile, "## Initial State: %e %e %e, id %d (left) | (right) %e %e %e, id %d.\n", 
      rhol, ul, pl, idl, rhor, ur, pr, idr);
  print(solFile, "## xi(x/t) | density | velocity | pressure | internal energy per mass | material id\n");

  for(auto it = sol1d.begin(); it != sol1d.end(); it++) 
    print(solFile,"% e    % e    % e    % e    % e    % d\n", (*it)[0], (*it)[1], (*it)[2], (*it)[3], 
	vf[(int)(*it)[4]]->GetInternalEnergyPerUnitMass((*it)[1], (*it)[3]), (int)(*it)[4]);

  fclose(solFile);
#endif

}

//-----------------------------------------------------

void
ExactRiemannSolverBase::FinalizeOneSidedSolution(double *dir, double *Vm, 
    double rhol, double ul, double pl, int idl, 
    double rhol2, double u2/*ustar*/, double p2,
    bool trans_rare, double Vrare_x0[3], /*inputs*/
    double *Vs, int &id, double *Vsm /*outputs*/)
{
  // find tangential velocity from input
  double utanl[3] = {Vm[1]-ul*dir[0], Vm[2]-ul*dir[1], Vm[3]-ul*dir[2]};

  // find material id at xi = x = 0
  if(u2>=0)
    id = idl;
  else
    id = INVALID_MATERIAL_ID;

#if PRINT_RIEMANN_SOLUTION == 1
  // the 2-wave
  sol1d.push_back(vector<double>{u2, rhol2, u2, p2, (double)idl});
#endif

  if(id != INVALID_MATERIAL_ID) {

    Vs[0] = Vs[1] = Vs[2] = Vs[3] = Vs[4] = 0.0;

    if(trans_rare) {
      Vs[0] = Vrare_x0[0];
      Vs[1] = Vrare_x0[1]*dir[0];
      Vs[2] = Vrare_x0[1]*dir[1];
      Vs[3] = Vrare_x0[1]*dir[2];
      Vs[4] = Vrare_x0[2];
    }
    else { 
      //find state variables at xi = x = 0

      bool is_star_state = false;

      if(pl >= p2) {//1-wave is rarefaction
	double el2 = vf[idl]->GetInternalEnergyPerUnitMass(rhol2, p2);
	double cl2 = vf[idl]->ComputeSoundSpeedSquare(rhol2, el2);

	if(rhol2<=0 || cl2<0) {
	  fprintf(stdout,"*** Error: Negative density or c^2 (square of sound speed) in ComputeRiemannSolution(l2)."
	      " rho = %e, p = %e, e = %e, c^2 = %e, id = %d.\n",
	      rhol2, pl, el2, cl2, idl);
	  exit(-1);
	} else
	  cl2 = sqrt(cl2);

	if(u2 - cl2 <= 0) //rarefaction tail speed
	  is_star_state = true;
      } 
      else {//1-wave is shock
	double us = (rhol2*u2 - rhol*ul)/(rhol2 - rhol); //shock speed
	if(us <= 0)      
	  is_star_state = true;
      }

      if(is_star_state) {
	Vs[0] = rhol2;
	Vs[1] = u2*dir[0];
	Vs[2] = u2*dir[1];
	Vs[3] = u2*dir[2];
	Vs[4] = p2;
      } else {
	Vs[0] = rhol;
	Vs[1] = ul*dir[0];
	Vs[2] = ul*dir[1];
	Vs[3] = ul*dir[2];
	Vs[4] = pl;
      }
    }

    // determine the tangential components of velocity -- upwinding
    for(int i=1; i<=3; i++)
      Vs[i] += utanl[i-1];
  }


  // determine Vsm, i.e. the star states on the minus side of the contact discontinuity
  Vsm[0] = rhol2;
  Vsm[1] = utanl[0] + u2*dir[0];
  Vsm[2] = utanl[1] + u2*dir[1];
  Vsm[3] = utanl[2] + u2*dir[2];
  Vsm[4] = p2;


#if PRINT_RIEMANN_SOLUTION == 1
  std::sort(sol1d.begin(), sol1d.end(), 
      [](vector<double> v1, vector<double> v2){return v1[0]<v2[0];});
  int last = sol1d.size()-1;
  double xi_span = sol1d[last][0] - sol1d[0][0];
  sol1d.insert(sol1d.begin(), vector<double>{sol1d[0][0]-xi_span, sol1d[0][1], sol1d[0][2], sol1d[0][3], sol1d[0][4]});

  FILE* solFile = fopen("RiemannSolution.txt", "w");
  print(solFile, "## One-Dimensional Riemann Problem.\n");
  print(solFile, "## Initial State: %e %e %e, id %d (left) | wall velocity: %e.\n", 
      rhol, ul, pl, idl, u2);
  print(solFile, "## xi(x/t) | density | velocity | pressure | internal energy per mass | material id\n");

  for(auto it = sol1d.begin(); it != sol1d.end(); it++) 
    print(solFile,"% e    % e    % e    % e    % e    % d\n", (*it)[0], (*it)[1], (*it)[2], (*it)[3], 
	vf[(int)(*it)[4]]->GetInternalEnergyPerUnitMass((*it)[1], (*it)[3]), (int)(*it)[4]);

  fclose(solFile);
#endif

}


//----------------------------------------------------------------------------------
//! find a bracketing interval [p0, p1] (f0*f1<=0)
  bool
ExactRiemannSolverBase::FindInitialInterval(double rhol, double ul, double pl, double el, double cl, int idl,
    double rhor, double ur, double pr, double er, double cr, int idr, /*inputs*/
    double &p0, double &rhol0, double &rhor0, double &ul0, double &ur0,
    double &p1, double &rhol1, double &rhor1, double &ul1, double &ur1/*outputs*/)
{
  /*convention: p0 < p1*/

  bool success = true;

  // Step 1: Find two feasible points (This step should never fail)
  success = FindInitialFeasiblePoints(rhol, ul, pl, el, cl, idl, rhor, ur, pr, er, cr, idr, /*inputs*/
      p0, rhol0, rhor0, ul0, ur0, p1, rhol1, rhor1, ul1, ur1/*outputs*/);

  if(!success) {//This should never happen (unless user's inputs have errors)!
    p0 = p1 = pressure_at_failure;
    return false;
  }

#if PRINT_RIEMANN_SOLUTION == 1
  fprintf(stdout, "Found two initial points: p0 = %e, f0 = %e, p1 = %e, f1 = %e.\n", p0, ul0-ur0, p1, ul1-ur1);
  fprintf(stdout, "Searching for a bracketing interval...\n");
#endif

  // Step 2: Starting from the two feasible points, find a bracketing interval
  //         This step may fail, which indicates a solution may not exist for arbitrary left & right states
  //         If this happens, return the point with smallest absolute value of "f"

  double p2, rhol2, rhor2, ul2, ur2;
  double f0, f1;
  int i;

  // These are variables corresponding to the smallest magnitude of "f" --- used only if a bracketing interval
  // cannot be found. This is to just to minimize the chance of code crashing...
  double fmin, p_fmin, rhol_fmin, rhor_fmin, ul_fmin, ur_fmin;
  if(fabs(ul0-ur0) < fabs(ul1-ur1)) {
    fmin = fabs(ul0 - ur0);
    p_fmin = p0;  rhol_fmin = rhol0;  rhor_fmin = rhor0;  ul_fmin = ul0;  ur_fmin = ur0;
  } else {
    fmin = fabs(ul1 - ur1);
    p_fmin = p1;  rhol_fmin = rhol1;  rhor_fmin = rhor1;  ul_fmin = ul1;  ur_fmin = ur1;
  }

  for(i=0; i<maxIts_bracket; i++) {

    f0 = ul0 - ur0;
    f1 = ul1 - ur1;

    if(f0*f1<=0.0)
      return true;

    // find a physical p2 that has the opposite sign
    if(fabs(f0-f1)>1e-9) {
      p2 = p1 - f1*(p1-p0)/(f1-f0); //the Secant method
      if(p2<p0)
	p2 -= 0.1*(p1-p0);
      else //p2 cannot be between p0 and p1, so p2>p1
	p2 += 0.1*(p1-p0);
    } else {//f0 == f1
      p2 = 1.1*p1; 
      //p2 = p1 + *(p1-p0); 
    }

    if(p2<min_pressure || i==int(maxIts_bracket/2) ) {//does not look right. reset to a small non-negative pressure
      p2 = 1.0e-8; 
    }

    success = ComputeRhoUStar(1, integrationPath1, rhol, ul, pl, p2, idl, rhol0, rhol1, rhol2, ul2);
    // compute the 3-wave only if the 1-wave is succeeded
    success = success && ComputeRhoUStar(3, integrationPath3, rhor, ur, pr,  p2, idr, rhor0, rhor1, rhor2, ur2);

    if(!success) {

#if PRINT_RIEMANN_SOLUTION == 1
      fprintf(stdout, "  -- p2 = %e (failed)\n", p2);
#endif
      //move closer to [p0, p1]
      for(int j=0; j<maxIts_bracket/2; j++) {
	if(p2<p0)     
	  p2 = p0 - 0.5*(p0-p2);
	else //p2>p1
	  p2 = p1 + 0.5*(p2-p1);

	success = ComputeRhoUStar(1, integrationPath1, rhol, ul, pl, p2, idl, rhol0, rhol1, rhol2, ul2);
	// compute the 3-wave only if the 1-wave is succeeded
	success = success && ComputeRhoUStar(3, integrationPath3, rhor, ur, pr,  p2, idr, rhor0, rhor1, rhor2, ur2);


	if(success)
	  break;
      }
    }

    if(!success)
      break;

    // check & update fmin (only used if the Riemann solver fails -- a "failsafe" feature)
    if(fabs(ul2-ur2)<fmin) {
      fmin = fabs(ul2-ur2);
      p_fmin = p2;  rhol_fmin = rhol2;  rhor_fmin = rhor2;  ul_fmin = ul2;  ur_fmin = ur2;
    }

    // update p0 or p1
    if(p2<p0) {
      p1 = p0;  rhol1 = rhol0;  rhor1 = rhor0;  ul1 = ul0;  ur1 = ur0;
      p0 = p2;  rhol0 = rhol2;  rhor0 = rhor2;  ul0 = ul2;  ur0 = ur2;
    } else {//p2>p1
      p0 = p1;  rhol0 = rhol1;  rhor0 = rhor1;  ul0 = ul1;  ur0 = ur1;
      p1 = p2;  rhol1 = rhol2;  rhor1 = rhor2;  ul1 = ul2;  ur1 = ur2;
    }

#if PRINT_RIEMANN_SOLUTION == 1
    fprintf(stdout, "  -- p0 = %e, f0 = %e, p1 = %e, f1 = %e (success)\n", p0, ul0-ur0, p1, ul1-ur1);
#endif

  }

  if(!success || i==maxIts_bracket) {
    if(verbose>=1) {
      cout << "Warning: Exact Riemann solver failed. (Unable to find a bracketing interval) " << endl;
      cout << "   left: " << std::setprecision(10) << rhol << ", " << std::setprecision(10) << ul << ", " 
	<< std::setprecision(10) << pl << " (" << idl << "); right: "
	<< std::setprecision(10) << rhor << ", " << std::setprecision(10) << ur << ", " 
	<< std::setprecision(10) << pr << " (" << idr << ")."
	<< " Residual (|ulstar-urstar|): " << std::setprecision(10) << fmin << endl;
    }
    if(fmin<failure_threshold*fabs(ul-ur)) {
      if(verbose>=1) 
	cout << "*** Best approximate solution: rhols = " << rhol_fmin << ", ps = " << p_fmin << ", us = ("
	  << ul_fmin << "(l) + " << ur_fmin << "(r))/2, rhors = " << rhor_fmin << "." << endl;
      p0    = p1    = p_fmin;
      rhol0 = rhol1 = rhol_fmin;
      rhor0 = rhor1 = rhor_fmin;
      ul0   = ul1   = ul_fmin; 
      ur0   = ur1   = ur_fmin;
    } else { //it could be that the Riemann problem has no solution!
      p2 = pressure_at_failure;
      success = ComputeRhoUStar(1, integrationPath1, rhol, ul, pl, p2, idl, rhol0, rhol1, rhol2, ul2);
      // compute the 3-wave only if the 1-wave is succeeded
      success = success && ComputeRhoUStar(3, integrationPath3, rhor, ur, pr,  p2, idr, rhor0, rhor1, rhor2, ur2);
      if(success) {
	if(verbose >= 1)
	  cout << "*** Prescribed solution: rhols = " << rhol2 << ", ps = " << p2 << ", us = ("
	    << ul2 << "(l) + " << ur2 << "(r))/2, rhors = " << rhor2 << "." << endl;
	p0    = p1    = p2;
	rhol0 = rhol1 = rhol2; 
	rhor0 = rhor1 = rhor2;
	ul0   = ul1   = ul2; 
	ur0   = ur1   = ur2;
      } else {
	if(verbose >= 1)
	  cout << "*** Best approximation: rhols = " << rhol_fmin << ", ps = " << p_fmin << ", us = ("
	    << ul_fmin << "(l) + " << ur_fmin << "(r))/2, rhors = " << rhor_fmin << "." << endl;
	p0    = p1    = p_fmin;
	rhol0 = rhol1 = rhol_fmin;
	rhor0 = rhor1 = rhor_fmin;
	ul0   = ul1   = ul_fmin; 
	ur0   = ur1   = ur_fmin;
      }
    }

    return false;
  }

  return true;
}

//----------------------------------------------------------------------------------
//! find a bracketing interval [p0, p1] (f0*f1<=0)
  bool
ExactRiemannSolverBase::FindInitialIntervalOneSided(double rhol, double ul, double pl, double el, double cl, int idl,
    double ustar, double &p0, double &rhol0, double &ul0, 
    double &p1, double &rhol1, double &ul1)
{

  assert(ul>ustar); //this function is only needed (and applicable) when there is a shock

  /*convention: p0 < p1*/

  bool success = true;

  // Step 1: Find two feasible points (This step should never fail)
  success = FindInitialFeasiblePointsOneSided(rhol, ul, pl, el, cl, idl, ustar, /*inputs*/
      p0, rhol0, ul0, p1, rhol1, ul1/*outputs*/);

  if(!success) {//This should never happen (unless user's inputs have errors)!
    p0 = p1 = pressure_at_failure;
    return false;
  }

  assert(p0>=pl);
  assert(p1>=pl);

#if PRINT_RIEMANN_SOLUTION == 1
  fprintf(stdout, "Found two initial points: p0 = %e, f0 = %e, p1 = %e, f1 = %e.\n", p0, ul0-ustar, p1, ul1-ustar);
  fprintf(stdout, "Searching for a bracketing interval...\n");
#endif

  // Step 2: Starting from the two feasible points, find a bracketing interval
  //         This step may fail, which indicates a solution may not exist for arbitrary left & right states
  //         If this happens, return the point with smallest absolute value of "f"

  double p2, rhol2, ul2;
  double f0, f1;
  int i;

  // These are variables corresponding to the smallest magnitude of "f" --- used only if a bracketing interval
  // cannot be found. This is to just to minimize the chance of code crashing...
  double fmin, p_fmin, rhol_fmin, ul_fmin;
  if(fabs(ul0-ustar) < fabs(ul1-ustar)) {
    fmin = fabs(ul0 - ustar);
    p_fmin = p0;  rhol_fmin = rhol0;  ul_fmin = ul0;
  } else {
    fmin = fabs(ul1 - ustar);
    p_fmin = p1;  rhol_fmin = rhol1;  ul_fmin = ul1;
  }

  for(i=0; i<maxIts_bracket; i++) {

    f0 = ul0 - ustar;
    f1 = ul1 - ustar;

    if(f0*f1<=0.0)
      return true;

    // find a physical p2 that has the opposite sign
    if(fabs(f0-f1)>1e-9) {
      p2 = p1 - f1*(p1-p0)/(f1-f0); //the Secant method
      if(p2<p0)
	p2 -= 0.1*(p1-p0);
      else //p2 cannot be between p0 and p1, so p2>p1
	p2 += 0.1*(p1-p0);
    } else {//f0 == f1
      p2 = 1.1*p1; 
      //p2 = p1 + *(p1-p0); 
    }


    if(p2<pl || i==int(maxIts_bracket/2) ) {//does not look right. reset to pl
      p2 = pl; 
    }
    success = ComputeRhoUStar(1, integrationPath1, rhol, ul, pl, p2, idl, rhol0, rhol1, rhol2, ul2);

    if(!success) {

#if PRINT_RIEMANN_SOLUTION == 1
      fprintf(stdout, "  -- p2 = %e (failed)\n", p2);
#endif
      //move closer to [p0, p1]
      for(int j=0; j<maxIts_bracket/2; j++) {
	if(p2<p0)     
	  p2 = p0 - 0.5*(p0-p2);
	else //p2>p1
	  p2 = p1 + 0.5*(p2-p1);

	if(p2<pl)
	  p2 = pl;
	success = ComputeRhoUStar(1, integrationPath1, rhol, ul, pl, p2, idl, rhol0, rhol1, rhol2, ul2);

	if(success)
	  break;
      }
    }

    if(!success)
      break;

    // check & update fmin (only used if the Riemann solver fails -- a "failsafe" feature)
    if(fabs(ul2-ustar)<fmin) {
      fmin = fabs(ul2-ustar);
      p_fmin = p2;  rhol_fmin = rhol2;  ul_fmin = ul2;
    }

    // update p0 or p1
    if(p2<p0) {
      p1 = p0;  rhol1 = rhol0;  ul1 = ul0;
      p0 = p2;  rhol0 = rhol2;  ul0 = ul2;
    } else {//p2>p1
      p0 = p1;  rhol0 = rhol1;  ul0 = ul1;
      p1 = p2;  rhol1 = rhol2;  ul1 = ul2;
    }

#if PRINT_RIEMANN_SOLUTION == 1
    fprintf(stdout, "  -- p0 = %e, f0 = %e, p1 = %e, f1 = %e (success)\n", p0, ul0-ustar, p1, ul1-ustar);
#endif

  }

  if(!success || i==maxIts_bracket) {
    if(verbose>=1) {
      cout << "Warning: Exact half-Riemann solver failed. (Unable to find a bracketing interval) " << endl;
      cout << "   left: " << std::setprecision(10) << rhol << ", " << std::setprecision(10) << ul << ", " 
	<< std::setprecision(10) << pl << " (" << idl << "); ustar: "
	<< std::setprecision(10) << ustar << "."
	<< " Residual (|ulstar-ustar|): " << std::setprecision(10) << fmin << endl;
    }
    if(fmin<failure_threshold*fabs(ul-ustar)) {
      if(verbose>=1) 
	cout << "*** Best approximate solution: rhols = " << rhol_fmin << ", ps = " << p_fmin << "." << endl;
      p0    = p1    = p_fmin;
      rhol0 = rhol1 = rhol_fmin;
      ul0   = ul1   = ul_fmin; 
    } else { //it could be that the Riemann problem has no solution!
      p2 = pressure_at_failure;
      success = ComputeRhoUStar(1, integrationPath1, rhol, ul, pl, p2, idl, rhol0, rhol1, rhol2, ul2);
      if(success) {
	if(verbose >= 1)
	  cout << "*** Prescribed solution: rhols = " << rhol2 << ", ps = " << p2 << "." << endl;
	p0    = p1    = p2;
	rhol0 = rhol1 = rhol2; 
	ul0   = ul1   = ul2; 
      } else {
	if(verbose >= 1)
	  cout << "*** Best approximation: rhols = " << rhol_fmin << ", ps = " << p_fmin << "." << endl;
	p0    = p1    = p_fmin;
	rhol0 = rhol1 = rhol_fmin;
	ul0   = ul1   = ul_fmin; 
      }
    }

    return false;
  }

  return true;
}

//----------------------------------------------------------------------------------

  bool
ExactRiemannSolverBase::FindInitialFeasiblePoints(double rhol, double ul, double pl, double el, double cl, int idl,
    double rhor, double ur, double pr, double er, double cr, int idr, /*inputs*/
    double &p0, double &rhol0, double &rhor0, double &ul0, double &ur0, 
    double &p1, double &rhol1, double &rhor1, double &ul1, double &ur1/*outputs*/)
{
  double dp;
  int found = 0;
  bool success = true;

  // Method 1: Use the acoustic theory (Eqs. (20)-(22) of Kamm) to find p0, p1
  found = FindInitialFeasiblePointsByAcousticTheory(rhol, ul, pl, el, cl, idl, 
      rhor, ur, pr, er, cr, idr, /*inputs*/
      p0, rhol0, rhor0, ul0, ur0, p1, rhol1, rhor1, ul1, ur1/*outputs*/);

  if(found==2)
    return true; //yeah

  if(found==1) //the first one (p0) is good --> only need to find p1
    goto myLabel;

  // Method 2: based on dp (fixed width search)

  // 2.1. find the first one (p0)
  dp = (pl!=pr) ? fabs(pl-pr) : 0.5*pl;
  for(int i=0; i<maxIts_bracket; i++) {
    p0 = std::min(pl,pr) + 0.01*(i+1)*(i+1)*std::min(dp, std::min(fabs(pl),fabs(pr)));

    if(p0<min_pressure)
      p0 = pressure_at_failure; 
    success = ComputeRhoUStar(1, integrationPath1, rhol, ul, pl, p0, idl, 
	rhol, (p0>pl) ? rhol*1.1 : rhol*0.9,
	rhol0, ul0);
    success = success && ComputeRhoUStar(3, integrationPath3, rhor, ur, pr, p0, idr, 
	rhor, (p0>pr) ? rhor*1.1 : rhor*0.9,
	rhor0, ur0);
    if(success)
      break;
  }
  if(!success) {//search in the opposite direction
    for(int i=0; i<maxIts_bracket; i++) {
      p0 = std::min(pl,pr) - 0.01*(i+1)*(i+1)*std::min(dp, std::min(fabs(pl),fabs(pr)));
      if(p0<min_pressure || i==(int)(maxIts_bracket/2)) //not right... set to a small pos. pressure
	p0 = pressure_at_failure; 
      success = ComputeRhoUStar(1, integrationPath1, rhol, ul, pl, p0, idl, 
	  rhol, (p0>pl) ? rhol*1.1 : rhol*0.9,
	  rhol0, ul0);
      success = success && ComputeRhoUStar(3, integrationPath3, rhor, ur, pr, p0, idr, 
	  rhor, (p0>pr) ? rhor*1.1 : rhor*0.9,
	  rhor0, ur0);
      if(success)
	break;
    } 
  }

  if(!success) {
    if(verbose>=1)
      fprintf(stdout,"Warning: Failed to find the first initial guess (p0) in the 1D Riemann solver (it = %d). "
	  "Left: %e %e %e (%d), Right: %e %e %e (%d).\n",
	  maxIts_bracket, rhol, ul, pl, idl, rhor, ur, pr, idr);
    return false;
  }

myLabel:
  // 2.2. find the second one (p1)
  dp = std::min(fabs(p0-pl), fabs(p0-pr));
  for(int i=0; i<maxIts_bracket; i++) {
    p1 = p0 + 0.01*(i+1)*(i+1)*dp;
    success = ComputeRhoUStar(1, integrationPath1, rhol, ul, pl, p1, idl, rhol, rhol0, rhol1, ul1);
    success = success && ComputeRhoUStar(3, integrationPath3, rhor, ur, pr, p1, idr, rhor, rhor0, rhor1, ur1);
    if(success)
      break;
  }
  if(!success) //search in the opposite direction
    for(int i=0; i<maxIts_bracket; i++) {
      p1 = p0 - 0.01*(i+1)*(i+1)*dp;
      if(p1<min_pressure)
	p1 = pressure_at_failure*1000.0; //so it is not the same as p0
      success = ComputeRhoUStar(1, integrationPath1, rhol, ul, pl, p1, idl, rhol, rhol0, rhol1, ul1);
      success = success && ComputeRhoUStar(3, integrationPath3, rhor, ur, pr, p1, idr, rhor, rhor0, rhor1, ur1);
      if(success)
	break;
    } 
  if(!success) {
    if(verbose>=1)
      fprintf(stdout,"Warning: Failed to find the second initial guess (p1) in the 1D Riemann solver (it = %d). "
	  "Left: %e %e %e (%d), Right: %e %e %e (%d).\n",
	  maxIts_bracket, rhol, ul, pl, idl, rhor, ur, pr, idr);
    return false;
  }

  // Make sure p0<p1
  if(p0>p1) {
    std::swap(p0,p1);
    std::swap(rhol0, rhol1);
    std::swap(rhor0, rhor1);
    std::swap(ul0, ul1);
    std::swap(ur0, ur1);
  } 

  return true;
}

//----------------------------------------------------------------------------------

int
ExactRiemannSolverBase::FindInitialFeasiblePointsByAcousticTheory(double rhol, double ul, double pl,
    [[maybe_unused]] double el, double cl, int idl,
    double rhor, double ur, double pr, [[maybe_unused]] double er, double cr, int idr, /*inputs*/
    double &p0, double &rhol0, double &rhor0, double &ul0, double &ur0, 
    double &p1, double &rhol1, double &rhor1, double &ul1, double &ur1/*outputs*/)
{
  int found = 0;
  bool success = true;

  // 1.1: Initialize p0 using acoustic theory ((20) of Kamm)
  double Cl = rhol*cl; //acoustic impedance
  double Cr = rhor*cr; //acoustic impedance
  p0 = (Cr*pl + Cl*pr + Cl*Cr*(ul - ur))/(Cl + Cr);

  success = ComputeRhoUStar(1, integrationPath1, rhol, ul, pl, p0, idl/*inputs*/,
      rhol, (p0>pl) ? rhol*1.1 : rhol*0.9/*initial guesses for Hugo. eq.*/,
      rhol0, ul0/*outputs*/);
  if(!success)
    return found;

  success = ComputeRhoUStar(3, integrationPath3, rhor, ur, pr, p0, idr/*inputs*/,
      rhor, (p0>pr) ? rhor*1.1 : rhor*0.9/*initial guesses for Hugo. eq.*/,
      rhor0, ur0/*outputs*/);
  if(!success)
    return found;

  found = 1; //found p0!

  // 1.2. Initialize p1 ((21)(22) of Kamm)
  double Clbar = (ul0 == ul) ? Cl : fabs(p0 - pl)/fabs(ul0 - ul);
  double Crbar = (ur0 == ur) ? Cr : fabs(p0 - pr)/fabs(ur0 - ur);
  p1 = (Crbar*pl + Clbar*pr + Clbar*Crbar*(ul - ur))/(Clbar + Crbar);
  double tmp = std::max(fabs(p0), fabs(p1));
  if(fabs(p1 - p0)/tmp<1.0e-8)
    p1 = p0 + 1.0e-8*tmp; //to avoid f0 = f1 (divide-by-zero)

  success = ComputeRhoUStar(1, integrationPath1, rhol, ul, pl, p1, idl/*inputs*/,
      rhol, rhol0/*initial guesses for Hugo. eq.*/,
      rhol1, ul1/*outputs*/);
  if(!success)
    return found;

  success = ComputeRhoUStar(3, integrationPath3, rhor, ur, pr, p1, idr/*inputs*/,
      rhor, rhor0/*initial guesses for Hugo. eq.*/,
      rhor1, ur1/*outputs*/);
  if(!success)
    return found;

  found = 2; //found p0 and p1!
  return found;
}

//----------------------------------------------------------------------------------

  int
ExactRiemannSolverBase::FindInitialFeasiblePointsOneSidedByAcousticTheory(double rhol, double ul, 
    double pl, [[maybe_unused]] double el, double cl, int idl, double ustar,
    double &p0, double &rhol0, double &ul0, double &p1, double &rhol1, double &ul1)
{

  assert(ul>ustar); //only needed in the case of a shock

  int found = 0;
  bool success = true;

  // 1.1: Initialize p0 using acoustic theory ((20) of Kamm)
  double Cl = rhol*cl; //acoustic impedance
  p0 = pl + Cl*(ul - ustar);

  success = ComputeRhoUStar(1, integrationPath1, rhol, ul, pl, p0, idl/*inputs*/,
      rhol, (p0>pl) ? rhol*1.1 : rhol*0.9/*initial guesses for Hugo. eq.*/,
      rhol0, ul0/*outputs*/);
  if(!success) {
    p0 = 1.5*pl;
    success = ComputeRhoUStar(1, integrationPath1, rhol, ul, pl, p0, idl/*inputs*/,
	rhol, (p0>pl) ? rhol*1.1 : rhol*0.9/*initial guesses for Hugo. eq.*/,
	rhol0, ul0/*outputs*/);
    if(!success)
      return found;
  }


  found = 1; //found p0!

  double Clbar = (ul0 == ul) ? Cl : fabs(p0 - pl)/fabs(ul0 - ul);
  p1 = pl + Clbar*(ul - ustar);
  double tmp = std::max(fabs(p0), fabs(p1));
  if(fabs(p1 - p0)/tmp<1.0e-8)
    p1 = p0 + 1.0e-8*tmp; //to avoid f0 = f1 (divide-by-zero)

  success = ComputeRhoUStar(1, integrationPath1, rhol, ul, pl, p1, idl/*inputs*/,
      rhol, rhol0/*initial guesses for Hugo. eq.*/,
      rhol1, ul1/*outputs*/);
  //  fprintf(stdout,"p1 = %e, success = %d.\n", p1, (int)success);
  if(!success)
    return found;

  found = 2; //found p0 and p1!
  return found;
}

//----------------------------------------------------------------------------------

//! Connect the left/right initial state with the left/right star state (the 1-wave or 3-wave)
  bool  //true: success  | false: failure
ExactRiemannSolverBase::ComputeRhoUStar(int wavenumber /*1 or 3*/,
    std::vector<std::vector<double>>& integrationPath /*3 by n, first index: 1-pressure, 2-density, 3-velocity*/,
    double rho, double u, double p, double ps, int id/*inputs*/,
    double rhos0, double rhos1/*initial guesses for Hugo. eq.*/,
    double &rhos, double &us/*outputs*/, 
    bool *trans_rare, double *Vrare_x0/*filled only if found tran rf*/)
{
  // default
  rhos = rho;
  us   = u;

  if(p > ps) {//rarefaction --- numerical integration

    // prepare for numerical integration
    double rhos_0 = rho, us_0 = u, ps_0 = p, xi_0; //start point of each step
    double rhos_1 = rho, us_1 = u, ps_1 = p, xi_1; //end point of each step
    double dp = (p-ps)/numSteps_rarefaction;  // initial step size
    double dp_min_adaption = dp / 2.5;
    double pressure_endpoint_tol = tol_rarefaction * std::max( 1.0, fabs(p) );

    double e = vf[id]->GetInternalEnergyPerUnitMass(rho,p);
    double c = vf[id]->ComputeSoundSpeedSquare(rho, e);
    if(rho<=0 || c<0) {
      fprintf(stdout,"Warning: Negative density or c^2 (square of sound speed) in ComputeRhoUStar." 
	  "rho = %e, p = %e, e = %e, c^2 = %e, id = %d.\n",
	  rho, p, e, c, id);
      return false; //failure
    } else
      c = sqrt(c);

    int index0 = 0; // index of new starting point
    // find the new starting point, and update dp accordingly
    if (integrationPath.size() > 1) { 
      for (int j = integrationPath.size()-1; j >= 0; j--) {
	if (integrationPath[j][0] > ps) {
	  index0 = j;
	  break;
	}
      }
      ps_0 = integrationPath[index0][0];
      rhos_0 = integrationPath[index0][1];
      us_0 = integrationPath[index0][2];
      dp_min_adaption = (ps_0 - ps) / numSteps_rarefaction / 2.5;
      if (index0 != (int)integrationPath.size()-1) { // new starting point is not the last on the trajectory
	dp = ps_0-ps; 
      } else { // dp from the last step 
	dp = std::min( integrationPath[index0-1][0]-integrationPath[index0][0], ps_0-ps );
	// dp = ps_0-ps;
      }
    }

    double xi = (wavenumber == 1) ? u - c : u + c; // xi = u -/+ c
    xi_0 = xi;

#if PRINT_RIEMANN_SOLUTION == 1
    sol1d.push_back(vector<double>{xi, rho, u, p, (double)id});
#endif

    //fprintf(stdout,"rho = %e, p = %e, ps = %e\n", rho, p, ps);
    // integration by Runge-Kutta 4
    bool done = false;
    double uErr = 0.;
    double rhoErr = 0.;
    int moreSteps = 3;

    for(int i=0; i<numSteps_rarefaction*moreSteps; i++) {
      // Check if we have reached the final pressure ps
      if(fabs(ps_1 - ps) <= pressure_endpoint_tol) {
	rhos = rhos_1;
	us   = us_1;

	if(vf[id]->CheckState(rhos,ps,true)) { //true: silence
#if PRINT_RIEMANN_SOLUTION == 1
	  cout << "Rarefaction solver reached a nonphysical state!" << endl;         
#endif
	  return false;
	} 
#if PRINT_RIEMANN_SOLUTION == 1
	cout << "  " << wavenumber << "-wave: rarefaction, integration completed in " << i << " steps" << endl;
	cout << "rhos_1, us_1, ps_1: " << rhos_1 << ", " << us_1 << ", " << ps_1 << "." << endl;
#endif
	done = true;
	break; //done!
      }

      bool success = Rarefaction_OneStepRK4(wavenumber/*1 or 3*/, id,
	  rhos_0, us_0, ps_0 /*start state*/, dp /*step size*/,
	  rhos_1, us_1, ps_1, xi_1 /*output: end state*/,
	  uErr, rhoErr /*output: absolute error in us*/);

      //      fprintf(stdout,"RK4 step: rhos_0 = %e, us_0 = %e, ps_0 = %e, drho = %e, rhos_1 = %e, us_1 = %e, ps_1 = %e | ps = %e | success = %d.\n",
      //              rhos_0, us_0, ps_0, drho, rhos_1, us_1, ps_1, ps, success);
      /*			  if (i > continueTolerance || dp < 1.0e-13) {
      //          fprintf(stdout,"*** Warning: step size halved %d times, force break the loop. id = %d, p = %e, ps = %e, dp = %e, ps_1 = %e, c = %e, rho = %e, path size = %ld.\n",
      //                       continueTimes, id, p, ps, dp, ps_1, c, rho, integrationPath.size());
      break;
      }
      */

      if(!success) {
	dp = dp/2.0;
	continue;
      }

      if (ps_1-ps < -pressure_endpoint_tol) {
	dp = ps_0 - ps;
	continue;
      }

      double tiny = 1.e-30; // a small number, used to prevent uErrScaled and rhoErrScaled from being zero
      double errBar = std::max(1.0e-8, std::min(1.0e-5, 100.0*tol_rarefaction) ); // one-step error of numerical integration 
      double uErrScaled = uErr / c + tiny;
      double rhoErrScaled = rhoErr / rho + tiny;
      if (isnan(rhoErrScaled)) {
	fprintf(stdout,"Warning: rhoErrScaled is nan. id = %d, p = %e, ps = %e, dp = %e, ps_0 = %e, ps_1 = %e, uErr = %e, uErrScaled = %e, rhos_0 = %e, rhos_1 = %e, rhoErr = %e, rhoErrScaled = %e, c = %e, rho = %e, path size = %ld.\n",
	    id, p, ps, dp, ps_0, ps_1, uErr, uErrScaled, rhos_0, rhos_1, rhoErr, rhoErrScaled, c, rho, integrationPath.size());
	exit(-1);
      }  

      if (isnan(uErrScaled)) {
	fprintf(stdout,"Warning: uErrScaled is nan. id = %d, p = %e, ps = %e, dp = %e, ps_0 = %e, ps_1 = %e, uErr = %e, uErrScaled = %e, rhoErr = %e, rhoErrScaled = %e, c = %e, rho = %e, path size = %ld.\n",
	    id, p, ps, dp, ps_0, ps_1, uErr, uErrScaled, rhoErr, rhoErrScaled, c, rho, integrationPath.size());
	exit(-1);
      }  

      double dpTemp = dp; // temporary dp used for step adaption
      double safety = 0.9; // safety factor for step adaption

      if (rhoErrScaled >= errBar && dp > dp_min_adaption) { 
	dpTemp = safety * dp * pow( fabs(errBar/rhoErrScaled) , 0.25 );
	dpTemp = std::max( std::max(dpTemp, 0.2*dp), dp_min_adaption); //don't decrease dp too much
	dp = std::min(dpTemp, ps_0-ps);
	continue;
      }

      if (uErrScaled >= errBar && dp > dp_min_adaption) {
	dpTemp = safety * dp * pow( fabs(errBar/uErrScaled) , 0.25 );
	dpTemp = std::max( std::max(dpTemp, 0.2*dp), dp_min_adaption); //don't decrease dp too much
	dp = std::min(dpTemp, ps_0-ps);
	continue;
      }

      if (ps_1 < integrationPath[integrationPath.size()-1][0]) { // store the new point if necessary
	std::vector<double> currentVect = {ps_1, rhos_1, us_1};
	integrationPath.push_back(currentVect);
      }

#if PRINT_RIEMANN_SOLUTION == 1
      sol1d.push_back(vector<double>{xi_1, rhos_1, us_1, ps_1, (double)id});
#endif

      if(trans_rare && Vrare_x0 && xi_0*xi_1<=0) {//transonic rarefaction, crossing x = xi = 0
	*trans_rare = true;
	double w0 = fabs(xi_1), w1 = fabs(xi_0);
	double ww = w0 + w1;
	w0 /= ww;
	w1 /= ww;
	Vrare_x0[0] = w0*rhos_0 + w1*rhos_1;
	Vrare_x0[1] = w0*us_0   + w1*us_1;
	Vrare_x0[2] = w0*ps_0   + w1*ps_1;

#if PRINT_RIEMANN_SOLUTION == 1
	sol1d.push_back(vector<double>{0.0, Vrare_x0[0], Vrare_x0[1], Vrare_x0[2], (double)id});
#endif

      }

      //fprintf(stdout,"drho = %e, rho: %e -> %e,  u: %e -> %e,  p: %e -> %e\n", drho, rhos_0, rhos_1, us_0, us_1, ps_0, ps_1);

      // If the solver gets here, it means it hasn't reached the final pressure ps.
      // Adjust step size, then update state
      //
      /*
	 if(i == moreSteps*numSteps_rarefaction-1) {
	 fprintf(stdout,"Warning: integrator used up all the steps specified. id = %d, p = %e, ps = %e, dp = %e, ps_1 = %e, uErrScaled = %e, rhoErrScaled = %e, c = %e, rho = %e, rhos_1 = %e, path size = %ld.\n",
	 id, p, ps, dp, ps_1, uErrScaled, rhoErrScaled, c, rho, rhos_1, integrationPath.size());
	 }
	 */        
      // increase dp is allowed 
      double dpTemp_rho = safety * dp * pow( fabs(errBar/rhoErrScaled) , 0.2 );
      double dpTemp_u = safety * dp * pow( fabs(errBar/uErrScaled) , 0.2 );
      dpTemp = std::min(dpTemp_rho, dpTemp_u);
      dpTemp = std::min(dpTemp, 10*dp); //don't increase dp too much
      dpTemp = std::max(dpTemp, dp_min_adaption);
      dp = std::min(dpTemp, ps_1-ps); //don't go beyond ps

      rhos_0 = rhos_1;
      us_0   = us_1;
      ps_0   = ps_1;
      xi_0   = xi_1;
    }

    if(!done) {
      if(vf[id]->CheckState(rhos_1,ps_1,true)) {
#if PRINT_RIEMANN_SOLUTION == 1
	cout << "  " << wavenumber << "-wave: rarefaction, solver failed (unphysical state: rhos = "
	  << rhos_1 << ", ps = " << ps_1 << "!)" << endl;
#endif
	return false; //failed
      } else {
#if PRINT_RIEMANN_SOLUTION == 1
	cout << "  " << wavenumber << "-wave: rarefaction, solver did not converge (final sol.: rhos_1 = "
	  << rhos_1 << ", ps_1 = " << ps_1 << "; inputs: rho = " << rho << ", p = " << p << ", ps = " << ps << ")" << endl;
#endif
	return false; 
      }
    }
  }

  else {// shock (p<=ps, rho<=rhos)

    HugoniotEquation hugo(vf[id],rho,p,ps);

    //find a bracketing interval
    double f0, f1;
    double drho = std::max(fabs(rhos0 - rhos1), 0.001*rhos0); 
    bool found_rhos0 = false, found_rhos1 = false;
    if(std::min(rhos0,rhos1)>=rho) {//both rhos0 and rhos1 are physically admissible
      f0 = hugo(rhos0);
      f1 = hugo(rhos1);
      if(f0*f1<=0) {
	/*found bracketing interval [rhos0, rhos1]*/
	if(rhos0>rhos1) {
	  std::swap(rhos0,rhos1);
	  std::swap(f0,f1);
	}
	found_rhos0 = found_rhos1 = true;
      } else {
	rhos0 = rhos1; //this is our starting point (presumably rhos1 is closer to sol'n)
	f0    = f1; 
      }
    } else { //at least, the smaller one among rhos0, rhos1 is non-physical
      if(rhos1>rhos0) {
	rhos0 = rhos1; 
	f0    = hugo(rhos0);
      }
      if(rhos0<rho) {//this one is also non-physical
	rhos0 = rho;
	f0    = hugo(rhos0);
	found_rhos0 = true;
      } else {
	/*rhos0 = rhos0;*/
	f0 = hugo(rhos0);
      }
    } 

    if(!found_rhos0 || !found_rhos1) {
      int i = 0;
      double factor = 1.5; 
      double tmp, ftmp;
      // before the search, rhos0 = rhos1 = an adimissible point > rho
      rhos1 = rhos0;
      f1    = f0;
      while(!found_rhos0) {
	if(++i>=maxIts_shock) {
	  //          cout << "*** Error: Unable to find a bracketing interval after " << maxIts_shock 
	  //               << " iterations (in the solution of the Hugoniot equation)." << endl;
	  return false;
	}
	tmp = rhos1;
	ftmp = f1;
	//move to the left (towards rho)
	rhos1 = rhos0;
	f1    = f0;
	rhos0 = rhos1 - factor*drho;
	if(rhos0<=rho) {
	  rhos0 = rho;
	  found_rhos0 = true;
	}
	f0 = hugo(rhos0);

	if(f0*f1<=0) {
	  found_rhos0 = found_rhos1 = true;
	} else {
	  //move to the right
	  rhos1 = tmp;
	  f1    = ftmp;
	  tmp   = rhos0; //don't forget the smallest point
	  ftmp  = f0;
	  rhos0 = rhos1;
	  f0    = f1;
	  rhos1 = rhos0 + factor*drho;
	  f1 = hugo(rhos1);
	  if(f0*f1<=0) {
	    found_rhos0 = found_rhos1 = true;
	  } else {
	    rhos0 = tmp;
	    f0    = ftmp;
	    drho  = rhos1 - rhos0; //update drho
	  }
	}
      }

      if(!found_rhos1) {//keep moving to the right
	i = 0;
	double factor = 2.5;
	while(!found_rhos1) {
	  if(++i>=maxIts_shock) {
	    //            cout << "*** Error: Unable to find a bracketing interval after " << maxIts_shock 
	    //                 << " iterations (in the solution of the Hugoniot equation (2))." << endl;
	    return false; //failure
	  }
	  rhos0 = rhos1;
	  f0    = f1;
	  rhos1 = rhos0 + factor*drho;
	  f1    = hugo(rhos1);
	  if(f0*f1<=0) {
	    found_rhos1 = true;
	  } else
	    drho = rhos1 - rhos0;
	}
      }

    }

    pair<double,double> sol;
    double loc_tol_shock = tol_shock*std::min(rhos0,rhos1);
#ifndef WITHOUT_BOOST
    //*******************************************************************
    // Calling boost function for root-finding
    // Warning: "maxit" is BOTH AN INPUT AND AN OUTPUT
    boost::uintmax_t maxit = maxIts_shock;
    if(f0==0.0)
      sol.first = sol.second = rhos0;
    else if(f1==0.0)
      sol.first = sol.second = rhos1;
    else {
      sol = toms748_solve(hugo, rhos0, rhos1, f0, f1,
	  [=](double r0, double r1){return r1-r0<std::min(loc_tol_shock,0.001*(rhos1-rhos0));}, 
	  maxit);
    }
    //*******************************************************************
#else
    //*******************************************************************
    // Using a hybrid (Brent) method for root-finding
    int maxit = 0;
    if(f0==0.0)
      sol.first = sol.second = rhos0;
    else if(f1==0.0)
      sol.first = sol.second = rhos1;
    else {
      double rhos2 = rhos1; //rhos2 is always the latest one
      double f2    = f1;
      int it;
      for(it = 0; it<maxIts_shock; it++) {
	drho = rhos1 - rhos0;
	rhos2 = rhos2 - f2*(rhos1 - rhos0)/(f1 - f0); //secant method
	if(rhos2 >= rhos1 || rhos2 <= rhos0) //discard and switch to bisection
	  rhos2 = 0.5*(rhos0+rhos1);
	f2 = hugo(rhos2);
	if(f2==0.0) {
	  sol.first = sol.second = rhos2;
	  break;
	}
	if(f2*f0<0) {
	  rhos1 = rhos2;
	  f1    = f2;
	} else {
	  rhos0 = rhos2;
	  f0    = f2;
	}
	if(rhos1-rhos0<loc_tol_shock) {
	  sol.first  = rhos0;
	  sol.second = rhos1;
	  break;
	}
      }
      if(it==maxIts_shock) {
	fprintf(stdout,"*** Error: Root-finding method failed to converge after %d iterations.\n", it);
	return false;
      }
      maxit = it;
    }
    //*******************************************************************
#endif


#if PRINT_RIEMANN_SOLUTION == 1
    cout << "  " << wavenumber << "-wave: shock, converged in " << maxit << " iterations. fun = " 
      << hugo(0.5*(sol.first+sol.second)) << "." << endl;
#endif

    rhos = 0.5*(sol.first+sol.second);

    double du = -(ps-p)*(1.0/rhos-1.0/rho);
    if(du<0) {
      //cout << "Warning: Violation of hyperbolicitiy when enforcing the Rankine-Hugoniot jump conditions (du = "
      //     << du << ")." << endl;
      return false;
    }

    if(vf[id]->CheckState(rhos,ps,true)) //true: silence
      return false; //nonphysical...

    us = (wavenumber==1) ? u - sqrt(du) : u + sqrt(du);


#if PRINT_RIEMANN_SOLUTION == 1
    double xi = (rhos*us - rho*u)/(rhos-rho);
    if(wavenumber==1) {
      sol1d.push_back(vector<double>{xi-0.0001*fabs(xi), rho, u, p, (double)id});
      sol1d.push_back(vector<double>{xi, rhos, us, ps, (double)id});
    } else {
      sol1d.push_back(vector<double>{xi, rhos, us, ps, (double)id});
      sol1d.push_back(vector<double>{xi+0.0001*fabs(xi), rho, u, p, (double)id});
    }
#endif

  }

  return true; //yeah!
}

//----------------------------------------------------------------------------------
//! Connect the left initial state with the left star state (the 1-wave) --- for one-sided Riemann problem
//! where the solution contains a rarefaction (not a shock).
bool  //true: success  | false: failure
ExactRiemannSolverBase::ComputeOneSidedRarefaction(double rho, double u, double p, [[maybe_unused]] double e,
    double c, int id, double us/*inputs*/,
    double &rhos, double &ps/*outputs*/, 
    bool *trans_rare, double *Vrare_x0/*filled only if found tran rf*/)
{

  assert(u<us);

  int wavenumber = 1;

  // default
  rhos = rho;
  ps   = p;

  //in this function, du, drho are all positive!

  double du_target = (us-u)/numSteps_rarefaction;
  double du_max = du_target*1.25;

  // initialize drho based on linear approximation & fixing c
  double dp_op1 = fabs(du_target*c*rho);
  double dp_op2 = (fabs(p)+1.0e-6)/(numSteps_rarefaction*2.5); //initial step size
  double dp = std::min(dp_op1, dp_op2);

  // prepare for numerical integration
  double rhos_0 = rho, us_0 = u, ps_0 = p, xi_0; //start point of each step
  double rhos_1 = rho, us_1 = u, ps_1 = p, xi_1; //end point of each step
  double du;

  double xi = u - c; 
  xi_0 = xi;

#if PRINT_RIEMANN_SOLUTION == 1
  sol1d.push_back(vector<double>{xi, rho, u, p, (double)id});
#endif

  // integration by Runge-Kutta 4
  // as we integrate, rho and p decreases, while u increases until reaching ustar.
  bool done = false;
  double loc_tol_rarefaction = tol_rarefaction*std::max(us-u, std::max(fabs(us), fabs(u)));
  double uErr(0), rhoErr(0);
  for(int i=0; i<numSteps_rarefaction*5; i++) {

    bool success = Rarefaction_OneStepRK4(wavenumber/*1 or 3*/, id,
	rhos_0, us_0, ps_0 /*start state*/, dp /*step size*/,
	rhos_1, us_1, ps_1, xi_1 /*output: end state*/,
	uErr, rhoErr); 

    if(!success) {
      dp /= 2.0;
      continue;
    }

    du = us_1 - us_0; //du and drho are positive in this function

    if(du > du_max) {
      dp = dp/du*du_target;
      continue;
    }
    if(us_1 - us > loc_tol_rarefaction) { //went over the hill...
      if(du!=0)
	dp = dp/du*(us - us_0);
      else
	dp /= 2.0;
      continue;
    } 

#if PRINT_RIEMANN_SOLUTION == 1
    sol1d.push_back(vector<double>{xi_1, rhos_1, us_1, ps_1, (double)id});
#endif

    if(trans_rare && Vrare_x0 && xi_0*xi_1<=0) {//transonic rarefaction, crossing x = xi = 0
      *trans_rare = true;
      double w0 = fabs(xi_1), w1 = fabs(xi_0);
      double ww = w0 + w1;
      w0 /= ww;
      w1 /= ww;
      Vrare_x0[0] = w0*rhos_0 + w1*rhos_1;
      Vrare_x0[1] = w0*us_0   + w1*us_1;
      Vrare_x0[2] = w0*ps_0   + w1*ps_1;

#if PRINT_RIEMANN_SOLUTION == 1
      sol1d.push_back(vector<double>{0.0, Vrare_x0[0], Vrare_x0[1], Vrare_x0[2], (double)id});
#endif

    }


    // Check if we have reached the final pressure ps
    if(fabs(us - us_1) <= loc_tol_rarefaction) {
      rhos = rhos_1;
      ps   = ps_1;

      if(vf[id]->CheckState(rhos,ps,true)) { //true: silence
#if PRINT_RIEMANN_SOLUTION == 1
	cout << "Rarefaction solver reached a nonphysical state!" << endl;         
#endif
	return false;
      }

#if PRINT_RIEMANN_SOLUTION == 1
      cout << "  " << wavenumber << "-wave: rarefaction, integration completed in " << i << " steps" << endl;
#endif
      done = true;
      break; //done!
    }

    // If the solver gets here, it means it hasn't reached the final velocity us
    // Adjust step size, then update state

    dp = std::min( dp/du*std::min(du_target,us-us_1), //don't go beyond us
	dp*5.0); //don't increase too much in one step

    rhos_0 = rhos_1;
    us_0   = us_1;
    ps_0   = ps_1;
    xi_0   = xi_1;

  }

  if(!done) {
    if(vf[id]->CheckState(rhos_1,ps_1,true)) {
#if PRINT_RIEMANN_SOLUTION == 1
      cout << "  " << wavenumber << "-wave: rarefaction, solver failed (unphysical state: rhos = "
	<< rhos_1 << ", ps = " << ps_1 << "!)" << endl;
#endif
      return false; //failed
    } else {
#if PRINT_RIEMANN_SOLUTION == 1
      cout << "  " << wavenumber << "-wave: rarefaction, solver did not converge (final sol.: rhos_1 = "
	<< rhos_1 << ", ps_1 = " << ps_1 << "; inputs: rho = " << rho << ", p = " << p << ", us = " << us << ")" << endl;
#endif
      return false; 
    }
  }

  return true; //yeah!
}

//----------------------------------------------------------------------------------
// Adaptive Runge-Kutta (aka. Runge-Kutta-Fehlberg). Ref. Section 25.5.2 of Chapra and Canale,
// Numerical Methods for Engineers (7th Edition)
  bool
ExactRiemannSolverBase::Rarefaction_OneStepRK4(int wavenumber/*1 or 3*/, int id,
    double rho_0, double u_0, double p_0 /*start state*/, 
    double dp /*step*/,
    double &rho, double &u, double &p, double &xi /*output*/,
    double & uErr, double & rhoErr /*output*/)
{
  dp = -dp; // dp is positive when passed in. It is actually negative if we follow Kamm's paper 
  // Equations (36 - 42)

  double e_0 = vf[id]->GetInternalEnergyPerUnitMass(rho_0, p_0);
  double c_0_square = vf[id]->ComputeSoundSpeedSquare(rho_0, e_0);

  if(rho_0<=0 || c_0_square<0) {
    //    fprintf(stdout,"*** Error: Negative density or c^2 (square of sound speed, %e) in Rarefaction_OneStepRK4(0)." 
    //            " rho = %e, p = %e, e = %e, id = %d.\n",
    //            c_0_square, rho_0, p_0, e_0, id);
    return false;
  } 

  double c_0 = sqrt(c_0_square);

  double rho_1 = rho_0 + 0.2*dp/c_0_square;
  double p_1 = p_0 + 0.2*dp;
  double e_1 = vf[id]->GetInternalEnergyPerUnitMass(rho_1, p_1);
  double c_1_square = vf[id]->ComputeSoundSpeedSquare(rho_1, e_1);

  if(rho_1<=0 || c_1_square<0) {
    //    fprintf(stdout,"*** Error: Negative density or c^2 (square of sound speed, %e) in Rarefaction_OneStepRK4(1)." 
    //                   " rho = %e, p = %e, e = %e, id = %d.\n",
    //                   c_1_square, rho_1, p_1, e_1, id);
    return false;
  } 

  //double c_1 = sqrt(c_1_square); //unused.
                             
  double rho_2 = rho_0 + 3./40.*dp/c_0_square + 9./40.*dp/c_1_square;
  double p_2 = p_1 + 3./10.*dp;
  double e_2 = vf[id]->GetInternalEnergyPerUnitMass(rho_2, p_2);
  double c_2_square = vf[id]->ComputeSoundSpeedSquare(rho_2, e_2);

  if(rho_2<=0 || c_2_square<0) {
    //    fprintf(stdout,"*** Error: Negative density or c^2 (square of sound speed, %e) in Rarefaction_OneStepRK4(2)." 
    //               " rho = %e, p = %e, e = %e, id = %d.\n",
    //                c_2_square, rho_2, p_2, e_2, id);
    return false;
  }

  double c_2 = sqrt(c_2_square); 

  double rho_3 = rho_0 + 3./10.*dp/c_0_square - 9./10.*dp/c_1_square + 6./5.*dp/c_2_square;
  double p_3 = p_0 + 3./5.*dp;
  double e_3 = vf[id]->GetInternalEnergyPerUnitMass(rho_3, p_3);
  double c_3_square = vf[id]->ComputeSoundSpeedSquare(rho_3, e_3);

  if(rho_3<=0 || c_3_square<0) {
    //    fprintf(stdout,"*** Error: Negative density or c^2 (square of sound speed, %e) in Rarefaction_OneStepRK4(3)." 
    //                " rho = %e, p = %e, e = %e, id = %d.\n",
    //                c_3_square, rho_3, p_3, e_3, id);
    return false;
  }  

  double c_3 = sqrt(c_3_square);

  double rho_4 = rho_0 - 11./54.*dp/c_0_square + 2.5*dp/c_1_square - 70./27.*dp/c_2_square + 35./27.*dp/c_3_square;
  double p_4 = p_0 + dp;
  double e_4 = vf[id]->GetInternalEnergyPerUnitMass(rho_4, p_4);
  double c_4_square = vf[id]->ComputeSoundSpeedSquare(rho_4, e_4);

  if (rho_4<=0 || c_4_square<0) {
    return false;
  }

  double c_4 = sqrt(c_4_square);

  double rho_5 = rho_0 + 1631./55296.*dp/c_0_square + 175./512.*dp/c_1_square + 575./13824.*dp/c_2_square + 44275./110592.*dp/c_3_square + 253./4096.*dp/c_4_square;
  double p_5 = p_0 + 7./8.*dp;
  double e_5 = vf[id]->GetInternalEnergyPerUnitMass(rho_5, p_5);
  double c_5_square = vf[id]->ComputeSoundSpeedSquare(rho_5, e_5);

  if (rho_5<=0 || c_5_square<0) {
    return false;
  }

  double c_5 = sqrt(c_5_square);

  // calculate the outputs
  //
  double drho = (37./378./c_0_square + 250./621./c_2_square + 125./594./c_3_square + 512./1771./c_5_square) * dp;
  double drho_err = (2825./27648./c_0_square + 18575./48384./c_2_square + 13525./55296./c_3_square + 277./14336./c_4_square + 0.25/c_5_square) * dp;
  double du = (37./378./c_0/rho_0 + 250./621./c_2/rho_2 + 125./594./c_3/rho_3 + 512./1771./c_5/rho_5) * dp;
  double du_err = (2825./27648./c_0/rho_0 + 18575./48384./c_2/rho_2 + 13525./55296./c_3/rho_3 + 277./14336./c_4/rho_4 + 0.25/c_5/rho_5) * dp;

  if (isnan(du) == 1 || isnan(du_err) ==1) {
    // fprintf(stdout, "*** Error: du or du_err is nan: du = %e, du_err = %e.\n", du, du_err);
    return false;
  }

  rhoErr = fabs(drho_err - drho);
  uErr = fabs(du_err - du);

  //  std::cout << "absolute uErr = " << uErr << std::endl;

  rho = rho_0 + drho;
  p = p_0 + dp;
  u = (wavenumber == 1) ? u_0 - du : u_0 + du; 

  double e = vf[id]->GetInternalEnergyPerUnitMass(rho, p);
  double c = vf[id]->ComputeSoundSpeedSquare(rho, e);

  //std::cout << std::setw(16) << p << std::setw(16) << rho << std::endl;

  if(rho<=0 || c<0) {
    //    fprintf(stdout,"*** Error: Negative density or c^2 (square of sound speed, %e) in Rarefaction_OneStepRK4(final)." 
    //               " rho = %e, p = %e, e = %e, id = %d.\n",
    //               c, rho, p, e, id);
    return false;
  } else

    c = sqrt(c);

  xi = (wavenumber == 1) ? u - c : u + c;

  return true;
}

//-----------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------

  void
ExactRiemannSolverBase::PrintStarRelations(double rhol, double ul, double pl, int idl,
    double rhor, double ur, double pr, int idr,
    double pmin, double pmax, double dp)
{

  vector<std::array<double,3> > left; //(p*, rhol*, ul*)
  vector<std::array<double,3> > right; //(p*, rhor*, ur*)

  double ps, rhols, rhors, uls,  urs;
  bool success;
  ps = pmin;

  while(true) {

    success = ComputeRhoUStar(1, integrationPath1, rhol, ul, pl, ps, idl/*inputs*/,
	rhol, (ps>pl) ? rhol*1.1 : rhol*0.9/*initial guesses for Hugo. eq.*/,
	rhols, uls/*outputs*/);
    if(success)
      left.push_back(std::array<double,3>{{ps,rhols,uls}});
    else 
      fprintf(stdout," -- ComputeRhoUStar(1) failed. left state: %e %e %e (%d), ps = %e.\n",
	  rhol, ul, pl, idl, ps);

    success = ComputeRhoUStar(3, integrationPath3, rhor, ur, pr, ps, idr/*inputs*/,
	rhor, (ps>pr) ? rhor*1.1 : rhor*0.9/*initial guesses for Hugo. eq.*/,
	rhors, urs/*outputs*/);
    if(success)
      right.push_back(std::array<double,3>{{ps,rhors,urs}});
    else
      fprintf(stdout," -- ComputeRhoUStar(3) failed. right state: %e %e %e (%d), ps = %e.\n",
	  rhor, ur, pr, idr, ps);

    if(ps>=pmax)
      break;

    ps = std::min(ps+dp, pmax);

  }


  FILE* file = fopen("LeftStarState.txt", "w");
  print(file, "## One-Dimensional Riemann Problem.\n");
  print(file, "## Initial State: %e %e %e, id %d (left) | (right) %e %e %e, id %d.\n", 
      rhol, ul, pl, idl, rhor, ur, pr, idr);
  print(file, "## pmin = %e, pmax = %e, dp = %e.\n", pmin, pmax, dp);

  for(auto it = left.begin(); it != left.end(); it++) 
    print(file, "%e    %e    %e\n", (*it)[0], (*it)[1], (*it)[2]);
  fclose(file);


  file = fopen("RightStarState.txt", "w");
  print(file, "## One-Dimensional Riemann Problem.\n");
  print(file, "## Initial State: %e %e %e, id %d (left) | (right) %e %e %e, id %d.\n", 
      rhol, ul, pl, idl, rhor, ur, pr, idr);
  print(file, "## pmin = %e, pmax = %e, dp = %e.\n", pmin, pmax, dp);

  for(auto it = right.begin(); it != right.end(); it++) 
    print(file, "%e    %e    %e\n", (*it)[0], (*it)[1], (*it)[2]);
  fclose(file);


}

//-----------------------------------------------------
/** Solves the one-sided Riemann problem. The valid material
 * is assumed to be on the "left", that is, "dir" points towards
 * the interface or wall
 * Return values:
 * 0: no errors
 * 1: riemann solver failed to find a bracketing interval
 */
  int
ExactRiemannSolverBase::ComputeOneSidedRiemannSolution(double *dir/*unit normal towards interface/wall*/,
    double *Vm, int idl /*left state*/,
    double *Ustar, /*interface/wall velocity (3D)*/
    double *Vs, int &id, /*solution at xi = 0 (i.e. x=0), id = -1 if invalid*/
    double *Vsm /*left 'star' solution*/)
{

  // Convert to a 1D problem (i.e. One-Dimensional Riemann)
  double rhol  = Vm[0];
  double ul    = Vm[1]*dir[0] + Vm[2]*dir[1] + Vm[3]*dir[2];
  double pl    = Vm[4];
  double ustar = Ustar[0]*dir[0] + Ustar[1]*dir[1] + Ustar[2]*dir[2];

  double el = vf[idl]->GetInternalEnergyPerUnitMass(rhol, pl);
  double cl = vf[idl]->ComputeSoundSpeedSquare(rhol, el);

  if(rhol<=0 || cl<0) {
    fprintf(stdout,"*** Error: Negative density or c^2 (square of sound speed) in ComputeOneSidedRiemannSolution." 
	" rho = %e, u = %e, p = %e, e = %e, c^2 = %e, ID = %d.\n",
	rhol, ul, pl, el, cl, idl);
    exit(-1);
  } else
    cl = sqrt(cl);

  integrationPath1.clear();
  std::vector<double> vectL{pl, rhol, ul};
  integrationPath1.push_back(vectL);


  // Declare variables in the "star region"
  double p0(DBL_MIN), ul0(0.0), rhol0(DBL_MIN);
  double p1(DBL_MIN), ul1(0.0), rhol1(DBL_MIN); //Secant Method ("k-1","k" in Kamm, (19))
  double p2(DBL_MIN), ul2(0.0), rhol2(DBL_MIN); // "k+1"
  double f0(0.0), f1(0.0), f2(0.0); //difference between ul and ur


  // -------------------------------
  // Now, Solve The One-Sided Riemann Problem
  // -------------------------------
  bool success = true;

  // monitor if the solution involves a transonic rarefaction. This is special as the solution 
  // at xi = x = 0 is within the rarefaction fan.
  bool trans_rare = false;
  double Vrare_x0[3]; //rho, u, and p at x = 0, in the case of a transonic rarefaction


  // A Trivial Case
  if(fabs(ul-ustar)<1.0e-20) {
    FinalizeOneSidedSolution(dir, Vm, rhol, ul, pl, idl, rhol, ul, pl, trans_rare, Vrare_x0, //inputs
	Vs, id, Vsm/*outputs*/);
    return 0;
  }



  if(ul < ustar) { //rarefaction

    success = ComputeOneSidedRarefaction(rhol, ul, pl, el, cl, idl, ustar, rhol2, p2, &trans_rare, Vrare_x0); 

    if(!success) {
      cout << "Warning: One-sided Riemann solver failed to complete (Returning a modified initial state)." << endl;

      double utanl[3] = {Vm[1]-ul*dir[0], Vm[2]-ul*dir[1], Vm[3]-ul*dir[2]};

      Vsm[0] = Vm[0];
      Vsm[1] = utanl[0] + ustar*dir[0];
      Vsm[2] = utanl[1] + ustar*dir[1];
      Vsm[3] = utanl[2] + ustar*dir[2];
      Vsm[4] = Vm[4];

      if(ustar >= 0) {
	id = idl;
	for(int i=0; i<5; i++)
	  Vs[i] = Vsm[i];
      } else {
	id = INVALID_MATERIAL_ID;
      }
      return 1;
    }

    // success!
    //
    FinalizeOneSidedSolution(dir, Vm, rhol, ul, pl, idl, rhol2, ustar, p2, trans_rare, Vrare_x0, /*inputs*/
	Vs, id, Vsm);
    return 0;

  }


  // -------------------------------
  // Now, ul > ustar --> SHOCK
  // -------------------------------
  // Step 1: Initialization
  //         (find initial interval [p0, p1])
  // -------------------------------
  success = FindInitialIntervalOneSided(rhol, ul, pl, el, cl, idl, ustar, /*inputs*/
      p0, rhol0, ul0, p1, rhol1, ul1/*outputs*/);
  /* our convention is that p0 < p1 */

  if(!success) { //failed to find a bracketing interval. Output the state corresponding smallest "f"

    // get sol1d, trans_rare and Vrare_x0
#if PRINT_RIEMANN_SOLUTION == 1
    sol1d.clear();
#endif
    success = ComputeRhoUStar(1, integrationPath1, rhol, ul, pl, p1, idl, rhol0, rhol0*1.1, rhol2, ul2,
	&trans_rare, Vrare_x0/*filled only if found a trans. rarefaction*/);

    if(!success) {
      cout << "Warning: One-sided Riemann solver failed to find an initial bracketing interval (Returning a modified initial state)." << endl;

      double utanl[3] = {Vm[1]-ul*dir[0], Vm[2]-ul*dir[1], Vm[3]-ul*dir[2]};

      Vsm[0] = Vm[0];
      Vsm[1] = utanl[0] + ustar*dir[0];
      Vsm[2] = utanl[1] + ustar*dir[1];
      Vsm[3] = utanl[2] + ustar*dir[2];
      Vsm[4] = Vm[4];

      if(ustar >= 0) {
	id = idl;
	for(int i=0; i<5; i++)
	  Vs[i] = Vsm[i];
      } else {
	id = INVALID_MATERIAL_ID;
      }
      return 1;
    }

    FinalizeOneSidedSolution(dir, Vm, rhol, ul, pl, idl, rhol2, ustar, p1, trans_rare, Vrare_x0, /*inputs*/
	Vs, id, Vsm);
    return 1;
  }

  assert(p0>=pl);
  assert(p1>=pl);

  f0 = ul0 - ustar;
  f1 = ul1 - ustar;

#if PRINT_RIEMANN_SOLUTION == 1
  cout << "Found initial interval: p0 = " << p0 << ", f0 = " << f0 << ", p1 = " << p1 << ", f1 = " << f1 << endl;
#endif

  // -------------------------------
  // Step 2: Main Loop (Secant Method, Safeguarded) 
  // -------------------------------
  double denom = 0; 

  int iter = 0;
  double err_p = 1.0, err_u = 1.0;

  //p2 (and f2) is always the latest one
  p2 = p1; 
  f2 = f1;

#if PRINT_RIEMANN_SOLUTION == 1
  sol1d.clear();
#endif

  for(iter=0; iter<maxIts_main; iter++) {

    // 2.1: Update p using the Brent method (safeguarded secant method)
    denom = f1 - f0;
    if(denom == 0) {
      cout << "Warning: Division-by-zero while using the secant method to solve the one-sided Riemann problem." << endl;
      cout << "         left state: " << rhol << ", " << ul << ", " << pl << ", " << idl << " | wall velocity (normal): " 
	<< ustar << endl;
      cout << "         dir = " << dir[0] << "," << dir[1] << "," << dir[2]
	<< ", f0 = " << f0 << ", f1 = " << f1 << endl;
      p2 = 0.5*(p0+p1);
    } else {
      p2 = p2 - f2*(p1-p0)/denom;  // update p2
      if(p2<=p0 || p2>=p1) //discard and switch to bisection
	p2 = 0.5*(p0+p1);
    }

    //fprintf(stdout,"iter = %d, p0 = %e, p1 = %e, p2 = %e, f0 = %e, f1 = %e.\n", iter, p0, p1, p2, f0, f1);


try_again:

    // 2.2: Calculate ul2 
    success = ComputeRhoUStar(1, integrationPath1, rhol, ul, pl, p2, idl/*inputs*/, 
	rhol0, rhol1/*initial guesses for Hugo. eq.*/,
	rhol2, ul2/*outputs*/, 
	&trans_rare, Vrare_x0/*filled only if found a trans. rarefaction*/);

    if(!success) {
      //fprintf(stdout,"*** Error: Exact Riemann solver failed. left: %e %e %e (%d) | right: %e %e %e (%d).\n", 
      //        rhol, ul, pl, idl, rhor, ur, pr, idr);
      p2 = 0.5*(p2+p0); //move closer to p0

      iter++;
      if(iter<maxIts_main)
	goto try_again;
      else
	break;
    }


    f2 = ul2 - ustar;


    // 2.3: Update for the next iteration
    if(f0*f2<0.0) {
      p1 = p2;
      f1 = f2;
      rhol1 = rhol2;
    } else {
      p0 = p2;
      f0 = f2;
      rhol0 = rhol2;
    }


    // 2.4: Check stopping criterion
    err_p = fabs(p1 - p0)/fabs(pl + 0.5*rhol*ul*ul);
    err_u = fabs(f2)/cl;

#if PRINT_RIEMANN_SOLUTION == 1
    cout << "Iter " << iter << ": p-interval = [" << p0 << ", " << p1 << "], p2 = " << p2 << ", err_p = " << err_p 
      << ", err_u = " << err_u << "." << endl;
#endif

    if( (err_p < tol_main && err_u < tol_main) || (err_p < tol_main*1e-3) || (err_u < tol_main*1e-3) )
      break; // converged

    trans_rare = false; //reset

#if PRINT_RIEMANN_SOLUTION == 1
    sol1d.clear();
#endif

  }


  // -------------------------------
  // Step 3: Find state at xi = x = 0 (for output)
  // -------------------------------
  FinalizeOneSidedSolution(dir, Vm, rhol, ul, pl, idl, rhol2, ustar, p2, trans_rare, Vrare_x0, /*inputs*/
      Vs, id, Vsm);

  if(iter == maxIts_main) {
    if(verbose>=1) {
      cout << "Warning: Exact Half-Riemann solver failed to converge. err_p = " << err_p
	<< ", err_u = " << err_u << "." << endl;
      cout << "    Vm = [" << Vm[0] << ", " << Vm[1] << ", " << Vm[2] << ", " << Vm[3] << ", " << Vm[4] << "] (" << idl
	<< "),  ustar = " << ustar << endl;
    }
    return 1;  //failed
  }


  //success!
  return 0;

}

//----------------------------------------------------------------------------------
/** Solves the one-dimensional Riemann problem. Extension of Kamm 2015 
 * to Two Materials. See KW's notes for details
 * Returns an integer error code
 * 0: no errors
 * 1: riemann solver failed to find a bracketing interval
 */
int
ExactRiemannSolverNonAdaptive::ComputeRiemannSolution(double *dir, 
    double *Vm, int idl /*"left" state*/, 
    double *Vp, int idr /*"right" state*/, 
    double *Vs, int &id /*solution at xi = 0 (i.e. x=0) */,
    double *Vsm /*left 'star' solution*/,
    double *Vsp /*right 'star' solution*/,
    double curvature)
{
  assert(curvature == 0.0); //the base class does not handle curvature!

  // Convert to a 1D problem (i.e. One-Dimensional Riemann)
  double rhol  = Vm[0];
  double ul    = Vm[1]*dir[0] + Vm[2]*dir[1] + Vm[3]*dir[2];
  double pl    = Vm[4];
  double rhor  = Vp[0];
  double ur    = Vp[1]*dir[0] + Vp[2]*dir[1] + Vp[3]*dir[2];
  double pr    = Vp[4];
  //fprintf(stdout,"1DRiemann: left = %e %e %e (%d) : right = %e %e %e (%d)\n", rhol, ul, pl, idl, rhor, ur, pr, idr);

  integrationPath1.clear();
  integrationPath3.clear();
  std::vector<double> vectL{pl, rhol, ul};
  std::vector<double> vectR{pr, rhor, ur};
  integrationPath1.push_back(vectL);
  integrationPath3.push_back(vectR); 

#if PRINT_RIEMANN_SOLUTION == 1
  std::cout << "Left State (rho, u, p): " << rhol << ", " << ul << ", " << pl << "." << std::endl;
  std::cout << "Right State (rho, u, p): " << rhor << ", " << ur << ", " << pr << "." << std::endl;
#endif

  double el = vf[idl]->GetInternalEnergyPerUnitMass(rhol, pl);
  double cl = vf[idl]->ComputeSoundSpeedSquare(rhol, el);

  if(rhol<=0 || cl<0) {
    fprintf(stdout,"*** Error: Negative density or c^2 (square of sound speed) in ComputeRiemannSolution(l)." 
	" rho = %e, u = %e, p = %e, e = %e, c^2 = %e, ID = %d.\n",
	rhol, ul, pl, el, cl, idl);
    exit(-1);
  } else
    cl = sqrt(cl);

  double er = vf[idr]->GetInternalEnergyPerUnitMass(rhor, pr);
  double cr = vf[idr]->ComputeSoundSpeedSquare(rhor, er);

  if(rhor<=0 || cr<0) {
    fprintf(stdout,"*** Error: Negative density or c^2 (square of sound speed) in ComputeRiemannSolution(r)."
	" rho = %e, u = %e, p = %e, e = %e, c^2 = %e, ID = %d.\n",
	rhor, ur, pr, er, cr, idr);
    exit(-1);
  } else
    cr = sqrt(cr);


  // Declare variables in the "star region"
  double p0(DBL_MIN), ul0(0.0), ur0(0.0), rhol0(DBL_MIN), rhor0(DBL_MIN);
  double p1(DBL_MIN), ul1(0.0), ur1(0.0), rhol1(DBL_MIN), rhor1(DBL_MIN); //Secand Method ("k-1","k" in Kamm, (19))
  double p2(DBL_MIN), ul2(0.0), ur2(0.0), rhol2(DBL_MIN), rhor2(DBL_MIN); // "k+1"
  double f0(0.0), f1(0.0), f2(0.0); //difference between ul and ur


  // -------------------------------
  // Now, Solve The Riemann Problem
  // -------------------------------
  bool success = true;

  // monitor if the solution involves a transonic rarefaction. This is special as the solution 
  // at xi = x = 0 is within the rarefaction fan.
  bool trans_rare = false;
  double Vrare_x0[3]; //rho, u, and p at x = 0, in the case of a transonic rarefaction


  // A Trivial Case
  if(ul == ur && pl == pr) {
    FinalizeSolution(dir, Vm, Vp, rhol, ul, pl, idl, rhor, ur, pr, idr, rhol, rhor, ul, pl, 
	trans_rare, Vrare_x0, //inputs
	Vs, id, Vsm, Vsp/*outputs*/);
    return 0;
  }

  // -------------------------------
  // Step 1: Initialization
  //         (find initial interval [p0, p1])
  // -------------------------------
  success = FindInitialInterval(rhol, ul, pl, el, cl, idl, rhor, ur, pr, er, cr, idr, /*inputs*/
      p0, rhol0, rhor0, ul0, ur0,
      p1, rhol1, rhor1, ul1, ur1/*outputs*/);
  /* our convention is that p0 < p1 */

  if(!success) { //failed to find a bracketing interval. Output the state corresponding smallest "f"

    // get sol1d, trans_rare and Vrare_x0
#if PRINT_RIEMANN_SOLUTION == 1
    sol1d.clear();
#endif
    success = ComputeRhoUStar(1, integrationPath1, rhol, ul, pl, p1, idl, rhol0, rhol0*1.1, rhol2, ul2,
	&trans_rare, Vrare_x0/*filled only if found a trans. rarefaction*/);
    success = success && ComputeRhoUStar(3, integrationPath3, rhor, ur, pr, p1, idr, rhor0, rhol0*1.1, rhor2, ur2,
	&trans_rare, Vrare_x0/*filled only if found a trans. rarefaction*/);

    if(!success) {
      cout << "Warning: Riemann solver failed to find an initial bracketing interval (Returning initial states as solution)." << endl;
      for(int i=0; i<5; i++) {
	Vsm[i] = Vm[i];
	Vsp[i] = Vp[i];
      }
      double u_avg = 0.5*(ul+ur);
      if(u_avg >= 0) {
	id = idl;
	for(int i=0; i<5; i++)
	  Vs[i] = Vm[i];
      } else {
	id = idr;
	for(int i=0; i<5; i++)
	  Vs[i] = Vp[i];
      }

      return 1;
    }

    FinalizeSolution(dir, Vm, Vp, rhol, ul, pl, idl, rhor, ur, pr, idr, rhol2, rhor2, 0.5*(ul2+ur2), p1,
	trans_rare, Vrare_x0, /*inputs*/
	Vs, id, Vsm, Vsp /*outputs*/);

    return 1;
  }

  f0 = ul0 - ur0;
  f1 = ul1 - ur1;

#if PRINT_RIEMANN_SOLUTION == 1
  cout << "Found initial interval: p0 = " << p0 << ", f0 = " << f0 << ", p1 = " << p1 << ", f1 = " << f1 << endl;
#endif

  // -------------------------------
  // Step 2: Main Loop (Secant Method, Safeguarded) 
  // -------------------------------
  double denom = 0; 

  int iter = 0;
  double err_p = 1.0, err_u = 1.0;

  //p2 (and f2) is always the latest one
  p2 = p1; 
  f2 = f1;

#if PRINT_RIEMANN_SOLUTION == 1
  sol1d.clear();
#endif

  for(iter=0; iter<maxIts_main; iter++) {

    // 2.1: Update p using the Brent method (safeguarded secant method)
    denom = f1 - f0;
    if(denom == 0) {
      cout << "Warning: Division-by-zero while using the secant method to solve the Riemann problem." << endl;
      cout << "         left state: " << rhol << ", " << ul << ", " << pl << ", " << idl << " | right: " 
	<< rhor << ", " << ur << ", " << pr << ", " << idr << endl;
      cout << "         dir = " << dir[0] << "," << dir[1] << "," << dir[2]
	<< ", f0 = " << f0 << ", f1 = " << f1 << endl;
      p2 = 0.5*(p0+p1);
    } else {
      p2 = p2 - f2*(p1-p0)/denom;  // update p2
      if(p2<=p0 || p2>=p1) //discard and switch to bisection
	p2 = 0.5*(p0+p1);
    }

    //fprintf(stdout,"iter = %d, p0 = %e, p1 = %e, p2 = %e, f0 = %e, f1 = %e.\n", iter, p0, p1, p2, f0, f1);


try_again:

    // 2.2: Calculate ul2, ur2 
    success = ComputeRhoUStar(1, integrationPath1, rhol, ul, pl, p2, idl/*inputs*/, 
	rhol0, rhol1/*initial guesses for Hugo. eq.*/,
	rhol2, ul2/*outputs*/, 
	&trans_rare, Vrare_x0/*filled only if found a trans. rarefaction*/);

    if(!success) {
      //fprintf(stdout,"*** Error: Exact Riemann solver failed. left: %e %e %e (%d) | right: %e %e %e (%d).\n", 
      //        rhol, ul, pl, idl, rhor, ur, pr, idr);
      p2 = 0.5*(p2+p0); //move closer to p0

      iter++;
      if(iter<maxIts_main)
	goto try_again;
      else
	break;
    }

    success = ComputeRhoUStar(3, integrationPath3, rhor, ur, pr,  p2, idr/*inputs*/, 
	rhor0, rhor1/*initial guesses for Hugo. erq.*/,
	rhor2, ur2/*outputs*/,
	&trans_rare, Vrare_x0/*filled only if found a trans. rarefaction*/);

    if(!success) {
      //fprintf(stdout,"*** Error: Exact Riemann solver failed (2). left: %e %e %e (%d) | right: %e %e %e (%d).\n", 
      //        rhol, ul, pl, idl, rhor, ur, pr, idr);
      p2 = 0.5*(p2+p1); //move closer to p1

      iter++;
      if(iter<maxIts_main)
	goto try_again;
      else
	break;
    }

    f2 = ul2 - ur2;


    // 2.3: Update for the next iteration
    if(f0*f2<0.0) {
      p1 = p2;
      f1 = f2;
      rhol1 = rhol2;
      rhor1 = rhor2;
    } else {
      p0 = p2;
      f0 = f2;
      rhol0 = rhol2;
      rhor0 = rhor2;
    }


    // 2.4: Check stopping criterion
    err_p = fabs(p1 - p0)/std::max(fabs(pl + 0.5*rhol*ul*ul), fabs(pr + 0.5*rhor*ur*ur));
    err_u = fabs(f2)/std::max(cl, cr);

#if PRINT_RIEMANN_SOLUTION == 1
    cout << "Iter " << iter << ": p-interval = [" << p0 << ", " << p1 << "], p2 = " << p2 << ", err_p = " << err_p 
      << ", err_u = " << err_u << "." << endl;
#endif

    if( (err_p < tol_main && err_u < tol_main) || (err_p < tol_main*1e-3) || (err_u < tol_main*1e-3) )
      break; // converged

    trans_rare = false; //reset

#if PRINT_RIEMANN_SOLUTION == 1
    sol1d.clear();
#endif

  }


  // -------------------------------
  // Step 3: Find state at xi = x = 0 (for output)
  // -------------------------------
  double u2 = 0.5*(ul2 + ur2);
  FinalizeSolution(dir, Vm, Vp, rhol, ul, pl, idl, rhor, ur, pr, idr, rhol2, rhor2, u2, p2, 
      trans_rare, Vrare_x0, //inputs
      Vs, id, Vsm, Vsp/*outputs*/);


  if(iter == maxIts_main) {
    if(verbose>=1) {
      cout << "Warning: Exact Riemann solver failed to converge. err_p = " << err_p
	<< ", err_u = " << err_u << "." << endl;
      cout << "    Vm = [" << Vm[0] << ", " << Vm[1] << ", " << Vm[2] << ", " << Vm[3] << ", " << Vm[4] << "] (" << idl
	<< "),  Vp = [" << Vp[0] << ", " << Vp[1] << ", " << Vp[2] << ", " << Vp[3] << ", " << Vp[4] << "] (" << idr
	<< ")" << endl;
    }

    return 1;
  }

#if PRINT_RIEMANN_SOLUTION == 1
  std::cout << "Star State: (rhols, rhors, us, ps): " << rhol2 << ", " << rhor2 << ", " << u2 << ", " << p2 << "." << std::endl;
#endif


  //success!
  return 0;

}

//----------------------------------------------------------------------------------

//! non-adptive version, Connect the left/right initial state with the left/right star state (the 1-wave or 3-wave)
  bool  //true: success  | false: failure
ExactRiemannSolverNonAdaptive::ComputeRhoUStar(int wavenumber /*1 or 3*/,
    std::vector<std::vector<double>>& integrationPath /*3 by n, first index: 1-pressure, 2-density, 3-velocity*/, double rho, double u, double p, double ps, int id/*inputs*/,
    double rhos0, double rhos1/*initial guesses for Hugo. eq.*/,
    double &rhos, double &us/*outputs*/, 
    bool *trans_rare, double *Vrare_x0/*filled only if found tran rf*/)
{
  // default
  rhos = rho;
  us   = u;


  if(p > ps) {//rarefaction --- numerical integration

    // prepare for numerical integration
    double rhos_0 = rho, us_0 = u, ps_0 = p, xi_0; //start point of each step
    double rhos_1 = rho, us_1 = u, ps_1 = p, xi_1; //end point of each step
    double dp_target = (p-ps)/numSteps_rarefaction;

    double e = vf[id]->GetInternalEnergyPerUnitMass(rho,p);
    double dpdrho = vf[id]->GetDpdrho(rho,e);
    double drho_op1 = (p-ps)/dpdrho/numSteps_rarefaction; //use EOS
    double drho_op2 = rho/(2.5*numSteps_rarefaction); //initial step size     

    double drho_target = std::min(drho_op1, drho_op2);
    double dp = std::min(dp_target, drho_target*dpdrho);

    double pressure_endpoint_tol = tol_rarefaction * std::max( 1.0, fabs(p) );

    double c = vf[id]->ComputeSoundSpeedSquare(rho, e);
    if(rho<=0 || c<0) {
      fprintf(stdout,"Warning: Negative density or c^2 (square of sound speed) in ComputeRhoUStar." 
	  "rho = %e, p = %e, e = %e, c^2 = %e, id = %d.\n",
	  rho, p, e, c, id);
      return false; //failure
    } else
      c = sqrt(c);

    /*    for (size_t j = 0; j < integrationPath1[0].size(); j++) {
	  std::cout << integrationPath1[0][j] << ", " << integrationPath1[1][j] << ", " << integrationPath1[2][j] << std::endl;
	  }
	  */

    int index0 = 0;
    if (integrationPath.size() > 1) {
      for (int j = integrationPath.size()-1; j >= 0; j--) {
	if (integrationPath[j][0] > ps) {
	  index0 = j;
	  break;
	}
      }
      ps_0 = integrationPath[index0][0];
      rhos_0 = integrationPath[index0][1];
      us_0 = integrationPath[index0][2];
      dp = std::min( (p-ps)/numSteps_rarefaction, ps_0 - ps );
    }

    double xi = (wavenumber == 1) ? u - c : u + c; // xi = u -/+ c
    xi_0 = xi;

#if PRINT_RIEMANN_SOLUTION == 1
    sol1d.push_back(vector<double>{xi, rho, u, p, (double)id});
#endif

    //fprintf(stdout,"rho = %e, p = %e, ps = %e\n", rho, p, ps);
    // integration by Runge-Kutta 4
    bool done = false;
    int moreSteps = 5;

    for(int i=0; i<numSteps_rarefaction*moreSteps; i++) {

      // Check if we have reached the final pressure ps
      if(fabs(ps_1 - ps) <= pressure_endpoint_tol) {
	rhos = rhos_1;
	us   = us_1;

	if(vf[id]->CheckState(rhos,ps,true)) { //true: silence
#if PRINT_RIEMANN_SOLUTION == 1
	  cout << "Rarefaction solver reached a nonphysical state!" << endl;         
#endif
	  return false;
	}

#if PRINT_RIEMANN_SOLUTION == 1
	cout << "  " << wavenumber << "-wave: rarefaction, integration completed in " << i << " steps" << endl;
	cout << "rhos_1, us_1, ps_1: " << rhos_1 << ", " << us_1 << ", " << ps_1 << "." << endl;
#endif
	done = true;

	break; //done!
      }


      bool success = Rarefaction_OneStepRK4(wavenumber/*1 or 3*/, id,
	  rhos_0, us_0, ps_0 /*start state*/, dp /*step size*/,
	  rhos_1, us_1, ps_1, xi_1 /*output: end state*/);

      //      fprintf(stdout,"RK4 step: rhos_0 = %e, us_0 = %e, ps_0 = %e, drho = %e, rhos_1 = %e, us_1 = %e, ps_1 = %e | ps = %e | success = %d.\n",
      //              rhos_0, us_0, ps_0, drho, rhos_1, us_1, ps_1, ps, success);
      if(!success) {
	dp = dp/2.0;
	if (dp == 0.0) {
	  break;
	}
	continue;
      }

      if (ps_1-ps < -pressure_endpoint_tol) {
	dp = ps_0 - ps;
	continue;
      }

      if (ps_1 < integrationPath[integrationPath.size()-1][0]) {
	std::vector<double> currentVect = {ps_1, rhos_1, us_1};
	integrationPath.push_back(currentVect);
      }  

#if PRINT_RIEMANN_SOLUTION == 1
      sol1d.push_back(vector<double>{xi_1, rhos_1, us_1, ps_1, (double)id});
#endif

      if(trans_rare && Vrare_x0 && xi_0*xi_1<=0) {//transonic rarefaction, crossing x = xi = 0
	*trans_rare = true;
	double w0 = fabs(xi_1), w1 = fabs(xi_0);
	double ww = w0 + w1;
	w0 /= ww;
	w1 /= ww;
	Vrare_x0[0] = w0*rhos_0 + w1*rhos_1;
	Vrare_x0[1] = w0*us_0   + w1*us_1;
	Vrare_x0[2] = w0*ps_0   + w1*ps_1;

#if PRINT_RIEMANN_SOLUTION == 1
	sol1d.push_back(vector<double>{0.0, Vrare_x0[0], Vrare_x0[1], Vrare_x0[2], (double)id});
#endif

      }

      //fprintf(stdout,"drho = %e, rho: %e -> %e,  u: %e -> %e,  p: %e -> %e\n", drho, rhos_0, rhos_1, us_0, us_1, ps_0, ps_1);

      // If the solver gets here, it means it hasn't reached the final pressure ps.
      // Adjust step size, then update state
      //

      dp = std::min( std::min(dp_target, ps_1-ps), //don't go beyond ps
	  4.0*dp ); //don't increase too much in one step

      rhos_0 = rhos_1;
      us_0   = us_1;
      ps_0   = ps_1;
      xi_0   = xi_1;
    }

    if(!done) {
      if(vf[id]->CheckState(rhos_1,ps_1,true)) {
#if PRINT_RIEMANN_SOLUTION == 1
	cout << "  " << wavenumber << "-wave: rarefaction, solver failed (unphysical state: rhos = "
	  << rhos_1 << ", ps = " << ps_1 << "!)" << endl;
#endif
	return false; //failed
      } else {
#if PRINT_RIEMANN_SOLUTION == 1
	cout << "  " << wavenumber << "-wave: rarefaction, solver did not converge (final sol.: rhos_1 = "
	  << rhos_1 << ", ps_1 = " << ps_1 << "; inputs: rho = " << rho << ", p = " << p << ", ps = " << ps << ")" << endl;
#endif
	return false; 
      }
    }
  }

  else {// shock (p<=ps, rho<=rhos)

    HugoniotEquation hugo(vf[id],rho,p,ps);

    //find a bracketing interval
    double f0, f1;
    double drho = std::max(fabs(rhos0 - rhos1), 0.001*rhos0); 
    bool found_rhos0 = false, found_rhos1 = false;
    if(std::min(rhos0,rhos1)>=rho) {//both rhos0 and rhos1 are physically admissible
      f0 = hugo(rhos0);
      f1 = hugo(rhos1);
      if(f0*f1<=0) {
	/*found bracketing interval [rhos0, rhos1]*/
	if(rhos0>rhos1) {
	  std::swap(rhos0,rhos1);
	  std::swap(f0,f1);
	}
	found_rhos0 = found_rhos1 = true;
      } else {
	rhos0 = rhos1; //this is our starting point (presumably rhos1 is closer to sol'n)
	f0    = f1; 
      }
    } else { //at least, the smaller one among rhos0, rhos1 is non-physical
      if(rhos1>rhos0) {
	rhos0 = rhos1; 
	f0    = hugo(rhos0);
      }
      if(rhos0<rho) {//this one is also non-physical
	rhos0 = rho;
	f0    = hugo(rhos0);
	found_rhos0 = true;
      } else {
	/*rhos0 = rhos0;*/
	f0 = hugo(rhos0);
      }
    } 

    if(!found_rhos0 || !found_rhos1) {
      int i = 0;
      double factor = 1.5; 
      double tmp, ftmp;
      // before the search, rhos0 = rhos1 = an adimissible point > rho
      rhos1 = rhos0;
      f1    = f0;
      while(!found_rhos0) {
	if(++i>=maxIts_shock) {
	  //          cout << "*** Error: Unable to find a bracketing interval after " << maxIts_shock 
	  //               << " iterations (in the solution of the Hugoniot equation)." << endl;
	  return false;
	}
	tmp = rhos1;
	ftmp = f1;
	//move to the left (towards rho)
	rhos1 = rhos0;
	f1    = f0;
	rhos0 = rhos1 - factor*drho;
	if(rhos0<=rho) {
	  rhos0 = rho;
	  found_rhos0 = true;
	}
	f0 = hugo(rhos0);

	if(f0*f1<=0) {
	  found_rhos0 = found_rhos1 = true;
	} else {
	  //move to the right
	  rhos1 = tmp;
	  f1    = ftmp;
	  tmp   = rhos0; //don't forget the smallest point
	  ftmp  = f0;
	  rhos0 = rhos1;
	  f0    = f1;
	  rhos1 = rhos0 + factor*drho;
	  f1 = hugo(rhos1);
	  if(f0*f1<=0) {
	    found_rhos0 = found_rhos1 = true;
	  } else {
	    rhos0 = tmp;
	    f0    = ftmp;
	    drho  = rhos1 - rhos0; //update drho
	  }
	}
      }

      if(!found_rhos1) {//keep moving to the right
	i = 0;
	double factor = 2.5;
	while(!found_rhos1) {
	  if(++i>=maxIts_shock) {
	    //            cout << "*** Error: Unable to find a bracketing interval after " << maxIts_shock 
	    //                 << " iterations (in the solution of the Hugoniot equation (2))." << endl;
	    return false; //failure
	  }
	  rhos0 = rhos1;
	  f0    = f1;
	  rhos1 = rhos0 + factor*drho;
	  f1    = hugo(rhos1);
	  if(f0*f1<=0) {
	    found_rhos1 = true;
	  } else
	    drho = rhos1 - rhos0;
	}
      }

    }

    pair<double,double> sol;
    double loc_tol_shock = tol_shock*std::min(rhos0,rhos1);
#ifndef WITHOUT_BOOST
    //*******************************************************************
    // Calling boost function for root-finding
    // Warning: "maxit" is BOTH AN INPUT AND AN OUTPUT
    boost::uintmax_t maxit = maxIts_shock;
    if(f0==0.0)
      sol.first = sol.second = rhos0;
    else if(f1==0.0)
      sol.first = sol.second = rhos1;
    else {
      sol = toms748_solve(hugo, rhos0, rhos1, f0, f1,
	  [=](double r0, double r1){return r1-r0<std::min(loc_tol_shock,0.001*(rhos1-rhos0));}, 
	  maxit);
    }
    //*******************************************************************
#else
    //*******************************************************************
    // Using a hybrid (Brent) method for root-finding
    int maxit = 0;
    if(f0==0.0)
      sol.first = sol.second = rhos0;
    else if(f1==0.0)
      sol.first = sol.second = rhos1;
    else {
      double rhos2 = rhos1; //rhos2 is always the latest one
      double f2    = f1;
      int it;
      for(it = 0; it<maxIts_shock; it++) {
	drho = rhos1 - rhos0;
	rhos2 = rhos2 - f2*(rhos1 - rhos0)/(f1 - f0); //secant method
	if(rhos2 >= rhos1 || rhos2 <= rhos0) //discard and switch to bisection
	  rhos2 = 0.5*(rhos0+rhos1);
	f2 = hugo(rhos2);
	if(f2==0.0) {
	  sol.first = sol.second = rhos2;
	  break;
	}
	if(f2*f0<0) {
	  rhos1 = rhos2;
	  f1    = f2;
	} else {
	  rhos0 = rhos2;
	  f0    = f2;
	}
	if(rhos1-rhos0<loc_tol_shock) {
	  sol.first  = rhos0;
	  sol.second = rhos1;
	  break;
	}
      }
      if(it==maxIts_shock) {
	fprintf(stdout,"*** Error: Root-finding method failed to converge after %d iterations.\n", it);
	return false;
      }
      maxit = it;
    }
    //*******************************************************************
#endif


#if PRINT_RIEMANN_SOLUTION == 1
    cout << "  " << wavenumber << "-wave: shock, converged in " << maxit << " iterations. fun = " 
      << hugo(0.5*(sol.first+sol.second)) << "." << endl;
#endif

    rhos = 0.5*(sol.first+sol.second);

    double du = -(ps-p)*(1.0/rhos-1.0/rho);
    if(du<0) {
      //cout << "Warning: Violation of hyperbolicitiy when enforcing the Rankine-Hugoniot jump conditions (du = "
      //     << du << ")." << endl;
      return false;
    }

    if(vf[id]->CheckState(rhos,ps,true)) //true: silence
      return false; //nonphysical...

    us = (wavenumber==1) ? u - sqrt(du) : u + sqrt(du);


#if PRINT_RIEMANN_SOLUTION == 1
    double xi = (rhos*us - rho*u)/(rhos-rho);
    if(wavenumber==1) {
      sol1d.push_back(vector<double>{xi-0.0001*fabs(xi), rho, u, p, (double)id});
      sol1d.push_back(vector<double>{xi, rhos, us, ps, (double)id});
    } else {
      sol1d.push_back(vector<double>{xi, rhos, us, ps, (double)id});
      sol1d.push_back(vector<double>{xi+0.0001*fabs(xi), rho, u, p, (double)id});
    }
#endif

  }

  return true; //yeah!
}

//----------------------------------------------------------------------------------
  bool
ExactRiemannSolverNonAdaptive::Rarefaction_OneStepRK4(int wavenumber/*1 or 3*/, int id,
    double rho_0, double u_0, double p_0 /*start state*/, 
    double dp /*step*/,
    double &rho, double &u, double &p, double &xi /*output*/)
{
  dp = -dp; // dp is positive when passed in. It is actually negative if we follow Kamm's paper 
  // Equations (36 - 42)

  double e_0 = vf[id]->GetInternalEnergyPerUnitMass(rho_0, p_0);
  double c_0_square = vf[id]->ComputeSoundSpeedSquare(rho_0, e_0);

  if(rho_0<=0 || c_0_square<0) {
    //    fprintf(stdout,"*** Error: Negative density or c^2 (square of sound speed, %e) in Rarefaction_OneStepRK4(0)." 
    //            " rho = %e, p = %e, e = %e, id = %d.\n",
    //            c_0_square, rho_0, p_0, e_0, id);
    return false;
  } 

  double c_0 = sqrt(c_0_square);

  double rho_1 = rho_0 + 0.5*dp/c_0_square;
  double p_1 = p_0 + 0.5*dp;
  double e_1 = vf[id]->GetInternalEnergyPerUnitMass(rho_1, p_1);
  double c_1_square = vf[id]->ComputeSoundSpeedSquare(rho_1, e_1);

  if(rho_1<=0 || c_1_square<0) {
    //    fprintf(stdout,"*** Error: Negative density or c^2 (square of sound speed, %e) in Rarefaction_OneStepRK4(1)." 
    //                   " rho = %e, p = %e, e = %e, id = %d.\n",
    //                   c_1_square, rho_1, p_1, e_1, id);
    return false;
  } 

  double c_1 = sqrt(c_1_square);

  double rho_2 = rho_0 + 0.5*dp/c_1_square;
  double p_2 = p_1;
  double e_2 = vf[id]->GetInternalEnergyPerUnitMass(rho_2, p_2);
  double c_2_square = vf[id]->ComputeSoundSpeedSquare(rho_2, e_2);

  if(rho_2<=0 || c_2_square<0) {
    //    fprintf(stdout,"*** Error: Negative density or c^2 (square of sound speed, %e) in Rarefaction_OneStepRK4(2)." 
    //               " rho = %e, p = %e, e = %e, id = %d.\n",
    //                c_2_square, rho_2, p_2, e_2, id);
    return false;
  }

  double c_2 = sqrt(c_2_square); 

  double rho_3 = rho_0 + dp/c_2_square;
  double p_3 = p_0 + dp;
  double e_3 = vf[id]->GetInternalEnergyPerUnitMass(rho_3, p_3);
  double c_3_square = vf[id]->ComputeSoundSpeedSquare(rho_3, e_3);

  if(rho_3<=0 || c_3_square<0) {
    //    fprintf(stdout,"*** Error: Negative density or c^2 (square of sound speed, %e) in Rarefaction_OneStepRK4(3)." 
    //                " rho = %e, p = %e, e = %e, id = %d.\n",
    //                c_3_square, rho_3, p_3, e_3, id);
    return false;
  }  

  double c_3 = sqrt(c_3_square);

  // calculate the outputs
  //
  double drho = 1.0/6.0*dp*(1./c_0_square + 2.0*(1./c_1_square+1./c_2_square) + 1./c_3_square);
  double du = 1.0/6.0*dp*( 1./c_0/rho_0 + 2.0*(1./c_1/rho_1 + 1./c_2/rho_2) + 1./c_3/rho_3);

  if (isnan(du) == 1) {
    return false;
  }

  rho = rho_0 + drho;
  p = p_0 + dp;
  u = (wavenumber == 1) ? u_0 - du : u_0 + du; 

  double e = vf[id]->GetInternalEnergyPerUnitMass(rho, p);
  double c = vf[id]->ComputeSoundSpeedSquare(rho, e);

  //std::cout << std::setw(16) << p << std::setw(16) << rho << std::endl;

  if(rho<=0 || c<0) {
    //    fprintf(stdout,"*** Error: Negative density or c^2 (square of sound speed, %e) in Rarefaction_OneStepRK4(final)." 
    //               " rho = %e, p = %e, e = %e, id = %d.\n",
    //               c, rho, p, e, id);
    return false;
  } else

    c = sqrt(c);

  xi = (wavenumber == 1) ? u - c : u + c;

  return true;
}

//----------------------------------------------------------------------------------
  bool
ExactRiemannSolverBase::FindInitialFeasiblePointsOneSided(double rhol, double ul, double pl, double el, 
    double cl, int idl, double ustar, double &p0, double &rhol0, double &ul0, 
    double &p1, double &rhol1, double &ul1)
{
  double dp;
  int found = 0;
  bool success = true;

  // Method 1: Use the acoustic theory (Eqs. (20)-(22) of Kamm) to find p0, p1
  found = FindInitialFeasiblePointsOneSidedByAcousticTheory(rhol, ul, pl, el, cl, idl, ustar,
      p0, rhol0, ul0, p1, rhol1, ul1/*outputs*/);

  if(found==2)
    return true; //yeah

  if(found==1) //the first one (p0) is good --> only need to find p1
    goto myLabel;

  // Method 2: based on dp (fixed width search)

  // 2.1. find the first one (p0)
  dp = ul>ustar ? 0.5*pl : -0.5*pl;
  for(int i=0; i<maxIts_bracket; i++) {
    p0 = pl + 0.01*(i+1)*(i+1)*dp;

    if(p0<min_pressure)
      p0 = pressure_at_failure; 
    success = ComputeRhoUStar(1, integrationPath1, rhol, ul, pl, p0, idl, 
	rhol, (p0>pl) ? rhol*1.1 : rhol*0.9,
	rhol0, ul0);
    if(success)
      break;
  }
  if(!success) {
    if(verbose>=1)
      fprintf(stdout,"Warning: Failed to find the first initial guess (p0) in the 1D half-Riemann solver (it = %d). "
	  "Left: %e %e %e (%d), ustar: %e.\n",
	  maxIts_bracket, rhol, ul, pl, idl, ustar);
    return false;
  }

myLabel:
  // 2.2. find the second one (p1)
  dp = p0-pl;
  for(int i=0; i<maxIts_bracket; i++) {
    p1 = p0 + 0.01*(i+1)*(i+1)*dp;
    success = ComputeRhoUStar(1, integrationPath1, rhol, ul, pl, p1, idl, rhol, rhol0, rhol1, ul1);
    if(success)
      break;
  }
  if(!success) {
    if(verbose>=1)
      fprintf(stdout,"Warning: Failed to find the second initial guess (p1) in the 1D half-Riemann solver (it = %d). "
	  "Left: %e %e %e (%d), ustar: %e.\n",
	  maxIts_bracket, rhol, ul, pl, idl, ustar);
    return false;
  }

  // Make sure p0<p1
  if(p0>p1) {
    std::swap(p0,p1);
    std::swap(rhol0, rhol1);
    std::swap(ul0, ul1);
  } 

  return true;
}

//----------------------------------------------------------------------------------

  double ExactRiemannSolverBase::GetSurfaceTensionCoefficient() {
    std::cout << "ExactRiemannSolverBase::GetSurfaceTensionCoefficient: This function in the base class should be overrided by the one in the derived class \"ExactRiemannSolverInterfaceJump\". Saying this message indicating something is wrong." << std::endl;
    exit(-1);
  }
