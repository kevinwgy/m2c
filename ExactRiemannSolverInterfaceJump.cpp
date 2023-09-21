/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<ExactRiemannSolverInterfaceJump.h>
#include<array>
#include<utility> //std::pair
#include<bits/stdc++.h> //std::swap
#include<iostream>

using std::pair;
using std::cout;
using std::endl;

extern int verbose;

#define INVALID_MATERIAL_ID -1

ExactRiemannSolverInterfaceJump::ExactRiemannSolverInterfaceJump(std::vector<VarFcnBase*> &vf_, ExactRiemannSolverData &iod_riemann_) : ExactRiemannSolverBase(vf_, iod_riemann_) {
  surface_tension_coefficient = iod_riemann.surface_tension_coefficient;
  surface_tension_materialid = iod_riemann.surface_tension_materialid;
  delta_p = 0.0; 
}

//-----------------------------------------------------
/** Solves the one-dimensional Riemann problem. Extension of Kamm 2015 
 * to Two Materials. See KW's notes for details
 * Returns an integer error code
 * 0: no errors
 * 1: riemann solver failed to find a bracketing interval
 */
int
ExactRiemannSolverInterfaceJump::ComputeRiemannSolution(double *dir, 
    double *Vm, int idl /*"left" state*/, 
    double *Vp, int idr /*"right" state*/, 
    double *Vs, int &id /*solution at xi = 0 (i.e. x=0) */,
    double *Vsm /*left 'star' solution*/,
    double *Vsp /*right 'star' solution*/,
    double curvature)
{

  ComputePressureJump(idr, curvature);

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
//  if(ul == ur && pl + delta_p == pr) {
//    FinalizeSolution(dir, Vm, Vp, rhol, ul, pl, idl, rhor, ur, pr, idr, rhol, rhor, ul, pl, 
//	trans_rare, Vrare_x0, //inputs
//	Vs, id, Vsm, Vsp/*outputs*/);
//    return 0;
//  }

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
    success = success && ComputeRhoUStar(3, integrationPath3, rhor, ur, pr, p1 + delta_p, idr, rhor0, rhol0*1.1, rhor2, ur2,
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

      if(verbose>=1)
	cout << "Warning: Riemann solver failed to find an initial bracketing interval." << endl; 
      return 1;
    }

    FinalizeSolution(dir, Vm, Vp, rhol, ul, pl, idl, rhor, ur, pr, idr, rhol2, rhor2, 0.5*(ul2+ur2), p1,
	trans_rare, Vrare_x0, /*inputs*/
	Vs, id, Vsm, Vsp /*outputs*/);

    if(verbose>=1)
      cout << "Warning: Riemann solver failed to find an initial bracketing interval." << endl;
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

    success = ComputeRhoUStar(3, integrationPath3, rhor, ur, pr,  p2 + delta_p, idr/*inputs*/, 
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
	<< ")" << "pressure jump = " << delta_p << ", curvature = " << curvature << endl;
      cout << "dir = [" << dir[0] << ", " << dir[1] << ", " << dir[2] << "]" << endl;
    }

    if(verbose>=1)
      cout << "Warning: Exact Riemann solver (adaptive) failed to converge." << endl;
    return 1;
  }

#if PRINT_RIEMANN_SOLUTION == 1
  std::cout << "Star State: (rhols, rhors, us, ps): " << rhol2 << ", " << rhor2 << ", " << u2 << ", " << p2 << "." << std::endl;
#endif


  //success!
  return 0;

}


//----------------------------------------------------------------------------------
//! find a bracketing interval [p0, p1] (f0*f1<=0)
  bool
ExactRiemannSolverInterfaceJump::FindInitialInterval(double rhol, double ul, double pl, double el, double cl, int idl,
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
    success = success && ComputeRhoUStar(3, integrationPath3, rhor, ur, pr,  p2 + delta_p, idr, rhor0, rhor1, rhor2, ur2);

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
	success = success && ComputeRhoUStar(3, integrationPath3, rhor, ur, pr,  p2 + delta_p, idr, rhor0, rhor1, rhor2, ur2);


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
      success = success && ComputeRhoUStar(3, integrationPath3, rhor, ur, pr,  p2 + delta_p, idr, rhor0, rhor1, rhor2, ur2);
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

  bool
ExactRiemannSolverInterfaceJump::FindInitialFeasiblePoints(double rhol, double ul, double pl, double el, double cl, int idl,
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
    success = success && ComputeRhoUStar(3, integrationPath3, rhor, ur, pr, p0 + delta_p, idr, 
	rhor, (p0+delta_p>pr) ? rhor*1.1 : rhor*0.9,
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
      success = success && ComputeRhoUStar(3, integrationPath3, rhor, ur, pr, p0 + delta_p, idr, 
	  rhor, (p0+delta_p>pr) ? rhor*1.1 : rhor*0.9,
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
    success = success && ComputeRhoUStar(3, integrationPath3, rhor, ur, pr, p1 + delta_p, idr, rhor, rhor0, rhor1, ur1);
    if(success)
      break;
  }
  if(!success) //search in the opposite direction
    for(int i=0; i<maxIts_bracket; i++) {
      p1 = p0 - 0.01*(i+1)*(i+1)*dp;
      if(p1<min_pressure)
	p1 = pressure_at_failure*1000.0; //so it is not the same as p0
      success = ComputeRhoUStar(1, integrationPath1, rhol, ul, pl, p1, idl, rhol, rhol0, rhol1, ul1);
      success = success && ComputeRhoUStar(3, integrationPath3, rhor, ur, pr, p1 + delta_p, idr, rhor, rhor0, rhor1, ur1);
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
ExactRiemannSolverInterfaceJump::FindInitialFeasiblePointsByAcousticTheory(double rhol, double ul, double pl,
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
  p0 = (Cr*pl + Cl*(pr-delta_p) + Cl*Cr*(ul - ur))/(Cl + Cr);

  success = ComputeRhoUStar(1, integrationPath1, rhol, ul, pl, p0, idl/*inputs*/,
      rhol, (p0>pl) ? rhol*1.1 : rhol*0.9/*initial guesses for Hugo. eq.*/,
      rhol0, ul0/*outputs*/);
  if(!success)
    return found;

  success = ComputeRhoUStar(3, integrationPath3, rhor, ur, pr, p0 + delta_p, idr/*inputs*/,
      rhor, (p0+delta_p>pr) ? rhor*1.1 : rhor*0.9/*initial guesses for Hugo. eq.*/,
      rhor0, ur0/*outputs*/);
  if(!success)
    return found;

  found = 1; //found p0!

  // 1.2. Initialize p1 ((21)(22) of Kamm)
  double Clbar = (ul0 == ul) ? Cl : fabs(p0 - pl)/fabs(ul0 - ul);
  double Crbar = (ur0 == ur) ? Cr : fabs(p0+delta_p - pr)/fabs(ur0 - ur);
  p1 = (Crbar*pl + Clbar*(pr-delta_p) + Clbar*Crbar*(ul - ur))/(Clbar + Crbar);
  double tmp = std::max(fabs(p0), fabs(p1));
  if(fabs(p1 - p0)/tmp<1.0e-8)
    p1 = p0 + 1.0e-8*tmp; //to avoid f0 = f1 (divide-by-zero)

  success = ComputeRhoUStar(1, integrationPath1, rhol, ul, pl, p1, idl/*inputs*/,
      rhol, rhol0/*initial guesses for Hugo. eq.*/,
      rhol1, ul1/*outputs*/);
  if(!success)
    return found;

  success = ComputeRhoUStar(3, integrationPath3, rhor, ur, pr, p1 + delta_p, idr/*inputs*/,
      rhor, rhor0/*initial guesses for Hugo. eq.*/,
      rhor1, ur1/*outputs*/);
  if(!success)
    return found;

  found = 2; //found p0 and p1!
  return found;
}


//-----------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------

  void
ExactRiemannSolverInterfaceJump::PrintStarRelations(double rhol, double ul, double pl, int idl,
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

    success = ComputeRhoUStar(3, integrationPath3, rhor, ur, pr, ps + delta_p, idr/*inputs*/,
	rhor, (ps+delta_p>pr) ? rhor*1.1 : rhor*0.9/*initial guesses for Hugo. eq.*/,
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

  void
ExactRiemannSolverInterfaceJump::FinalizeSolution(double *dir, double *Vm, double *Vp,
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
  sol1d.push_back(vector<double>{u2, rhor2, u2, p2+delta_p, (double)idr});
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
	double er2 = vf[idr]->GetInternalEnergyPerUnitMass(rhor2, p2+delta_p);
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
  Vsp[4] = p2 + delta_p;



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

//--------------------------------------------------------------
double ExactRiemannSolverInterfaceJump::GetSurfaceTensionCoefficient() { return surface_tension_coefficient;}











