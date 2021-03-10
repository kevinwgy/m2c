#include<ExactRiemannSolverBase.h>
#include<assert.h>
#include<utility> //std::pair
#include<bits/stdc++.h> //std::swap
#include<boost/math/tools/roots.hpp>
using namespace boost::math::tools;
using std::pair;
using std::cout;
using std::endl;

//-----------------------------------------------------

ExactRiemannSolverBase::ExactRiemannSolverBase(std::vector<VarFcnBase*> &vf_, 
                            ExactRiemannSolverData &iod_riemann) : vf(vf_)
{
  maxIts_main          = iod_riemann.maxIts_main;
  maxIts_shock         = iod_riemann.maxIts_shock;
  numSteps_rarefaction = iod_riemann.numSteps_rarefaction;
  tol_main             = iod_riemann.tol_main;
  tol_shock            = iod_riemann.tol_shock;
  tol_rarefaction      = iod_riemann.tol_rarefaction;

  // TODO: some of these parameters should goto IoData
  maxIts_main = 200;
  maxIts_shock = 200;   
  numSteps_rarefaction = 100;
  tol_main = 1.0e-5;
  tol_shock = 1.0e-6;
  tol_rarefaction = 1.0e-5;
}

//-----------------------------------------------------
/** Solves the one-dimensional Riemann problem. Extension of Kamm 2015 
 * to Two Materials. See KW's notes for details
 */
void
ExactRiemannSolverBase::ComputeRiemannSolution(int dir/*0~x,1~y,2~z*/, 
                            double *Vm, int idl /*"left" state*/, 
                            double *Vp, int idr /*"right" state*/, 
                            double *V, int &id /*solution at xi = 0 (i.e. x=0) */)
{
  // Convert to a 1D problem (i.e. One-Dimensional Riemann)
  double rhol  = Vm[0];
  double ul    = Vm[dir+1];
  double pl    = Vm[4];
  double rhor  = Vp[0];
  double ur    = Vp[dir+1];
  double pr    = Vp[4];

  double el = vf[idl]->GetInternalEnergyPerUnitMass(rhol, pl);
  double cl = vf[idl]->ComputeSoundSpeed(rhol, el);
  double er = vf[idr]->GetInternalEnergyPerUnitMass(rhor, pr);
  double cr = vf[idr]->ComputeSoundSpeed(rhor, er);

  // Declare vairables in the "star region"
  double p0, ul0, ur0, rhol0, rhor0;
  double p1, ul1, ur1, rhol1, rhor1; //Secand Method ("k-1","k" in Kamm, (19))
  double p2, ul2, ur2, rhol2, rhor2; // "k+1"

  // -------------------------------
  // Now, Solve The Riemann Problem
  // -------------------------------

  // -------------------------------
  // Step 1: Initialization
  // -------------------------------
 
  // 1.1: Initialize p0 using acoustic theory ((20) of Kamm) 
  double Cl = rhol*cl; //acoustic impedance
  double Cr = rhor*cr; //acoustic impedance
  p0 = (Cr*pl + Cl*pr + Cl*Cr*(ul - ur))/(Cl + Cr);

  // 1.2: Calculate ul0, ur0
  ComputeRhoUStar(1, rhol, ul, pl, p0, idl/*inputs*/, 
                  rhol, (p0>pl) ? rhol*1.1 : rhol*0.9/*initial guesses for Hugo. eq.*/,
                  rhol0, ul0/*outputs*/, NULL, NULL, true);
  ComputeRhoUStar(3, rhor, ur, pr, p0, idr/*inputs*/, 
                  rhor, (p0>pr) ? rhor*1.1 : rhor*0.9/*initial guesses for Hugo. eq.*/,
                  rhor0, ur0/*outputs*/, NULL, NULL, true);
  double f0 = ul0 - ur0; // Eq. (18) in Kamm

  // 1.3: Initialize p1 ((21)(22) of Kamm) 
  double Clbar = (ul0 == ul) ? Cl : fabs(p0 - pl)/fabs(ul0 - ul);
  double Crbar = (ur0 == ur) ? Cr : fabs(p0 - pr)/fabs(ur0 - ur);
  p1 = (Crbar*pl + Clbar*pr + Clbar*Crbar*(ul - ur))/(Clbar + Crbar); 

  // 1.4 Calculate ul1, ur1 
  ComputeRhoUStar(1, rhol, ul, pl, p1, idl/*inputs*/, 
                  rhol, rhol0/*initial guesses for Hugo. eq.*/,
                  rhol1, ul1/*outputs*/, NULL, NULL, true);
  ComputeRhoUStar(3, rhor, ur, pr, p1, idr/*inputs*/, 
                  rhor, rhor0/*initial guesses for Hugo. eq.*/,
                  rhor1, ur1/*outputs*/, NULL, NULL, true);
  double f1 = ul1 - ur1; // Eq. (18) in Kamm

  // -------------------------------
  // Step 2: Main Loop (Secant Method) 
  // -------------------------------
  double denom = 0, f2 = 0;

  // monitor if the solution involves a transonic rarefaction. This is special as the solution 
  // at xi = x = 0 is within the rarefaction fan.
  bool trans_rare = false;
  double Vrare_x0[3]; //rho, u, and p at x = 0, in the case of a transonic rarefaction

  int iter = 0;
  double err_p = 1.0, err_u = 1.0;

  p2 = p1; //to avoid compiler warning

  for(iter=0; iter<maxIts_main; iter++) {

    // 2.1: Update p using the Secant method
    denom = f1 - f0;
    if(denom == 0) {
      cout << "ERROR: Division-by-zero while using the secant method to solve the Riemann problem." << endl;
      exit_mpi();
    }
    p2 = p1 - f1*(p1-p0)/denom;

    // 2.2: Calculate ul2, ur2 
    ComputeRhoUStar(1, rhol, ul, pl, p2, idl/*inputs*/, 
                    rhol0, rhol1/*initial guesses for Hugo. eq.*/,
                    rhol2, ul2/*outputs*/, 
                    &trans_rare, Vrare_x0/*filled only if found a trans. rarefaction*/);
    ComputeRhoUStar(3, rhor, ur, pr,  p2, idr/*inputs*/, 
                    rhor0, rhor1/*initial guesses for Hugo. erq.*/,
                    rhor2, ur2/*outputs*/,
                    &trans_rare, Vrare_x0/*filled only if found a trans. rarefaction*/);
    f2 = ul2 - ur2;
    
    // 2.3: Check stopping criterion
    err_p = fabs(p2 - p1)/std::max(fabs(p1), fabs(p2));
    err_u = fabs(f2)/std::max(fabs(ul2), fabs(ur2));

#if PRINT_RIEMANN_SOLUTION == 1
    cout << "Iter " << iter << ": err_p = " << err_p << ", err_u = " << err_u << "." << endl;
#endif

    if( err_p < tol_main && err_u < tol_main )
      break; // converged

    // 2.4: Update for the next iteration
    p0 = p1;
    p1 = p2;
    f0 = f1;
    f1 = f2;
    rhol0 = rhol1;
    rhol1 = rhol2;
    rhor0 = rhor1;
    rhor1 = rhor2;
    trans_rare = false;

#if PRINT_RIEMANN_SOLUTION == 1
    sol1d.clear();
#endif

  }

  if(iter == maxIts_main) {
    cout << "ERROR: Exact Riemann solver failed to converge. err_p = " << err_p
         << ", err_u = " << err_u << "." << endl;
    exit_mpi();
  }
  

  // -------------------------------
  // Step 3: Find state at xi = x = 0 (for output)
  // -------------------------------
  double u2 = 0.5*(ul2 + ur2);

  // find material id at xi = x = 0
  if(u2>=0)
    id = idl;
  else
    id = idr;

#if PRINT_RIEMANN_SOLUTION == 1
  sol1d.push_back(vector<double>{u2 - std::max(1e-6, 0.001*fabs(u2)), rhol2, u2, p2});
  sol1d.push_back(vector<double>{u2, rhor2, u2, p2});
#endif

  V[0] = V[1] = V[2] = V[3] = V[4] = 0.0;

  if(trans_rare) {
    V[0]     = Vrare_x0[0];
    V[dir+1] = Vrare_x0[1];
    V[4]     = Vrare_x0[2];
  }
  else { 
    //find state variables at xi = x = 0

    if(u2>=0) { //either Vl or Vlstar --- check the 1-wave

      bool is_star_state = false;

      if(pl >= p2) {//1-wave is rarefaction
        double el2 = vf[idl]->GetInternalEnergyPerUnitMass(rhol2, p2);
        double cl2 = vf[idl]->ComputeSoundSpeed(rhol2, el2);
        if(u2 - cl2 <= 0) //rarefaction tail speed
          is_star_state = true;
      } 
      else {//1-wave is shock
        double us = (rhol2*u2 - rhol*ul)/(rhol2 - rhol); //shock speed
        if(us <= 0)      
          is_star_state = true;
      }

      if(is_star_state) {
        V[0]     = rhol2;
        V[dir+1] = u2;
        V[4]     = p2;
      } else {
        V[0]     = rhol;
        V[dir+1] = ul;
        V[4]     = pl;
      }

    } else { //either Vr or Vrstar --- check the 3-wave

      bool is_star_state = false;

      if(pr >= p2) {//3-wave is rarefaction
        double er2 = vf[idr]->GetInternalEnergyPerUnitMass(rhor2, p2);
        double cr2 = vf[idr]->ComputeSoundSpeed(rhor2, er2);
        if(u2 - cr2 >= 0)
          is_star_state = true;
      }
      else {//3-wave is shock
        double us = (rhor2*u2 - rhor*ur)/(rhor2 - rhor); //shock speed
        if(us >= 0)
          is_star_state = true;
      }

      if(is_star_state) {
        V[0]     = rhor2;
        V[dir+1] = u2;
        V[4]     = p2;
      } else {
        V[0]     = rhor;
        V[dir+1] = ur;
        V[4]     = pr;
      }

    }
  }

  // determine the tangential components of velocity -- upwinding
  int k;
  if(u2>0) {
    for(int i=1; i<=2; i++) {
      k = (dir + i)%3 + 1;
      V[k] = Vm[k];
    }
  } else if(u2<0) {
    for(int i=1; i<=2; i++) {
      k = (dir + i)%3 + 1;
      V[k] = Vp[k];
    }
  } else {//u2 == 0
    for(int i=1; i<=2; i++) {
      k = (dir + i)%3 + 1;
      V[k] = 0.5*(Vm[k]+Vp[k]);
    }
  }


#if PRINT_RIEMANN_SOLUTION == 1
  std::sort(sol1d.begin(), sol1d.end(), 
            [](vector<double> v1, vector<double> v2){return v1[0]<v2[0];});
  int last = sol1d.size()-1;
  double xi_span = sol1d[last][0] - sol1d[0][0];
  sol1d.insert(sol1d.begin(), vector<double>{sol1d[0][0]-xi_span, sol1d[0][1], sol1d[0][2], sol1d[0][3]});
  last++;
  sol1d.push_back(vector<double>{sol1d[last][0]+xi_span, sol1d[last][1], sol1d[last][2], sol1d[last][3]});

  FILE* solFile = fopen("RiemannSolution.txt", "w");
  print(solFile, "## One-Dimensional Riemann Problem.\n");
  print(solFile, "## Initial State: %e %e %e, id %d (left) | (right) %e %e %e, id %d.\n", 
        rhol, ul, pl, idl, rhor, ur, pr, idr);

  for(auto it = sol1d.begin(); it != sol1d.end(); it++) 
    print(solFile,"%e    %e    %e    %e\n", (*it)[0], (*it)[1], (*it)[2], (*it)[3]);

  fclose(solFile);
#endif

}

//----------------------------------------------------------------------------------

//! Connect the left/right initial state with the left/right star state (the 1-wave or 3-wave)
void
ExactRiemannSolverBase::ComputeRhoUStar(int wavenumber /*1 or 3*/,
                                        double rho, double u, double p, double ps, int id/*inputs*/,
                                        double rhos0, double rhos1/*initial guesses for Hugo. eq.*/,
                                        double &rhos, double &us/*outputs*/, 
                                        bool *trans_rare, double *Vrare_x0/*filled only if found tran rf*/,
                                        bool rough_estimate)
{
  // default
  rhos = rho;
  us   = u;

  int    N_rare   = rough_estimate ? 10   : numSteps_rarefaction;
  double tol_rare = rough_estimate ? 1e-2 : tol_rarefaction;
  double tol_sh   = rough_estimate ? 1e-2 : tol_shock;

  if(p > ps) {//rarefaction --- numerical integration

    double dp_max = 1.25*(p-ps)/N_rare;
    double dp_target = dp_max/1.25;

    // initialize drho
    double drho = rho/(N_rare*2.5); //initial step size

    // prepare for numerical integration
    double rhos_0 = rho, us_0 = u, ps_0 = p, xi_0; //start point of each step
    double rhos_1 = rho, us_1 = u, ps_1 = p, xi_1; //end point of each step
    double dp;

    double e = vf[id]->GetInternalEnergyPerUnitMass(rho, p);
    double c = vf[id]->ComputeSoundSpeed(rho, e);
    double xi = (wavenumber == 1) ? u - c : u + c; // xi = u -/+ c
    xi_0 = xi;

#if PRINT_RIEMANN_SOLUTION == 1
    sol1d.push_back(vector<double>{xi, rho, u, p});
#endif

    // integration by Runge-Kutta 4
    for(int i=0; i<numSteps_rarefaction*5; i++) {

      Rarefaction_OneStepRK4(wavenumber/*1 or 3*/, id,
                             rhos_0, us_0, ps_0 /*start state*/, drho /*step size*/,
                             rhos_1, us_1, ps_1, xi_1 /*output: end state*/); 

      dp = ps_0 - ps_1;

      //Check if it went too far. If so, rewind and reduce step size
      if(dp > dp_max) {
        drho = drho/dp*dp_target;
        continue;
      }
      if(ps_1 - ps < -tol_rare) {
        drho = drho/dp*(ps_0 - ps);
        continue;
      } 

#if PRINT_RIEMANN_SOLUTION == 1
    sol1d.push_back(vector<double>{xi_1, rhos_1, us_1, ps_1});
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
    sol1d.push_back(vector<double>{0.0, Vrare_x0[0], Vrare_x0[1], Vrare_x0[2]});
#endif

      }

      // Check if we have reached the final pressure ps
      if(fabs(ps_1 - ps) <= tol_rare) {
        rhos = rhos_1;
        us   = us_1;
      
#if PRINT_RIEMANN_SOLUTION == 1
        cout << "  " << wavenumber << "-wave: rarefaction, integration completed in " << i << " steps" << endl;
#endif

        break; //done!
      }

      // If the solver gets here, it means it hasn't reached the final pressure ps.
      // Adjust step size, then update state
      //
      drho = std::min( (rhos_0-rhos_1)/dp*std::min(dp_target,ps_1-ps), //don't go beyond ps
                       drho*4.0); //don't increase too much in one step
       
      rhos_0 = rhos_1;
      us_0   = us_1;
      ps_0   = ps_1;
      xi_0   = xi_1;

    }
  }

  else {// shock (p<=ps, rho<=rhos)

    HugoniotEquation hugo(vf[id],rho,p,ps);

    //find a bracketing interval
    double f0, f1;
    double drho = fabs(rhos0 - rhos1); 
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
          cout << "Error: Unable to find a bracketing interval after " << maxIts_shock 
               << " iterations (in the solution of the Hugoniot equation)." << endl;
          exit_mpi();
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
            cout << "Error: Unable to find a bracketing interval after " << maxIts_shock 
                 << " iterations (in the solution of the Hugoniot equation (2))." << endl;
            exit_mpi();
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

    //*******************************************************************
    // Calling boost function for root-finding
    // Warning: "maxit" is BOTH AN INPUT AND AN OUTPUT
    boost::uintmax_t maxit = maxIts_shock;
    pair<double,double> sol = toms748_solve(hugo, rhos0, rhos1, f0, f1,
                                            [=](double r0, double r1){return r1-r0<tol_sh;}, 
                                            maxit);
    //*******************************************************************
    
#if PRINT_RIEMANN_SOLUTION == 1
    cout << "  " << wavenumber << "-wave: shock, converged in " << maxit << " iterations. fun = " 
         << hugo(0.5*(sol.first+sol.second)) << "." << endl;
#endif

    rhos = 0.5*(sol.first+sol.second);

    double du = -(ps-p)*(1.0/rhos-1.0/rho);
    if(du<0) {
      cout << "Error: Violation of hyperbolicitiy when enforcing the Rankine-Hugoniot jump conditions (du = "
           << du << ")." << endl;
      exit_mpi();
    }
    us = (wavenumber==1) ? u - sqrt(du) : u + sqrt(du);
    

#if PRINT_RIEMANN_SOLUTION == 1
    double xi = (rhos*us - rho*u)/(rhos-rho);
    if(wavenumber==1) {
      sol1d.push_back(vector<double>{xi-0.0001*fabs(xi), rho, u, p});
      sol1d.push_back(vector<double>{xi, rhos, us, ps});
    } else {
      sol1d.push_back(vector<double>{xi, rhos, us, ps});
      sol1d.push_back(vector<double>{xi+0.0001*fabs(xi), rho, u, p});
    }
#endif

  }
}

//----------------------------------------------------------------------------------

void 
ExactRiemannSolverBase::Rarefaction_OneStepRK4(int wavenumber/*1 or 3*/, int id,
                            double rho_0, double u_0, double p_0 /*start state*/, 
                            double drho /*step size*/,
                            double &rho, double &u, double &p, double &xi /*output*/)
{
  drho = -drho; //drho is positive when passed in. It is actually negative if we follow Kamm's paper
                // Equations (36 - 42)
                
  double e_0 = vf[id]->GetInternalEnergyPerUnitMass(rho_0, p_0);
  double c_0 = vf[id]->ComputeSoundSpeed(rho_0, e_0);
  double c_0_square = c_0*c_0;

  double p_1 = p_0 + 0.5*drho*c_0_square;
  double rho_1 = rho_0 + 0.5*drho;
  double e_1 = vf[id]->GetInternalEnergyPerUnitMass(rho_1, p_1);
  double c_1 = vf[id]->ComputeSoundSpeed(rho_1, e_1);
  double c_1_square = c_1*c_1;

  double p_2 = p_0 + 0.5*drho*c_1_square;
  double rho_2 = rho_1;
  double e_2 = vf[id]->GetInternalEnergyPerUnitMass(rho_2, p_2);
  double c_2 = vf[id]->ComputeSoundSpeed(rho_2, e_2);
  double c_2_square = c_2*c_2;
 
  double p_3 = p_0 + drho*c_2_square;
  double rho_3 = rho_0 + drho;
  double e_3 = vf[id]->GetInternalEnergyPerUnitMass(rho_3, p_3);
  double c_3 = vf[id]->ComputeSoundSpeed(rho_3, e_3);
  double c_3_square = c_3*c_3;

  // now, calculate the outputs
  //
  p = p_0 + 1.0/6.0*drho*(c_0_square + 2.0*(c_1_square+c_2_square) + c_3_square);

  double du = 1.0/6.0*drho*( c_0/rho_0 + 2.0*(c_1/rho_1 + c_2/rho_2) + c_3/rho_3);
  u = (wavenumber == 1) ? u_0 - du : u_0 + du;

  rho = rho_0 + drho;

  double e = vf[id]->GetInternalEnergyPerUnitMass(rho, p);
  double c = vf[id]->ComputeSoundSpeed(rho, e);
  xi = (wavenumber == 1) ? u - c : u + c;

}

//----------------------------------------------------------------------------------

//----------------------------------------------------------------------------------
