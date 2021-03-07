#include<ExactRiemannSolverBase.h>

//-----------------------------------------------------

ExactRiemannSolverBase::ExactRiemannSolverBase(FluxFcnBase &ff_, 
                            std::vector<VarFcnBase*> &vf_, IoData &iod)
                      : ff(ff_), vf(vf_)
{
  // TODO: some of these parameters should goto IoData
  maxIts_main = 500;
  maxIts_shock = 100;   
  numSteps_rarefaction = 100;
  tol_main = 1.0e-3;
  tol_shock = 1.0e-3;
  tol_rarefaction = 1.0e-3;
}

//-----------------------------------------------------
/** Solves the one-dimensional Riemann problem. Extension of Kamm 2015 
 * to Two Materials. See KW's notes for details
 */
virtual void
ExactRiemannSolverBase::ComputeRiemannSolution(int dir/*0~x,1~y,2~z*/, 
                            double *Vm, int idl, /*"left" state*/, 
                            double *Vp, int idr, /*"right" state*/, 
                            double *V /*solution at xi = 0 (i.e. x=0) */)
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
  double p1, ul1, ur1, rhol1, rhor1; //two states in Secand Method ("k-1","k" in Kamm, (19))
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
  ComputeRhoUStarLeft(rhol, ul, pl, p0, idl/*inputs*/, 
                      rhol, (p0>pl) ? rhol*1.01 : rhol*0.99/*initial guesses for Hugo. eq.*/,
                      rhol0, ul0/*outputs*/);
  ComputeRhoUStarRight(rhor, ur, pr, p0, idr/*inputs*/, 
                      rhor, (p0>pr) ? rhor*1.01 : rhor*0.99/*initial guesses for Hugo. eq.*/,
                      rhor0, ur0/*outputs*/);
  double f0 = ul0 - ur0; // Eq. (18) in Kamm

  // 1.3: Initialize p1 ((21)(22) of Kamm) 
  double Clbar = (ul0 == ul) ? Cl : fabs(p0 - pl)/fabs(ul0 - ul);
  double Crbar = (ur0 == ur) ? Cr : fabs(p0 - pr)/fabs(ur0 - ur);
  p1 = (Crbar*pl + Clbar*pr + Clbar*Crbar*(ul - ur))/(Clbar + Crbar); 

  // 1.4 Calculate ul1, ur1 
  ComputeRhoUStarLeft(rhol, ul, pl, p1, idl/*inputs*/, 
                      rhol, rhol0/*initial guesses for Hugo. eq.*/,
                      rhol1, ul1/*outputs*/);
  ComputeRhoUStarRight(rhor, ur, pr, p1, idr/*inputs*/, 
                      rhor, rhor0/*initial guesses for Hugo. eq.*/,
                      rhor1, ur1/*outputs*/);
  double f1 = ul1 - ur1; // Eq. (18) in Kamm

  // -------------------------------
  // Step 2: Main Loop (Secant Method) 
  // -------------------------------
  double denom, f2;

  // monitor if the solution involves a transonic rarefaction. This is special as the solution 
  // at xi = x = 0 is within the rarefaction fan.
  bool trans_rare = false;
  double Vrare_x0[3]; //rho, u, and p at x = 0, in the case of a transonic rarefaction

  int iter = 0;
  double err_p = 1.0, err_u = 1.0;

  for(iter=0; iter<maxIts_main; iter++) {

    // 2.1: Update p using the Secant method
    denom = f1 - f0;
    if(denom == 0) {
      cout << "ERROR: Division-by-zero while using the secant method to solve the Riemann problem." << endl;
      exit_mpi();
    }
    p2 = p1 - f1*(p1-p0)/denom;

    // 2.2: Calculate ul2, ur2 
    ComputeRhoUStarLeft(rhol, ul, pl, p2, idl/*inputs*/, 
                        rhol0, rhol1/*initial guesses for Hugo. eq.*/,
                        rhol2, ul2/*outputs*/, 
                        &trans_rare, Vrare_x0/*filled only if found a trans. rarefaction*/);
    ComputeRhoUStarRight(rhor, ur, pr,  p2, idr/*inputs*/, 
                         rhor0, rhor1/*initial guesses for Hugo. erq.*/,
                         rhor2, ur2/*outputs*/,
                         &trans_rare, Vrare_x0/*filled only if found a trans. rarefaction*/);
    f2 = ul2 - ur2;
    
    // 2.3: Check stopping criterion
    err_p = fabs(p2 - p1)/max(fabs(p1), fabs(p2));
    err_u = fabs(f2)/max(fabs(ul2), fabs(ur2));
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
        double us = (rhor2*u2 - rhor*ur)/(rhor2 - rhor);
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
  if(u2>0) 
    for(int i=1; i<=2; i++) {
      k = (dir + i)%3 + 1;
      V[k] = Vm[k];
    }
  else if(u2<0)
    for(int i=1; i<=2; i++) {
      k = (dir + i)%3 + 1;
      V[k] = Vp[k];
    }
  else //u2 == 0
    for(int i=1; i<=2; i++) {
      k = (dir + i)%3 + 1;
      V[k] = 0.5*(Vm[k]+Vp[k]);
    }

}

//----------------------------------------------------------------------------------

//! Connect the left initial state with the left star state (the 1-wave)
virtual void
ExactRiemannSolverBase::ComputeRhoUStarLeft(double rhol, double ul, double pl, double ps, int idl/*inputs*/,
                                            double rhols0, double rhols1/*initial guesses for Hugo. eq.*/,
                                            double &rhols, double uls/*outputs*/, 
                                            bool *trans_rare, double *Vrare_x0/*filled only if found tran rf*/)
{
  if(pl >= ps) {//rarefaction --- numerical integration

    double dp_max = (pl-ps)/numSteps_rarefaction;

    // initialize drho
    double drho = rhol/(numSteps_rarefaction*2); //initial step size

    // The "basic" 4th-order Runge-Kutta method for integration
    double p1, p2, p3, p4, c1, c2, c3, c4, u1, u2, u3, u4;

    double ps_0 = pl, ps_1 = pl;

    while(ps_1 - ps > tol_rarefaction) {

      ps_1 = RarefactionIntegration_OneStepRK4(rhos_0, us_0, ps_0, drho); //XXX I AM HERE






    }
  }



}

//----------------------------------------------------------------------------------

//! Connect the right initial state with the right star state (the 3-wave)
virtual void
ExactRiemannSolverBase::ComputeRhoUStarRight(double rhor, double ur, double pr, double ps, int idr/*inputs*/,
                                             double rhors0, double rhors1/*initial guesses for Hugo. eq.*/,
                                             double &rhors, double urs/*outputs*/,
                                             bool *trans_rare, double *Vrare_x0/*filled only if found tran rf*/)
{




}

//----------------------------------------------------------------------------------

//----------------------------------------------------------------------------------
