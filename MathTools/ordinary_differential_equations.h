/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#pragma once
#include<cmath> //std::isfinite
#include<vector>
#include<algorithm> //std:min, std::max

namespace MathTools {

/****************************************************************************
 * Solve first-order ODEs of the form dU/dt = F(U,t) using Runge-Kutta methods.
 * The solvers implemented here allow the user to pass in a "state checker" 
 * that checks whether the intermediate states are "meaningful" for the user's 
 * application (e.g., density cannot be negative). This checker can even make 
 * changes to the intermediate states.
 ***************************************************************************/

/****************************************************************************
 * Solve first-order ODEs of the form dU/dt = F(U,t) using 4-th order Runge-Kutta
 *
 * Reference: This is the most common 4-th order Runge-Kutta method.
 *
 * Note: 
 *   - "fun" must be of type "void fun(double*, double, double*)." It should implement
 *     F(U,t), with the first and second input argments being "U" (pointer) and "t". 
 *     The third input argment is supposed to be the pointer to the output "F(U,t)"
 *   - "check_state" must be of type "bool check_state(double*)." It should take the pointer
 *     to the state variable (U), and return true or false. True means the state variable
 *     FAILED THE CHECK. False means it has passed the check. If needed, this function can 
 *     also change U, since a pointer is passed to it.
 *
 * size  : size of U
 * t0, U0: Initial condition, i.e. U(t0) = U0
 * dt    : step size
 * tf, U : final time (t, input) and the solution (U, output). 
 *         If tf is smaller than t0, dt must be negative.
 * Uall  : (optional output) pointer to a 1D double vector that stores the entire trajectory.
 * last_dt : (optional output) pointer to a double value that stores the last time step
 *
 * Returns 0 if successful; -1 if there is input error; 1 if the solution process produced error
 ***************************************************************************/

template<typename Functor1, typename Functor2>
int runge_kutta_4(Functor1 fun, int size, double t0, double* U0, double dt, double tf,
                  double* U, //U: computed solution at tf
                  Functor2 check_state, //state checker (If not needed, set it to "[](double*){return false;}"
                  std::vector<double>* Uall = NULL, // U at all time steps (optional output)
                  double* last_dt = NULL) // last time step (may not be dt!)
{
  if(tf==t0) { //trivial
    for(int j=0; j<size; j++)
      U[j] = U0[j];
    if(Uall) {
      Uall->clear();
      for(int j=0; j<size; j++)
        Uall->push_back(U[j]);
    }
    if(last_dt) *last_dt = 0.0;
    return 0;
  }

  if(dt==0.0)
    return -1;

  int N = floor((tf-t0)/dt);
  if(!std::isfinite(N) || N<0)
    return -1;

  if(Uall) {
    Uall->clear();
    Uall->reserve(size*(N+2));
  }

  double t = t0;
  for(int j=0; j<size; j++) {
    U[j] = U0[j];
    if(!std::isfinite(U[j])) return -1;
    if(Uall) Uall->push_back(U[j]);
  }
  if(check_state(U)) return -1;


  double* Utmp = new double[size];
  double* k1 = new double[size];
  double* k2 = new double[size];
  double* k3 = new double[size];
  double* k4 = new double[size];
  
  double one_sixth = 1.0/6.0;
  double one_third = 1.0/3.0;
  double half_dt   = 0.5*dt;

  bool finished = false;
  while(!finished) {
    for(int i=0; i<N; i++) {

      //1. Calculate k1
      fun(U, t, k1); 
    
      //2. Calculate k2
      for(int j=0; j<size; j++) {
        Utmp[j] = U[j] + half_dt*k1[j];
        if(!std::isfinite(Utmp[j]))
          goto END_OF_LOOP_RK4;
      }
      if(check_state(Utmp)) 
        goto END_OF_LOOP_RK4;

      fun(Utmp, t+half_dt, k2);


      //3. Calculate k3
      for(int j=0; j<size; j++) {
        Utmp[j] = U[j] + half_dt*k2[j]; 
        if(!std::isfinite(Utmp[j]))
          goto END_OF_LOOP_RK4;
      }
      if(check_state(Utmp)) 
        goto END_OF_LOOP_RK4;

      fun(Utmp, t+half_dt, k3);


      //4. Calculate k4
      for(int j=0; j<size; j++) {
        Utmp[j] = U[j] + dt*k3[j]; 
        if(!std::isfinite(Utmp[j]))
          goto END_OF_LOOP_RK4;
      }
      if(check_state(Utmp)) 
        goto END_OF_LOOP_RK4;

      fun(Utmp, t+dt, k4);


      for(int j=0; j<size; j++) {
        U[j] += dt*(one_sixth*(k1[j]+k4[j]) + one_third*(k2[j]+k3[j]));
        if(!std::isfinite(U[j]))
          goto END_OF_LOOP_RK4;
        if(Uall) Uall->push_back(U[j]);
      }
      if(check_state(U)) 
        goto END_OF_LOOP_RK4;

      t += dt;
    }

    if((dt>=0 && t-tf>=-1e-9*dt) ||
       (dt<0  && t-tf<=-1e-9*dt)) { 
      finished = true;
      if(last_dt) *last_dt = dt;
    } else { //run one more step.
      N = 1;
      dt = tf-t;
    }
  }

END_OF_LOOP_RK4:

  delete[] k1; delete[] k2; delete[] k3; delete[] k4; delete[] Utmp;

  if(!finished) //either got infinite number or failed "check_state"
    return 1;

  return 0;  
}

//------------------------------------------------------------------------

/****************************************************************************
 * Solve first-order ODEs of the form dU/dt = F(U,t) using the adative 4-th / 5-th
 * order Runge-Kutta-Fehlberg method, with coefficients given by Cash and Karp, 1990
 *
 * Note: 
 *   - "fun" must be of type "void fun(double*, double, double*)." It should implement
 *     F(U,t), with the first and second input argments being "U" (Pointer) and "t". 
 *     The third input argment is supposed to be the pointer to the output "F(U,t)"
 *   - "check_state" must be of type "bool check_state(double*)." It should take the pointer
 *     to the state variable (U), and return true or false. True means the state variable
 *     FAILED THE CHECK. False means it has passed the check. If needed, this function can 
 *     also change U, since a pointer is passed to it.
 *
 * size  : size of U
 * t0, U0: Initial condition, i.e. U(t0) = U0
 * dt0   : INITIAL step size
 * tol   : Error tolerance for the relative error in each component of U
 * tf, U : final time (t, input) and the solution (U, output). 
 *         If tf is smaller than t0, dt0 must be negative.
 * (tall,Uall) : (optional output) pointer to a 1D double vector that stores the entire trajectory.
 *
 * Returns 0 if successful; -1 if there is input error; 1 if the solution process produced error,
 *         2 if max. iteration number is reached
 ***************************************************************************/


template<typename Functor1, typename Functor2>
int runge_kutta_45(Functor1 fun, int size, double t0, double* U0, double dt0, double tf,
                   double* U, //U: computed solution at tf
                   Functor2 check_state, //state checker (If not needed, set it to "[](double*){return false;}"
                   double tol = 1.0e-8, int Nmax = 1e7, //Max. number of time steps
                   std::vector<double>* tall = NULL, // t at all time steps (optional output)
                   std::vector<double>* Uall = NULL) // U at all time steps (optional output)
{

  if(tf==t0) { //trivial
    for(int j=0; j<size; j++) 
      U[j] = U0[j];
    if(tall && Uall) {
      tall->assign(1,t0);
      Uall->clear();
      for(int j=0; j<size; j++)
        Uall->push_back(U[j]);
    }
    return 0;
  }

  if(dt0==0.0 || tol<=0.0)
    return -1;

  int N = floor((tf-t0)/dt0);
  if(!std::isfinite(N) || N<=0.0)
    return -1;

  if(tall && Uall) {
    tall->clear();
    tall->reserve(std::min((int)1e5,N+2));
    Uall->clear();
    Uall->reserve(size*std::min((int)1e5,N+2));
  }


  // Coefficients
  double c2 = 0.2, c3 = 0.3, c4 = 0.6, c5 = 1.0, c6 = 0.875;
  double p2_1 = 0.2, p3_1 = 3.0/40.0, p3_2 = 9.0/40.0, p4_1 = 0.3, p4_2 = -0.9, p4_3 = 1.2;
  double p5_1 = -11.0/54.0, p5_2 = 2.5, p5_3 = -70.0/27.0, p5_4 = 35.0/27.0;
  double p6_1 = 1631.0/55296.0, p6_2 = 175.0/512.0, p6_3 = 575.0/13824.0, p6_4 = 44275.0/110592.0, p6_5 = 253.0/4096.0;
  double o4_1 = 2825.0/27648.0, o4_3 = 18575.0/48384.0, o4_4 = 13525.0/55296.0, o4_5 = 277.0/14336.0, o4_6 = 0.25;
  double o5_1 = 37.0/378.0, o5_3 = 250.0/621.0, o5_4 = 125.0/594.0, o5_6 = 512.0/1771.0;


  // Set initial time and solution
  double t = t0;
  for(int j=0; j<size; j++) {
    U[j] = U0[j];
    if(!std::isfinite(U[j])) return -1;
  }
  if(tall && Uall) {
    tall->push_back(t);
    for(int j=0; j<size; j++)
      if(Uall) Uall->push_back(U[j]);
  }
  if(check_state(U)) return -1;


  double* Utmp = new double[size];
  double* Utmp2 = new double[size];
  double* k1 = new double[size];
  double* k2 = new double[size];
  double* k3 = new double[size];
  double* k4 = new double[size];
  double* k5 = new double[size];
  double* k6 = new double[size];

  double dt = dt0; //initial step size
  if(t+dt>tf) dt = tf-t;
  double safety = 0.9, alpha_1 = -0.2, alpha_2 = -0.25, err0 = 1.0e-30; //params in adaptation
  double max_increase = 10.0, max_decrease = 0.2;
  double ratio, rel_error;

  int i;
  for(i=0; i<Nmax; i++) {

    //1. Calculate k1
    fun(U, t, k1); 
    

    //2. Calculate k2
    for(int j=0; j<size; j++) {
      Utmp[j] = U[j] + dt*p2_1*k1[j];
      if(!std::isfinite(Utmp[j]))
        goto END_OF_LOOP_RKF45;
    }
    if(check_state(Utmp))
      goto END_OF_LOOP_RKF45;

    fun(Utmp, t+c2*dt, k2);


    //3. Calculate k3
    for(int j=0; j<size; j++) {
      Utmp[j] = U[j] + dt*(p3_1*k1[j] + p3_2*k2[j]);
      if(!std::isfinite(Utmp[j]))
        goto END_OF_LOOP_RKF45;
    }
    if(check_state(Utmp))
      goto END_OF_LOOP_RKF45;

    fun(Utmp, t+c3*dt, k3);


    //4. Calculate k4
    for(int j=0; j<size; j++) {
      Utmp[j] = U[j] + dt*(p4_1*k1[j] + p4_2*k2[j] + p4_3*k3[j]);
      if(!std::isfinite(Utmp[j]))
        goto END_OF_LOOP_RKF45;
    }
    if(check_state(Utmp))
      goto END_OF_LOOP_RKF45;

    fun(Utmp, t+c4*dt, k4);


    //5. Calculate k5
    for(int j=0; j<size; j++) {
      Utmp[j] = U[j] + dt*(p5_1*k1[j] + p5_2*k2[j] + p5_3*k3[j] + p5_4*k4[j]);
      if(!std::isfinite(Utmp[j]))
        goto END_OF_LOOP_RKF45;
    }
    if(check_state(Utmp))
      goto END_OF_LOOP_RKF45;

    fun(Utmp, t+c5*dt, k5);


    //6. Calculate k6
    for(int j=0; j<size; j++) {
      Utmp[j] = U[j] + dt*(p6_1*k1[j] + p6_2*k2[j] + p6_3*k3[j] + p6_4*k4[j] + p6_5*k5[j]);
      if(!std::isfinite(Utmp[j]))
        goto END_OF_LOOP_RKF45;
    }
    if(check_state(Utmp))
      goto END_OF_LOOP_RKF45;

    fun(Utmp, t+c6*dt, k6);

    //fprintf(stderr,"%e : %e %e %e %e %e %e\n", t, k1[0], k2[0], k3[0], k4[0], k5[0], k6[0]);

    //4th and 5th order approximations
    for(int j=0; j<size; j++) {
      Utmp[j]  = U[j] + dt*(o4_1*k1[j] + o4_3*k3[j] + o4_4*k4[j] + o4_5*k5[j] + o4_6*k6[j]); //4th order
      Utmp2[j] = U[j] + dt*(o5_1*k1[j] + o5_3*k3[j] + o5_4*k4[j] + o5_6*k6[j]); //5th order
      if(!std::isfinite(Utmp[j]) || !std::isfinite(Utmp2[j]))
        goto END_OF_LOOP_RKF45;
    }
    if(check_state(Utmp) || check_state(Utmp2))
      goto END_OF_LOOP_RKF45;


    //Check relative error to determine (1) whether we should repeat the step and (2) the new dt
    ratio = max_increase;
    for(int j=0; j<size; j++) {
      rel_error = std::max(fabs(Utmp[j]), fabs(Utmp2[j]));
      if(rel_error!=0.0)
        rel_error = fabs((Utmp[j] - Utmp2[j])/rel_error + err0)/tol;
      else
        rel_error = err0/tol;

      if(rel_error>1.0) // Should decrease dt
        ratio = std::min(ratio, safety*pow(rel_error, alpha_2));
      else // this component (j) wants to increase dt
        ratio = std::min(ratio, safety*pow(rel_error, alpha_1));
    }
    ratio = std::max(ratio, max_decrease); //cannot decrease by more than max_decrease


    if(ratio>=1.0) { //continue to the next time step
      t  += dt;
      for(int j=0; j<size; j++)
        U[j] = Utmp2[j];
      if(tall && Uall) {
        for(int j=0; j<size; j++)
          Uall->push_back(U[j]);
        tall->push_back(t);
      }

      if((dt>=0 && t-tf>=-1e-9*dt) || //marching forward
         (dt<0  && t-tf<=-1e-9*dt)) {  //marching backward
        delete[] k1; delete[] k2; delete[] k3; delete[] k4; delete[] k5; delete[] k6; 
        delete[] Utmp; delete[] Utmp2;
        return 0; //DONE! 
      }

      dt *= ratio; 
      if((dt>=0 && t+dt-tf>=1e-9*dt) ||//next step should be at tf exactly
         (dt<0  && t+dt-tf<=1e-9*dt)) 
        dt = tf-t;
    }
    // otherwise, repeat the same step with a smaller dt

    dt *= ratio; 

  }


END_OF_LOOP_RKF45 :

  delete[] k1; delete[] k2; delete[] k3; delete[] k4; delete[] k5; delete[] k6; 
  delete[] Utmp; delete[] Utmp2;

  if(i>=Nmax)
    return 2; //it has exhaused Nmax iterations w/o reaching tf.

  return 1; //it must have failed an infinite check or a user-specified state check.

}

//------------------------------------------------------------------------






//------------------------------------------------------------------------

} //end of namespace
