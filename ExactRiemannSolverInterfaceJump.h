#ifndef _EXACT_RIEMANN_SOLVER_INTERFACE_JUMP_H_
#define _EXACT_RIEMANN_SOLVER_INTERFACE_JUMP_H_

#include <ExactRiemannSolverBase.h>


/****************************************************************************************
 * A derived class that include the interface pressure jump due to surface tension
 * **************************************************************************************/

class ExactRiemannSolverInterfaceJump: public ExactRiemannSolverBase {

  protected:
    double delta_p; // pressure jump across the contact discontinuity
    double surface_tension_coefficient; // In general, it depends on the material property on both sides of the interface and temperature;
                                        // currently, we assume it to be a constant specified from the input file
    double interface_curvature;
                                        
  public:
    ExactRiemannSolverInterfaceJump(std::vector<VarFcnBase*> &vf_, ExactRiemannSolverData &iod_riemann_) : ExactRiemannSolverBase(vf_, iod_riemann_) {
       surface_tension_coefficient = iod_riemann.surface_tension_coefficient;
       interface_curvature = 0.5; // hard-coded for specific case at moment, will be implemented in a better way 
       delta_p = surface_tension_coefficient * interface_curvature; 
    };

    int ComputeRiemannSolution(double *dir /*unit normal*/, double *Vm, int idm /*"left" state*/, 
	double *Vp, int idp /*"right" state*/, 
	double *Vs, int &id /*solution at xi = 0 (i.e. x=0) */,
	double *Vsm /*left 'star' solution*/,
	double *Vsp /*right 'star' solution*/);

    bool FindInitialInterval(double rhol, double ul, double pl, double el, double cl, int idl,
	double rhor, double ur, double pr, double er, double cr, int idr, /*inputs*/
	double &p0, double &rhol0, double &rhor0, double &ul0, double &ur0,
	double &p1, double &rhol1, double &rhor1, double &ul1, double &ur1/*outputs*/);

    bool FindInitialFeasiblePoints(double rhol, double ul, double pl, double el, double cl, int idl,
	double rhor, double ur, double pr, double er, double cr, int idr, /*inputs*/
	double &p0, double &rhol0, double &rhor0, double &ul0, double &ur0, 
	double &p1, double &rhol1, double &rhor1, double &ul1, double &ur1/*outputs*/);

    int FindInitialFeasiblePointsByAcousticTheory(double rhol, double ul, double pl,
	[[maybe_unused]] double el, double cl, int idl,
	double rhor, double ur, double pr, [[maybe_unused]] double er, double cr, int idr, /*inputs*/
	double &p0, double &rhol0, double &rhor0, double &ul0, double &ur0, 
	double &p1, double &rhol1, double &rhor1, double &ul1, double &ur1/*outputs*/);

    void PrintStarRelations(double rhol, double ul, double pl, int idl,
	double rhor, double ur, double pr, int idr,
	double pmin, double pmax, double dp);

    void FinalizeSolution(double *dir, double *Vm, double *Vp,
	double rhol, double ul, double pl, int idl, 
	double rhor, double ur, double pr, int idr, 
	double rhol2, double rhor2, double u2, double p2,
	bool trans_rare, double Vrare_x0[3], /*inputs*/
	double *Vs, int &id, double *Vsm, double *Vsp /*outputs*/);


};  

#endif

