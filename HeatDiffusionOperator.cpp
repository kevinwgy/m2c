#include <HeatDiffusionOperator.h>
#include <Vector5D.h>
#include <Utils.h>

//--------------------------------------------------------------------------

HeatDiffusionOperator::HeatDiffusionOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, 
                                             EquationsData &iod_eqs_, vector<VarFcnBase*>& varFcn_,
                                             SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_,
                                             InterpolatorBase &interpolator, GradientCalculatorBase &grad_)
                     : iod_eqs(iod_eqs_), coordinates(coordinates_), delta_xyz(delta_xyz_),
                       varFcn(varFcn_), interpolator(interpolator), grad(grad_),
                       T(comm_, &(dm_all_.ghosted1_1dof)),
                       dTdx_i_minus_half(comm_, &(dm_all_.ghosted1_1dof)),
                       dTdy_j_minus_half(comm_, &(dm_all_.ghosted1_1dof)),
                       dTdz_k_minus_half(comm_, &(dm_all_.ghosted1_1dof))
{

  // Get i0, j0, etc.
  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);

  // Set internal variables to 0 (including ghost layer)
  T.SetConstantValue(0.0, true);
  dTdx_i_minus_half.SetConstantValue(0.0, true);
  dTdy_j_minus_half.SetConstantValue(0.0, true);
  dTdz_k_minus_half.SetConstantValue(0.0, true);

  // TODO

}

//--------------------------------------------------------------------------

HeatDiffusionOperator::~HeatDiffusionOperator()
{
}

//--------------------------------------------------------------------------
// Destroy internal "SpaceVariables"
void
HeatDiffusionOperator::Destroy()
{
  T.Destroy();
  dTdx_i_minus_half.Destroy();
  dTdy_j_minus_half.Destroy();
  dTdz_k_minus_half.Destroy();
}

//--------------------------------------------------------------------------
// Add diffusion fluxes on the left hand side of the N-S equations
void
HeatDiffusionOperator::AddDiffusionFluxes(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &R)
{

  //TODO

}

//--------------------------------------------------------------------------



