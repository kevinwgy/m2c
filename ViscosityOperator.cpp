#include <ViscosityOperator.h>
#include <Utils.h>
#include <bits/stdc++.h>

//--------------------------------------------------------------------------

ViscosityOperator::ViscosityOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, 
                                     ViscosityModelData &iod_viscosity_,
                                     SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_)
                  : iod_viscosity(iod_viscosity_), coordinates(coordinates_), delta_xyz(delta_xyz_),
                    G(comm_, &(dm_all_.ghosted1_5dof)),
                    V_i_minus_half(comm_, &(dm_all_.ghosted1_3dof)),
                    V_j_minus_half(comm_, &(dm_all_.ghosted1_3dof)),
                    V_k_minus_half(comm_, &(dm_all_.ghosted1_3dof)),
                    DU(comm_, &(dm_all_.ghosted1_3dof)),
                    DV(comm_, &(dm_all_.ghosted1_3dof)),
                    DW(comm_, &(dm_all_.ghosted1_3dof)),
                    Var(comm_, &(dm_all_.ghosted1_1dof)),
                    interpolator(NULL)
{
  //set internal variables to 0 (including ghost layer)
  G.SetConstantValue(0.0, true);
  V_i_minus_half.SetConstantValue(0.0, true);
  V_j_minus_half.SetConstantValue(0.0, true);
  V_k_minus_half.SetConstantValue(0.0, true);
  DU.SetConstantValue(0.0, true);
  DV.SetConstantValue(0.0, true);
  DW.SetConstantValue(0.0, true);
  Var.SetConstantValue(0.0, true);
}

//--------------------------------------------------------------------------

ViscosityOperator::~ViscosityOperator()
{

}

//--------------------------------------------------------------------------

void ViscosityOperator::AddDiffusionFluxes(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &R)
{
  //1. Calculate the x, y, and z velocities at cell interfaces by interpolation 
  assert(interpolator!=NULL);
  interpolator->InterpolateAtCellInterfaces(0/*x-dir*/, V, std::vector<int>{1,2,3}, 
                                            V_i_minus_half, std::vector<int>{0,1,2});
  interpolator->InterpolateAtCellInterfaces(1/*y-dir*/, V, std::vector<int>{1,2,3}, 
                                            V_j_minus_half, std::vector<int>{0,1,2});
  interpolator->InterpolateAtCellInterfaces(2/*z-dir*/, V, std::vector<int>{1,2,3}, 
                                            V_k_minus_half, std::vector<int>{0,1,2});

  //2. Calculate DU, DV, DW at cell interfaces;
  Var.AXPlusBY(0.0, 1.0, V, std::vector<int>{0}/*Var index*/, std::vector<int>{1}/*V index*/, true/*including ghost*/);
  differentiator.CalculateFirstDerivativesAtCellInterfaces(Var, DU);

  Var.AXPlusBY(0.0, 1.0, V, std::vector<int>{0}/*Var index*/, std::vector<int>{2}/*V index*/, true/*including ghost*/);
  differentiator.CalculateFirstDerivativesAtCellInterfaces(Var, DV);

  Var.AXPlusBY(0.0, 1.0, V, std::vector<int>{0}/*Var index*/, std::vector<int>{3}/*V index*/, true/*including ghost*/);
  differentiator.CalculateFirstDerivativesAtCellInterfaces(Var, DW);

  //3. Loop through cells (in the physical domain) and calculate viscous fluxes
  XXX 
}

//--------------------------------------------------------------------------






