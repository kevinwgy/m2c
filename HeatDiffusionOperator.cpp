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

  // Creat a HeatDiffusionFcn for each material.
  // For problems involving multiple materials, it could happen that the user specified heat diffusion model
  // for only some of the materials, but not all. In this case, we still create a dummy (i.e. base)
  // HeatDiffusionFcn for those materials without heat diffusion. If all the materials are not specified heat
  // diffusion, this constructor should not be called in the first place.
  //
  for(auto it = iod_eqs.materials.dataMap.begin(); it != iod_eqs.materials.dataMap.end(); it++)
    heatdiffFcn.push_back(NULL);//allocate space
  for(auto it = iod_eqs.materials.dataMap.begin(); it != iod_eqs.materials.dataMap.end(); it++){
    int matid = it->first;
    if(matid < 0 || matid >= heatdiffFcn.size()) {
      print_error("*** Error: Detected error in the specification of material indices (id = %d).\n", matid);
    }
    switch (it->second->heat_diffusion.type) {
      case HeatDiffusionModelData::NONE :
        heatdiffFcn[matid] = new HeatDiffuFcnBase(*varFcn[matid]);
        break;
      case HeatDiffusionModelData::CONSTANT :
        heatdiffFcn[matid] = new HeatDiffuFcnConstant(it->second->heat_diffusion,*varFcn[matid]);
    }
  }
}

//--------------------------------------------------------------------------

HeatDiffusionOperator::~HeatDiffusionOperator()
{
  for(int i =0; i<heatdiffFcn.size();i++)
    delete heatdiffFcn[i];
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

  //1. Calculate the temperature at node
  Vec5D*** v    = (Vec5D***)V.GetDataPointer();
  double*** id  = (double***)ID.GetDataPointer();
  double*** Te  = (double***)T.GetDataPointer();

  int myid = 0;
  double e;
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++) {

        myid = id[k][j][i]; 
        e = varFcn[myid]->GetInternalEnergyPerUnitMass(v[k][j][i][0], v[k][j][i][4]);
        Te[k][j][i] = varFcn[myid]->GetTemperature(v[k][j][i][0], e);
  }
  
  T.RestoreDataPointerToLocalVector(); 
  
  //2. Calculate temperature derivatives at cell interfaces
  std::vector<int> Ti0{0};
  grad.CalculateFirstDerivativeAtCellInterfaces(0, 0, T, Ti0, dTdx_i_minus_half, Ti0);
  grad.CalculateFirstDerivativeAtCellInterfaces(1, 1, T, Ti0, dTdy_j_minus_half, Ti0);
  grad.CalculateFirstDerivativeAtCellInterfaces(2, 2, T, Ti0, dTdz_k_minus_half, Ti0);

  //3. loop through cell interfaces and calculate heat diffusion fluxes
  double*** dTdx_i = (double***)dTdx_i_minus_half.GetDataPointer();
  double*** dTdy_j = (double***)dTdy_j_minus_half.GetDataPointer();
  double*** dTdz_k = (double***)dTdz_k_minus_half.GetDataPointer();
  Vec3D*** dxyz = (Vec3D***)delta_xyz.GetDataPointer();
  Vec5D*** res  = (Vec5D***)R.GetDataPointer();

  double dx = 0.0, dy = 0.0, dz = 0.0;
  Vec5D flux;
  flux[0] = flux[1] = flux[2] = flux[3] = 0.0;
  int neighid = 0;
  double myk = 0.0;
  double neighk = 0.0;
  double denom = 0.0;
  for(int k=k0; k<kkmax; k++)
    for(int j=j0; j<jjmax; j++)
      for(int i=i0; i<iimax; i++) {

        myid = id[k][j][i];
        dx   = dxyz[k][j][i][0];
        dy   = dxyz[k][j][i][1];
        dz   = dxyz[k][j][i][2];

        //*****************************************
        // calculate flux function F_{i-1/2,j,k}
        //*****************************************
        if(k!=kkmax-1 && j!=jjmax-1) {
          neighid = id[k][j][i-1];
          myk = heatdiffFcn[myid]->conduct;
          neighk = heatdiffFcn[neighid]->conduct;
          denom = myk + neighk;
          flux[4] = (denom == 0) ? 0.0 : myk*neighk/denom*dTdx_i[k][j][i];
          flux *= dy*dz;
          res[k][j][i]   += flux;
          res[k][j][i-1] -= flux;
        }
        //Debug
        if(j == 300)
          fprintf(stderr,"- Node (%d, %d, %d), heat diffusivity is %e \n, diffusivity of its neighbor is %e, diffusion flux at i-1/2 is: %e.\n", i,j,k,myk,neighk,flux[4]);
        
        //*****************************************
        //calculate flux function G_{i,j-1/2,k}
        //*****************************************
        if(i!=iimax-1 && k!=kkmax-1) {
          neighid = id[k][j-1][i];
          myk = heatdiffFcn[myid]->conduct;
          neighk = heatdiffFcn[neighid]->conduct;
          flux[4] = (denom == 0) ? 0.0 : myk*neighk/denom*dTdy_j[k][j][i];
          flux *= dx*dz;
          res[k][j][i]   += flux;
          res[k][j-1][i] -= flux;
        }
        //Debug
        if(j == 300)
          fprintf(stderr,"diffusivity of its neighbor is %e, heat diffusion flux at j-1/2 is: %e.\n", neighk, flux[4]);


        //*****************************************
        //calculate flux function H_{i,j,k-1/2}
        //*****************************************   
        if(i!=iimax-1 && j!=jjmax-1) {
          neighid = id[k-1][j][i];
          myk = heatdiffFcn[myid]->conduct;
          neighk = heatdiffFcn[neighid]->conduct;
          flux[4] = (denom == 0) ? 0.0 : myk*neighk/denom*dTdz_k[k][j][i];
          flux *= dx*dy;
          res[k][j][i]   += flux;
          res[k-1][j][i] -= flux;
        }
        //debug
        if(j == 300)
          fprintf(stderr,"diffusivity of its neighbor is %e, heat diffusion flux at k-1/2 is: %e.\n", neighk, flux[4]);
  }

  dTdx_i_minus_half.RestoreDataPointerToLocalVector();
  dTdy_j_minus_half.RestoreDataPointerToLocalVector();
  dTdz_k_minus_half.RestoreDataPointerToLocalVector();
  delta_xyz.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();
  V.RestoreDataPointerToLocalVector();

  R.RestoreDataPointerToLocalVector(); //NOTE: although R has been updated, there is no need of
                                       //cross-subdomain communications. So, no need to
                                       //update the global vec.

}


//--------------------------------------------------------------------------



