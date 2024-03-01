/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <HeatDiffusionOperator.h>
#include <EmbeddedBoundaryDataSet.h>
#include <Vector5D.h>
#include <Utils.h>
using std::unique_ptr;

extern int INACTIVE_MATERIAL_ID;

//--------------------------------------------------------------------------

HeatDiffusionOperator::HeatDiffusionOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, MeshData &iod_mesh_, 
                                             EquationsData &iod_eqs_, vector<VarFcnBase*>& varFcn_,
                                             SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_, SpaceVariable3D &volume_,
                                             InterpolatorBase &interpolator, GradientCalculatorBase &grad_)
                     : iod_mesh(iod_mesh_), iod_eqs(iod_eqs_), coordinates(coordinates_), 
                       delta_xyz(delta_xyz_), volume(volume_),
                       varFcn(varFcn_), interpolator(interpolator), grad(grad_),
                       T(comm_, &(dm_all_.ghosted1_1dof)),
                       dTdx_i_minus_half(comm_, &(dm_all_.ghosted1_1dof)),
                       dTdy_j_minus_half(comm_, &(dm_all_.ghosted1_1dof)),
                       dTdz_k_minus_half(comm_, &(dm_all_.ghosted1_1dof)),
                       dTdr(comm_, &(dm_all_.ghosted1_1dof))
{

  // Get i0, j0, etc.
  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);

  // Set internal variables to 0 (including ghost layer)
  T.SetConstantValue(0.0, true);
  dTdx_i_minus_half.SetConstantValue(0.0, true);
  dTdy_j_minus_half.SetConstantValue(0.0, true);
  dTdz_k_minus_half.SetConstantValue(0.0, true);
  dTdr.SetConstantValue(0.0,true);

  // Creat a HeatDiffusionFcn for each material.
  // For problems involving multiple materials, it could happen that the user specified heat diffusion model
  // for only some of the materials, but not all. In this case, we still create a dummy (i.e. base)
  // HeatDiffusionFcn for those materials without heat diffusion. If heat diffusion is not specified for
  // any material, this constructor should not be called in the first place.
  //
  heatdiffFcn.resize(varFcn.size(), NULL);
  for(auto it = iod_eqs.materials.dataMap.begin(); it != iod_eqs.materials.dataMap.end(); it++){
    int matid = it->first;
    if(matid < 0 || matid >= (int)heatdiffFcn.size()) {
      print_error("*** Error: Detected error in the specification of material indices (id = %d).\n", matid);
    }
    switch (it->second->heat_diffusion.type) {
      case HeatDiffusionModelData::NONE :
        heatdiffFcn[matid] = new HeatDiffusionFcnBase();
        break;
      case HeatDiffusionModelData::CONSTANT :
        heatdiffFcn[matid] = new HeatDiffusionFcnConstant(it->second->heat_diffusion);
    }
  }
  for(auto&& diff : heatdiffFcn)
    if(!diff)
     diff = new HeatDiffusionFcnBase();


}

//--------------------------------------------------------------------------

HeatDiffusionOperator::~HeatDiffusionOperator()
{
  for(int i =0; i<(int)heatdiffFcn.size();i++)
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
  dTdr.Destroy();
}

//--------------------------------------------------------------------------
// Add diffusion fluxes on the left hand side of the N-S equations
void
HeatDiffusionOperator::AddDiffusionFluxes(SpaceVariable3D &V, SpaceVariable3D &ID, 
                                          [[maybe_unused]] vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS,
                                          SpaceVariable3D &R)
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
  double flux(0.0);
  int neighid = 0;
  double myk = 0.0; //conductivity
  double neighk = 0.0; //neighbor's conductivity
  double denom = 0.0;

  for(int k=k0; k<kkmax; k++)
    for(int j=j0; j<jjmax; j++)
      for(int i=i0; i<iimax; i++) {

        myid = (int)id[k][j][i];

        if(myid==INACTIVE_MATERIAL_ID)
          continue; 

        dx   = dxyz[k][j][i][0];
        dy   = dxyz[k][j][i][1];
        dz   = dxyz[k][j][i][2];

        myk = heatdiffFcn[myid]->GetConductivity();


        //*****************************************
        // calculate flux function F_{i-1/2,j,k}
        //*****************************************
        if(k!=kkmax-1 && j!=jjmax-1) {
          neighid = id[k][j][i-1];
          neighk = heatdiffFcn[neighid]->GetConductivity();
          denom = myk + neighk;
          flux = (denom == 0) ? 0.0 : 2.0*myk*neighk/denom*dTdx_i[k][j][i]; //If both are inactive, flux is 0
          flux *= dy*dz;
          res[k][j][i][4]   += flux;
          res[k][j][i-1][4] -= flux;
        }
        
        //*****************************************
        //calculate flux function G_{i,j-1/2,k}
        //*****************************************
        if(i!=iimax-1 && k!=kkmax-1) {
          neighid = id[k][j-1][i];
          neighk = heatdiffFcn[neighid]->GetConductivity();
          flux = (denom == 0) ? 0.0 : 2.0*myk*neighk/denom*dTdy_j[k][j][i];
          flux *= dx*dz;
          res[k][j][i][4]   += flux;
          res[k][j-1][i][4] -= flux;
        }

        //*****************************************
        //calculate flux function H_{i,j,k-1/2}
        //*****************************************   
        if(i!=iimax-1 && j!=jjmax-1) {
          neighid = id[k-1][j][i];
          neighk = heatdiffFcn[neighid]->GetConductivity();
          flux = (denom == 0) ? 0.0 : 2.0*myk*neighk/denom*dTdz_k[k][j][i];
          flux *= dx*dy;
          res[k][j][i][4]   += flux;
          res[k-1][j][i][4] -= flux;
        }
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
// The symmetry terms are placed on the left-hand-side of the N-S equations,
// and *added* to residual R (which is assumed to be on the left-hand-side)
void
HeatDiffusionOperator::AddSymmetryDiffusionTerms(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &R)
{
  if(iod_mesh.type == MeshData::SPHERICAL)
    AddSphericalSymmetryDiffusionTerms(V, ID, R);
  else if(iod_mesh.type == MeshData::CYLINDRICAL)
    AddCylindricalSymmetryDiffusionTerms(V, ID, R);
}

//--------------------------------------------------------------------------

void
HeatDiffusionOperator::AddSphericalSymmetryDiffusionTerms(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &R)
{
  //Get mesh data
  double***    vol = (double***)volume.GetDataPointer();
  Vec3D***  coords = (Vec3D***)coordinates.GetDataPointer();

  //Get data
  Vec5D***  v  = (Vec5D***) V.GetDataPointer();
  double*** id = (double***) ID.GetDataPointer();
  Vec5D***   r = (Vec5D***) R.GetDataPointer();
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

  //Get \partial(T)/partial(r)
  std::vector<int> Ti0{0};
  grad.CalculateFirstDerivativeAtNodes(0, T, Ti0, dTdr, Ti0);

  //Loop through the interior of the subdomain
  double*** dTdx = (double***)dTdr.GetDataPointer();
  double radial; //radial coord.
  double coeff;
  double myk = 0.0;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        radial = coords[k][j][i][0];
        assert(radial>0);
        coeff = vol[k][j][i]*2.0/radial;
 
        myid = id[k][j][i];
        myk = heatdiffFcn[myid]->GetConductivity();

        r[k][j][i][4] -= coeff*myk*dTdx[k][j][i]; // vol*(2*k*partialT)/(r*partialr)
      }

  //Restore data
  R.RestoreDataPointerToLocalVector(); //although data has been udpated, no need to communicate.
                                       //one subdomain does not need the residual info of another
                                       //subdomain
  volume.RestoreDataPointerToLocalVector();
  coordinates.RestoreDataPointerToLocalVector();
  V.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();
  dTdr.RestoreDataPointerToLocalVector();

}

//--------------------------------------------------------------------------

void
HeatDiffusionOperator::AddCylindricalSymmetryDiffusionTerms(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &R)
{

  //Get data
  Vec5D***  v  = (Vec5D***) V.GetDataPointer();
  double*** id = (double***) ID.GetDataPointer();
  Vec5D***   r = (Vec5D***) R.GetDataPointer();
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

  //Get \partial(T)/partial(r)
  std::vector<int> Ti0{0};
  grad.CalculateFirstDerivativeAtNodes(1, T, Ti0, dTdr, Ti0);

  //Get mesh data
  double***    vol = (double***)volume.GetDataPointer();
  Vec3D***  coords = (Vec3D***)coordinates.GetDataPointer();

  //Loop through the interior of the subdomain
  double*** dTdy = (double***)dTdr.GetDataPointer();
  double radial; //radial coord.
  double coeff;
  double myk = 0.0;
  
  //fprintf(stderr,"i0-imax: %d->%d, j0-jmax: %d->%d, k0-kmax: %d->%d", i0, imax, j0, jmax, k0, kmax);
  //fprintf(stderr,"ii0-iimax: %d->%d, jj0-jjmax: %d->%d, kk0-kkmax: %d->%d", ii0, iimax, jj0, jjmax, kk0, kkmax);

  
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        
        radial = coords[k][j][i][1];
        assert(radial>0);
        
        coeff = vol[k][j][i]/radial;

        //fprintf(stderr,"Check at (%d,%d,%d) coeff = %e\n", i,j,k,coeff);

        myid = id[k][j][i];
        myk = heatdiffFcn[myid]->GetConductivity();

        r[k][j][i][4] -= coeff*myk*dTdy[k][j][i]; // vol*(2*k*partialT)/(r*partialr)
        //fprintf(stderr,"Check at (%d,%d,%d) R[4] = %e\n", i,j,k,r[k][j][i][4]);
        

//        if((imax-i == 1 && jmax-jmax==1)&& (kmax-k == 1))
//          fprintf(stderr,"Get the end point of subdomain!\n");
      }

  //Restore data
  R.RestoreDataPointerToLocalVector(); //although data has been udpated, no need to communicate.
                                       //one subdomain does not need the residual info of another
                                       //subdomain
  volume.RestoreDataPointerToLocalVector();
  coordinates.RestoreDataPointerToLocalVector();
  V.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();
  dTdr.RestoreDataPointerToLocalVector();
}
