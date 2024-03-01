/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <SymmetryOperator.h>
#include <Vector5D.h>
#include <cassert>

//--------------------------------------------------------------------------

SymmetryOperator::SymmetryOperator(MPI_Comm &comm_, [[maybe_unused]] DataManagers3D &dm_all_, MeshData &iod_mesh_,
                                   vector<VarFcnBase*> &varFcn_, SpaceVariable3D &coordinates_,
                                   SpaceVariable3D &delta_xyz_, SpaceVariable3D &volume_)
                 : comm(comm_), iod_mesh(iod_mesh_), varFcn(varFcn_), 
                   coordinates(coordinates_), delta_xyz(delta_xyz_), volume(volume_)
{
  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);
}

//--------------------------------------------------------------------------

SymmetryOperator::~SymmetryOperator()
{ }

//--------------------------------------------------------------------------

void
SymmetryOperator::Destroy()
{ }

//--------------------------------------------------------------------------
// The symmetry terms are placed on the left-hand-side of the N-S equations,
// and *added* to residual R (which is assumed to be on the left-hand-side)
void
SymmetryOperator::AddSymmetryTerms(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &R)
{
  if(iod_mesh.type == MeshData::SPHERICAL)
    AddSphericalSymmetryTerms(V, ID, R);
  else if(iod_mesh.type == MeshData::CYLINDRICAL)
    AddCylindricalSymmetryTerms(V, ID, R);
}

//--------------------------------------------------------------------------
// The symmetry terms are placed on the left-hand-side of the N-S equations,
// and *added* to residual R (which is assumed to be on the left-hand-side)
void 
SymmetryOperator::AddSphericalSymmetryTerms(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &R)
{
  //Get mesh data
  double***    vol = (double***)volume.GetDataPointer();
  Vec3D***  coords = (Vec3D***)coordinates.GetDataPointer();

  //Get data
  Vec5D***  v  = (Vec5D***) V.GetDataPointer();
  double*** id = (double***) ID.GetDataPointer();
  Vec5D***   r = (Vec5D***) R.GetDataPointer();

  //Loop through the interior of the subdomain
  double radial; //radial coord.
  double H; //total enthalpy per unit mass
  double coeff;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        radial = coords[k][j][i][0];
        assert(radial>0);
        coeff = vol[k][j][i]*2.0/radial;

        r[k][j][i][0] += coeff*v[k][j][i][0]*v[k][j][i][1]; // vol*(2*rho*u_r)/r
        r[k][j][i][1] += coeff*v[k][j][i][0]*v[k][j][i][1]*v[k][j][i][1]; // vol*(2*rho*u_r*u_r)/r

        H = varFcn[id[k][j][i]]->ComputeTotalEnthalpyPerUnitMass(v[k][j][i]);
        r[k][j][i][4] += coeff*v[k][j][i][0]*H*v[k][j][i][1]; // vol*(2*rho*H*u_r)/r

      }

  //Restore data
  R.RestoreDataPointerToLocalVector(); //although data has been udpated, no need to communicate.
                                       //one subdomain does not need the residual info of another
                                       //subdomain

  volume.RestoreDataPointerToLocalVector();
  coordinates.RestoreDataPointerToLocalVector();
  V.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();

}

//--------------------------------------------------------------------------
// The symmetry terms are placed on the left-hand-side of the N-S equations,
// and *added* to residual R (which is assumed to be on the left-hand-side)
void 
SymmetryOperator::AddCylindricalSymmetryTerms(SpaceVariable3D &V, SpaceVariable3D &ID, SpaceVariable3D &R)
{
  //Get mesh data
  double***    vol = (double***)volume.GetDataPointer();
  Vec3D***  coords = (Vec3D***)coordinates.GetDataPointer();

  //Get data
  Vec5D***  v  = (Vec5D***) V.GetDataPointer();
  double*** id = (double***) ID.GetDataPointer();
  Vec5D***   r = (Vec5D***) R.GetDataPointer();

  //Loop through the interior of the subdomain
  double radial; //radial coord.
  double ur; //radial velocity
  double w; //axial velocity
  double rho; //density
  double H; //total enthalpy per unit mass
  double coeff;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        ur  = v[k][j][i][2];
        w   = v[k][j][i][1];
        rho = v[k][j][i][0];
        H   = varFcn[id[k][j][i]]->ComputeTotalEnthalpyPerUnitMass(v[k][j][i]);

        radial = coords[k][j][i][1];
        assert(radial>0);

        coeff = vol[k][j][i]*rho*ur/radial;

        r[k][j][i][0] += coeff; 
        r[k][j][i][1] += coeff*w;
        r[k][j][i][2] += coeff*ur;
        r[k][j][i][4] += coeff*H;

      }

  //Restore data
  R.RestoreDataPointerToLocalVector(); //although data has been udpated, no need to communicate.
                                       //one subdomain does not need the residual info of another
                                       //subdomain

  volume.RestoreDataPointerToLocalVector();
  coordinates.RestoreDataPointerToLocalVector();
  V.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();

}

//--------------------------------------------------------------------------





