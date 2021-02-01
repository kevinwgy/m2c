#include<Reconstructor.h>
#include<Vector3D.h>
#include<Vector5D.h>
#include<Utils.h>
using std::round;


//--------------------------------------------------------------------------

Reconstructor::Reconstructor(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_, 
                             SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_)
                   : iod(iod_), delta_xyz(delta_xyz_),
                     CoeffA(comm_, &(dm_all_.ghosted1_3dof)), 
                     CoeffB(comm_, &(dm_all_.ghosted1_3dof)), 
                     CoeffK(comm_, &(dm_all_.ghosted1_3dof))
{

}

//--------------------------------------------------------------------------

Reconstructor::~Reconstructor()
{

}

//--------------------------------------------------------------------------
/** Compute AB and K */
void Reconstructor::Setup()
{
  //! Note: This function should be called after coordinates and delta_xyz have been computed
  
  //! Get domain info
  int i0, j0, k0, imax, jmax, kmax;
  delta_xyz.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);

  Vec3D*** dxyz = (Vec3D***)delta_xyz.GetDataPointer();
  Vec3D*** A    = (Vec3D***)CoeffA.GetDataPointer();
  Vec3D*** B    = (Vec3D***)CoeffB.GetDataPointer();
  Vec3D*** K    = (Vec3D***)CoeffK.GetDataPointer();

  for(int k=k0; k<kmax; k++) {
    for(int j=j0; j<jmax; j++) {
      for(int i=i0; i<imax; i++) {

        A[k][j][i][0] = (dxyz[k][j][i-1][0] + dxyz[k][j][i][0]) / (dxyz[k][j][i][0] + dxyz[k][j][i+1][0]);
        A[k][j][i][1] = (dxyz[k][j-1][i][1] + dxyz[k][j][i][1]) / (dxyz[k][j][i][1] + dxyz[k][j+1][i][1]);
        A[k][j][i][2] = (dxyz[k-1][j][i][2] + dxyz[k][j][i][2]) / (dxyz[k][j][i][2] + dxyz[k+1][j][i][2]);

        B[k][j][i][0] = 2.0*dxyz[k][j][i][0] / (dxyz[k][j][i][0] + dxyz[k][j][i+1][0]);
        B[k][j][i][1] = 2.0*dxyz[k][j][i][1] / (dxyz[k][j][i][1] + dxyz[k][j+1][i][1]);
        B[k][j][i][2] = 2.0*dxyz[k][j][i][2] / (dxyz[k][j][i][2] + dxyz[k+1][j][i][2]);

        K[k][j][i][0] = CalculateSlopeLimiterCoefficientK(A[k][j][i][0],B[k][j][i][0]);    
        K[k][j][i][1] = CalculateSlopeLimiterCoefficientK(A[k][j][i][1],B[k][j][i][1]);    
        K[k][j][i][2] = CalculateSlopeLimiterCoefficientK(A[k][j][i][2],B[k][j][i][2]);    
      }
    }
  }

  delta_xyz.RestoreDataPointerToLocalVector(); //!< no changes to vector
  CoeffA.RestoreDataPointerAndInsert();
  CoeffB.RestoreDataPointerAndInsert();
  CoeffK.RestoreDataPointerAndInsert();
}

//--------------------------------------------------------------------------

int Reconstructor::CalculateSlopeLimiterCoefficientK(double A, double B)
{
  if(iod.schemes.ns.limiter == SchemeData::VANALBADA || 
     iod.schemes.ns.limiter == SchemeData::MODIFIED_VANALBADA) {

    int k;
    double rhs;
    for(k=2; k<1000; k++) {
      rhs = 2.0 / (1.0 + pow((k-1.0)/k, k-1)/k) * std::min(A,1.0);
      if(B <= rhs) break;
    }
    if(k>=999) {
      print_error("Error: Cannot setup slope limiter parameter k for A = %e, B = %e.\n", A, B);
      exit_mpi();
    }
    return k;
  }

  //Otherwise, k is not needed.
  return 0;
}

//--------------------------------------------------------------------------

void Reconstructor::Destroy()
{
  CoeffA.Destroy();
  CoeffB.Destroy();
  CoeffK.Destroy();
}

//--------------------------------------------------------------------------

void Reconstructor::Reconstruct(SpaceVariable3D &U, SpaceVariable3D &Ul, SpaceVariable3D &Ur,
           SpaceVariable3D &Ub, SpaceVariable3D &Ut, SpaceVariable3D &Uk, SpaceVariable3D &Uf)
{

  //! Get mesh info
  int i0, j0, k0, imax, jmax, kmax, ii0, jj0, kk0, iimax, jjmax, kkmax, NX, NY, NZ;
  delta_xyz.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  delta_xyz.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);
  delta_xyz.GetGlobalSize(&NX, &NY, &NZ);

  //! Extract "natural" vectors
  Vec5D*** u    = (Vec5D***) U.GetDataPointer(); 
  Vec5D*** ul   = (Vec5D***) Ul.GetDataPointer(); 
  Vec5D*** ur   = (Vec5D***) Ur.GetDataPointer(); 
  Vec5D*** ub   = (Vec5D***) Ub.GetDataPointer(); 
  Vec5D*** ut   = (Vec5D***) Ut.GetDataPointer(); 
  Vec5D*** uk   = (Vec5D***) Uk.GetDataPointer(); 
  Vec5D*** uf   = (Vec5D***) Uf.GetDataPointer(); 
  Vec3D*** A    = (Vec3D***)CoeffA.GetDataPointer();
  Vec3D*** B    = (Vec3D***)CoeffB.GetDataPointer();
  Vec3D*** K    = (Vec3D***)CoeffK.GetDataPointer();
  Vec3D*** dxyz = (Vec3D***)delta_xyz.GetDataPointer();

  //! Number of DOF per cell
  int nDOF = U.NumDOF();

  /***************************************************************
   *  Loop through all the (real and ghost) cells.
   *  Calculate slope limiter --> slope --> face values
   *  Rules:
   *    - Ouside the physical domain (i.e. the external ghost layer): constant reconstruction
   *    - Within the first layer of cells inside the physical domain: constant reconstruction in
   *      the normal direction (normal to the physical domain boundary
   *    - Within the internal ghost layer (i.e. shared with some other subdomains): do not compute
   *    - Otherwise: Normal computation based on user input.
   ***************************************************************/
  double phi[3]; //!< slope limiter
  double sigma[3]; //!< slope
  double theta[3]; //!< input of slope limiter function

  double a[3], b[3];
  double alpha = iod.schemes.ns.generalized_minmod_coeff; //!< only needed for gen. minmod
  int kk[3]; //!< only for Van Albada

  for(int k=kk0; k<kkmax; k++) {
    for(int j=jj0; j<jjmax; j++) {
      for(int i=ii0; i<iimax; i++) {

        if(i==ii0 || i==iimax-1 || j==jj0 || j==jjmax-1 || k==kk0 || k==kkmax-1) {// ghost layer
          if(i == -1 || i==NX || j==-1 || j==NY || k==-1 || k==NZ) { //!< outside the physical domain
            for(int dof=0; dof<nDOF; dof++) {
              ul[k][j][i][dof] = u[k][j][i][dof];
              ur[k][j][i][dof] = u[k][j][i][dof];
              ub[k][j][i][dof] = u[k][j][i][dof];
              ut[k][j][i][dof] = u[k][j][i][dof];
              uk[k][j][i][dof] = u[k][j][i][dof];
              uf[k][j][i][dof] = u[k][j][i][dof];
            }
          } else { //!< inside the physical domain (overlapped with another subdomain)
            /*do nothing*/
          }
          continue;
        }

        //! In the real part of the subdomain 
        for(int dof=0; dof<nDOF; dof++) {

          //! get constant coefficients
          a[0] = A[k][j][i][0];
          a[1] = A[k][j][i][1];
          a[2] = A[k][j][i][2];
          b[0] = B[k][j][i][0];
          b[1] = B[k][j][i][1];
          b[2] = B[k][j][i][2];
          kk[0] = round(K[k][j][i][0]);
          kk[1] = round(K[k][j][i][1]);
          kk[2] = round(K[k][j][i][2]);

          //! calculate theta: input argument of slope limiter function
          theta[0] = (u[k][j][i][dof] - u[k][j][i-1][dof]) / (u[k][j][i+1][dof] - u[k][j][i][dof]);
          theta[1] = (u[k][j][i][dof] - u[k][j-1][i][dof]) / (u[k][j+1][i][dof] - u[k][j][i][dof]);
          theta[2] = (u[k][j][i][dof] - u[k-1][j][i][dof]) / (u[k+1][j][i][dof] - u[k][j][i][dof]);

          //! calculate slope limiter phi within cell (i,j) 
          phi[0] = phi[1] = phi[2] = 0.0; //!< zero slope
          if(iod.schemes.ns.reconstruction == SchemeData::LINEAR) {
            switch (iod.schemes.ns.limiter) {
              case SchemeData::GENERALIZED_MINMOD :
                phi[0] = GeneralizedMinMod(a[0], b[0], alpha, theta[0]);
                phi[1] = GeneralizedMinMod(a[1], b[1], alpha, theta[1]);
                phi[2] = GeneralizedMinMod(a[2], b[2], alpha, theta[2]);
                break;
              case SchemeData::VANALBADA :
                phi[0] = VanAlbada(a[0], b[0], kk[0], theta[0]);
                phi[1] = VanAlbada(a[1], b[1], kk[1], theta[1]);
                phi[2] = VanAlbada(a[2], b[2], kk[2], theta[2]);
                break;
              case SchemeData::MODIFIED_VANALBADA :
                phi[0] = ModifiedVanAlbada(a[0], b[0], kk[0], theta[0]);
                phi[1] = ModifiedVanAlbada(a[1], b[1], kk[1], theta[1]);
                phi[2] = ModifiedVanAlbada(a[2], b[2], kk[2], theta[2]);
                break;
              case SchemeData::NONE :
                phi[0] = 1.0;
                phi[1] = 1.0;
                phi[2] = 1.0;
                break;
            }
          }
  
          //! calculate the slope sigma (times half cell width) within cell (i,j)
          sigma[0] = 0.5*phi[0]*(u[k][j][i+1][dof]-u[k][j][i][dof]);
          sigma[1] = 0.5*phi[1]*(u[k][j+1][i][dof]-u[k][j][i][dof]);
          sigma[2] = 0.5*phi[2]*(u[k+1][j][i][dof]-u[k][j][i][dof]);

          //! calculate face values
          ul[k][j][i][dof] = u[k][j][i][dof] - sigma[0];
          ur[k][j][i][dof] = u[k][j][i][dof] + sigma[0];
          ub[k][j][i][dof] = u[k][j][i][dof] - sigma[1];
          ut[k][j][i][dof] = u[k][j][i][dof] + sigma[1];
          uk[k][j][i][dof] = u[k][j][i][dof] - sigma[2];
          uf[k][j][i][dof] = u[k][j][i][dof] + sigma[2];

          //! For first-layer cells, switch back to constant reconstruction in the normal dir.
          if(i==0 || i==NX-1) {
            ul[k][j][i][dof] = u[k][j][i][dof];
            ur[k][j][i][dof] = u[k][j][i][dof];
          }
          if(j==0 || j==NY-1) {
            ub[k][j][i][dof] = u[k][j][i][dof];
            ut[k][j][i][dof] = u[k][j][i][dof];
          }
          if(k==0 || k==NZ-1) {
            uk[k][j][i][dof] = u[k][j][i][dof];
            uf[k][j][i][dof] = u[k][j][i][dof];
          }

        }

      }
    }
  }

  //! Restore vectors
  CoeffA.RestoreDataPointerToLocalVector(); //!< no changes to vector
  CoeffB.RestoreDataPointerToLocalVector(); //!< no changes to vector
  CoeffK.RestoreDataPointerToLocalVector(); //!< no changes to vector
  U.RestoreDataPointerToLocalVector(); //!< no changes to vector
  delta_xyz.RestoreDataPointerToLocalVector(); //!< no changes to vector
  Ul.RestoreDataPointerAndInsert();
  Ur.RestoreDataPointerAndInsert();
  Ub.RestoreDataPointerAndInsert();
  Ut.RestoreDataPointerAndInsert();
  Uk.RestoreDataPointerAndInsert();
  Uf.RestoreDataPointerAndInsert();

}

//--------------------------------------------------------------------------














