#include<Reconstructor.h>
#include<Vector2D.h>
#include<Vector5D.h>
#include<Utils.h>
using std::round;


//--------------------------------------------------------------------------

Reconstructor::Reconstructor(MPI_Comm &comm_, DataManagers2D &dm_all_, IoData &iod_, 
                             SpaceVariable2D &coordinates_, SpaceVariable2D &delta_xy_)
                   : iod(iod_), delta_xy(delta_xy_),
                     CoeffA(comm_, &(dm_all_.ghosted1_2dof)), 
                     CoeffB(comm_, &(dm_all_.ghosted1_2dof)), 
                     CoeffK(comm_, &(dm_all_.ghosted1_2dof))
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
  //! Note: This function should be called after coordinates and delta_xy have been computed
  
  //! Get domain info
  int i0, j0, imax, jmax;
  delta_xy.GetCornerIndices(&i0, &j0, &imax, &jmax);

  Vec2D** dxy = (Vec2D**)delta_xy.GetDataPointer();
  Vec2D** A   = (Vec2D**)CoeffA.GetDataPointer();
  Vec2D** B   = (Vec2D**)CoeffB.GetDataPointer();
  Vec2D** K   = (Vec2D**)CoeffK.GetDataPointer();

  for(int j=j0; j<jmax; j++) {
    for(int i=i0; i<imax; i++) {
      A[j][i][0] = (dxy[j][i-1][0] + dxy[j][i][0]) / (dxy[j][i][0] + dxy[j][i+1][0]);
      A[j][i][1] = (dxy[j-1][i][1] + dxy[j][i][1]) / (dxy[j][i][1] + dxy[j+1][i][1]);
      B[j][i][0] = 2.0*dxy[j][i][0] / (dxy[j][i][0] + dxy[j][i+1][0]);
      B[j][i][1] = 2.0*dxy[j][i][1] / (dxy[j][i][1] + dxy[j+1][i][1]);
      K[j][i][0] = CalculateSlopeLimiterCoefficientK(A[j][i][0],B[j][i][0]);    
      K[j][i][1] = CalculateSlopeLimiterCoefficientK(A[j][i][1],B[j][i][1]);    
    }
  }

  delta_xy.RestoreDataPointerToLocalVector(); //!< no changes to vector
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

void Reconstructor::Reconstruct(SpaceVariable2D &U, SpaceVariable2D &Ul, SpaceVariable2D &Ur,
           SpaceVariable2D &Ub, SpaceVariable2D &Ut)
{

  //! Get mesh info
  int i0, j0, imax, jmax, ii0, jj0, iimax, jjmax, NX, NY;
  delta_xy.GetCornerIndices(&i0, &j0, &imax, &jmax);
  delta_xy.GetGhostedCornerIndices(&ii0, &jj0, &iimax, &jjmax);
  delta_xy.GetGlobalSize(&NX, &NY);

  //! Extract "natural" vectors
  Vec5D** u  = (Vec5D**) U.GetDataPointer(); 
  Vec5D** ul = (Vec5D**) Ul.GetDataPointer(); 
  Vec5D** ur = (Vec5D**) Ur.GetDataPointer(); 
  Vec5D** ub = (Vec5D**) Ub.GetDataPointer(); 
  Vec5D** ut = (Vec5D**) Ut.GetDataPointer(); 
  Vec2D** A   = (Vec2D**)CoeffA.GetDataPointer();
  Vec2D** B   = (Vec2D**)CoeffB.GetDataPointer();
  Vec2D** K   = (Vec2D**)CoeffK.GetDataPointer();
  Vec2D** dxy = (Vec2D**)delta_xy.GetDataPointer();

  //! Number of DOF per cell
  int nDOF = U.NumDOF();

  /***************************************************************
   *  Loop through all the (real and ghost) cells.
   *  Calculate slope limiter --> slope --> face values
   *  Rules:
   *    - Ouside the physical domain (i.e. the external ghost layer): constant reconstruction
   *    - Within the first layer of cells inside the physical domain: constant reconstruction in
   *      the normal direction (normal to the physical boundary
   *    - Within the internal ghost layer (i.e. shared with some other subdomains): do not compute
   *    - Otherwise: Normal computation based on user input.
   ***************************************************************/
  double phi[2]; //!< slope limiter
  double sigma[2]; //!< slope
  double theta[2]; //!< input argument of slope limiter function

  double a[2], b[2];
  double alpha = iod.schemes.ns.generalized_minmod_coeff; //!< only needed for gen. minmod
  int kk[2]; //!< only for Van Albada

  for(int j=jj0; j<jjmax; j++) {
    for(int i=ii0; i<iimax; i++) {

      if(i==ii0 || i==iimax-1 || j==jj0 || j==jjmax-1) {// ghost layer
        if(i == -1 || i==NX-1 || j==-1 || j==NY-1) { //!< outside the physical domain
          for(int dof=0; dof<nDOF; dof++) {
            ul[j][i][dof] = u[j][i][dof];
            ur[j][i][dof] = u[j][i][dof];
            ub[j][i][dof] = u[j][i][dof];
            ut[j][i][dof] = u[j][i][dof];
          }
        } else { //!< inside the physical domain (overlapped with another subdomain)
          /*do nothing*/
        }
        goto label;
      }

      // In the real part of the subdomain 
      for(int dof=0; dof<nDOF; dof++) {

        //! get constant coefficients
        a[0] = A[j][i][0];
        a[1] = A[j][i][1];
        b[0] = B[j][i][0];
        b[1] = B[j][i][1];
        kk[0] = round(K[j][i][0]);
        kk[1] = round(K[j][i][1]);

        //! calculate theta: input argument of slope limiter function
        theta[0] = (u[j][i][dof] - u[j][i-1][dof]) / (u[j][i+1][dof] - u[j][i][dof]);
        theta[1] = (u[j][i][dof] - u[j-1][i][dof]) / (u[j+1][i][dof] - u[j][i][dof]);

        //! calculate slope limiter phi within cell (i,j) 
        phi[0] = phi[1] = 0.0; //!< zero slope
        if(iod.schemes.ns.reconstruction == SchemeData::LINEAR) {
          switch (iod.schemes.ns.limiter) {
            case SchemeData::GENERALIZED_MINMOD :
              phi[0] = GeneralizedMinMod(a[0], b[0], alpha, theta[0]);
              phi[1] = GeneralizedMinMod(a[1], b[1], alpha, theta[1]);
              break;
            case SchemeData::VANALBADA :
              phi[0] = VanAlbada(a[0], b[0], kk[0], theta[0]);
              phi[1] = VanAlbada(a[1], b[1], kk[1], theta[1]);
              break;
            case SchemeData::MODIFIED_VANALBADA :
              phi[0] = ModifiedVanAlbada(a[0], b[0], kk[0], theta[0]);
              phi[1] = ModifiedVanAlbada(a[1], b[1], kk[1], theta[1]);
              break;
            case SchemeData::NONE :
              phi[0] = 1.0;
              phi[1] = 1.0;
              break;
          }
        }
  
        //! calculate the slope sigma (times half cell width) within cell (i,j)
        sigma[0] = 0.5*phi[0]*(u[j][i+1][dof]-u[j][i][dof]);
        sigma[1] = 0.5*phi[1]*(u[j+1][i][dof]-u[j][i][dof]);

        //! calculate face values
        ul[j][i][dof] = u[j][i][dof] - sigma[0];
        ur[j][i][dof] = u[j][i][dof] + sigma[0];
        ub[j][i][dof] = u[j][i][dof] - sigma[1];
        ut[j][i][dof] = u[j][i][dof] + sigma[1];

        //! For first-layer cells, switch back to constant reconstruction in the normal dir.
        if(i==0 || i==NX-2) {
          ul[j][i][dof] = u[j][i][dof];
          ur[j][i][dof] = u[j][i][dof];
        }
        if(j==0 || j==NY-2) {
          ub[j][i][dof] = u[j][i][dof] - sigma[1];
          ut[j][i][dof] = u[j][i][dof] + sigma[1];
        }

      }
    }
  }

label:
  //! Restore vectors
  CoeffA.RestoreDataPointerToLocalVector(); //!< no changes to vector
  CoeffB.RestoreDataPointerToLocalVector(); //!< no changes to vector
  CoeffK.RestoreDataPointerToLocalVector(); //!< no changes to vector
  U.RestoreDataPointerToLocalVector(); //!< no changes to vector
  delta_xy.RestoreDataPointerToLocalVector(); //!< no changes to vector
  Ul.RestoreDataPointerAndInsert();
  Ur.RestoreDataPointerAndInsert();
  Ub.RestoreDataPointerAndInsert();
  Ut.RestoreDataPointerAndInsert();

}

//--------------------------------------------------------------------------














