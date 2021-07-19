#include<Reconstructor.h>
#include<Vector3D.h>
#include<Utils.h>
using std::round;


//--------------------------------------------------------------------------

Reconstructor::Reconstructor(MPI_Comm &comm_, DataManagers3D &dm_all_, ReconstructionData &iod_rec_, 
                             SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_,
                             vector<VarFcnBase*>* vf_, FluxFcnBase* ff_)
                   : iod_rec(iod_rec_), delta_xyz(delta_xyz_), varFcn(vf_), fluxFcn(ff_),
                     CoeffA(comm_, &(dm_all_.ghosted1_3dof)), 
                     CoeffB(comm_, &(dm_all_.ghosted1_3dof)), 
                     CoeffK(comm_, &(dm_all_.ghosted1_3dof)),
                     U(comm_, &(dm_all_.ghosted1_5dof)),
                     ghost_nodes_inner(NULL), ghost_nodes_outer(NULL)
{
  if(iod_rec.varType != ReconstructionData::PRIMITIVE && (!varFcn || !fluxFcn)) {
    print_error("*** Error: Reconstructor needs to know VarFcn and FluxFcn. (Software bug)\n");
    exit_mpi();
  }
}

//--------------------------------------------------------------------------

Reconstructor::~Reconstructor()
{

}

//--------------------------------------------------------------------------
/** Compute AB and K */
void Reconstructor::Setup(vector<GhostPoint> *inner, vector<GhostPoint> *outer)
{
  //! Note: This function should be called after coordinates and delta_xyz have been computed
  
  //! Store pointers to ghost nodes
  ghost_nodes_inner = inner;
  ghost_nodes_outer = outer;

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

        //fprintf(stderr,"(%d,%d,%d): A = %e %e %e, B = %e %e %e, k = %e %e %e\n", i,j,k, A[k][j][i][0], A[k][j][i][1], A[k][j][i][2], B[k][j][i][0], B[k][j][i][1], B[k][j][i][2], K[k][j][i][0], K[k][j][i][1], K[k][j][i][2]);
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
  if(iod_rec.limiter == ReconstructionData::VANALBADA){
    int k;
    double rhs;
    for(k=2; k<1000; k++) {
      rhs = 2.0 / (1.0 + pow((k-1.0)/k, k-1)/k) * std::min(A,1.0);
      if(B <= rhs) break;
    }
    if(k>=999) {
      print_error("*** Error: Cannot setup slope limiter parameter k for A = %e, B = %e.\n", A, B);
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
  U.Destroy();
}

//--------------------------------------------------------------------------
/** Linear reconstruction in 3D
 *  The input and output variables are assumed to be primitive variables. If IoData specifies a
 *  different variable to be reconstructed (e.g., primitive or characterstic), a conversion is
 *  performed within this function. */
void Reconstructor::Reconstruct(SpaceVariable3D &V, SpaceVariable3D &Vl, SpaceVariable3D &Vr,
           SpaceVariable3D &Vb, SpaceVariable3D &Vt, SpaceVariable3D &Vk, SpaceVariable3D &Vf,
           SpaceVariable3D *ID)
{

  //! Constant reconstruction is trivial.
  if(iod_rec.type == ReconstructionData::CONSTANT) {
    Vl.AXPlusBY(0.0, 1.0, V, true);
    Vr.AXPlusBY(0.0, 1.0, V, true);
    Vb.AXPlusBY(0.0, 1.0, V, true);
    Vt.AXPlusBY(0.0, 1.0, V, true);
    Vk.AXPlusBY(0.0, 1.0, V, true);
    Vf.AXPlusBY(0.0, 1.0, V, true);
    return;
  }


  //! Linear reconstruction from now on

  //! Check for obvious errors
  if(iod_rec.varType != ReconstructionData::PRIMITIVE && ID == NULL) {
    print_error("*** Error: Reconstructor::Reconstruct needs materiald ID, which is not provided.\n");
    exit_mpi();
  }
  if(iod_rec.varType != ReconstructionData::PRIMITIVE && V.NumDOF() != 5) {
    print_error("*** Error: In Reconstructor::Reconstruct, expect NumDOF = 5, but got %d.\n", V.NumDOF());
    exit_mpi();
  }
  //! Get mesh info
  int i0, j0, k0, imax, jmax, kmax, ii0, jj0, kk0, iimax, jjmax, kkmax, NX, NY, NZ;
  delta_xyz.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  delta_xyz.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);
  delta_xyz.GetGlobalSize(&NX, &NY, &NZ);

  //! Extract "natural" vectors
  double*** v    = (double***) V.GetDataPointer(); 
  double*** id   = ID ? (double***) ID->GetDataPointer() : NULL;
  double*** vl   = (double***) Vl.GetDataPointer(); 
  double*** vr   = (double***) Vr.GetDataPointer(); 
  double*** vb   = (double***) Vb.GetDataPointer(); 
  double*** vt   = (double***) Vt.GetDataPointer(); 
  double*** vk   = (double***) Vk.GetDataPointer(); 
  double*** vf   = (double***) Vf.GetDataPointer(); 
  Vec3D*** A    = (Vec3D***)CoeffA.GetDataPointer();
  Vec3D*** B    = (Vec3D***)CoeffB.GetDataPointer();
  Vec3D*** K    = (Vec3D***)CoeffK.GetDataPointer();
  Vec3D*** dxyz = (Vec3D***)delta_xyz.GetDataPointer();

  //! Number of DOF per cell
  int nDOF = V.NumDOF();

  //! Convert V to U (conservative) if needed
  double*** u = (double***) U.GetDataPointer(); 
  if(iod_rec.varType == ReconstructionData::CONSERVATIVE ||
     iod_rec.varType == ReconstructionData::CONSERVATIVE_CHARACTERISTIC) {
    for(int k=kk0; k<kkmax; k++)
      for(int j=jj0; j<jjmax; j++)
        for(int i=ii0; i<iimax; i++)
          (*varFcn)[id[k][j][i]]->PrimitiveToConservative(&v[k][j][i*nDOF], &u[k][j][i*nDOF]);          
  }

  /***************************************************************
   *  Loop through all the real cells.
   *  Calculate slope limiter --> slope --> face values
   ***************************************************************/
  double dql[nDOF], dqr[nDOF], dqb[nDOF], dqt[nDOF], dqk[nDOF], dqf[nDOF];
  double sigmax[nDOF], sigmay[nDOF], sigmaz[nDOF]; 

  double a[3], b[3];
  double alpha = iod_rec.generalized_minmod_coeff; //!< only needed for gen. minmod
  int kay[3]; //!< only for Van Albada

  //----------------------------------------------------------------
  // Step 1: Reconstruction within the interior of each subdomain
  //----------------------------------------------------------------
  int vType;
  for(int k=k0; k<kmax; k++) {
    for(int j=j0; j<jmax; j++) {
      for(int i=i0; i<imax; i++) {
        
        vType = iod_rec.varType;

RETRY: 
        //---------------------------
        // Step 1.1. Calculate dq
        //---------------------------
        if(vType == ReconstructionData::PRIMITIVE ||
           vType == ReconstructionData::PRIMITIVE_CHARACTERISTIC) {

          for(int dof=0; dof<nDOF; dof++) {
            dql[dof] = v[k][j][i*nDOF+dof]     - v[k][j][(i-1)*nDOF+dof];
            dqr[dof] = v[k][j][(i+1)*nDOF+dof] - v[k][j][i*nDOF+dof];
            dqb[dof] = v[k][j][i*nDOF+dof]     - v[k][j-1][i*nDOF+dof];
            dqt[dof] = v[k][j+1][i*nDOF+dof]   - v[k][j][i*nDOF+dof];
            dqk[dof] = v[k][j][i*nDOF+dof]     - v[k-1][j][i*nDOF+dof];
            dqf[dof] = v[k+1][j][i*nDOF+dof]   - v[k][j][i*nDOF+dof];
          }

          if(vType == ReconstructionData::PRIMITIVE_CHARACTERISTIC) {
            // convert to characteristic
            int myid = (int)id[k][j][i];
            double dw[5];
            fluxFcn->PrimitiveToPrimitiveCharacteristic(0, &v[k][j][i*nDOF], dql, myid, dw);
            copyarray(dw, dql, 5);
            fluxFcn->PrimitiveToPrimitiveCharacteristic(0, &v[k][j][i*nDOF], dqr, myid, dw);
            copyarray(dw, dqr, 5);
            fluxFcn->PrimitiveToPrimitiveCharacteristic(1, &v[k][j][i*nDOF], dqb, myid, dw);
            copyarray(dw, dqb, 5);
            fluxFcn->PrimitiveToPrimitiveCharacteristic(1, &v[k][j][i*nDOF], dqt, myid, dw);
            copyarray(dw, dqt, 5);
            fluxFcn->PrimitiveToPrimitiveCharacteristic(2, &v[k][j][i*nDOF], dqk, myid, dw);
            copyarray(dw, dqk, 5);
            fluxFcn->PrimitiveToPrimitiveCharacteristic(2, &v[k][j][i*nDOF], dqf, myid, dw);
            copyarray(dw, dqf, 5);
          }

        } 
        else if (vType == ReconstructionData::CONSERVATIVE ||
                 vType == ReconstructionData::CONSERVATIVE_CHARACTERISTIC) { 

          for(int dof=0; dof<nDOF; dof++) {
            dql[dof] = u[k][j][i*nDOF+dof]     - u[k][j][(i-1)*nDOF+dof];
            dqr[dof] = u[k][j][(i+1)*nDOF+dof] - u[k][j][i*nDOF+dof];
            dqb[dof] = u[k][j][i*nDOF+dof]     - u[k][j-1][i*nDOF+dof];
            dqt[dof] = u[k][j+1][i*nDOF+dof]   - u[k][j][i*nDOF+dof];
            dqk[dof] = u[k][j][i*nDOF+dof]     - u[k-1][j][i*nDOF+dof];
            dqf[dof] = u[k+1][j][i*nDOF+dof]   - u[k][j][i*nDOF+dof];
          }

          if(vType == ReconstructionData::CONSERVATIVE_CHARACTERISTIC) {
            // convert to characteristic
            int myid = (int)id[k][j][i];
            double dw[5];
            fluxFcn->ConservativeToConservativeCharacteristic(0, &v[k][j][i*nDOF], dql, myid, dw);
            copyarray(dw, dql, 5);
            fluxFcn->ConservativeToConservativeCharacteristic(0, &v[k][j][i*nDOF], dqr, myid, dw);
            copyarray(dw, dqr, 5);
            fluxFcn->ConservativeToConservativeCharacteristic(1, &v[k][j][i*nDOF], dqb, myid, dw);
            copyarray(dw, dqb, 5);
            fluxFcn->ConservativeToConservativeCharacteristic(1, &v[k][j][i*nDOF], dqt, myid, dw);
            copyarray(dw, dqt, 5);
            fluxFcn->ConservativeToConservativeCharacteristic(2, &v[k][j][i*nDOF], dqk, myid, dw);
            copyarray(dw, dqk, 5);
            fluxFcn->ConservativeToConservativeCharacteristic(2, &v[k][j][i*nDOF], dqf, myid, dw);
            copyarray(dw, dqf, 5);
          }

        }

        //------------------------------------------------------------
        // Step 2.2. Calculate slope, component by component
        //------------------------------------------------------------

        for(int dof=0; dof<nDOF; dof++) {
          
          sigmax[dof] = sigmay[dof] = sigmaz[dof] = 0.0;

          //! get constant coefficients
          a[0] = A[k][j][i][0];
          a[1] = A[k][j][i][1];
          a[2] = A[k][j][i][2];
          b[0] = B[k][j][i][0];
          b[1] = B[k][j][i][1];
          b[2] = B[k][j][i][2];

          //! calculate theta: input argument of slope limiter function
          switch (iod_rec.limiter) {
            case ReconstructionData::GENERALIZED_MINMOD :
              sigmax[dof] = GeneralizedMinMod(a[0], b[0], alpha, dql[dof], dqr[dof]);
              sigmay[dof] = GeneralizedMinMod(a[1], b[1], alpha, dqb[dof], dqt[dof]);
              sigmaz[dof] = GeneralizedMinMod(a[2], b[2], alpha, dqk[dof], dqf[dof]);
              break;
            case ReconstructionData::VANALBADA :
              kay[0] = round(K[k][j][i][0]);
              kay[1] = round(K[k][j][i][1]);
              kay[2] = round(K[k][j][i][2]);
              sigmax[dof] = VanAlbada(a[0], b[0], kay[0], dql[dof], dqr[dof]);
              sigmay[dof] = VanAlbada(a[1], b[1], kay[1], dqb[dof], dqt[dof]);
              sigmaz[dof] = VanAlbada(a[2], b[2], kay[2], dqk[dof], dqf[dof]);
              break;
            case ReconstructionData::NONE :
              sigmax[dof] = 0.5*(dql[dof]+dqr[dof]);
              sigmay[dof] = 0.5*(dqb[dof]+dqt[dof]);
              sigmaz[dof] = 0.5*(dqk[dof]+dqf[dof]);
              break;
          }

        }

        //------------------------------------------------------------
        // Step 2.3. Convert back to differences in primitive or conservative variables
        //------------------------------------------------------------
        if(vType == ReconstructionData::PRIMITIVE_CHARACTERISTIC) {
          int myid = (int)id[k][j][i];
          double dw[5];
          fluxFcn->PrimitiveCharacteristicToPrimitive(0, &v[k][j][i*nDOF], sigmax, myid, dw);
          copyarray(dw, sigmax, 5);
          fluxFcn->PrimitiveCharacteristicToPrimitive(1, &v[k][j][i*nDOF], sigmay, myid, dw);
          copyarray(dw, sigmay, 5);
          fluxFcn->PrimitiveCharacteristicToPrimitive(2, &v[k][j][i*nDOF], sigmaz, myid, dw);
          copyarray(dw, sigmaz, 5);
        }
        else if (vType == ReconstructionData::CONSERVATIVE_CHARACTERISTIC) {
          int myid = (int)id[k][j][i];
          double dw[5];
          fluxFcn->ConservativeCharacteristicToConservative(0, &v[k][j][i*nDOF], sigmax, myid, dw);
          copyarray(dw, sigmax, 5);
          fluxFcn->ConservativeCharacteristicToConservative(1, &v[k][j][i*nDOF], sigmay, myid, dw);
          copyarray(dw, sigmay, 5);
          fluxFcn->ConservativeCharacteristicToConservative(2, &v[k][j][i*nDOF], sigmaz, myid, dw);
          copyarray(dw, sigmaz, 5);
        }

        //------------------------------------------------------------------------------------
        // Step 2.4. (Optional) Switch back to constant reconstruction near material interface
        //------------------------------------------------------------------------------------
        if(ID && iod_rec.slopeNearInterface == ReconstructionData::ZERO) {
          if(id[k][j][i] != id[k][j][i-1] || id[k][j][i] != id[k][j][i+1])
            setValue(sigmax, 0.0, nDOF);
          if(id[k][j][i] != id[k][j-1][i] || id[k][j][i] != id[k][j+1][i])
            setValue(sigmay, 0.0, nDOF);
          if(id[k][j][i] != id[k-1][j][i] || id[k][j][i] != id[k+1][j][i])
            setValue(sigmaz, 0.0, nDOF);
        }

        //------------------------------------------------------------
        // Step 2.5. Calculate interface values
        //------------------------------------------------------------
        if(vType == ReconstructionData::PRIMITIVE || 
           vType == ReconstructionData::PRIMITIVE_CHARACTERISTIC) {
          for(int dof=0; dof<nDOF; dof++) {
            vl[k][j][i*nDOF+dof] = v[k][j][i*nDOF+dof] - 0.5*sigmax[dof];
            vr[k][j][i*nDOF+dof] = v[k][j][i*nDOF+dof] + 0.5*sigmax[dof];
            vb[k][j][i*nDOF+dof] = v[k][j][i*nDOF+dof] - 0.5*sigmay[dof];
            vt[k][j][i*nDOF+dof] = v[k][j][i*nDOF+dof] + 0.5*sigmay[dof];
            vk[k][j][i*nDOF+dof] = v[k][j][i*nDOF+dof] - 0.5*sigmaz[dof];
            vf[k][j][i*nDOF+dof] = v[k][j][i*nDOF+dof] + 0.5*sigmaz[dof];
          }
        }
        else if(vType == ReconstructionData::CONSERVATIVE ||
                vType == ReconstructionData::CONSERVATIVE_CHARACTERISTIC) {

          int myid = (int)id[k][j][i];
          double u2[nDOF]; //temporary variable

          for(int dof=0; dof<nDOF; dof++)
            u2[dof] = u[k][j][i*nDOF+dof] - 0.5*sigmax[dof];
          (*varFcn)[myid]->ConservativeToPrimitive(u2, &vl[k][j][i*nDOF]);

          for(int dof=0; dof<nDOF; dof++)
            u2[dof] = u[k][j][i*nDOF+dof] + 0.5*sigmax[dof];
          (*varFcn)[myid]->ConservativeToPrimitive(u2, &vr[k][j][i*nDOF]);

          for(int dof=0; dof<nDOF; dof++)
            u2[dof] = u[k][j][i*nDOF+dof] - 0.5*sigmay[dof];
          (*varFcn)[myid]->ConservativeToPrimitive(u2, &vb[k][j][i*nDOF]);

          for(int dof=0; dof<nDOF; dof++)
            u2[dof] = u[k][j][i*nDOF+dof] + 0.5*sigmay[dof];
          (*varFcn)[myid]->ConservativeToPrimitive(u2, &vt[k][j][i*nDOF]);

          for(int dof=0; dof<nDOF; dof++)
            u2[dof] = u[k][j][i*nDOF+dof] - 0.5*sigmaz[dof];
          (*varFcn)[myid]->ConservativeToPrimitive(u2, &vk[k][j][i*nDOF]);

          for(int dof=0; dof<nDOF; dof++)
            u2[dof] = u[k][j][i*nDOF+dof] + 0.5*sigmaz[dof];
          (*varFcn)[myid]->ConservativeToPrimitive(u2, &vf[k][j][i*nDOF]);

        }

        //------------------------------------------------------------
        // Step 2.6. Check reconstructed values
        //------------------------------------------------------------
        if(ID && nDOF==5) { //reconstructing the fluid state variables
          int myid = (int)id[k][j][i];
          if((*varFcn)[myid]->CheckState(&vl[k][j][i*nDOF])) {
            fprintf(stderr,"Warning: Found nonphysical reconstructed state (vType = %d). Retrying...\n", vType);
            if(vType == ReconstructionData::PRIMITIVE) // constant rec...
              copyarray(&v[k][j][i*nDOF], &vl[k][j][i*nDOF], 5);
            else {
              vType = ReconstructionData::PRIMITIVE;
              goto RETRY;
            }
          }
          else if((*varFcn)[myid]->CheckState(&vr[k][j][i*nDOF])) {
            fprintf(stderr,"Warning: Found nonphysical reconstructed state (vType = %d). Retrying...\n", vType);
            if(vType == ReconstructionData::PRIMITIVE) // constant rec...
              copyarray(&v[k][j][i*nDOF], &vr[k][j][i*nDOF], 5);
            else {
              vType = ReconstructionData::PRIMITIVE;
              goto RETRY;
            }
          }
          else if((*varFcn)[myid]->CheckState(&vb[k][j][i*nDOF])) {
            fprintf(stderr,"Warning: Found nonphysical reconstructed state (vType = %d). Retrying...\n", vType);
            if(vType == ReconstructionData::PRIMITIVE) // constant rec...
              copyarray(&v[k][j][i*nDOF], &vb[k][j][i*nDOF], 5);
            else {
              vType = ReconstructionData::PRIMITIVE;
              goto RETRY;
            }
          }
          else if((*varFcn)[myid]->CheckState(&vt[k][j][i*nDOF])) {
            fprintf(stderr,"Warning: Found nonphysical reconstructed state (vType = %d). Retrying...\n", vType);
            if(vType == ReconstructionData::PRIMITIVE) // constant rec...
              copyarray(&v[k][j][i*nDOF], &vt[k][j][i*nDOF], 5);
            else {
              vType = ReconstructionData::PRIMITIVE;
              goto RETRY;
            }
          }
          else if((*varFcn)[myid]->CheckState(&vk[k][j][i*nDOF])) {
            fprintf(stderr,"Warning: Found nonphysical reconstructed state (vType = %d). Retrying...\n", vType);
            if(vType == ReconstructionData::PRIMITIVE) // constant rec...
              copyarray(&v[k][j][i*nDOF], &vk[k][j][i*nDOF], 5);
            else {
              vType = ReconstructionData::PRIMITIVE;
              goto RETRY;
            }
          }
          else if((*varFcn)[myid]->CheckState(&vf[k][j][i*nDOF])) {
            fprintf(stderr,"Warning: Found nonphysical reconstructed state (vType = %d). Retrying...\n", vType);
            if(vType == ReconstructionData::PRIMITIVE) // constant rec...
              copyarray(&v[k][j][i*nDOF], &vf[k][j][i*nDOF], 5);
            else {
              vType = ReconstructionData::PRIMITIVE;
              goto RETRY;
            }
          }
        }

      }
    }
  }

  
  //----------------------------------------------------------------
  // Step 2: Exchange info with neighbors
  //----------------------------------------------------------------
  Vl.RestoreDataPointerAndInsert();
  Vr.RestoreDataPointerAndInsert();
  Vb.RestoreDataPointerAndInsert();
  Vt.RestoreDataPointerAndInsert();
  Vk.RestoreDataPointerAndInsert();
  Vf.RestoreDataPointerAndInsert();


  //----------------------------------------------------------------
  // Step 3: Update ghost layer outside the physical domain
  //----------------------------------------------------------------
  vl = (double***) Vl.GetDataPointer(); 
  vr = (double***) Vr.GetDataPointer(); 
  vb = (double***) Vb.GetDataPointer(); 
  vt = (double***) Vt.GetDataPointer(); 
  vk = (double***) Vk.GetDataPointer(); 
  vf = (double***) Vf.GetDataPointer(); 
  int i,j,k,ii,jj,kk;

  for(auto gp = ghost_nodes_outer->begin(); gp != ghost_nodes_outer->end(); gp++) {

    if(gp->type_projection != GhostPoint::FACE)
      continue; //skip edges and corners (not needed)

    i  = gp->ijk[0];
    j  = gp->ijk[1];
    k  = gp->ijk[2];
    ii = gp->image_ijk[0];
    jj = gp->image_ijk[1];
    kk = gp->image_ijk[2];
    
    // check boundary condition
    switch (gp->bcType) {

      case MeshData::INLET :
      case MeshData::OUTLET :
        //constant reconstruction (Dirichlet b.c.)
      
        if     (i<0)   copyarray(&v[kk][jj][ii*nDOF], &vr[k][j][i*nDOF], nDOF);
        else if(i>=NX) copyarray(&v[kk][jj][ii*nDOF], &vl[k][j][i*nDOF], nDOF);
        else if(j<0)   copyarray(&v[kk][jj][ii*nDOF], &vt[k][j][i*nDOF], nDOF);
        else if(j>=NY) copyarray(&v[kk][jj][ii*nDOF], &vb[k][j][i*nDOF], nDOF);
        else if(k<0)   copyarray(&v[kk][jj][ii*nDOF], &vf[k][j][i*nDOF], nDOF);
        else if(k>=NZ) copyarray(&v[kk][jj][ii*nDOF], &vk[k][j][i*nDOF], nDOF);

        break;

      case MeshData::SYMMETRY :
      case MeshData::WALL :
        //constant or linear reconstruction, matching the image

        if     (i<0)   copyarray_flip(&vl[kk][jj][ii*nDOF], &vr[k][j][i*nDOF], nDOF, 1);
        else if(i>=NX) copyarray_flip(&vr[kk][jj][ii*nDOF], &vl[k][j][i*nDOF], nDOF, 1);
        else if(j<0)   copyarray_flip(&vb[kk][jj][ii*nDOF], &vt[k][j][i*nDOF], nDOF, 2);
        else if(j>=NY) copyarray_flip(&vt[kk][jj][ii*nDOF], &vb[k][j][i*nDOF], nDOF, 2);
        else if(k<0)   copyarray_flip(&vk[kk][jj][ii*nDOF], &vf[k][j][i*nDOF], nDOF, 3);
        else if(k>=NZ) copyarray_flip(&vf[kk][jj][ii*nDOF], &vk[k][j][i*nDOF], nDOF, 3);

        break;

      default :
        fprintf(stderr,"*** Error: Cannot perform reconstruction for b.c. %d.\n",
                gp->bcType);
    }
  }
  Vl.RestoreDataPointerToLocalVector(); //no need to communicate
  Vr.RestoreDataPointerToLocalVector(); //no need to communicate
  Vb.RestoreDataPointerToLocalVector(); //no need to communicate
  Vt.RestoreDataPointerToLocalVector(); //no need to communicate
  Vk.RestoreDataPointerToLocalVector(); //no need to communicate
  Vf.RestoreDataPointerToLocalVector(); //no need to communicate

  //! Restore vectors
  CoeffA.RestoreDataPointerToLocalVector(); //!< no changes to vector
  CoeffB.RestoreDataPointerToLocalVector(); //!< no changes to vector
  CoeffK.RestoreDataPointerToLocalVector(); //!< no changes to vector
  V.RestoreDataPointerToLocalVector(); //!< no changes to vector
  delta_xyz.RestoreDataPointerToLocalVector(); //!< no changes to vector

  if(ID) ID->RestoreDataPointerToLocalVector(); //!< no changes to vector

  U.RestoreDataPointerToLocalVector(); //!< internal variable

}

//--------------------------------------------------------------------------

void Reconstructor::ReconstructIn1D(int dir/*0~x,1~y,2~z*/, SpaceVariable3D &U, 
                                    SpaceVariable3D &Um, SpaceVariable3D &Up, SpaceVariable3D *Slope)
{
  if(iod_rec.varType != ReconstructionData::PRIMITIVE) {
    print_error("*** Error: Calling Reconstructor::ReconstructIn1D to reconstruct a 'non-primitive' variable.\n");
    exit_mpi();
  }

  //! Get mesh info
  int i0, j0, k0, imax, jmax, kmax, ii0, jj0, kk0, iimax, jjmax, kkmax, NX, NY, NZ;
  delta_xyz.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  delta_xyz.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);
  delta_xyz.GetGlobalSize(&NX, &NY, &NZ);

  //! Extract "natural" vectors
  double*** u    = (double***) U.GetDataPointer(); 
  double*** um   = (double***) Um.GetDataPointer(); 
  double*** up   = (double***) Up.GetDataPointer(); 
  Vec3D*** A    = (Vec3D***)CoeffA.GetDataPointer();
  Vec3D*** B    = (Vec3D***)CoeffB.GetDataPointer();
  Vec3D*** K    = (Vec3D***)CoeffK.GetDataPointer();
  Vec3D*** dxyz = (Vec3D***)delta_xyz.GetDataPointer();

  double***slope = (Slope) ? (double***)Slope->GetDataPointer() : NULL;
    
  //! Number of DOF per cell
  int nDOF = U.NumDOF();
  if(nDOF != 1) {
    print_error("*** Error: ReconstructIn1D only works with nDOF = 1 at the moment. "
                "Detected nDOF = %d.\n", nDOF);
    exit_mpi();
  }

  /***************************************************************
   *  Loop through all the real cells.
   *  Calculate slope limiter --> slope --> face values
   ***************************************************************/
  double sigma; //!< slope
  double dq0 = 0.0, dq1 = 0.0;

  double a, b;
  double alpha = iod_rec.generalized_minmod_coeff; //!< only needed for gen. minmod
  int kay; //!< only for Van Albada

  for(int k=k0; k<kmax; k++) {
    for(int j=j0; j<jmax; j++) {
      for(int i=i0; i<imax; i++) {

        //! In the real part of the subdomain 
        for(int dof=0; dof<nDOF; dof++) {

          //! calculate slope limiter phi within cell (i,j) 
          sigma = 0.0;

          if(iod_rec.type == ReconstructionData::LINEAR) {

            //! get constant coefficients
            a = A[k][j][i][dir];
            b = B[k][j][i][dir];

            //! calculate theta: input argument of slope limiter function
            switch (dir) {
              case 0:
                dq0 = (u[k][j][i*nDOF+dof]     - u[k][j][(i-1)*nDOF+dof]);
                dq1 = (u[k][j][(i+1)*nDOF+dof] - u[k][j][i*nDOF+dof]);
                break;
              case 1:
                dq0 = (u[k][j][i*nDOF+dof]     - u[k][j-1][i*nDOF+dof]);
                dq1 = (u[k][j+1][i*nDOF+dof]   - u[k][j][i*nDOF+dof]);
                break;
              case 2:
                dq0 = (u[k][j][i*nDOF+dof]     - u[k-1][j][i*nDOF+dof]);
                dq1 = (u[k+1][j][i*nDOF+dof]   - u[k][j][i*nDOF+dof]);
                break;
              default:
                print_error("*** Error: dir(%d) not recognized.\n", dir);
                exit_mpi();
            }

            switch (iod_rec.limiter) {
              case ReconstructionData::GENERALIZED_MINMOD :
                sigma = GeneralizedMinMod(a, b, alpha, dq0, dq1);
                break;
              case ReconstructionData::VANALBADA :
                kay = round(K[k][j][i][dir]);
                sigma = VanAlbada(a, b, kay, dq0, dq1);
                break;
              case ReconstructionData::NONE :
                sigma = 0.5*(dq0+dq1);
                break;
            }
          }
  
          //! calculate face values
          um[k][j][i*nDOF+dof] = u[k][j][i*nDOF+dof] - 0.5*sigma;
          up[k][j][i*nDOF+dof] = u[k][j][i*nDOF+dof] + 0.5*sigma;

          if(slope) slope[k][j][i*nDOF+dof] = sigma/dxyz[k][j][i][dir];

        }

      }
    }
  }

  Um.RestoreDataPointerAndInsert();
  Up.RestoreDataPointerAndInsert();

  if(Slope)
    Slope->RestoreDataPointerAndInsert(); 


  // Now, go over the ghost layer
  um   = (double***) Um.GetDataPointer(); 
  up   = (double***) Up.GetDataPointer(); 
  slope = (Slope) ? (double***)Slope->GetDataPointer() : NULL;

  int i,j,k,ii,jj,kk;

  for(auto gp = ghost_nodes_outer->begin(); gp != ghost_nodes_outer->end(); gp++) {

    if(gp->type_projection != GhostPoint::FACE)
      continue; //skip edges and corners (not needed)

    i  = gp->ijk[0];
    j  = gp->ijk[1];
    k  = gp->ijk[2];
    ii = gp->image_ijk[0];
    jj = gp->image_ijk[1];
    kk = gp->image_ijk[2];
    
    // check boundary condition
    switch (gp->bcType) {

      case MeshData::INLET :
      case MeshData::OUTLET :
        //constant reconstruction (Dirichlet b.c.)
      
        if(dir==0) {
          if     (i<0)   copyarray(&u[kk][jj][ii*nDOF], &up[k][j][i*nDOF], nDOF);
          else if(i>=NX) copyarray(&u[kk][jj][ii*nDOF], &um[k][j][i*nDOF], nDOF);
        }
        else if(dir==1) {
          if     (j<0)   copyarray(&u[kk][jj][ii*nDOF], &up[k][j][i*nDOF], nDOF);
          else if(j>=NY) copyarray(&u[kk][jj][ii*nDOF], &um[k][j][i*nDOF], nDOF);
        }
        else if(dir==2) {
          if(k<0)        copyarray(&u[kk][jj][ii*nDOF], &up[k][j][i*nDOF], nDOF);
          else if(k>=NZ) copyarray(&u[kk][jj][ii*nDOF], &um[k][j][i*nDOF], nDOF);
        }
        if(slope) setValue(&slope[k][j][i*nDOF], 0.0, nDOF);
        break;

      case MeshData::SYMMETRY :
      case MeshData::WALL :
        //constant or linear reconstruction, matching the image

        //relying on nDOF = 1!!
        if(dir==0) {
          if     (i<0)   um[k][j][i] = -up[kk][jj][ii];
          else if(i>=NX) up[k][j][i] = -um[kk][jj][ii]; 
          slope[k][j][i] = -slope[kk][jj][ii];
        }
        if(dir==1) {
          if     (j<0)   um[k][j][i] = -up[kk][jj][ii];
          else if(j>=NY) up[k][j][i] = -um[kk][jj][ii]; 
          slope[k][j][i] = -slope[kk][jj][ii];
        }
        if(dir==2) {
          if     (k<0)   um[k][j][i] = -up[kk][jj][ii];
          else if(k>=NZ) up[k][j][i] = -um[kk][jj][ii]; 
          slope[k][j][i] = -slope[kk][jj][ii];
        }

        break;

      default :
        fprintf(stderr,"*** Error: Cannot perform reconstruction for b.c. %d.\n",
                gp->bcType);
    }
  }
  Um.RestoreDataPointerToLocalVector(); //no need to communicate
  Up.RestoreDataPointerToLocalVector(); //no need to communicate


  //! Restore vectors
  CoeffA.RestoreDataPointerToLocalVector(); //!< no changes to vector
  CoeffB.RestoreDataPointerToLocalVector(); //!< no changes to vector
  CoeffK.RestoreDataPointerToLocalVector(); //!< no changes to vector
  U.RestoreDataPointerToLocalVector(); //!< no changes to vector
  delta_xyz.RestoreDataPointerToLocalVector(); //!< no changes to vector
}

//--------------------------------------------------------------------------














