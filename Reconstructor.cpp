/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <Reconstructor.h>
#include <DistancePointToSpheroid.h>
#include <DistancePointToParallelepiped.h>
#include <memory> //std::unique_ptr
using std::round;

extern int INACTIVE_MATERIAL_ID;

//--------------------------------------------------------------------------

Reconstructor::Reconstructor(MPI_Comm &comm_, DataManagers3D &dm_all_, ReconstructionData &iod_rec_, 
                             SpaceVariable3D &coordinates_, SpaceVariable3D &delta_xyz_,
                             vector<VarFcnBase*>* vf_, FluxFcnBase* ff_)
                   : comm(comm_), iod_rec(iod_rec_), coordinates(coordinates_),
                     delta_xyz(delta_xyz_), varFcn(vf_), fluxFcn(ff_),
                     CoeffA(comm_, &(dm_all_.ghosted1_3dof)), 
                     CoeffB(comm_, &(dm_all_.ghosted1_3dof)), 
                     CoeffK(comm_, &(dm_all_.ghosted1_3dof)),
                     U(comm_, &(dm_all_.ghosted1_5dof)),
                     ghost_nodes_inner(NULL), ghost_nodes_outer(NULL),
                     FixedByUser(NULL)
{
  if(iod_rec.varType != ReconstructionData::PRIMITIVE && (!varFcn || !fluxFcn)) {
    print_error(comm, "*** Error: Reconstructor needs to know VarFcn and FluxFcn. (Software bug)\n");
    exit_mpi();
  }

  if(iod_rec.fixes.sphereMap.dataMap.empty() &&
     iod_rec.fixes.parallelepipedMap.dataMap.empty() &&
     iod_rec.fixes.spheroidMap.dataMap.empty() &&
     iod_rec.fixes.cylindersphereMap.dataMap.empty() &&
     iod_rec.fixes.cylinderconeMap.dataMap.empty())
    FixedByUser = NULL;
  else
    FixedByUser = new SpaceVariable3D(comm_, &(dm_all_.ghosted1_1dof));

}

//--------------------------------------------------------------------------

Reconstructor::~Reconstructor()
{
  if(FixedByUser) delete FixedByUser;
}

//--------------------------------------------------------------------------
/** Compute AB and K, and Figure out which nodes are inside fixes */
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

        //fprintf(stdout,"(%d,%d,%d): A = %e %e %e, B = %e %e %e, k = %e %e %e\n", i,j,k, A[k][j][i][0], A[k][j][i][1], A[k][j][i][2], B[k][j][i][0], B[k][j][i][1], B[k][j][i][2], K[k][j][i][0], K[k][j][i][1], K[k][j][i][2]);
      }
    }
  }

  delta_xyz.RestoreDataPointerToLocalVector(); //!< no changes to vector
  CoeffA.RestoreDataPointerAndInsert();
  CoeffB.RestoreDataPointerAndInsert();
  CoeffK.RestoreDataPointerAndInsert();

  if(FixedByUser)
    TagNodesFixedByUser();

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
      print_error(comm, "*** Error: Cannot setup slope limiter parameter k for A = %e, B = %e.\n", A, B);
      exit_mpi();
    }
    return k;
  }

  //Otherwise, k is not needed.
  return 0;
}

//--------------------------------------------------------------------------

void Reconstructor::TagNodesFixedByUser()
{
  if(!FixedByUser)
    return;

  FixedByUser->SetConstantValue(0.0, true); //default value: not fixed

  // Get domain info
  int i0, j0, k0, imax, jmax, kmax;
  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);

  // get data
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
  double*** fixed = (double***)FixedByUser->GetDataPointer();

  // spheres
  for(auto it=iod_rec.fixes.sphereMap.dataMap.begin();
          it!=iod_rec.fixes.sphereMap.dataMap.end(); it++) {

    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);

    print(comm, "- Applying constant reconstruction within sphere %d:\n", it->first);
    print(comm, "  o center: %e %e %e;  radius: %e.\n", x0[0], x0[1], x0[2], it->second->radius);

    if(it->second->side != SphereData::INTERIOR) {
      print_error(comm, "*** Error: Only supports Side = Interior at the moment.\n");
      exit_mpi();
    }

    double dist;
    // loop through the subdomain interior (i.e. No tags applied to ghost nodes outside physical domain)
    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++) {
          dist = (coords[k][j][i]-x0).norm() - it->second->radius;
          if (dist<0)
            fixed[k][j][i] = 1;
        } 
  }

  // parallelepipeds
  for(auto it=iod_rec.fixes.parallelepipedMap.dataMap.begin(); 
          it!=iod_rec.fixes.parallelepipedMap.dataMap.end(); it++) {
  
    Vec3D x0(it->second->x0, it->second->y0, it->second->z0);
    Vec3D oa(it->second->ax, it->second->ay, it->second->az); oa -= x0;
    Vec3D ob(it->second->bx, it->second->by, it->second->bz); ob -= x0;
    Vec3D oc(it->second->cx, it->second->cy, it->second->cz); oc -= x0;

    print(comm, "- Applying constant reconstruction within parallelepiped %d:\n", it->first);
    print(comm, "  o origin: %e %e %e;  edge 1: %e %e %e.\n", x0[0], x0[1], x0[2], oa[0], oa[1], oa[2]);
    print(comm, "  o edge 2: %e %e %e;  edge 3: %e %e %e.\n", ob[0], ob[1], ob[2], oc[0], oc[1], oc[2]);

    if(oa.norm()==0 || ob.norm()==0 || oc.norm()==0 || (oa^ob)*oc<=0.0) {
      print_error(comm, "*** Error: Detected error in a user-specified parallelepiped. "
                  "Overlapping vertices or violation of right-hand rule.\n");
      exit_mpi();
    }

    if(it->second->side != ParallelepipedData::INTERIOR) {
      print_error(comm, "*** Error: Only supports Side = Interior at the moment.\n");
      exit_mpi();
    }

    GeoTools::DistanceFromPointToParallelepiped distCal(x0, oa, ob, oc);

    double dist;
    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++) {
          dist = distCal.Calculate(coords[k][j][i]); //>0 outside the spheroid
          if (dist<0)
            fixed[k][j][i] = 1;
        }
  }   


  // spheroids
  for(auto it=iod_rec.fixes.spheroidMap.dataMap.begin(); 
          it!=iod_rec.fixes.spheroidMap.dataMap.end(); it++) {

    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
    Vec3D axis(it->second->axis_x, it->second->axis_y, it->second->axis_z);

    print(comm, "- Applying constant reconstruction within spheroid %d:\n", it->first);
    print(comm, "  o center: %e %e %e;  axis: %e %e %e.\n", x0[0], x0[1], x0[2], axis[0], axis[1], axis[2]);
    print(comm, "  o semi-length: %e;  radius: %e.\n", it->second->semi_length, it->second->radius);

    if(it->second->side != SpheroidData::INTERIOR) {
      print_error(comm, "*** Error: Only supports Side = Interior at the moment.\n");
      exit_mpi();
    }

    GeoTools::DistanceFromPointToSpheroid distCal(x0, axis, it->second->semi_length, it->second->radius);

    double dist;
    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++) {
          dist = distCal.Calculate(coords[k][j][i]); //>0 outside the spheroid
          if (dist<0)
            fixed[k][j][i] = 1;
        }
  }

  // cylinder-cone
  for(auto it=iod_rec.fixes.cylinderconeMap.dataMap.begin(); 
          it!=iod_rec.fixes.cylinderconeMap.dataMap.end(); it++) {

    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
    Vec3D dir(it->second->nx, it->second->ny, it->second->nz);
    dir /= dir.norm();

    double L = it->second->L; //cylinder height
    double R = it->second->r; //cylinder radius
    double tan_alpha = tan(it->second->opening_angle_degrees/180.0*acos(-1.0));//opening angle
    double Hmax = R/tan_alpha;
    double H = min(it->second->cone_height, Hmax); //cone's height

    print(comm, "- Applying constant reconstruction within cylinder-cone %d:\n", it->first);
    print(comm, "  o base center: %e %e %e;  axis: %e %e %e.\n", 
          x0[0], x0[1], x0[2], dir[0], dir[1], dir[2]);
    print(comm, "  o cylinder height %e;  radius: %e;  openning angle: %e.\n", 
          L, R, it->second->opening_angle_degrees);

    if(it->second->side != CylinderConeData::INTERIOR) {
      print_error(comm, "*** Error: Only supports Side = Interior at the moment.\n");
      exit_mpi();
    }

    double x, r;
    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++) {
          x = (coords[k][j][i]-x0)*dir;
          r = (coords[k][j][i] - x0 - x*dir).norm();
          if( (x>0 && x<L && r<R) || (x>=L && x<L+H && r<(L+Hmax-x)*tan_alpha) ) //inside
            fixed[k][j][i] = 1;
        }
  }


  // cylinder with spherical caps
  for(auto it=iod_rec.fixes.cylindersphereMap.dataMap.begin(); 
          it!=iod_rec.fixes.cylindersphereMap.dataMap.end(); it++) {

    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z); //base center
    Vec3D dir(it->second->nx, it->second->ny, it->second->nz);
    dir /= dir.norm();

    double L = it->second->L; //cylinder height
    double R = it->second->r; //cylinder and sphere's radius
    double Lhalf = 0.5*L;
    bool front_cap = (it->second->front_cap == CylinderSphereData::On);
    bool back_cap = (it->second->back_cap == CylinderSphereData::On);

    x0 += Lhalf*dir; //now, x0 becomes the center of the cylinder

    print(comm, "- Applying constant reconstruction within cylinder-sphere %d:\n", it->first);
    print(comm, "  o cylinder center: %e %e %e;  axis: %e %e %e.\n", 
          x0[0], x0[1], x0[2], dir[0], dir[1], dir[2]);
    print(comm, "  o cylinder height %e;  radius: %e.\n", L, R);

    if(it->second->side != CylinderSphereData::INTERIOR) {
      print_error(comm, "*** Error: Only supports Side = Interior at the moment.\n");
      exit_mpi();
    }

    Vec3D xf = x0 + Lhalf*dir;
    Vec3D xb = x0 - Lhalf*dir;
    double x, r;
    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++) {
          x = (coords[k][j][i]-x0)*dir;
          r = (coords[k][j][i] - x0 - x*dir).norm();
          if(x>-Lhalf && x<Lhalf && r<R)
            fixed[k][j][i] = 1;
          else if(front_cap && x>=Lhalf && (coords[k][j][i]-xf).norm() < R)
            fixed[k][j][i] = 1;
          else if(back_cap && x<=-Lhalf && (coords[k][j][i]-xb).norm() < R)
            fixed[k][j][i] = 1;
        }
  }


  coordinates.RestoreDataPointerToLocalVector();

  FixedByUser->RestoreDataPointerAndInsert();

}

//--------------------------------------------------------------------------

void Reconstructor::Destroy()
{
  CoeffA.Destroy();
  CoeffB.Destroy();
  CoeffK.Destroy();
  U.Destroy();

  if(FixedByUser)
    FixedByUser->Destroy();
}

//--------------------------------------------------------------------------
/** Linear reconstruction in 3D
 *  The input and output variables are assumed to be primitive variables. If IoData specifies a
 *  different variable to be reconstructed (e.g., primitive or characterstic), a conversion is
 *  performed within this function. */
void Reconstructor::Reconstruct(SpaceVariable3D &V, SpaceVariable3D &Vl, SpaceVariable3D &Vr,
           SpaceVariable3D &Vb, SpaceVariable3D &Vt, SpaceVariable3D &Vk, SpaceVariable3D &Vf,
           SpaceVariable3D *ID, vector<std::unique_ptr<EmbeddedBoundaryDataSet> > *EBDS,
           SpaceVariable3D *Selected, bool do_nothing_if_not_selected)
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
    print_error(comm, "*** Error: Reconstructor::Reconstruct needs materiald ID, which is not provided.\n");
    exit_mpi();
  }
  if(iod_rec.varType != ReconstructionData::PRIMITIVE && V.NumDOF() != 5) {
    print_error(comm, "*** Error: In Reconstructor::Reconstruct, expect NumDOF = 5, but got %d.\n", V.NumDOF());
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
  //Vec3D*** dxyz = (Vec3D***)delta_xyz.GetDataPointer();

  //user-specified "Fixes"
  double*** fixed = FixedByUser ? FixedByUser->GetDataPointer() : NULL;

  //selected cells specified
  double*** sel = Selected ? Selected->GetDataPointer() : NULL;

  //May need intersections
  vector<Vec3D***> xf;
  if(EBDS && iod_rec.slopeNearInterface == ReconstructionData::ZERO)
    for(auto it = EBDS->begin(); it != EBDS->end(); it++) 
      xf.push_back((Vec3D***)(*it)->XForward_ptr->GetDataPointer());

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
        
        //---------------------------
        // If node is not selected, skip
        //---------------------------
        if(sel && do_nothing_if_not_selected && !sel[k][j][i])
          continue;

        //---------------------------
        // If node is 'fixed', trivial (const. rec.)
        //---------------------------
        if((fixed && fixed[k][j][i]) ||
           (sel && !sel[k][j][i])) {
          for(int dof=0; dof<nDOF; dof++) {
            vl[k][j][i*nDOF+dof] = v[k][j][i*nDOF+dof];
            vr[k][j][i*nDOF+dof] = v[k][j][i*nDOF+dof];
            vb[k][j][i*nDOF+dof] = v[k][j][i*nDOF+dof];
            vt[k][j][i*nDOF+dof] = v[k][j][i*nDOF+dof];
            vk[k][j][i*nDOF+dof] = v[k][j][i*nDOF+dof];
            vf[k][j][i*nDOF+dof] = v[k][j][i*nDOF+dof];
          }
          continue;
        }

        //---------------------------
        // If node is occluded, nothing needs to be done
        //---------------------------
        if(id && ((int)id[k][j][i] == INACTIVE_MATERIAL_ID)) {
          for(int dof=0; dof<nDOF; dof++) {
            vl[k][j][i*nDOF+dof] = v[k][j][i*nDOF+dof];
            vr[k][j][i*nDOF+dof] = v[k][j][i*nDOF+dof];
            vb[k][j][i*nDOF+dof] = v[k][j][i*nDOF+dof];
            vt[k][j][i*nDOF+dof] = v[k][j][i*nDOF+dof];
            vk[k][j][i*nDOF+dof] = v[k][j][i*nDOF+dof];
            vf[k][j][i*nDOF+dof] = v[k][j][i*nDOF+dof];
          }
          continue;
        }

        // Get variable type
        vType = iod_rec.varType;

RETRY: 
        //---------------------------
        // Step 1.1. Calculate dq
        //---------------------------
        if(vType == ReconstructionData::PRIMITIVE ||
           vType == ReconstructionData::PRIMITIVE_CHARACTERISTIC) {

          for(int dof=0; dof<nDOF; dof++) {
            if(id && ((int)id[k][j][i-1] == INACTIVE_MATERIAL_ID))
              dql[dof] = 0.0;
            else
              dql[dof] = v[k][j][i*nDOF+dof]     - v[k][j][(i-1)*nDOF+dof];

            if(id && ((int)id[k][j][i+1] == INACTIVE_MATERIAL_ID))
              dqr[dof] = 0.0;
            else
              dqr[dof] = v[k][j][(i+1)*nDOF+dof] - v[k][j][i*nDOF+dof];

            if(id && ((int)id[k][j-1][i] == INACTIVE_MATERIAL_ID))
              dqb[dof] = 0.0;
            else
              dqb[dof] = v[k][j][i*nDOF+dof]     - v[k][j-1][i*nDOF+dof];

            if(id && ((int)id[k][j+1][i] == INACTIVE_MATERIAL_ID))
              dqt[dof] = 0.0;
            else
              dqt[dof] = v[k][j+1][i*nDOF+dof]   - v[k][j][i*nDOF+dof];

            if(id && ((int)id[k-1][j][i] == INACTIVE_MATERIAL_ID))
              dqk[dof] = 0.0;
            else
              dqk[dof] = v[k][j][i*nDOF+dof]     - v[k-1][j][i*nDOF+dof];

            if(id && ((int)id[k+1][j][i] == INACTIVE_MATERIAL_ID))
              dqf[dof] = 0.0;
            else
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
            if(id && ((int)id[k][j][i-1] == INACTIVE_MATERIAL_ID))
              dql[dof] = 0.0;
            else
              dql[dof] = u[k][j][i*nDOF+dof]     - u[k][j][(i-1)*nDOF+dof];

            if(id && ((int)id[k][j][i+1] == INACTIVE_MATERIAL_ID))
              dqr[dof] = 0.0;
            else
              dqr[dof] = u[k][j][(i+1)*nDOF+dof] - u[k][j][i*nDOF+dof];

            if(id && ((int)id[k][j-1][i] == INACTIVE_MATERIAL_ID))
              dqb[dof] = 0.0;
            else
              dqb[dof] = u[k][j][i*nDOF+dof]     - u[k][j-1][i*nDOF+dof];

            if(id && ((int)id[k][j+1][i] == INACTIVE_MATERIAL_ID))
              dqt[dof] = 0.0;
            else
              dqt[dof] = u[k][j+1][i*nDOF+dof]   - u[k][j][i*nDOF+dof];

            if(id && ((int)id[k-1][j][i] == INACTIVE_MATERIAL_ID))
              dqk[dof] = 0.0;
            else
              dqk[dof] = u[k][j][i*nDOF+dof]     - u[k-1][j][i*nDOF+dof];

            if(id && ((int)id[k+1][j][i] == INACTIVE_MATERIAL_ID))
              dqf[dof] = 0.0;
            else
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
        } else if (vType == ReconstructionData::CONSERVATIVE_CHARACTERISTIC) {
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
        // Step 2.4. (Optional) Switch back to constant reconstruction near material interface (not FS interface)
        //------------------------------------------------------------------------------------
        if(iod_rec.slopeNearInterface == ReconstructionData::ZERO) {
          if(ID) {
            if(id[k][j][i] != id[k][j][i-1] || id[k][j][i] != id[k][j][i+1])
              setValue(sigmax, 0.0, nDOF);
            if(id[k][j][i] != id[k][j-1][i] || id[k][j][i] != id[k][j+1][i])
              setValue(sigmay, 0.0, nDOF);
            if(id[k][j][i] != id[k-1][j][i] || id[k][j][i] != id[k+1][j][i])
              setValue(sigmaz, 0.0, nDOF);
          }

          if(xf.size()>0) {
            for(int s=0; s<(int)xf.size(); s++) {
              if(xf[s][k][j][i][0]>=0 || (i+1<NX && xf[s][k][j][i+1][0]>=0))
                setValue(sigmax, 0.0, nDOF);
              if(xf[s][k][j][i][1]>=0 || (j+1<NY && xf[s][k][j+1][i][1]>=0))
                setValue(sigmay, 0.0, nDOF);
              if(xf[s][k][j][i][2]>=0 || (k+1<NZ && xf[s][k+1][j][i][2]>=0))
                setValue(sigmaz, 0.0, nDOF);
            }
          }
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
        extern int verbose;
        if(ID && nDOF==5) { //reconstructing the fluid state variables
          int myid = (int)id[k][j][i];
          if((*varFcn)[myid]->CheckState(&vl[k][j][i*nDOF])) {
            if(verbose>1) fprintf(stdout,"Warning: Found nonphysical reconstructed state (vType = %d). Retrying...\n", vType);
            if(vType == ReconstructionData::PRIMITIVE) // constant rec...
              copyarray(&v[k][j][i*nDOF], &vl[k][j][i*nDOF], 5);
            else {
              vType = ReconstructionData::PRIMITIVE;
              goto RETRY;
            }
          }
          else if((*varFcn)[myid]->CheckState(&vr[k][j][i*nDOF])) {
            if(verbose>1) fprintf(stdout,"Warning: Found nonphysical reconstructed state (vType = %d). Retrying...\n", vType);
            if(vType == ReconstructionData::PRIMITIVE) // constant rec...
              copyarray(&v[k][j][i*nDOF], &vr[k][j][i*nDOF], 5);
            else {
              vType = ReconstructionData::PRIMITIVE;
              goto RETRY;
            }
          }
          else if((*varFcn)[myid]->CheckState(&vb[k][j][i*nDOF])) {
            if(verbose>1) fprintf(stdout,"Warning: Found nonphysical reconstructed state (vType = %d). Retrying...\n", vType);
            if(vType == ReconstructionData::PRIMITIVE) // constant rec...
              copyarray(&v[k][j][i*nDOF], &vb[k][j][i*nDOF], 5);
            else {
              vType = ReconstructionData::PRIMITIVE;
              goto RETRY;
            }
          }
          else if((*varFcn)[myid]->CheckState(&vt[k][j][i*nDOF])) {
            if(verbose>1) fprintf(stdout,"Warning: Found nonphysical reconstructed state (vType = %d). Retrying...\n", vType);
            if(vType == ReconstructionData::PRIMITIVE) // constant rec...
              copyarray(&v[k][j][i*nDOF], &vt[k][j][i*nDOF], 5);
            else {
              vType = ReconstructionData::PRIMITIVE;
              goto RETRY;
            }
          }
          else if((*varFcn)[myid]->CheckState(&vk[k][j][i*nDOF])) {
            if(verbose>1) fprintf(stdout,"Warning: Found nonphysical reconstructed state (vType = %d). Retrying...\n", vType);
            if(vType == ReconstructionData::PRIMITIVE) // constant rec...
              copyarray(&v[k][j][i*nDOF], &vk[k][j][i*nDOF], 5);
            else {
              vType = ReconstructionData::PRIMITIVE;
              goto RETRY;
            }
          }
          else if((*varFcn)[myid]->CheckState(&vf[k][j][i*nDOF])) {
            if(verbose>1) fprintf(stdout,"Warning: Found nonphysical reconstructed state (vType = %d). Retrying...\n", vType);
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
    
    if(sel && do_nothing_if_not_selected && !sel[k][j][i])
      continue;

    // check boundary condition
    switch (gp->bcType) {

      case MeshData::INLET :
      case MeshData::INLET2 :
        //constant reconstruction (Dirichlet b.c.)
      
        if     (i<0)   copyarray(&v[k][j][i*nDOF], &vr[k][j][i*nDOF], nDOF);
        else if(i>=NX) copyarray(&v[k][j][i*nDOF], &vl[k][j][i*nDOF], nDOF);
        else if(j<0)   copyarray(&v[k][j][i*nDOF], &vt[k][j][i*nDOF], nDOF);
        else if(j>=NY) copyarray(&v[k][j][i*nDOF], &vb[k][j][i*nDOF], nDOF);
        else if(k<0)   copyarray(&v[k][j][i*nDOF], &vf[k][j][i*nDOF], nDOF);
        else if(k>=NZ) copyarray(&v[k][j][i*nDOF], &vk[k][j][i*nDOF], nDOF);

        break;

      case MeshData::OVERSET :
        //constant reconstruction (Dirichlet b.c.)
      
        if     (i<0)   copyarray(&v[k][j][i*nDOF], &vr[k][j][i*nDOF], nDOF);
        else if(i>=NX) copyarray(&v[k][j][i*nDOF], &vl[k][j][i*nDOF], nDOF);
        else if(j<0)   copyarray(&v[k][j][i*nDOF], &vt[k][j][i*nDOF], nDOF);
        else if(j>=NY) copyarray(&v[k][j][i*nDOF], &vb[k][j][i*nDOF], nDOF);
        else if(k<0)   copyarray(&v[k][j][i*nDOF], &vf[k][j][i*nDOF], nDOF);
        else if(k>=NZ) copyarray(&v[k][j][i*nDOF], &vk[k][j][i*nDOF], nDOF);

        break;

      case MeshData::OUTLET :
        //constant or linear reconstruction, matching the image
      
        if     (i<0)   copyarray(&vl[kk][jj][ii*nDOF], &vr[k][j][i*nDOF], nDOF);
        else if(i>=NX) copyarray(&vr[kk][jj][ii*nDOF], &vl[k][j][i*nDOF], nDOF);
        else if(j<0)   copyarray(&vb[kk][jj][ii*nDOF], &vt[k][j][i*nDOF], nDOF);
        else if(j>=NY) copyarray(&vt[kk][jj][ii*nDOF], &vb[k][j][i*nDOF], nDOF);
        else if(k<0)   copyarray(&vk[kk][jj][ii*nDOF], &vf[k][j][i*nDOF], nDOF);
        else if(k>=NZ) copyarray(&vf[kk][jj][ii*nDOF], &vk[k][j][i*nDOF], nDOF);

        break;

      case MeshData::SYMMETRY :
      case MeshData::SLIPWALL :
        //constant or linear reconstruction, matching the image

        if     (i<0)   copyarray_flip(&vl[kk][jj][ii*nDOF], &vr[k][j][i*nDOF], nDOF, 1);
        else if(i>=NX) copyarray_flip(&vr[kk][jj][ii*nDOF], &vl[k][j][i*nDOF], nDOF, 1);
        else if(j<0)   copyarray_flip(&vb[kk][jj][ii*nDOF], &vt[k][j][i*nDOF], nDOF, 2);
        else if(j>=NY) copyarray_flip(&vt[kk][jj][ii*nDOF], &vb[k][j][i*nDOF], nDOF, 2);
        else if(k<0)   copyarray_flip(&vk[kk][jj][ii*nDOF], &vf[k][j][i*nDOF], nDOF, 3);
        else if(k>=NZ) copyarray_flip(&vf[kk][jj][ii*nDOF], &vk[k][j][i*nDOF], nDOF, 3);

        break;

      case MeshData::STICKWALL :
        if     (i<0)   copyarray_flip(&vl[kk][jj][ii*nDOF], &vr[k][j][i*nDOF], nDOF, 1, 3);
        else if(i>=NX) copyarray_flip(&vr[kk][jj][ii*nDOF], &vl[k][j][i*nDOF], nDOF, 1, 3);
        else if(j<0)   copyarray_flip(&vb[kk][jj][ii*nDOF], &vt[k][j][i*nDOF], nDOF, 1, 3);
        else if(j>=NY) copyarray_flip(&vt[kk][jj][ii*nDOF], &vb[k][j][i*nDOF], nDOF, 1, 3);
        else if(k<0)   copyarray_flip(&vk[kk][jj][ii*nDOF], &vf[k][j][i*nDOF], nDOF, 1, 3);
        else if(k>=NZ) copyarray_flip(&vf[kk][jj][ii*nDOF], &vk[k][j][i*nDOF], nDOF, 1, 3);

        break;

      default :
        fprintf(stdout,"\033[0;31m*** Error: Cannot perform reconstruction for b.c. %d.\033[0m\n",
                gp->bcType);
    }
  }
  //NOTE: Should not communicate. Otherwise the ghost layer will be corrupated.
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
  //delta_xyz.RestoreDataPointerToLocalVector(); //!< no changes to vector

  if(FixedByUser) FixedByUser->RestoreDataPointerToLocalVector();

  if(ID) ID->RestoreDataPointerToLocalVector(); //!< no changes to vector

  if(Selected) Selected->RestoreDataPointerToLocalVector(); //!< no changes to vector

  if(xf.size()>0) {
    for(auto it = EBDS->begin(); it != EBDS->end(); it++) 
      (*it)->XForward_ptr->RestoreDataPointerToLocalVector();
  }

  U.RestoreDataPointerToLocalVector(); //!< internal variable

}

//--------------------------------------------------------------------------

void Reconstructor::ReconstructIn1D(int dir/*0~x,1~y,2~z*/, SpaceVariable3D &U, 
                                    SpaceVariable3D &Um, SpaceVariable3D &Up, SpaceVariable3D *Slope,
                                    SpaceVariable3D *Selected)
{
  if(iod_rec.varType != ReconstructionData::PRIMITIVE) {
    print_error(comm, "*** Error: Calling Reconstructor::ReconstructIn1D to reconstruct a 'non-primitive' variable.\n");
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

  double***sel = (Selected) ? (double***)Selected->GetDataPointer() : NULL;
    
  //! Number of DOF per cell
  int nDOF = U.NumDOF();
  if(nDOF != 1) {
    print_error(comm, "*** Error: ReconstructIn1D only works with nDOF = 1 at the moment. "
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

        //! if not selected, skip it
        if(sel && !sel[k][j][i])
          continue;

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
                print_error(comm, "*** Error: dir(%d) not recognized.\n", dir);
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
    
    // if not selected, skip it
    if(sel && !sel[k][j][i])
      continue;

    // check boundary condition
    switch (gp->bcType) {

      case MeshData::INLET :
      case MeshData::INLET2 :
        //constant reconstruction (Dirichlet b.c.)
      
        if(dir==0) {
          if     (i<0)   copyarray(&u[k][j][i*nDOF], &up[k][j][i*nDOF], nDOF);
          else if(i>=NX) copyarray(&u[k][j][i*nDOF], &um[k][j][i*nDOF], nDOF);
        }
        else if(dir==1) {
          if     (j<0)   copyarray(&u[k][j][i*nDOF], &up[k][j][i*nDOF], nDOF);
          else if(j>=NY) copyarray(&u[k][j][i*nDOF], &um[k][j][i*nDOF], nDOF);
        }
        else if(dir==2) {
          if(k<0)        copyarray(&u[k][j][i*nDOF], &up[k][j][i*nDOF], nDOF);
          else if(k>=NZ) copyarray(&u[k][j][i*nDOF], &um[k][j][i*nDOF], nDOF);
        }
        if(slope) setValue(&slope[k][j][i*nDOF], 0.0, nDOF);
        break;

      case MeshData::OVERSET :
        //constant reconstruction (Dirichlet b.c.)
      
        if(dir==0) {
          if     (i<0)   copyarray(&u[k][j][i*nDOF], &up[k][j][i*nDOF], nDOF);
          else if(i>=NX) copyarray(&u[k][j][i*nDOF], &um[k][j][i*nDOF], nDOF);
        }
        else if(dir==1) {
          if     (j<0)   copyarray(&u[k][j][i*nDOF], &up[k][j][i*nDOF], nDOF);
          else if(j>=NY) copyarray(&u[k][j][i*nDOF], &um[k][j][i*nDOF], nDOF);
        }
        else if(dir==2) {
          if(k<0)        copyarray(&u[k][j][i*nDOF], &up[k][j][i*nDOF], nDOF);
          else if(k>=NZ) copyarray(&u[k][j][i*nDOF], &um[k][j][i*nDOF], nDOF);
        }
        if(slope) setValue(&slope[k][j][i*nDOF], 0.0, nDOF);
        break;

      case MeshData::OUTLET :

        //relying on nDOF = 1!!
        if(dir==0) {
          if     (i<0)   up[k][j][i] = um[kk][jj][ii];
          else if(i>=NX) um[k][j][i] = up[kk][jj][ii]; 
          slope[k][j][i] = slope[kk][jj][ii];
        }
        if(dir==1) {
          if     (j<0)   up[k][j][i] = um[kk][jj][ii];
          else if(j>=NY) um[k][j][i] = up[kk][jj][ii]; 
          slope[k][j][i] = slope[kk][jj][ii];
        }
        if(dir==2) {
          if     (k<0)   up[k][j][i] = um[kk][jj][ii];
          else if(k>=NZ) um[k][j][i] = up[kk][jj][ii]; 
          slope[k][j][i] = slope[kk][jj][ii];
        }

        break;

      case MeshData::SYMMETRY :
      case MeshData::SLIPWALL :
      case MeshData::STICKWALL : //TODO: Correct??
        //constant or linear reconstruction, matching the image

        //relying on nDOF = 1!!
        if(dir==0) {
          if     (i<0)   up[k][j][i] = -um[kk][jj][ii];
          else if(i>=NX) um[k][j][i] = -up[kk][jj][ii]; 
          slope[k][j][i] = -slope[kk][jj][ii];
        }
        if(dir==1) {
          if     (j<0)   up[k][j][i] = -um[kk][jj][ii];
          else if(j>=NY) um[k][j][i] = -up[kk][jj][ii]; 
          slope[k][j][i] = -slope[kk][jj][ii];
        }
        if(dir==2) {
          if     (k<0)   up[k][j][i] = -um[kk][jj][ii];
          else if(k>=NZ) um[k][j][i] = -up[kk][jj][ii]; 
          slope[k][j][i] = -slope[kk][jj][ii];
        }

        break;

      default :
        fprintf(stdout,"\033[0;31m*** Error: Cannot perform reconstruction for b.c. %d.\033[0m\n",
                gp->bcType);
    }
  }

  //Should not communicate.
  Um.RestoreDataPointerToLocalVector(); //no need to communicate
  Up.RestoreDataPointerToLocalVector(); //no need to communicate
  if(Slope)
    Slope->RestoreDataPointerToLocalVector(); 

  //! Restore vectors
  CoeffA.RestoreDataPointerToLocalVector(); //!< no changes to vector
  CoeffB.RestoreDataPointerToLocalVector(); //!< no changes to vector
  CoeffK.RestoreDataPointerToLocalVector(); //!< no changes to vector
  U.RestoreDataPointerToLocalVector(); //!< no changes to vector
  delta_xyz.RestoreDataPointerToLocalVector(); //!< no changes to vector
  if(Selected)
    Selected->RestoreDataPointerToLocalVector();

}

//--------------------------------------------------------------------------














