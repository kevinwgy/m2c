#include<MultiPhaseOperator.h>
#include<SpaceOperator.h>
#include<LevelSetOperator.h>
#include<Vector5D.h>
#include<RiemannSolutions.h>
#include<algorithm>//find
using std::cout;
using std::endl;
extern int verbose;
//-----------------------------------------------------

MultiPhaseOperator::MultiPhaseOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_,
                                       vector<VarFcnBase*> &varFcn_, SpaceOperator &spo, vector<LevelSetOperator*> &lso)
                  : comm(comm_), iod(iod_), varFcn(varFcn_),
                    coordinates(spo.GetMeshCoordinates()),
                    delta_xyz(spo.GetMeshDeltaXYZ()),
                    Tag(comm_, &(dm_all_.ghosted1_1dof))
{
  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);

  for(int i=0; i<lso.size(); i++)
    ls2matid[i] = lso[i]->GetMaterialID();
}

//-----------------------------------------------------

MultiPhaseOperator::~MultiPhaseOperator()
{ }

//-----------------------------------------------------

void
MultiPhaseOperator::Destroy()
{
  Tag.Destroy();
}

//-----------------------------------------------------

void 
MultiPhaseOperator::UpdateMaterialID(vector<SpaceVariable3D*> &Phi, SpaceVariable3D &ID)
{
  // reset tag to 0
  Tag.SetConstantValue(0, true/*workOnGhost*/);
  ID.SetConstantValue(0, true/*workOnGhost*/);
  int overlap = 0;

  double*** tag = (double***)Tag.GetDataPointer();
  double*** id  = (double***)ID.GetDataPointer();

  for(int ls = 0; ls<Phi.size(); ls++) {//loop through all the level set functions
  
    double*** phi = (double***)Phi[ls]->GetDataPointer();
    int matid = ls2matid[ls];

    for(int k=kk0; k<kkmax; k++)
      for(int j=jj0; j<jjmax; j++)
        for(int i=ii0; i<iimax; i++) {
          if(phi[k][j][i]<0) {
            if(id[k][j][i] != 0) {
              overlap++;
              tag[k][j][i] = 1;
            } 
            id[k][j][i] = matid; 
          }

        }

    Phi[ls]->RestoreDataPointerToLocalVector(); //no changes made
  }

  MPI_Allreduce(MPI_IN_PLACE, &overlap, 1, MPI_INT, MPI_SUM, comm);


  if(overlap) {
    print_error("*** Error: Found overlapping material interfaces. Number of overlapped cells: %d.\n", overlap);
    exit_mpi();
  } 


  if(overlap) 
    Tag.RestoreDataPointerAndInsert();
  else
    Tag.RestoreDataPointerToLocalVector();
}

//-----------------------------------------------------
// Section 4.2.4 of Arthur Rallu's thesis
void
MultiPhaseOperator::UpdateStateVariablesAfterInterfaceMotion(SpaceVariable3D &IDn, 
                        SpaceVariable3D &ID, SpaceVariable3D &V, RiemannSolutions &riemann_solutions)
{
  switch (iod.multiphase.phasechange_type) {

    case MultiPhaseData::RIEMANN_SOLUTION :
      UpdateStateVariablesByRiemannSolutions(IDn, ID, V, riemann_solutions);
      break;

    case MultiPhaseData::EXTRAPOLATION :
      UpdateStateVariablesByExtrapolation(IDn, ID, V);
      break;

    default :
      print_error("*** Error: Specified method for phase-change update (%d) has not been implemented.\n", 
                  (int)iod.multiphase.phasechange_type);
  }

}

//-----------------------------------------------------

void
MultiPhaseOperator::UpdateStateVariablesByRiemannSolutions(SpaceVariable3D &IDn, 
                        SpaceVariable3D &ID, SpaceVariable3D &V, RiemannSolutions &riemann_solutions)
{
  // extract info
  double*** idn = (double***)IDn.GetDataPointer();
  double*** id  = (double***)ID.GetDataPointer();
  Vec5D***  v   = (Vec5D***) V.GetDataPointer();

  // create a vector that temporarily stores unresolved nodes (which will be resolved separately)
  vector<Int3> unresolved;

  // work inside the real domain
  int counter = 0;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        if(id[k][j][i] == idn[k][j][i]) //id remains the same. Skip
          continue;

        counter = LocalUpdateByRiemannSolutions(i, j, k, id[k][j][i], v[k][j][i-1], v[k][j][i+1], 
                      v[k][j-1][i], v[k][j+1][i], v[k-1][j][i], v[k+1][j][i], riemann_solutions,
                      v[k][j][i], true);
        if(counter==0)
          counter = LocalUpdateByRiemannSolutions(i, j, k, id[k][j][i], v[k][j][i-1], v[k][j][i+1], 
              v[k][j-1][i], v[k][j+1][i], v[k-1][j][i], v[k+1][j][i], riemann_solutions,
              v[k][j][i], false);

        if(counter==0) //add it to unresolved nodes...
          unresolved.push_back(Int3(k,j,i)); //note the order: k,j,i

      }  


  V.RestoreDataPointerAndInsert(); //insert data & communicate with neighbor subd's
  ID.RestoreDataPointerToLocalVector();
  IDn.RestoreDataPointerToLocalVector();

  // Fix the unresolved nodes (if any)
  int nUnresolved = unresolved.size();
  MPI_Allreduce(MPI_IN_PLACE, &nUnresolved, 1, MPI_INT, MPI_SUM, comm);
  if(nUnresolved) //some of the subdomains have unresolved nodes
    FixUnresolvedNodes(unresolved, IDn, ID, V); 

} 

//-----------------------------------------------------

int
MultiPhaseOperator::LocalUpdateByRiemannSolutions(int i, int j, int k, int id, Vec5D &vl, Vec5D &vr, 
                        Vec5D &vb, Vec5D &vt, Vec5D &vk, Vec5D &vf, RiemannSolutions &riemann_solutions, 
                        Vec5D &v, bool upwind)
{
  int counter = 0;
  double weight = 0, sum_weight = 0;
  Int3 ind(k,j,i);

  // left
  auto it = riemann_solutions.left.find(ind);
  if(it != riemann_solutions.left.end()) {
    if(it->second.second/*ID*/ == id && (!upwind || vl[1] > 0)) {
      Vec3D v1(vl[1], vl[2], vl[3]);
      weight = upwind ? vl[1]/v1.norm() : 1.0;
      sum_weight += weight;
      if(counter==0) 
        v = weight*it->second.first; /*riemann solution*/
      else 
        v += weight*it->second.first; /*riemann solution*/
      counter++;
    }
  }

  // right
  it = riemann_solutions.right.find(ind);
  if(it != riemann_solutions.right.end()) {
    if(it->second.second/*ID*/ == id && (!upwind || vr[1] < 0)) {
      Vec3D v1(vr[1], vr[2], vr[3]);
      weight = upwind ? -vr[1]/v1.norm() : 1.0;
      sum_weight += weight;
      if(counter==0)
        v = weight*it->second.first; /*riemann solution*/
      else
        v += weight*it->second.first; /*riemann solution*/
      counter++;
    }
  }

  // bottom
  it = riemann_solutions.bottom.find(ind);
  if(it != riemann_solutions.bottom.end()) {
    if(it->second.second/*ID*/ == id && (!upwind || vb[2] > 0)) {
      Vec3D v1(vb[1], vb[2], vb[3]);
      weight = upwind ? vb[2]/v1.norm() : 1.0;
      sum_weight += weight;
      if(counter==0)
        v = weight*it->second.first; /*riemann solution*/
      else
        v += weight*it->second.first; /*riemann solution*/
      counter++;
    }
  }

  // top
  it = riemann_solutions.top.find(ind);
  if(it != riemann_solutions.top.end()) {
    if(it->second.second/*ID*/ == id && (!upwind || vt[2] < 0)) {
      Vec3D v1(vt[1], vt[2], vt[3]);
      weight = upwind ? -vt[2]/v1.norm() : 1.0;
      sum_weight += weight;
      if(counter==0)
        v = weight*it->second.first; /*riemann solution*/
      else
        v += weight*it->second.first; /*riemann solution*/
      counter++;
    }
  }

  // back
  it = riemann_solutions.back.find(ind);
  if(it != riemann_solutions.back.end()) {
    if(it->second.second/*ID*/ == id && (!upwind || vk[3] > 0)) {
      Vec3D v1(vk[1], vk[2], vk[3]);
      weight = upwind ? vk[3]/v1.norm() : 1.0;
      sum_weight += weight;
      if(counter==0)
        v = weight*it->second.first; /*riemann solution*/
      else
        v += weight*it->second.first; /*riemann solution*/
      counter++;
    }
  }

  // front
  it = riemann_solutions.front.find(ind);
  if(it != riemann_solutions.front.end()) {
    if(it->second.second/*ID*/ == id && (!upwind || vf[3] < 0)) {
      Vec3D v1(vf[1], vf[2], vf[3]);
      weight = upwind ? -vf[3]/v1.norm() : 1.0;
      sum_weight += weight;
      if(counter==0)
        v = weight*it->second.first; /*riemann solution*/
      else
        v += weight*it->second.first; /*riemann solution*/
      counter++;
    }
  }

  if(sum_weight > 0.0)
    v /= sum_weight;
  else if(upwind) {
    if(verbose>1) 
      fprintf(stderr,"Warning: Unable to update phase change at (%d,%d,%d) by Riemann solutions w/ upwinding. Retrying.\n", i,j,k);
  } else {
    if(verbose>1) 
      fprintf(stderr,"Warning: Unable to update phase change at (%d,%d,%d) by Riemann solutions. Retrying.\n", i,j,k);
  }

  return counter;
}

//-----------------------------------------------------

void
MultiPhaseOperator::UpdateStateVariablesByExtrapolation(SpaceVariable3D &IDn, 
                        SpaceVariable3D &ID, SpaceVariable3D &V)
{
  // extract info
  double*** idn = (double***)IDn.GetDataPointer();
  double*** id  = (double***)ID.GetDataPointer();
  Vec5D***  v   = (Vec5D***) V.GetDataPointer();

  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();

  double weight, sum_weight;
  Vec3D v1, x1x0;
  double v1norm;

  // create a vector that temporarily stores unresolved nodes (which will be resolved separately)
  vector<Int3> unresolved;

  // work inside the real domain
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        if(id[k][j][i] == idn[k][j][i]) //id remains the same. Skip
          continue;

        // coordinates of this node
        Vec3D& x0(coords[k][j][i]);

        sum_weight = 0.0;

        bool reset = false; //whether v[k][j][i] has been reset to 0

        //go over the neighboring nodes 
        for(int neighk = k-1; neighk <= k+1; neighk++)         
          for(int neighj = j-1; neighj <= j+1; neighj++)
            for(int neighi = i-1; neighi <= i+1; neighi++) {

              if(id[neighk][neighj][neighi] != id[k][j][i])
                continue; //this neighbor has a different ID. Skip it.

              if(id[neighk][neighj][neighi] != idn[neighk][neighj][neighi])
                continue; //this neighbor also changed ID. Skip it. (Also skipping node [k][j][i])

              if(ID.OutsidePhysicalDomain(neighi, neighj, neighk))
                continue; //this neighbor is outside the physical domain. Skip.

              // coordinates and velocity at the neighbor node
              Vec3D& x1(coords[neighk][neighj][neighi]);
              v1[0] = v[neighk][neighj][neighi][1];
              v1[1] = v[neighk][neighj][neighi][2];
              v1[2] = v[neighk][neighj][neighi][3];

              // compute weight
              v1norm = v1.norm();
              if(v1norm != 0)
                v1 /= v1norm;
              x1x0 = x0 - x1; 
              x1x0 /= x1x0.norm();

              weight = max(0.0, x1x0*v1);

              // add weighted s.v. at neighbor node
              if(weight>0) {
                sum_weight += weight;
                if(reset)
                  v[k][j][i] += weight*v[neighk][neighj][neighi];
                else {
                  v[k][j][i] = weight*v[neighk][neighj][neighi];
                  reset = true;
                }
              }
            }

        if(sum_weight==0) {
          if(verbose>1) 
            fprintf(stderr,"Warning: Unable to update phase change at (%d,%d,%d)(%e,%e,%e) "
                    "by extrapolation w/ upwinding.\n", i,j,k, x0[0],x0[1],x0[2]);
          unresolved.push_back(Int3(k,j,i)); //note the order: k,j,i          
        } else
          v[k][j][i] /= sum_weight; 
      }

  V.RestoreDataPointerAndInsert(); //insert data & communicate with neighbor subd's
  ID.RestoreDataPointerToLocalVector();
  IDn.RestoreDataPointerToLocalVector();
  coordinates.RestoreDataPointerToLocalVector();

  // Fix the unresolved nodes (if any)
  int nUnresolved = unresolved.size();
  MPI_Allreduce(MPI_IN_PLACE, &nUnresolved, 1, MPI_INT, MPI_SUM, comm);
  if(nUnresolved) //some of the subdomains have unresolved nodes
    FixUnresolvedNodes(unresolved, IDn, ID, V); 

}

//-----------------------------------------------------

void MultiPhaseOperator::FixUnresolvedNodes(vector<Int3> &unresolved, SpaceVariable3D &IDn, SpaceVariable3D &ID,
                                            SpaceVariable3D &V)
{
  // Note: all the processor cores will enter this function even if only one or a few have unresolved nodes

  // extract info
  double*** idn = (double***)IDn.GetDataPointer();
  double*** id  = (double***)ID.GetDataPointer();
  Vec5D***  v   = (Vec5D***) V.GetDataPointer();

  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();


  // loop through unresolved nodes
  int i,j,k;
  double weight, sum_weight = 0.0; 
  Vec3D v1, x1x0;
  double v1norm;

  bool reset = false; //whether v[k][j][i] has been reset

  for(auto it = unresolved.begin(); it != unresolved.end(); it++) { 

    k = (*it)[0];
    j = (*it)[1];
    i = (*it)[2];

    Vec3D& x0(coords[k][j][i]);
    sum_weight = 0.0;

    //go over the neighboring nodes 
    for(int neighk = k-1; neighk <= k+1; neighk++)         
      for(int neighj = j-1; neighj <= j+1; neighj++)
        for(int neighi = i-1; neighi <= i+1; neighi++) {

          if(ID.OutsidePhysicalDomain(neighi, neighj, neighk))
            continue; //this neighbor is outside the physical domain. Skip.

          if(id[neighk][neighj][neighi] != id[k][j][i])
            continue; //this neighbor has a different ID. Skip it.

          if(neighk==k && neighj==j && neighi==i) 
            continue; //the same node. Skip

          bool neighbor_unresolved = false;
          for(auto it2 = unresolved.begin(); it2 != unresolved.end(); it2++) {
            if(neighk==(*it2)[0] && neighj==(*it2)[1] && neighi==(*it2)[2]) {
              neighbor_unresolved = true; 
              break; 
            }
          }
          if(neighbor_unresolved)
            continue; //this neighbor is also unresolved. Skip

          // coordinates and velocity at the neighbor node
          Vec3D& x1(coords[neighk][neighj][neighi]);
          v1[0] = v[neighk][neighj][neighi][1];
          v1[1] = v[neighk][neighj][neighi][2];
          v1[2] = v[neighk][neighj][neighi][3];

          // compute weight
          v1norm = v1.norm();
          if(v1norm != 0)
            v1 /= v1norm;
          x1x0 = x0 - x1; 
          x1x0 /= x1x0.norm();

          weight = max(0.0, x1x0*v1);

          // add weighted s.v. at neighbor node
          if(weight>0) {
            sum_weight += weight;
            if(reset)
              v[k][j][i] += weight*v[neighk][neighj][neighi];
            else {
              v[k][j][i] = weight*v[neighk][neighj][neighi];
              reset = true;
            }
          }
        }


    if(sum_weight>0) {
      v[k][j][i] /= sum_weight; //Done!
      if(verbose>1) fprintf(stderr,"*** (%d,%d,%d): Updated state variables by extrapolation w/ upwinding. (2nd attempt)\n",
                          i,j,k);
      continue;
    }


    // Our last resort: keep the pressure and velocity (both normal & TANGENTIAL) at the current node, 
    // find a valid density nearby. (In this case, the solution may be different for different domain
    // partitions --- but this should rarely happen.)
          
    //go over the neighboring nodes & interpolate velocity and pressure
    int max_layer = 10;
    double density = 0.0;
    for(int layer = 1; layer <= max_layer; layer++) {

      for(int neighk = k-layer; neighk <= k+layer; neighk++)         
        for(int neighj = j-layer; neighj <= j+layer; neighj++)
          for(int neighi = i-layer; neighi <= i+layer; neighi++) {

            if(ID.OutsidePhysicalDomain(neighi, neighj, neighk))
              continue; //this neighbor is outside the physical domain. Skip.
  
            if(!ID.IsHere(neighi,neighj,neighk,true/*include_ghost*/))
              continue; //this neighbor is outside the current subdomain

            if(id[neighk][neighj][neighi] != id[k][j][i])
              continue; //this neighbor has a different ID. Skip it.

            if(neighk==k && neighj==j && neighi==i) 
              continue; //the same node. Skip

            bool neighbor_unresolved = false;
            for(auto it2 = unresolved.begin(); it2 != unresolved.end(); it2++) {
              if(neighk==(*it2)[0] && neighj==(*it2)[1] && neighi==(*it2)[2]) {
                neighbor_unresolved = true; 
                break; 
              }
            }
            if(neighbor_unresolved)
              continue; //this neighbor is also unresolved. Skip


            Vec3D& x1(coords[neighk][neighj][neighi]);
            double dist = (x1-x0).norm();

            sum_weight += dist;
            density += dist*v[neighk][neighj][neighi][0];
          }

      if(sum_weight>0) {
        v[k][j][i][0] = density/sum_weight;
        if(verbose>1) fprintf(stderr,"*** (%d,%d,%d): Updated density by interpolation w/ stencil width = %d: %e %e %e %e %e\n",
                            i,j,k, layer, v[k][j][i][0], v[k][j][i][1], v[k][j][i][2], v[k][j][i][3], v[k][j][i][4]);
        break; //done with this node
      }

    }

    //Very unlikely. This means there is no neighbors within max_layer that have the same ID!!
    if(sum_weight==0) { 
      fprintf(stderr,"\033[0;35mWarning: Updating phase change at (%d,%d,%d)(%e,%e,%e) with pre-specified density (%e). No valid "
                     "neighbors within %d layers.\033[0m\n", 
                     i,j,k, coords[k][j][i][0], coords[k][j][i][1],
                     coords[k][j][i][2], varFcn[id[k][j][i]]->failsafe_density, max_layer);
      v[k][j][i][0] = varFcn[id[k][j][i]]->failsafe_density;
    }
  }

  V.RestoreDataPointerAndInsert(); //insert data & communicate with neighbor subd's

  ID.RestoreDataPointerToLocalVector();
  IDn.RestoreDataPointerToLocalVector();
  coordinates.RestoreDataPointerToLocalVector();

}

//-----------------------------------------------------

