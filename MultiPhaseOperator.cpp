#include<MultiPhaseOperator.h>
#include<SpaceOperator.h>
#include<LevelSetOperator.h>
#include<Vector5D.h>
#include<RiemannSolutions.h>
#include<algorithm>//find
#include<set>
using std::cout;
using std::endl;
using std::set;
using std::tuple;
using std::get;
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

  // Initialize phase/material transition functions (if specified by user)
  if(iod.eqs.transitions.dataMap.size()) {

    // start with creating an empty vector for each material ID
    int numMaterials = varFcn.size();
    trans.resize(numMaterials);

    // fill in the user-specified transitions
    for(auto it = iod.eqs.transitions.dataMap.begin(); it != iod.eqs.transitions.dataMap.end(); it++) {
      int i = it->second->from_id;
      int j = it->second->to_id;
      if(i<0 || i>=varFcn.size() || j<0 || j>=varFcn.size() || i==j) {
        print_error("*** Error: Detected input error in Material/Phase Transition [%d] (%d -> %d).\n", it->first,
                    i, j); 
        exit_mpi();
      }
      if(true) //no other choices at the moment
        trans[i].push_back(new PhaseTransitionBase(*it->second, *varFcn[i], *varFcn[j]));
      else {
        print_error("*** Error: Unknown phase transition type (%d).\n");
        exit_mpi();
      }
    }

    // make sure the needed level set functions are available
    for(auto it = iod.eqs.transitions.dataMap.begin(); it != iod.eqs.transitions.dataMap.end(); it++) {
      int i = it->second->from_id;
      int j = it->second->to_id;
      if(i!=0) {
        bool found = false;
        for(int ls = 0; ls < lso.size(); ls++) {
          if(ls2matid[ls] == i) {
            found = true;
            break;
          }
        }
        if(!found) {
          print_error("*** Error: Phase transitions involve material ID %d, but a level set solver is not specified.\n",
                      i);
          exit_mpi();
        }
      }
      if(j!=0) {
        bool found = false;
        for(int ls = 0; ls < lso.size(); ls++) {
          if(ls2matid[ls] == j) {
            found = true;
            break;
          }
        }
        if(!found) {
          print_error("*** Error: Phase transitions involve material ID %d, but a level set solver is not specified.\n",
                      j);
          exit_mpi();
        }
      }
    }

  } else
    trans.resize(0);

}

//-----------------------------------------------------

MultiPhaseOperator::~MultiPhaseOperator()
{ }

//-----------------------------------------------------

void
MultiPhaseOperator::Destroy()
{
  Tag.Destroy();

  for(int i=0; i<trans.size(); i++)
    for(auto it = trans[i].begin(); it != trans[i].end(); it++)
      delete *it;
}

//-----------------------------------------------------

void 
MultiPhaseOperator::UpdateMaterialID(vector<SpaceVariable3D*> &Phi, SpaceVariable3D &ID)
{

#ifdef LEVELSET_TEST
  return; //testing the level set solver w/o solving the N-S / Euler equations
#endif

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

//-------------------------------------------------------------------------
// This function checks for physical phase transitions based on varFcn. If
// found, the levelset (phi), the material ID, and possibly also the state
// variables V will be updated. The function returns the total number of
// nodes undergoing phase transitions. If it is non-zero, all the level set
// functions should be reinitialized (done outside of MultiPhaseOperator)
int 
MultiPhaseOperator::UpdatePhaseTransitions(vector<SpaceVariable3D*> &Phi, SpaceVariable3D &ID, 
                                           SpaceVariable3D &V, vector<int> &phi_updated, 
                                           vector<Int3> *new_useful_nodes)
{
  if(trans.size()==0)
    return 0; //nothing to do

  //---------------------------------------------
  // Step 1: Check for phase transitions; Update
  //         ID and V
  //---------------------------------------------
  double*** id  = ID.GetDataPointer();
  Vec5D***  v   = (Vec5D***) V.GetDataPointer();

  int counter = 0;
  vector<tuple<Int3,int,int> >  changed; //(i,j,k), old id, new id -- incl. ghosts inside physical domain

  set<int> affected_ids;

  //DEBUG!!
/*
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        //if(coords[k][j][i].norm()<0.02) {
        if(fabs(coords[k][j][i][0]-0.01)<0.005 && fabs(coords[k][j][i][1]-0.01)<0.005) {
          fprintf(stderr,"Changing state at (%d,%d,%d) (%e, %e, %e).\n", i,j,k, 
                  coords[k][j][i][0], coords[k][j][i][1], coords[k][j][i][2]);
          v[k][j][i][0] = 8.9e-4;
          v[k][j][i][4] = 1.2e10;
          int myid = id[k][j][i]; 
          double e = varFcn[myid]->GetInternalEnergyPerUnitMass(v[k][j][i][0], v[k][j][i][4]);
          double T = varFcn[myid]->GetTemperature(v[k][j][i][0], e);
          double h = varFcn[myid]->ComputeEnthalpyPerUnitMass(v[k][j][i][0], v[k][j][i][4]);
          fprintf(stderr," ...  e = %e, h = %e, T = %e.\n", e, h, T);
        }
      }
  coordinates.RestoreDataPointerToLocalVector();
  V.RestoreDataPointerAndInsert();
  v = (Vec5D***) V.GetDataPointer();
*/

  int myid;
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++) {
        myid = (int)id[k][j][i];
        // skip nodes outside the physical domain
        if(coordinates.OutsidePhysicalDomain(i,j,k))
          continue;

        for(auto it = trans[myid].begin(); it != trans[myid].end(); it++) {
          if((*it)->Transition(v[k][j][i])) { //detected phase transition

            // register the node
            changed.push_back(std::make_tuple(Int3(i,j,k), myid, (*it)->toID));

            // register the involved material ids
            affected_ids.insert(myid);
            affected_ids.insert((*it)->toID);

            // update id
            id[k][j][i] = (*it)->toID;

            // update p (rho and e fixed)
            double e = varFcn[myid]->GetInternalEnergyPerUnitMass(v[k][j][i][0], v[k][j][i][4]);
            v[k][j][i][4] = varFcn[(*it)->toID]->GetPressure(v[k][j][i][0], e);

            counter++;
            break;
          }
        }
      }

  MPI_Allreduce(MPI_IN_PLACE, &counter, 1, MPI_INT, MPI_SUM, comm);

  if(counter>0) {
    ID.RestoreDataPointerAndInsert();
    V.RestoreDataPointerAndInsert();
  } else {
    ID.RestoreDataPointerToLocalVector();
    V.RestoreDataPointerToLocalVector();
    return 0;
  }


  //---------------------------------------------
  // Step 2: Figure out which level set functions
  //         need to be updated
  //---------------------------------------------
  for(int ls = 0; ls < Phi.size(); ls++)
    phi_updated[ls] = (affected_ids.find(ls2matid[ls]) != affected_ids.end());
  MPI_Allreduce(MPI_IN_PLACE, (int*)phi_updated.data(), Phi.size(), MPI_INT, MPI_MAX, comm);

    
  //---------------------------------------------
  // Step 3: Update Phi
  //---------------------------------------------
  UpdatePhiAfterPhaseTransitions(Phi, ID, changed, phi_updated, new_useful_nodes);


  if(verbose>=1)
    print("- Detected phase/material transitions at %d node(s).\n", counter);

  return counter;

}

//-----------------------------------------------------

void
MultiPhaseOperator::UpdatePhiAfterPhaseTransitions(vector<SpaceVariable3D*> &Phi, SpaceVariable3D &ID, 
                                                   vector<tuple<Int3,int,int> > &changed, 
                                                   vector<int> &phi_updated, vector<Int3> *new_useful_nodes)
{
  // This function will provide correct value of phi (up to dx error) ONLY for first layer nodes.
  // Reinitialization is needed to find the value of phi elsewhere
  // Note that "changed" should include ghost nodes inside the physical domain

  int NX, NY, NZ;
  coordinates.GetGlobalSize(&NX, &NY, &NZ);

  Vec3D*** dxyz = (Vec3D***)delta_xyz.GetDataPointer();
  double*** id  = ID.GetDataPointer();

  for(int ls = 0; ls<Phi.size(); ls++) {//loop through all the level set functions
  
    if(phi_updated[ls] == 0)
      continue; //this level set function is not involved

    double*** phi = Phi[ls]->GetDataPointer();
    int matid = ls2matid[ls];

    int i,j,k;
    for(auto it = changed.begin(); it != changed.end(); it++) {
 
      if(matid != get<1>(*it) && matid != get<2>(*it))
        continue;

      i = get<0>(*it)[0];
      j = get<0>(*it)[1];
      k = get<0>(*it)[2];

      //first, push new nodes into the vector of new_useful_nodes
      new_useful_nodes[ls].push_back(get<0>(*it));
      if(i-1>=ii0)  new_useful_nodes[ls].push_back(Int3(i-1,j,k));
      if(i+1<iimax) new_useful_nodes[ls].push_back(Int3(i+1,j,k));
      if(j-1>=jj0)  new_useful_nodes[ls].push_back(Int3(i,j-1,k));
      if(j+1<jjmax) new_useful_nodes[ls].push_back(Int3(i,j+1,k));
      if(k-1>=kk0)  new_useful_nodes[ls].push_back(Int3(i,j,k-1));
      if(k+1<kkmax) new_useful_nodes[ls].push_back(Int3(i,j,k+1));


      if(matid == get<1>(*it)) { //this node is moving outside of the subdomain

        // update phi at this node
        phi[k][j][i] = 0.5*std::min(dxyz[k][j][i][0], std::min(dxyz[k][j][i][1],dxyz[k][j][i][2]));        

        // update phi at neighbors that are in the physical domain, and have opposite sign 
        if(i-1>=ii0  && i-1>=0 && phi[k][j][i-1]<=0)
          phi[k][j][i-1] = std::max(phi[k][j][i-1], -0.5*dxyz[k][j][i-1][0]);
        if(i+1<iimax && i+1<NX && phi[k][j][i+1]<=0)
          phi[k][j][i+1] = std::max(phi[k][j][i+1], -0.5*dxyz[k][j][i+1][0]);
        if(j-1>=jj0  && j-1>=0 && phi[k][j-1][i]<=0) 
          phi[k][j-1][i] = std::max(phi[k][j-1][i], -0.5*dxyz[k][j-1][i][1]);
        if(j+1<jjmax && j+1<NY && phi[k][j+1][i]<=0)
          phi[k][j+1][i] = std::max(phi[k][j+1][i], -0.5*dxyz[k][j+1][i][1]);
        if(k-1>=kk0  && k-1>=0 && phi[k-1][j][i]<=0) 
          phi[k-1][j][i] = std::max(phi[k-1][j][i], -0.5*dxyz[k-1][j][i][2]);
        if(k+1<kkmax && k+1<NZ && phi[k+1][j][i]<=0)
          phi[k+1][j][i] = std::max(phi[k+1][j][i], -0.5*dxyz[k+1][j][i][2]);
      }
      else if(matid == get<2>(*it)) { //this node is moving inside the subdomain

        // update phi at this node
        phi[k][j][i] = -0.5*std::min(dxyz[k][j][i][0], std::min(dxyz[k][j][i][1],dxyz[k][j][i][2]));        

        // update phi at neighbors that are in the physical domain, and have opposite sign 
        if(i-1>=ii0  && i-1>=0 && phi[k][j][i-1]>=0)
          phi[k][j][i-1] = std::min(phi[k][j][i-1], 0.5*dxyz[k][j][i-1][0]);
        if(i+1<iimax && i+1<NX && phi[k][j][i+1]>=0)
          phi[k][j][i+1] = std::min(phi[k][j][i+1], 0.5*dxyz[k][j][i+1][0]);
        if(j-1>=jj0  && j-1>=0 && phi[k][j-1][i]>=0) 
          phi[k][j-1][i] = std::min(phi[k][j-1][i], 0.5*dxyz[k][j-1][i][1]);
        if(j+1<jjmax && j+1<NY && phi[k][j+1][i]>=0)
          phi[k][j+1][i] = std::min(phi[k][j+1][i], 0.5*dxyz[k][j+1][i][1]);
        if(k-1>=kk0  && k-1>=0 && phi[k-1][j][i]>=0) 
          phi[k-1][j][i] = std::min(phi[k-1][j][i], 0.5*dxyz[k-1][j][i][2]);
        if(k+1<kkmax && k+1<NZ && phi[k+1][j][i]>=0)
          phi[k+1][j][i] = std::min(phi[k+1][j][i], 0.5*dxyz[k+1][j][i][2]);
      }

    }

    Phi[ls]->RestoreDataPointerAndInsert();
  }

  // For debug only
  for(int ls = 0; ls<Phi.size(); ls++) {//loop through all the level set functions
  
    if(phi_updated[ls] == 0)
      continue; //this level set function is not involved

    double*** phi = Phi[ls]->GetDataPointer();
    int matid = ls2matid[ls];

    int i,j,k;
    for(auto it = changed.begin(); it != changed.end(); it++) {
      if(matid == get<1>(*it) || matid == get<2>(*it)) { 
 
        i = get<0>(*it)[0];
        j = get<0>(*it)[1];
        k = get<0>(*it)[2];

        if(phi[k][j][i]<0)
          assert(id[k][j][i] == matid);
        else
          assert(id[k][j][i] != matid);

        if(i-1>=ii0 && i-1>=0) {
          if(phi[k][j][i-1]<0)
            assert(id[k][j][i-1] == matid);
          else
            assert(id[k][j][i-1] != matid);
        }

        if(i+1<iimax && i+1<NX) {
          if(phi[k][j][i+1]<0)
            assert(id[k][j][i+1] == matid);
          else
            assert(id[k][j][i+1] != matid);
        }

        if(j-1>=jj0 && j-1>=0) {
          if(phi[k][j-1][i]<0)
            assert(id[k][j-1][i] == matid);
          else
            assert(id[k][j-1][i] != matid);
        }

        if(j+1<jjmax && j+1<NY) {
          if(phi[k][j+1][i]<0)
            assert(id[k][j+1][i] == matid);
          else
            assert(id[k][j+1][i] != matid);
        }

        if(k-1>=kk0 && k-1>=0) {
          if(phi[k-1][j][i]<0)
            assert(id[k-1][j][i] == matid);
          else
            assert(id[k-1][j][i] != matid);
        }

        if(k+1<kkmax && k+1<NZ) {
          if(phi[k+1][j][i]<0) 
            assert(id[k+1][j][i] == matid);
          else
            assert(id[k+1][j][i] != matid);
        }

      }
    }
    Phi[ls]->RestoreDataPointerToLocalVector();
  }

  delta_xyz.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();

}

//-----------------------------------------------------

