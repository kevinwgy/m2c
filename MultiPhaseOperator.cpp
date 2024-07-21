/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<MultiPhaseOperator.h>
#include<SpaceOperator.h>
#include<LevelSetOperator.h>
#include<Vector5D.h>
#include<RiemannSolutions.h>
#include<Intersector.h>
#include<algorithm>//find
#include<set>
using std::cout;
using std::endl;
using std::set;
using std::tuple;
using std::get;
using std::unique_ptr;
extern int verbose;
extern int INACTIVE_MATERIAL_ID;
//-----------------------------------------------------

MultiPhaseOperator::MultiPhaseOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_,
                                       vector<VarFcnBase*> &varFcn_, GlobalMeshInfo &global_mesh_,
                                       SpaceOperator &spo, vector<LevelSetOperator*> &lso)
                  : comm(comm_), iod(iod_), varFcn(varFcn_), global_mesh(global_mesh_),
                    coordinates(spo.GetMeshCoordinates()),
                    delta_xyz(spo.GetMeshDeltaXYZ()),
                    ghost_nodes_inner(spo.GetPointerToInnerGhostNodes()), 
                    ghost_nodes_outer(spo.GetPointerToOuterGhostNodes()),
                    Tag(comm_, &(dm_all_.ghosted1_1dof)),
                    Lambda(comm_, &(dm_all_.ghosted1_1dof)),
                    volume(spo.GetMeshCellVolumes())
{

  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  coordinates.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);

  for(int i=0; i<(int)lso.size(); i++)
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
      if(i<0 || i>=(int)varFcn.size() || j<0 || j>=(int)varFcn.size() || i==j) {
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
        for(int ls = 0; ls < (int)lso.size(); ls++) {
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
        for(int ls = 0; ls < (int)lso.size(); ls++) {
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

  Lambda.Destroy();

  for(int i=0; i<(int)trans.size(); i++)
    for(auto it = trans[i].begin(); it != trans[i].end(); it++)
      delete *it;
}

//-----------------------------------------------------

void 
MultiPhaseOperator::UpdateMaterialIDByLevelSet(vector<SpaceVariable3D*> &Phi0, vector<SpaceVariable3D*> &Phi,
                                               vector<Intersector*> *intersector, SpaceVariable3D &ID)
{

#ifdef LEVELSET_TEST
  return; //testing the level set solver w/o solving the N-S / Euler equations
#endif


  double*** tag = (double***)Tag.GetDataPointer();
  double*** id  = (double***)ID.GetDataPointer();

  int ls_size = Phi.size();
  vector<double***> phi(ls_size, NULL);
  for(int ls=0; ls<ls_size; ls++) 
    phi[ls] = Phi[ls]->GetDataPointer();
  
  vector<double***> phi0(ls_size, NULL);
  for(int ls=0; ls<ls_size; ls++) 
    phi0[ls] = Phi0[ls]->GetDataPointer();
  

  int myls(-1);
  int total_swept_nodes = 0;
  set<std::pair<int,int> > swept; //pairs "ls" with -1 = phi[ls]<0 / 0 = phi[ls]=0 / 1 = phi[ls] > 0
  set<Int3> remaining_nodes;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        tag[k][j][i] = 0;

        if(id[k][j][i] == INACTIVE_MATERIAL_ID)
          continue; //if occluded, id should not be changed

        swept.clear();
        for(int ls=0; ls<ls_size; ls++) {
          if(phi[ls][k][j][i]*phi0[ls][k][j][i]>0.0)
            continue; //no change
          else {
            if     (phi[ls][k][j][i]<0)  swept.insert(std::make_pair(ls, -1));
            else if(phi[ls][k][j][i]==0) swept.insert(std::make_pair(ls, 0));
            else                         swept.insert(std::make_pair(ls, 1));
          }
        }

        if(swept.empty())
          continue; // no sign change: Should not update ID, particularly important if there are also embedded surfaces
                      // that also influence material ID

        total_swept_nodes++;

        myls = -1;
        for(auto&& ls_status : swept) {
          if(ls_status.second == -1) {
            if(myls!=-1) {
              fprintf(stdout,"\033[0;31m*** Error: Node (%d,%d,%d) belongs to two material subdomains. "
                             "phi[%d(matid:%d)] = %e, phi[%d(matid:%d)] = %e.\033[0m\n", i,j,k, myls, ls2matid[myls],
                             phi[myls][k][j][i], ls_status.first, ls2matid[ls_status.first],
                             phi[ls_status.first][k][j][i]); 
              exit(-1);
            }
            myls = ls_status.first;
          } 
        }

     //   fprintf(stdout,"(%d,%d,%d)(%e,%e,%e) is swept, phi0: %e -> %e, phi1: %e -> %e. myls = %d\n", i,j,k,
     //           global_mesh.GetX(i), global_mesh.GetY(j), global_mesh.GetZ(k),
     //           phi0[0][k][j][i], phi[0][k][j][i], phi0[1][k][j][i], phi[1][k][j][i], myls);

        if(myls != -1)
          id[k][j][i] = ls2matid[myls];
        else {// the node does not belong to any subdomain tracked by level set(s) ==> tag it
          tag[k][j][i] = 1;
          remaining_nodes.insert(Int3(i,j,k));
        }
      }

  int total_remaining = remaining_nodes.size();
  int progress(0);
  MPI_Allreduce(MPI_IN_PLACE, &total_remaining, 1, MPI_INT, MPI_SUM, comm);
  if(total_remaining>0)
    Tag.RestoreDataPointerAndInsert();
  else
    Tag.RestoreDataPointerToLocalVector();

  MPI_Allreduce(MPI_IN_PLACE, &total_swept_nodes, 1, MPI_INT, MPI_SUM, comm);
  if(total_swept_nodes>0)
    ID.RestoreDataPointerAndInsert();
  else
    ID.RestoreDataPointerToLocalVector();


  // Fix "remaining" nodes, which are swept nodes that do not belong to any subdomains tracked by level sets
  set<int> ids_tracked_by_levelsets;
  for(auto&& ls2id : ls2matid)
    ids_tracked_by_levelsets.insert(ls2id.second);

  for(int iter = 0; iter < 100; iter++) {// usually one iteration should be enough

    if(total_remaining==0)
      break; //done

    set<int> matid_cands;
     
    tag = Tag.GetDataPointer();
    id  = ID.GetDataPointer();

    progress = 0;

    set<Int3> tmp;
    tmp = remaining_nodes;
 
    int i,j,k, neighborid;
    for(auto&& ijk : tmp) {
      i = ijk[0];
      j = ijk[1];
      k = ijk[2];

      // check neighbors
      Vec3D x0(global_mesh.GetX(i), global_mesh.GetY(j), global_mesh.GetZ(k));
      matid_cands.clear();
      for(int neighk = k-1; neighk <= k+1; neighk++)         
        for(int neighj = j-1; neighj <= j+1; neighj++)
          for(int neighi = i-1; neighi <= i+1; neighi++) {
            if(ID.OutsidePhysicalDomain(neighi, neighj, neighk))
              continue; //this neighbor is outside the physical domain. Skip.
            if(tag[neighk][neighj][neighi]==1)
              continue; //this neighbor is undecided too

            neighborid = id[neighk][neighj][neighi];

            //fprintf(stdout,"(%d,%d,%d): neighbor (%d %d %d), id = %d.\n", i,j,k, neighi, neighj, neighk, neighborid);

            if(neighborid==INACTIVE_MATERIAL_ID)
              continue; //this neighbor is occluded
            if(ids_tracked_by_levelsets.find(neighborid) != ids_tracked_by_levelsets.end())
              continue; //this neighbor is inside a subdomain tracked by a level set --> not my material

            if(intersector) {
              Vec3D x1(global_mesh.GetX(neighi), global_mesh.GetY(neighj), global_mesh.GetZ(neighk));
              bool connected = true;
              for(auto&& xter : *intersector) {
                if(xter->Intersects(x0,x1)) {
                  connected = false;
                  break; //this neighbor is blocked to current node by an embedded surface
                }
              }
              if(!connected)
                continue;
            }

            matid_cands.insert(neighborid);
          }

      if(matid_cands.empty())
        continue; //nothing I can do
      else {
        id[k][j][i] = *(matid_cands.begin());  //take the first one, if multiple
        tag[k][j][i] = 0;
        remaining_nodes.erase(remaining_nodes.find(ijk));
        //fprintf(stdout,"  o Giving node (%d,%d,%d) ID %d.\n", i,j,k, (int)id[k][j][i]);
        progress++;

        if(matid_cands.size()>1)
          fprintf(stdout,"\033[0;35mWarning: Detected a node with multiple possible material IDs (%d,%d,%d), "
                  "setting to %d.\n\033[0m\n", i,j,k, (int)id[k][j][i]);
      }
    }


    MPI_Allreduce(MPI_IN_PLACE, &progress, 1, MPI_INT, MPI_SUM, comm);
    if(progress>0) {
      Tag.RestoreDataPointerAndInsert();
      ID.RestoreDataPointerAndInsert();
      total_remaining = remaining_nodes.size();
      MPI_Allreduce(MPI_IN_PLACE, &total_remaining, 1, MPI_INT, MPI_SUM, comm);
    } 
    else {
      // Likely some orphans. Material ID should be 0 (best guess)
      for(auto&& ijk : tmp)
        id[ijk[2]][ijk[1]][ijk[0]] = 0;

      Tag.RestoreDataPointerToLocalVector();
      ID.RestoreDataPointerToLocalVector();
      if(verbose>=1)
        print_warning("Warning: Found %d orphan nodes swept by interfaces tracked by level sets. Set ID = 0.\n",
                      total_remaining);
      break;
    }

  }


  // update external ghosts
  UpdateMaterialIDAtGhostNodes(ID);


  for(int ls=0; ls<ls_size; ls++) 
    Phi[ls]->RestoreDataPointerToLocalVector(); //no changes made

  for(int ls=0; ls<ls_size; ls++) 
    Phi0[ls]->RestoreDataPointerToLocalVector(); //no changes made

}

//-----------------------------------------------------

int
MultiPhaseOperator::ResolveConflictsWithEmbeddedSurfaces(vector<SpaceVariable3D*> &Phi, 
                        SpaceVariable3D &IDn, SpaceVariable3D &ID, 
                        vector<unique_ptr<EmbeddedBoundaryDataSet> > *EBDS, vector<Intersector*> *intersector)
{

  if(!EBDS || !intersector || Phi.empty())
    return 0; //nothing to do

  set<Int3> corrected;

/*
  IDn.StoreMeshCoordinates(coordinates);
  IDn.WriteToVTRFile("IDn.vtr", "IDn");
  ID.StoreMeshCoordinates(coordinates);
  ID.WriteToVTRFile("ID.vtr", "ID");
  print("I am here.\n");
  exit_mpi();
*/

  double*** idn = (double***)IDn.GetDataPointer();
  double*** id  = (double***)ID.GetDataPointer();

  vector<double***> phi; //may remain empty if there are no orphans
  vector<double***> phi_ebm; //..
  for(int ls=0; ls<(int)Phi.size(); ls++) 
    phi.push_back(Phi[ls]->GetDataPointer());
  for(auto&& ebds : *EBDS)
    phi_ebm.push_back(ebds->Phi_ptr->GetDataPointer()); //unsigned distance, always>=0


  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        if(idn[k][j][i] == id[k][j][i])
          continue;

        // check on nodes who have changed ID (due to level sets)
        
        if(IsOrphanAcrossEmbeddedSurfaces(i,j,k,idn,id,intersector)) {

          // update phi first (id will be switched back later)
          // first, we determine the absolute value of phi (checking both level sets and dist to surf)
          double newphi = DBL_MAX;
          for(int ls=0; ls<(int)phi.size(); ls++)
            newphi = std::min(newphi, fabs(phi[ls][k][j][i]));
          for(auto&& phis : phi_ebm)
            newphi = std::min(newphi, phis[k][j][i]);
          if(newphi==0)
            newphi = 1.0e-20; //avoid exactly 0

          for(int ls=0; ls<(int)phi.size(); ls++) {
            if(ls2matid[ls] == id[k][j][i]) { //this is the ID we are swtiching away from
              assert(phi[ls][k][j][i]<0.0); 
              phi[ls][k][j][i] = newphi; //make it positive
            }
            else if(ls2matid[ls] == idn[k][j][i]) {//this is the ID we are swtiching to
              assert(phi[ls][k][j][i]>=0.0);
              phi[ls][k][j][i] = -newphi; //make it negative
            }
            else 
              assert(phi[ls][k][j][i]>=0.0);
          }
 
          // switching back to the previous ID
          id[k][j][i] = idn[k][j][i]; 

          corrected.insert(Int3(i,j,k)); 
          
        }
      }

  int boundary_corrected = 0;
  // If some external ghost nodes need to be updated, update ID here (Phi and V updated outside anyway)
  for(auto&& g : *ghost_nodes_outer) {
    if(g.type_projection != GhostPoint::FACE)
      continue;
    auto it = corrected.find(g.image_ijk);
    if(it != corrected.end()) {
      id[g.ijk[2]][g.ijk[1]][g.ijk[0]] = id[(*it)[2]][(*it)[1]][(*it)[0]];
      boundary_corrected++;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &boundary_corrected, 1, MPI_INT, MPI_SUM, comm);

  int total_corrected = corrected.size();
  MPI_Allreduce(MPI_IN_PLACE, &total_corrected, 1, MPI_INT, MPI_SUM, comm);
  

  // Restore
  IDn.RestoreDataPointerToLocalVector();
  for(int surf=0; surf<(int)phi_ebm.size(); surf++)
    (*EBDS)[surf]->Phi_ptr->RestoreDataPointerToLocalVector();

  if(total_corrected>0) {
    for(int ls=0; ls<(int)phi.size(); ls++) //phi may be empty, which is fine here.
      Phi[ls]->RestoreDataPointerAndInsert();
    ID.RestoreDataPointerAndInsert();
  }  
  else {
    for(int ls=0; ls<(int)phi.size(); ls++) //phi may be empty, which is fine here.
      Phi[ls]->RestoreDataPointerToLocalVector();
    ID.RestoreDataPointerToLocalVector();
  }
    
/*
  MPI_Barrier(comm);
  print("Good 6. total_corrected = %d\n", total_corrected);
*/
  return total_corrected;

}

//-----------------------------------------------------

int 
MultiPhaseOperator::CheckLevelSetOverlapping(vector<SpaceVariable3D*> &Phi)
{

  int ls_size = Phi.size();

  if(ls_size<=1)
    return 0; //cannot have overlapping...

  vector<double***> phi(ls_size, NULL);
  for(int ls=0; ls<ls_size; ls++) 
    phi[ls] = Phi[ls]->GetDataPointer();
  
  int overlap = 0;
  bool inside = false;
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++) {

        inside = false;
        for(int ls = 0; ls<ls_size; ls++) {//loop through all the level set functions
          if(phi[ls][k][j][i]<0) {
            if(!inside) {
              inside = true;
            } else {
              overlap++;
              break;
            }
          }
        }

      }


  for(int ls=0; ls<ls_size; ls++) 
    Phi[ls]->RestoreDataPointerToLocalVector(); //no changes made

  MPI_Allreduce(MPI_IN_PLACE, &overlap, 1, MPI_INT, MPI_SUM, comm);

  return overlap;
}

//----------------------------------------------------- 
int
MultiPhaseOperator::UpdateCellsSweptByEmbeddedSurfaces(SpaceVariable3D &V, SpaceVariable3D &ID,
                                                       vector<SpaceVariable3D*> &Phi,
                                                       unique_ptr<vector<unique_ptr<EmbeddedBoundaryDataSet> > > EBDS,
                                                       vector<Intersector*> *intersector)
{

  // --------------------------------------------------------
  // Step 1: Extract useful information from EBDS
  // --------------------------------------------------------
  assert(EBDS);
  vector<std::set<Int3> > swept; //NOT pointers. Instead, it actually stores the "sets"
  vector<double***> phi_ebm;
  for(auto&& ebds : *EBDS) {
    phi_ebm.push_back(ebds->Phi_ptr->GetDataPointer());

    swept.push_back(*ebds->swept_ptr);
    // remove cells that are internal ghosts
    for(auto it = swept.back().begin(); it != swept.back().end();) {
      if(!ID.IsHere((*it)[0],(*it)[1],(*it)[2],false))
        it = swept.back().erase(it); //returns the iterator to the next element 
      else
        it++;
    }
  }

  // Clean up tag
  double*** tag = Tag.GetDataPointer();
  for(int k=kk0; k<kkmax; k++)
    for(int j=jj0; j<jjmax; j++)
      for(int i=ii0; i<iimax; i++) 
        tag[k][j][i] = 0;



  // --------------------------------------------------------
  // Step 2: Update occluded cells 
  // --------------------------------------------------------
  int has_occ(0);
  for(auto&& ebds : *EBDS) {
    if(ebds->occluded_ptr->size()>0) {
      has_occ = 1;
      break;
    }
    if(ebds->imposed_occluded_ptr->size()>0) {
      has_occ = 1;
      break;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &has_occ, 1, MPI_INT, MPI_MAX, comm);

  if(has_occ) {
    double*** id = ID.GetDataPointer();
    Vec5D***  v  = (Vec5D***) V.GetDataPointer();
    int ls_size = Phi.size();
    vector<double***> phi(ls_size, NULL);
    for(int ls=0; ls<ls_size; ls++) 
      phi[ls] = Phi[ls]->GetDataPointer();

    int i,j,k;
    int phi_updated = 0;
    for(int surf=0; surf<(int)EBDS->size(); surf++) {
      for(auto&& ijk : *(*EBDS)[surf]->occluded_ptr) {
        i = ijk[0];
        j = ijk[1];
        k = ijk[2];

        if(!ID.IsHere(i,j,k,false))
          continue; //this is an internal ghost of this subdomain

        id[k][j][i]   = INACTIVE_MATERIAL_ID;
        v[k][j][i][0] = iod.eqs.dummy_state.density;
        v[k][j][i][1] = iod.eqs.dummy_state.velocity_x;
	v[k][j][i][2] = iod.eqs.dummy_state.velocity_y;
        v[k][j][i][3] = iod.eqs.dummy_state.velocity_z;
        v[k][j][i][4] = iod.eqs.dummy_state.pressure;

        // Fix any issues in phi: phi must be non-negative
        for(int ls=0; ls<ls_size; ls++) {
          if(phi[ls][k][j][i]<0) {
            phi[ls][k][j][i] = phi_ebm[surf][k][j][i];
            if(phi_updated == 0)
              phi_updated = 1;
          }
        }

        tag[k][j][i] = -1;
      }

      for(auto&& ijk : *(*EBDS)[surf]->imposed_occluded_ptr) {
        i = ijk[0];
        j = ijk[1];
        k = ijk[2];

        if(!ID.IsHere(i,j,k,false))
          continue; //this is an internal ghost of this subdomain

        id[k][j][i]   = INACTIVE_MATERIAL_ID;
        v[k][j][i][0] = iod.eqs.dummy_state.density;
        v[k][j][i][1] = iod.eqs.dummy_state.velocity_x;
	v[k][j][i][2] = iod.eqs.dummy_state.velocity_y;
        v[k][j][i][3] = iod.eqs.dummy_state.velocity_z;
        v[k][j][i][4] = iod.eqs.dummy_state.pressure;

        // Fix any issues in phi: phi must be non-negative
        for(int ls=0; ls<ls_size; ls++) {
          if(phi[ls][k][j][i]<0) {
            phi[ls][k][j][i] = phi_ebm[surf][k][j][i];
            if(phi_updated == 0)
              phi_updated = 1;
          }
        }

        tag[k][j][i] = -2;
      }
    }  
    ID.RestoreDataPointerAndInsert();
    V.RestoreDataPointerAndInsert();
    MPI_Allreduce(MPI_IN_PLACE, &phi_updated, 1, MPI_INT, MPI_MAX, comm);
    if(phi_updated>0) {
      for(int ls=0; ls<ls_size; ls++) 
        Phi[ls]->RestoreDataPointerAndInsert();
    } else {
      for(int ls=0; ls<ls_size; ls++) 
        Phi[ls]->RestoreDataPointerToLocalVector();
    }
  }


  // --------------------------------------------------------
  // Step 3: Tag swept nodes (including occluded nodes, tagged above)
  // --------------------------------------------------------
  for(auto&& sub : swept)
    for(auto&& ijk : sub) {
      if(tag[ijk[2]][ijk[1]][ijk[0]]>0)
        fprintf(stdout,"\033[0;35mWarning: Node (%d %d %d) is swept by multiple (>1) embedded surfaces.\n\033[0m",
                ijk[0], ijk[1], ijk[2]);

      else if(tag[ijk[2]][ijk[1]][ijk[0]]==0)
        tag[ijk[2]][ijk[1]][ijk[0]] = 1; 
    }


  // --------------------------------------------------------
  // Step 4: Find out external ghosts that need to be updated at the end.
  // --------------------------------------------------------
  std::list<int> ghosts;
  for(int i=0; i<(int)ghost_nodes_outer->size(); i++) {
    GhostPoint& g((*ghost_nodes_outer)[i]);
    if(g.type_projection != GhostPoint::FACE)
      continue;
    int im_i(g.image_ijk[0]), im_j(g.image_ijk[1]), im_k(g.image_ijk[2]);
    if(tag[im_k][im_j][im_i] > 0)
      ghosts.push_back(i);
  }


  // --------------------------------------------------------
  // Step 5: Untag occluded nodes --- they have been resolved in Step 2.
  // --------------------------------------------------------
  if(has_occ) {
    for(auto&& ebds : *EBDS) {
      for(auto&& ijk : *ebds->occluded_ptr)
        tag[ijk[2]][ijk[1]][ijk[0]] = 0;
      for(auto&& ijk : *ebds->imposed_occluded_ptr)
        tag[ijk[2]][ijk[1]][ijk[0]] = 0;
    }  
  }
  Tag.RestoreDataPointerAndInsert();
  

  // --------------------------------------------------------
  // Step 6: Updating swept nodes, iteratively. (Usually one iteration should be enough.)
  // --------------------------------------------------------
  int total_count;
  int iter, max_iteration = 30;
  bool last_resort = false; //will be true if (1) progress = 0 or (2) max_iteration is reached

  for(iter=0; iter<max_iteration+2; iter++) {

    if(iter>=max_iteration)
      last_resort = true;

    // 6.1: Count number of cells within domain interior that need to be updated
    total_count = 0;
    for(auto&& sub : swept)
      total_count += sub.size();
    MPI_Allreduce(MPI_IN_PLACE, &total_count, 1, MPI_INT, MPI_SUM, comm);

    if(total_count==0)
      break;

    int progress = 0;

    // 6.2: Get data
    tag           = Tag.GetDataPointer();
    double*** id  = ID.GetDataPointer();
    Vec5D***  v   = (Vec5D***) V.GetDataPointer();
    int ls_size = Phi.size();
    vector<double***> phi(ls_size, NULL);
    for(int ls=0; ls<ls_size; ls++) 
      phi[ls] = Phi[ls]->GetDataPointer();
 

    // 6.3: Loop through unresolved swept nodes. 
    //      Note: One node should not be swept by multiple surfaces!
    int i,j,k,i2,j2,k2;
    vector<std::pair<Int3,bool> > neighbors; //useful neighbors, paired w/ connected or not
    for(int surf=0; surf<(int)swept.size(); surf++) {
      for(auto it = swept[surf].begin(); it != swept[surf].end();) {
        i = (*it)[0];
        j = (*it)[1];
        k = (*it)[2];

        if(tag[k][j][i]==0) {
          //It is possible (but uncommon) that a swept node is also imposed-occluded.
          //It should have already been updated in Step 2.
          it = swept[surf].erase(it);
          progress++;
          continue;
        }

/*
        bool gotyou = false;
        if(i==9 && j==45 && k==0)
          gotyou = true;
*/

        // Preparation: Collect useful neighbors and intersection info
        FindNeighborsForUpdatingSweptNode(i,j,k,tag,id,intersector,
                                          *(*EBDS)[surf]->occluded_ptr,
                                          *(*EBDS)[surf]->imposed_occluded_ptr,
                                          neighbors);
 
          
        // Determine its id: Check the ID of "valid" neighbors
        set<int> tmp_id_s;
        for(auto&& nei : neighbors)
          if(nei.second) //connected
            tmp_id_s.insert(id[nei.first[2]][nei.first[1]][nei.first[0]]);

        int myid(-1);
        if(tmp_id_s.size()==1)
          myid = *tmp_id_s.begin(); //update id[k][j][i] = myid later
        else {
          if(!last_resort) {
            it++;
            continue; //not ready to update this node.
          }
          else { //last resort / fail-safe
            vector<int> id_votes(varFcn.size(),0);
            for(auto&& nei : neighbors)
              if(nei.second) //connected
                id_votes[id[nei.first[2]][nei.first[1]][nei.first[0]]]++;
            int myvote = 0;
            for(int s=0; s<(int)id_votes.size(); s++)
              if(id_votes[s]>myvote) {
                myid = s;
                myvote = id_votes[s];
              }
            if(myid==-1) {
              fprintf(stdout,"\033[0;31m*** Error: Unable to determine the material ID at (%d,%d,%d).\n\033[0m",
                      i,j,k);
              exit(-1);
            }
          }
        }

        // Update state variables
        double sum_weight(0.0);
        double weight(0.0);
        Vec5D  vsum(0.0); 
        for(auto&& nei : neighbors) {
          if(nei.second) {//connected
            i2 = nei.first[0];
            j2 = nei.first[1];
            k2 = nei.first[2];
            if(myid == id[k2][j2][i2]) {
              double dist = (Vec3D(global_mesh.GetX(i) ,global_mesh.GetY(j) ,global_mesh.GetZ(k)) -
                             Vec3D(global_mesh.GetX(i2),global_mesh.GetY(j2),global_mesh.GetZ(k2))).norm();
              weight        = 1.0/dist; //weighted by 1/dist
              vsum         += weight*v[k2][j2][i2]; 
              sum_weight += weight;
            }
          }
        }
        if(sum_weight==0) {
          if(!last_resort) {
            it++;
            continue;  //not ready to update this node
          }
          else {
            fprintf(stdout,"\033[0;31m*** Error: Unable to determine the state at (%d,%d,%d), id:%d.\n\033[0m",
                    i,j,k, myid);
            exit(-1);
          }
        }

        v[k][j][i]  = vsum/sum_weight;
        id[k][j][i] = myid;


        // Fix any issues in phi: If the material subdomain corresponding to this ID is tracked by a phi, make sure the phi
        // value is negative; if it is not tracked by phi, the phi value should be non-negative.
        for(int ls=0; ls<ls_size; ls++) {
          if(id[k][j][i] == ls2matid[ls]) {          
            if(phi[ls][k][j][i]>=0) 
              phi[ls][k][j][i] = phi_ebm[surf][k][j][i]==0 ? -1.0e-20 : -phi_ebm[surf][k][j][i];
          }
          else {
            if(phi[ls][k][j][i]<0)
              phi[ls][k][j][i] = phi_ebm[surf][k][j][i];
          }
        }

        // Remove this swept node from workload
        it = swept[surf].erase(it);
        tag[k][j][i] = 0;
        progress++;
      }
    }


    MPI_Allreduce(MPI_IN_PLACE, &progress, 1, MPI_INT, MPI_SUM, comm);

    if(progress>0) {
      Tag.RestoreDataPointerAndInsert();
      ID.RestoreDataPointerAndInsert();
      V.RestoreDataPointerAndInsert();
      for(int ls=0; ls<ls_size; ls++) 
        Phi[ls]->RestoreDataPointerAndInsert();
    } 
    else {
      Tag.RestoreDataPointerToLocalVector();
      ID.RestoreDataPointerToLocalVector();
      V.RestoreDataPointerToLocalVector();
      for(int ls=0; ls<ls_size; ls++) 
        Phi[ls]->RestoreDataPointerToLocalVector();
      last_resort = true;
    }
  }

  // print warning and error messages
  if(iter>10) {
    if(!last_resort && verbose>=1) 
      print("- Performed a large number (%d) of iterations to update nodes/cells swept by embedded surfaces.\n",
            iter);
  }
  if(last_resort && total_count==0 && verbose>=1)
    print_warning("Warning: Activated the fail-safe procedure to update nodes swept by embedded surfaces.\n");
  if(total_count>0) {
    print_error("*** Error: Unable to update %d cells/nodes swept by embedded surfaces.\n", total_count);
    exit_mpi();
  }


  // Step 7: Apply boundary conditions to ID
  int total_ghosts = ghosts.size();
  MPI_Allreduce(MPI_IN_PLACE, &total_ghosts, 1, MPI_INT, MPI_SUM, comm);
  if(total_ghosts>0) {
    double*** id  = ID.GetDataPointer();
    for(auto&& gid : ghosts) {
      GhostPoint &g((*ghost_nodes_outer)[gid]);
      int i(g.ijk[0]), j(g.ijk[1]), k(g.ijk[2]);
      int im_i(g.image_ijk[0]), im_j(g.image_ijk[1]), im_k(g.image_ijk[2]);

      id[k][j][i] = id[im_k][im_j][im_i];
    }
    ID.RestoreDataPointerAndInsert();
  }


  // Step 8: Clean up
  for(auto&& ebds : *EBDS)
    ebds->Phi_ptr->RestoreDataPointerToLocalVector();


  return total_ghosts;
}

//-----------------------------------------------------

void 
MultiPhaseOperator::UpdateMaterialIDAtGhostNodes(SpaceVariable3D &ID)
{

#ifdef LEVELSET_TEST
  return; //testing the level set solver w/o solving the N-S / Euler equations
#endif

  double*** id  = (double***)ID.GetDataPointer();

  for(auto it = ghost_nodes_outer->begin(); it != ghost_nodes_outer->end();  it++) {

    if(it->type_projection != GhostPoint::FACE)
      continue; //corner (i.e. edge or vertex) nodes are not populated

    int i(it->ijk[0]), j(it->ijk[1]), k(it->ijk[2]);
    int im_i(it->image_ijk[0]), im_j(it->image_ijk[1]), im_k(it->image_ijk[2]);

    id[k][j][i] = id[im_k][im_j][im_i];

  }

  ID.RestoreDataPointerAndInsert();

}

//-----------------------------------------------------
// Section 4.2.4 of Arthur Rallu's thesis
int 
MultiPhaseOperator::UpdateStateVariablesAfterInterfaceMotion(SpaceVariable3D &IDn, 
                        SpaceVariable3D &ID, SpaceVariable3D &V, RiemannSolutions &riemann_solutions,
                        vector<Intersector*> *intersector, vector<Int3> &unresolved)
{

#ifdef LEVELSET_TEST
  return 0; //testing the level set solver w/o solving the N-S / Euler equations
#endif

  int unresolved_cells = 0;

  switch (iod.multiphase.phasechange_type) {

    case MultiPhaseData::RIEMANN_SOLUTION :
      unresolved_cells = UpdateStateVariablesByRiemannSolutions(IDn, ID, V, riemann_solutions, intersector, unresolved);
      break;

    case MultiPhaseData::EXTRAPOLATION :
      unresolved_cells = UpdateStateVariablesByExtrapolation(IDn, ID, V, intersector, unresolved);
      break;

    default :
      print_error("*** Error: Specified method for phase-change update (%d) has not been implemented.\n", 
                  (int)iod.multiphase.phasechange_type);
  }


  return unresolved_cells;
}

//-----------------------------------------------------

int
MultiPhaseOperator::UpdateStateVariablesByRiemannSolutions(SpaceVariable3D &IDn, 
                        SpaceVariable3D &ID, SpaceVariable3D &V, RiemannSolutions &riemann_solutions,
                        vector<Intersector*> *intersector, vector<Int3> &still_unresolved)
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

        if(id[k][j][i] == INACTIVE_MATERIAL_ID)
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
  int nStillUnresolved = 0;
  MPI_Allreduce(MPI_IN_PLACE, &nUnresolved, 1, MPI_INT, MPI_SUM, comm);
  if(nUnresolved) //some of the subdomains have unresolved nodes
    nStillUnresolved = FixUnresolvedNodes(unresolved, IDn, ID, V, intersector, still_unresolved,
                                          iod.multiphase.apply_failsafe_density==MultiPhaseData::On); 

  return nStillUnresolved;

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
      fprintf(stdout,"\033[0;35mWarning: Unable to update phase change at (%d,%d,%d) by Riemann solutions"
              " w/ upwinding. Retrying.\033[0m\n", i,j,k);
  } else {
    if(verbose>1) 
      fprintf(stdout,"\033[0;35mWarning: Unable to update phase change at (%d,%d,%d) by Riemann solutions."
              " Retrying.\033[0m\n", i,j,k);
  }

  return counter;
}

//-----------------------------------------------------

int
MultiPhaseOperator::UpdateStateVariablesByExtrapolation(SpaceVariable3D &IDn, 
                        SpaceVariable3D &ID, SpaceVariable3D &V, vector<Intersector*> *intersector, 
                        vector<Int3> &still_unresolved)
{
  // extract info
  double*** idn = (double***)IDn.GetDataPointer();
  double*** id  = (double***)ID.GetDataPointer();
  Vec5D***  v   = (Vec5D***) V.GetDataPointer();

  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();

  double weight, sum_weight;
  Vec5D vsum;
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

        if(id[k][j][i] == INACTIVE_MATERIAL_ID)
          continue;

        // coordinates of this node
        Vec3D& x0(coords[k][j][i]);

        sum_weight = 0.0;
        vsum       = 0.0;

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

              // coordinates of this neighbor
              Vec3D& x1(coords[neighk][neighj][neighi]);

              if(intersector) {
                bool connected = true;
                for(auto&& xter : *intersector) {
                  if(xter->Intersects(x0,x1)) {
                    connected = false;
                    break; //this neighbor is blocked to current node by an embedded surface
                  }
                }
                if(!connected)
                  continue;
              }

              if(iod.multiphase.phasechange_dir == MultiPhaseData::ALL) 
                weight = 1.0/((x0-x1).norm());
              else {//Upwind
                // velocity at the neighbor node
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
              }

              // add weighted s.v. at neighbor node
              if(weight>0) {
                sum_weight += weight;
                vsum       += weight*v[neighk][neighj][neighi];
              }
            }

        if(sum_weight==0) {
          if(verbose>1) {
            if(iod.multiphase.phasechange_dir == MultiPhaseData::ALL) 
              fprintf(stdout,"\033[0;35mWarning: Unable to update phase change at (%d,%d,%d)(%e,%e,%e) "
                      "by extrapolation.\n\033[0m", i,j,k, x0[0],x0[1],x0[2]);
            else
              fprintf(stdout,"\033[0;35mWarning: Unable to update phase change at (%d,%d,%d)(%e,%e,%e) "
                      "by extrapolation w/ upwinding.\n\033[0m", i,j,k, x0[0],x0[1],x0[2]);
          }
          unresolved.push_back(Int3(k,j,i)); //note the order: k,j,i          
        } else
          v[k][j][i] = vsum/sum_weight; 
      }

  V.RestoreDataPointerAndInsert(); //insert data & communicate with neighbor subd's
  ID.RestoreDataPointerToLocalVector();
  IDn.RestoreDataPointerToLocalVector();
  coordinates.RestoreDataPointerToLocalVector();


  // Fix the unresolved nodes (if any)
  int nUnresolved = unresolved.size();
  int nStillUnresolved = 0;
  MPI_Allreduce(MPI_IN_PLACE, &nUnresolved, 1, MPI_INT, MPI_SUM, comm);
  if(nUnresolved) //some of the subdomains have unresolved nodes
    nStillUnresolved = FixUnresolvedNodes(unresolved, IDn, ID, V, intersector, still_unresolved,
                                          iod.multiphase.apply_failsafe_density==MultiPhaseData::On); 

  return nStillUnresolved;

}

//-----------------------------------------------------

int
MultiPhaseOperator::FixUnresolvedNodes(vector<Int3> &unresolved, SpaceVariable3D &IDn, SpaceVariable3D &ID,
                                       SpaceVariable3D &V, vector<Intersector*> *intersector,
                                       vector<Int3> &still_unresolved,
                                       bool apply_failsafe_density) 
{

  // Note: all the processor cores will enter this function even if only one or a few have unresolved nodes

  still_unresolved.clear();

  // extract info
  double*** idn = (double***)IDn.GetDataPointer();
  double*** id  = (double***)ID.GetDataPointer();
  Vec5D***  v   = (Vec5D***) V.GetDataPointer();

  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();


  // loop through unresolved nodes
  int i,j,k;
  double weight, sum_weight = 0.0; 
  double sum_weight2 = 0.0; //similar, but regardless of "upwinding"
  Vec5D vtmp;
  Vec3D v1, x1x0;
  double v1norm, x1x0norm;

  bool reset = false; //whether v[k][j][i] has been reset

  int nStillUnresolved = 0;

  for(auto it = unresolved.begin(); it != unresolved.end(); it++) { 

    k = (*it)[0]; //order is k,j,i
    j = (*it)[1];
    i = (*it)[2];

    assert(idn[k][j][i] != id[k][j][i]);

    Vec3D& x0(coords[k][j][i]);

    sum_weight = 0.0;
    sum_weight2 = 0.0;
    vtmp = 0.0;

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

          if(intersector) {
            bool connected = true; 
            for(auto&& xter : *intersector) {
              if(xter->Intersects(x0,x1)) {
                connected = false;
                break; //this neighbor is blocked to current node by an embedded surface
              }
            }
            if(!connected)
              continue;
          }


          v1[0] = v[neighk][neighj][neighi][1];
          v1[1] = v[neighk][neighj][neighi][2];
          v1[2] = v[neighk][neighj][neighi][3];

          // compute weight
          v1norm = v1.norm();
          if(v1norm != 0)
            v1 /= v1norm;

          x1x0 = x0 - x1; 
          x1x0norm = x1x0.norm();
          x1x0 /= x1x0norm;

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

          // add to sum_weight2, regardless of upwinding
          vtmp += 1.0/x1x0norm*v[neighk][neighj][neighi];
          sum_weight2 += 1.0/x1x0norm;

        }


    if(iod.multiphase.phasechange_dir == MultiPhaseData::UPWIND) {
      if(sum_weight>0) {
        v[k][j][i] /= sum_weight; //Done!
        if(verbose>1) fprintf(stdout,"*** (%d,%d,%d): Updated state variables by extrapolation w/ upwinding. (2nd attempt)\n",
                              i,j,k);
        continue;
      }
    }

    // if still unresolved, try to apply an averaging w/o     
    if(sum_weight2>0) {
      v[k][j][i] = vtmp/sum_weight2; //Done!
      if(verbose>1) fprintf(stdout,"*** (%d,%d,%d): Updated state variables by extrapolation w/o enforcing upwinding."
                            " (2nd attempt)\n", i,j,k);
      continue;
    }

    // Our last resort: keep the pressure and velocity (both normal & TANGENTIAL) at the current node, 
    // find a valid density nearby. (In this case, the solution may be different for different domain
    // partitions --- but this should rarely happen. Also, we are no longer rigorously checking intersections...)
          
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
              continue; //this neighbor is outside the current subdomain (TODO: Hence, different subdomain partitions
                        //may affect the results. This can be fixed in future)

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

            if(intersector) {
              bool connected = true;
              for(auto&& xter : *intersector) {
                if(xter->Intersects(x0,x1)) {
                  connected = false;
                  break; //this neighbor is blocked to current node by an embedded surface
                }
              }
              if(!connected)
                continue;
            }

            double dist = (x1-x0).norm();

            sum_weight += 1.0/dist;
            density += (1.0/dist)*v[neighk][neighj][neighi][0];
          }

      if(sum_weight>0) {
        v[k][j][i][0] = density/sum_weight;
        if(verbose>1) fprintf(stdout,"*** (%d,%d,%d): Updated density by interpolation w/ stencil width = %d: %e %e %e %e %e\n",
                            i,j,k, layer, v[k][j][i][0], v[k][j][i][1], v[k][j][i][2], v[k][j][i][3], v[k][j][i][4]);
        break; //done with this node
      }

    }

    //Very unlikely. This means there is no neighbors within max_layer that are valid.
    if(sum_weight==0) { 

      //Treatment 1 (default): Store the (still) unresolved cells, and update Phi. Then, update material
      //    ID and state variables again.
      //Treatment 2: apply a constant density.

      if(id[k][j][i]!=0 || apply_failsafe_density) {
        fprintf(stdout,"\033[0;35mWarning: Updating phase change at (%d,%d,%d)(%e,%e,%e) with pre-specified density (%e). "
                       "Id:%d->%d. No valid neighbors within %d layers.\033[0m\n", 
                       i,j,k, coords[k][j][i][0], coords[k][j][i][1],
                       coords[k][j][i][2], varFcn[id[k][j][i]]->failsafe_density, (int)idn[k][j][i], (int)id[k][j][i], 
                       max_layer);
        v[k][j][i][0] = varFcn[id[k][j][i]]->failsafe_density;
      } else { //id[k][j][i] = 0 && not applying failsafe
        fprintf(stdout,"\033[0;35mWarning: Updating phase change at (%d,%d,%d)(%e,%e,%e). "
                       "Id:%d->%d. No valid neighbors within %d layers. Trying to correct the level set functions\033[0m\n", 
                       i,j,k, coords[k][j][i][0], coords[k][j][i][1],
                       coords[k][j][i][2], (int)idn[k][j][i], (int)id[k][j][i], max_layer);
        still_unresolved.push_back(Int3(k,j,i)); //order: k,j,i
        nStillUnresolved++;
      }
    }
  }


  if(!apply_failsafe_density) //check if there are cells still unresolved. If yes, it would prompt local update of level sets
    MPI_Allreduce(MPI_IN_PLACE, &nStillUnresolved, 1, MPI_INT, MPI_SUM, comm);


  V.RestoreDataPointerAndInsert(); //insert data & communicate with neighbor subd's

  ID.RestoreDataPointerToLocalVector();
  IDn.RestoreDataPointerToLocalVector();
  coordinates.RestoreDataPointerToLocalVector();


  return nStillUnresolved; //if the failsafe strategy has been applied, this will always be 0

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


  int NX, NY, NZ;
  coordinates.GetGlobalSize(&NX, &NY, &NZ);

  //---------------------------------------------
  // Step 1: Check for phase transitions; Update
  //         ID and V
  //---------------------------------------------
  double*** id  = ID.GetDataPointer();
  Vec5D***  v   = (Vec5D***) V.GetDataPointer();
  double*** lam = Lambda.GetDataPointer();

  int counter = 0;
  vector<tuple<Int3,int,int> >  changed; //(i,j,k), old id, new id -- incl. ghosts inside physical domain

  set<int> affected_ids;

  //DEBUG!!
/*
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {
        if(coords[k][j][i].norm()<0.05) {
          //v[k][j][i][0] = 8.9e-4;
          //v[k][j][i][4] = 1.2e10;
          int myid = id[k][j][i]; 
          double p0 = v[k][j][i][4];
          double e0 = varFcn[myid]->GetInternalEnergyPerUnitMass(v[k][j][i][0], p0);
          double T0 = varFcn[myid]->GetTemperature(v[k][j][i][0], e0);
          double T  = T0 + 0.4;
          double e  = varFcn[myid]->GetInternalEnergyPerUnitMassFromTemperature(v[k][j][i][0], T);
          double p  = varFcn[myid]->GetPressure(v[k][j][i][0], e);
          v[k][j][i][4] = p;
          if(coords[k][j][i].norm()<0.015)
            fprintf(stdout,"Changing state at (%d,%d,%d) p: %e->%e, T: %e->%e, h: %e->%e\n", i,j,k, 
                    p0, p, T0, T, e0 + p0/v[k][j][i][0], e + p/v[k][j][i][0]);
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

//        if(lam[k][j][i]>0)
//          fprintf(stdout,"lam[%d][%d][%d] = %e.\n", k,j,i, lam[k][j][i]);

        for(auto it = trans[myid].begin(); it != trans[myid].end(); it++) {

          //for debug only
          double rho0 = v[k][j][i][0];
          double p0 = v[k][j][i][4]; 
          double e0 = varFcn[myid]->GetInternalEnergyPerUnitMass(rho0, p0);
          double T0 = varFcn[myid]->GetTemperature(rho0, e0);

          if((*it)->Transition(v[k][j][i], lam[k][j][i])) { //NOTE: v[k][j][i][4] (p) and lam may be modified even if the 
                                                            //      return value is FALSE
            // detected phase transition

            // register the node
            changed.push_back(std::make_tuple(Int3(i,j,k), myid, (*it)->toID));

            // register the involved material ids
            affected_ids.insert(myid);
            affected_ids.insert((*it)->toID);

            // update id
            id[k][j][i] = (*it)->toID;

            // ------------------------------------------------------------------------
            // print to screen (for debug only)
            double rho1 = v[k][j][i][0];
            double p1 = v[k][j][i][4];
            double e1 = varFcn[(*it)->toID]->GetInternalEnergyPerUnitMass(rho1, p1);
            double T1 = varFcn[(*it)->toID]->GetTemperature(rho1, e1);
            fprintf(stdout,"Detected phase transition at (%d,%d,%d)(%d->%d). rho: %e->%e, p: %e->%e, T: %e->%e, h: %e->%e.\n", 
                    i,j,k, myid, (*it)->toID, rho0, rho1, p0, p1, T0, T1, e0+p0/rho0, e1+p1/rho1);
            // ------------------------------------------------------------------------


            // if node is next to a symmetry or wall boundary, update the ID of the ghost node (V will be updated by spo)
            if(i==0 && (iod.mesh.bc_x0==MeshData::SLIPWALL || iod.mesh.bc_x0==MeshData::STICKWALL ||
                        iod.mesh.bc_x0==MeshData::SYMMETRY))  
              id[k][j][i-1] = id[k][j][i];
            if(i==NX-1 && (iod.mesh.bc_xmax==MeshData::SLIPWALL || iod.mesh.bc_xmax==MeshData::STICKWALL || 
                           iod.mesh.bc_xmax==MeshData::SYMMETRY))
              id[k][j][i+1] = id[k][j][i];
         
            if(j==0 && (iod.mesh.bc_y0==MeshData::SLIPWALL || iod.mesh.bc_y0==MeshData::STICKWALL || 
                        iod.mesh.bc_y0==MeshData::SYMMETRY))  
              id[k][j-1][i] = id[k][j][i];
            if(j==NY-1 && (iod.mesh.bc_ymax==MeshData::SLIPWALL || iod.mesh.bc_ymax==MeshData::STICKWALL || 
                           iod.mesh.bc_ymax==MeshData::SYMMETRY))
              id[k][j+1][i] = id[k][j][i];
         
            if(k==0 && (iod.mesh.bc_z0==MeshData::SLIPWALL || iod.mesh.bc_z0==MeshData::STICKWALL || 
                        iod.mesh.bc_z0==MeshData::SYMMETRY))  
              id[k-1][j][i] = id[k][j][i];
            if(k==NZ-1 && (iod.mesh.bc_zmax==MeshData::SLIPWALL || iod.mesh.bc_zmax==MeshData::STICKWALL ||
                           iod.mesh.bc_zmax==MeshData::SYMMETRY))
              id[k+1][j][i] = id[k][j][i];
         

            counter++;
            break;
          }
        }
      }

  MPI_Allreduce(MPI_IN_PLACE, &counter, 1, MPI_INT, MPI_SUM, comm);

  Lambda.RestoreDataPointerAndInsert();

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
  for(int ls = 0; ls < (int)Phi.size(); ls++)
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


//-------------------------------------------------------------------------

int
MultiPhaseOperator::ResolveConflictsInLevelSets(int time_step, vector<SpaceVariable3D*> &Phi)
{

  int ls_size = Phi.size();

  if(ls_size==0)
    return 0; //nothing to do


  // whether perform some optional "improvements"
  bool optional_improvements = false;
  if(iod.multiphase.levelset_correction_frequency>0 &&
     time_step % iod.multiphase.levelset_correction_frequency == 0)
    optional_improvements = true;


  int resolved_conflicts = 0;

  vector<double***> phi(ls_size, NULL);

  for(int ls=0; ls<ls_size; ls++)
    phi[ls] = Phi[ls]->GetDataPointer();


  // ------------------------------------------
  // PART I: Find & resolve cells that are
  //         covered by more than one subdomain
  // ------------------------------------------

  if(ls_size>=2) {
    // loop through all the nodes
    vector<int> boundaries;
    vector<int> owner, inter;
    for(int k=kk0; k<kkmax; k++)
      for(int j=jj0; j<jjmax; j++)
        for(int i=ii0; i<iimax; i++) {

          boundaries.clear();

          //----------------------------------------------------------
          // Find subdomains that have this node next to its boundary
          //----------------------------------------------------------
          for(int ls=0; ls<ls_size; ls++) {
            if((i-1>=ii0  && phi[ls][k][j][i]*phi[ls][k][j][i-1]<=0) ||
               (i+1<iimax && phi[ls][k][j][i]*phi[ls][k][j][i+1]<=0) ||
               (j-1>=jj0  && phi[ls][k][j][i]*phi[ls][k][j-1][i]<=0) ||
               (j+1<jjmax && phi[ls][k][j][i]*phi[ls][k][j+1][i]<=0) ||
               (k-1>=kk0  && phi[ls][k][j][i]*phi[ls][k-1][j][i]<=0) ||
               (k+1<kkmax && phi[ls][k][j][i]*phi[ls][k+1][j][i]<=0))
              boundaries.push_back(ls);
          }

          if(boundaries.size()<=1) //nothing to worry about
            continue;

          owner.clear();
          inter.clear();
          for(auto&& ls : boundaries) {
            if(phi[ls][k][j][i]<0) 
              owner.push_back(ls);
            else if(phi[ls][k][j][i]==0)
              inter.push_back(ls);
          }

          if(owner.size()==0 && inter.size()==0)
            continue; //this node does not belong to any of the subdomains.

          if(owner.size()==1) {//great

            if(!optional_improvements)
              continue; // do nothing

            double new_phi = 0.0;
            for(auto&& ls : boundaries)
              new_phi += fabs(phi[ls][k][j][i]); 
            new_phi /= boundaries.size();
            for(auto&& ls : boundaries) {
              if(phi[ls][k][j][i]<0)
                phi[ls][k][j][i] = -new_phi;
              else
                phi[ls][k][j][i] = new_phi;

            }
          }
          else if(owner.size()==0) {// inter.size() is not 0

            if(!optional_improvements)
              continue; // do nothing

            double new_phi = 0.0;
            for(auto&& ls : boundaries)
              new_phi += fabs(phi[ls][k][j][i]); 
            new_phi /= boundaries.size();
            for(auto&& ls : boundaries) {
              if(ls==inter[0]) //you are the owner
                phi[ls][k][j][i] = -new_phi;
              else
                phi[ls][k][j][i] = new_phi;

            }
          }
          else {//owner.size()>1

            //1. find a unique owner
            int new_owner = owner[0];
            double max_phi = fabs(phi[owner[0]][k][j][i]);
            for(int ind=1; ind<(int)owner.size(); ind++) {
              if(fabs(phi[owner[ind]][k][j][i])>max_phi) {
                new_owner = owner[ind];
                max_phi = fabs(phi[owner[ind]][k][j][i]);
              }
            }

            //2. get a new phi (abs. value)
            double new_phi = 0.0;
            for(auto it = owner.begin(); it != owner.end(); it++)
              new_phi += fabs(phi[*it][k][j][i]); 
            new_phi /= owner.size();

            //3. update all the involved level set functions
            for(auto it = owner.begin(); it != owner.end(); it++) {
              if(*it==new_owner) //you are the owner
                phi[*it][k][j][i] = -new_phi;
              else
                phi[*it][k][j][i] = new_phi;

            }

            resolved_conflicts++;
          }
        }
  }


  // ------------------------------------------
  // PART II: (Optional) Find & resolve cells 
  //          that are trapped between material
  //          interfaces
  // ------------------------------------------

  if(optional_improvements) {

    int NX, NY, NZ;
    coordinates.GetGlobalSize(&NX, &NY, &NZ);
    
    // loop through the domain interior
    for(int k=k0; k<kmax; k++)
      for(int j=j0; j<jmax; j++)
        for(int i=i0; i<imax; i++) {

          bool background_cell = true; 
          for(int ls=0; ls<ls_size; ls++) {
            if(phi[ls][k][j][i]<0) {
              background_cell = false; 
              break; //not a background cell (ID != 0)
            }
          }
          if(!background_cell)
            continue;

          int qi = 0; //like "Qi" in Weiqi/"Go"
          // check neighbors to see if this is an isolated cell
          // left neighbor
          if(i-1>=0) {
            bool connected = true;
            for(int ls=0; ls<ls_size; ls++)
              if(phi[ls][k][j][i-1]<0) {
                connected = false;
                break;
              } 
            if(connected) {
              qi++;
              if(qi>=2)
                continue;
            }
          }
  
          // right neighbor
          if(i+1<NX) {
            bool connected = true;
            for(int ls=0; ls<ls_size; ls++)
              if(phi[ls][k][j][i+1]<0) {
                connected = false;
                break;
              } 
            if(connected) {
              qi++;
              if(qi>=2)
                continue;
            }
          }
 
          // bottom neighbor
          if(j-1>=0) {
            bool connected = true;
            for(int ls=0; ls<ls_size; ls++)
              if(phi[ls][k][j-1][i]<0) {
                connected = false;
                break;
              } 
            if(connected) {
              qi++;
              if(qi>=2)
                continue;
            }
          }
  
          // top neighbor
          if(j+1<NY) {
            bool connected = true;
            for(int ls=0; ls<ls_size; ls++)
              if(phi[ls][k][j+1][i]<0) {
                connected = false;
                break;
              } 
            if(connected) {
              qi++;
              if(qi>=2)
                continue;
            }
          }
 
          // back neighbor
          if(k-1>=0) {
            bool connected = true;
            for(int ls=0; ls<ls_size; ls++)
              if(phi[ls][k-1][j][i]<0) {
                connected = false;
                break;
              } 
            if(connected) {
              qi++;
              if(qi>=2)
                continue;
            }
          }
  
          // front neighbor
          if(k+1<NZ) {
            bool connected = true;
            for(int ls=0; ls<ls_size; ls++)
              if(phi[ls][k+1][j][i]<0) {
                connected = false;
                break;
              } 
            if(connected) {
              qi++;
              if(qi>=2)
                continue;
            }
          }
 
          // qi has to be 0 or 1 now
 
          // ---------------------------------------
          // This is an isolated background cell.
          // ---------------------------------------
          double min_phi = DBL_MAX; 
          int new_owner = -1;
          for(int ls=0; ls<ls_size; ls++) {
            if(phi[ls][k][j][i]<min_phi) {
              new_owner = ls;
              min_phi   = phi[ls][k][j][i];
            }
          }
          assert(min_phi>=0);
          phi[new_owner][k][j][i] = -min_phi;
            
          resolved_conflicts++;
        }

  }


  MPI_Allreduce(MPI_IN_PLACE, &resolved_conflicts, 1, MPI_INT, MPI_SUM, comm);


  for(int ls=0; ls<ls_size; ls++) {
    if(resolved_conflicts>0)
      Phi[ls]->RestoreDataPointerAndInsert();
    else
      Phi[ls]->RestoreDataPointerToLocalVector();
  }


  return resolved_conflicts;

}

//-------------------------------------------------------------------------

void
MultiPhaseOperator::UpdateLevelSetsInUnresolvedCells(vector<SpaceVariable3D*> &Phi, vector<Int3> &unresolved)
{

  //Note that "unresolved" may be non-empty only for some of the processor cores / subdomains
  //What we do here is very similar to "Part II" in function ResolveConflictsInLevelSets
  //TODO:Currently, this only handles unresolved *background* cells (i.e. ID = 0)

  int ls_size = Phi.size();

  if(ls_size==0)
    return; //nothing to do

  vector<double***> phi(ls_size, NULL);

  for(int ls=0; ls<ls_size; ls++)
    phi[ls] = Phi[ls]->GetDataPointer();


  // ---------------------------------------------------------------
  // Resolve *background* cells that are trapped between material interfaces
  // ---------------------------------------------------------------

  int NX, NY, NZ, i, j, k;
  coordinates.GetGlobalSize(&NX, &NY, &NZ);

  for(auto it = unresolved.begin(); it != unresolved.end(); it++) {

    k = (*it)[0]; //order should be k,j,i
    j = (*it)[1];
    i = (*it)[2];

    bool background_cell = true; 
    for(int ls=0; ls<ls_size; ls++) {
      if(phi[ls][k][j][i]<0) {
        background_cell = false; 
        break; //not a background cell (ID != 0)
      }
    }
    if(!background_cell)
      continue;

    // ---------------------------------------
    // This is an "isolated" background cell.
    // ---------------------------------------
    double min_phi = DBL_MAX; 
    int new_owner = -1;
    for(int ls=0; ls<ls_size; ls++) {
      if(phi[ls][k][j][i]<min_phi) {
        new_owner = ls;
        min_phi   = phi[ls][k][j][i];
      }
    }
    assert(min_phi>=0);
    phi[new_owner][k][j][i] = -min_phi;
      
  }


  for(int ls=0; ls<ls_size; ls++)
    Phi[ls]->RestoreDataPointerAndInsert();

}

//-------------------------------------------------------------------------

void
MultiPhaseOperator::AddLambdaToEnthalpyAfterInterfaceMotion(SpaceVariable3D &IDn, SpaceVariable3D &ID, 
                                                            SpaceVariable3D &V)
{

  if(trans.size()==0 || iod.multiphase.latent_heat_transfer!=MultiPhaseData::On)
    return; //nothing to do

  double*** idn = IDn.GetDataPointer();
  double*** id  = ID.GetDataPointer();
  Vec5D***  v   = (Vec5D***) V.GetDataPointer();
  double*** lam = Lambda.GetDataPointer();

  int myidn, myid;
  int counter = 0;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        myidn = (int)idn[k][j][i];
        myid  = (int)id[k][j][i];

        if(myidn == myid) //id remains the same. Skip
          continue;

        if(lam[k][j][i]<=0.0) //no latent heat here
          continue;

        for(auto it = trans[myidn].begin(); it != trans[myidn].end(); it++) {

          // check if this is the phase transition from "myidn" to "myid"
          if((*it)->ToID() != myid)
            continue;

          //---------------------------------------------------------------------
          // Now, do the actual work: Add lam to enthalpy
          double rho = v[k][j][i][0];
          //double p   = v[k][j][i][4];
          double e   = varFcn[myid]->GetInternalEnergyPerUnitMass(v[k][j][i][0], v[k][j][i][4]);
          //double h   = e + p/rho + lam[k][j][i]; //adding lam
          e += lam[k][j][i]; //Correction of phase-change model (see Xuning JFM)
          lam[k][j][i] = 0.0;
          // update p to account for the increase of enthalpy (rho is fixed)
          v[k][j][i][4] = varFcn[myid]->GetPressure(rho, e);

          counter++;
          //---------------------------------------------------------------------
          
        }
      }

  MPI_Allreduce(MPI_IN_PLACE, &counter, 1, MPI_INT, MPI_SUM, comm);

  IDn.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();
  if(counter>0) {  
    Lambda.RestoreDataPointerAndInsert();
    V.RestoreDataPointerAndInsert();
  } else { //nothing is changed...
    Lambda.RestoreDataPointerToLocalVector();
    V.RestoreDataPointerToLocalVector();
  }

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

  for(int ls = 0; ls<(int)Phi.size(); ls++) {//loop through all the level set functions
  
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
  for(int ls = 0; ls<(int)Phi.size(); ls++) {//loop through all the level set functions
  
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

void
MultiPhaseOperator::FindNeighborsForUpdatingSweptNode(int i, int j, int k, double*** tag, [[maybe_unused]] double*** id,
                                                      vector<Intersector*> *intersector, set<Int3> &occluded,
                                                      set<Int3> &imposed_occluded,
                                                      vector<std::pair<Int3,bool> > &neighbors)
{
  neighbors.clear();

  Vec3D X0(global_mesh.GetX(i), global_mesh.GetY(j), global_mesh.GetZ(k)), X1;

  // Loop through neighbors
  for(int k2=k-1; k2<=k+1; k2++)
    for(int j2=j-1; j2<=j+1; j2++)
      for(int i2=i-1; i2<=i+1; i2++) {

        if(k2==k && j2==j && i2==i) //itself
          continue; 

        if(!coordinates.IsHereOrInternalGhost(i2,j2,k2)) //outside physical domain
          continue;

        if(tag[k2][j2][i2] != 0) //swept, and unresolved yet
          continue;

        Int3 ijk2(i2,j2,k2);
        if(occluded.find(ijk2) != occluded.end())
          continue; //ignore occluded neighbors
        if(imposed_occluded.find(ijk2) != imposed_occluded.end())
          continue; //ignore occluded neighbors

        bool connected = true;
        X1 = Vec3D(global_mesh.GetX(i2), global_mesh.GetY(j2), global_mesh.GetZ(k2));
        for(auto&& xter : *intersector) {
          if(xter->Intersects(X0, X1)) {
            connected = false;
            break;
          }
        }

        neighbors.push_back(std::make_pair(ijk2, connected));
      }
}

//-----------------------------------------------------

bool
MultiPhaseOperator::IsOrphanAcrossEmbeddedSurfaces(int i, int j, int k, double*** idn, double*** id,
                                                   vector<Intersector*> *intersector)
{

  Vec3D X0(global_mesh.GetX(i), global_mesh.GetY(j), global_mesh.GetZ(k)), X1;

  int myid = id[k][j][i];

  // Loop through neighbors
  for(int k2=k-1; k2<=k+1; k2++)
    for(int j2=j-1; j2<=j+1; j2++)
      for(int i2=i-1; i2<=i+1; i2++) {

        if(k2==k && j2==j && i2==i) //itself
          continue; 

        if(!coordinates.IsHereOrInternalGhost(i2,j2,k2)) //outside physical domain
          continue;

        if(id[k2][j2][i2] != idn[k2][j2][i2]) //this node just changed id
          continue;

        if(id[k2][j2][i2] != myid)
          continue;

        bool connected = true;
        X1 = Vec3D(global_mesh.GetX(i2), global_mesh.GetY(j2), global_mesh.GetZ(k2));
        for(auto&& xter : *intersector) {
          if(xter->Intersects(X0, X1)) {
            connected = false;
            break;
          }
        }

        if(connected)
          return false; //has a valid neighbor

      }

  return true;

}

//-----------------------------------------------------


