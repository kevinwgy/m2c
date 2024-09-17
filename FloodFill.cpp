/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<FloodFill.h>
#include<Utils.h>

using std::vector;
using std::set;
using std::pair;
using std::make_pair;

//-----------------------------------------------------------------------------------

FloodFill::FloodFill(MPI_Comm &comm_, DataManagers3D &dms_, vector<GhostPoint> &ghost_nodes_inner_,
            vector<GhostPoint> &ghost_nodes_outer_)
         : comm(comm_), TMP(comm_, &(dms_.ghosted1_1dof)),
           ghost_nodes_inner(ghost_nodes_inner_), ghost_nodes_outer(ghost_nodes_outer_)
{
  TMP.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  TMP.GetGhostedCornerIndices(&ii0, &jj0, &kk0, &iimax, &jjmax, &kkmax);
  TMP.GetInternalGhostedCornerIndices(&ii0_in, &jj0_in, &kk0_in, &iimax_in, &jjmax_in, &kkmax_in);
  TMP.GetGlobalSize(&NX, &NY, &NZ);
}

//-----------------------------------------------------------------------------------

FloodFill::~FloodFill()
{ }

//-----------------------------------------------------------------------------------

void
FloodFill::Destroy()
{
  TMP.Destroy();
}

//-----------------------------------------------------------------------------------

int
FloodFill::FillBasedOnEdgeObstructions(SpaceVariable3D& Obs, int non_obstruction_flag,
                                       set<Int3>& occluded_nodes, SpaceVariable3D& Color)
{
  const int UNDECIDED = -1, IN_QUEUE = -2;

  //Note: Only fills nodes within the physical domain.

  int mpi_rank(-1), mpi_size(0);
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);

  assert(Obs.NumDOF() == 3);

  Vec3D***     ob = (Vec3D***) Obs.GetDataPointer();
  double*** color = Color.GetDataPointer();

/*
  if(mpi_rank==1) {
  fprintf(stdout,"ii0 = %d, jj0= %d, kk0 = %d, iimax_in = %d, jjmax_in = %d, kkmax_in = %d.\n", ii0, jj0, kk0, iimax_in, jjmax_in, kkmax_in);
  for(int k=kk0; k<kkmax_in; k++)
    for(int j=jj0; j<jjmax_in; j++)
      for(int i=ii0; i<iimax_in; i++) {
        if(ob[k][j][i][0] != non_obstruction_flag)
          fprintf(stdout,"[%d][%d][%d] --/-- [%d][%d][%d].\n", k,j,i-1, k,j,i);
        if(ob[k][j][i][1] != non_obstruction_flag)
          fprintf(stdout,"[%d][%d][%d] --/-- [%d][%d][%d].\n", k,j-1,i, k,j,i);
        if(ob[k][j][i][2] != non_obstruction_flag)
          fprintf(stdout,"[%d][%d][%d] --/-- [%d][%d][%d].\n", k-1,j,i, k,j,i);
      }
  }
*/

  //---------------------------------
  // Part I. Fill the subdomain
  //---------------------------------
  // 1. Preparation
  int mycolor(0), nColored(0); 
  for(int k=kk0_in; k<kkmax_in; k++)
    for(int j=jj0_in; j<jjmax_in; j++)
      for(int i=ii0_in; i<iimax_in; i++)
        color[k][j][i] = UNDECIDED; 

  for(auto it = occluded_nodes.begin(); it != occluded_nodes.end(); it++) {
    color[(*it)[2]][(*it)[1]][(*it)[0]] = 0;  //0: occluded
    nColored++;
  }
  
  // 2. Flood fill
  std::queue<Int3> Q;
  Int3 seed(ii0_in, jj0_in, kk0_in);
  while(nColored < Color.NumNodesIncludingInternalGhosts()) {

    mycolor++; //color starts at 1

    // find an undecided node as seed
    for(int k=kk0_in; k<kkmax_in; k++)
      for(int j=jj0_in; j<jjmax_in; j++)
        for(int i=ii0_in; i<iimax_in; i++) {
          if(color[k][j][i] != UNDECIDED)
            continue; 
          seed[0] = i;
          seed[1] = j;
          seed[2] = k;
          goto FOUND_SEED;
        }
    FOUND_SEED:
    //assert(seed != Int3(iimax_in-1,jjmax_in-1,kkmax_in-1)); //everyone else are decided. this is the last one...
    Q.push(seed);

    int i,j,k;
    while(!Q.empty()) {

      Int3& head(Q.front());
      i = head[0];
      j = head[1];
      k = head[2];

      // color the "head"
      color[k][j][i] = mycolor;
      nColored++;
      Q.pop();
          
      // go over the neighbors of head
      if(i-1>=ii0_in && color[k][j][i-1] == UNDECIDED && ob[k][j][i][0] == non_obstruction_flag) {
        Q.push(Int3(i-1,j,k)); color[k][j][i-1] = IN_QUEUE;}
      if(j-1>=jj0_in && color[k][j-1][i] == UNDECIDED && ob[k][j][i][1] == non_obstruction_flag) {
        Q.push(Int3(i,j-1,k)); color[k][j-1][i] = IN_QUEUE;}
      if(k-1>=kk0_in && color[k-1][j][i] == UNDECIDED && ob[k][j][i][2] == non_obstruction_flag) {
        Q.push(Int3(i,j,k-1)); color[k-1][j][i] = IN_QUEUE;}
      if(i+1<iimax_in && color[k][j][i+1] == UNDECIDED && ob[k][j][i+1][0] == non_obstruction_flag) {
        Q.push(Int3(i+1,j,k)); color[k][j][i+1] = IN_QUEUE;}
      if(j+1<jjmax_in && color[k][j+1][i] == UNDECIDED && ob[k][j+1][i][1] == non_obstruction_flag) {
        Q.push(Int3(i,j+1,k)); color[k][j+1][i] = IN_QUEUE;}
      if(k+1<kkmax_in && color[k+1][j][i] == UNDECIDED && ob[k+1][j][i][2] == non_obstruction_flag) {
        Q.push(Int3(i,j,k+1)); color[k+1][j][i] = IN_QUEUE;}

    }
  }

  Obs.RestoreDataPointerToLocalVector();

  // store local colors at inner ghost nodes (used in Part II)
  vector<int> ghost_nodes_inner_color(ghost_nodes_inner.size(), -1);
  for(int i=0; i<(int)ghost_nodes_inner.size(); i++) {
    Int3 &ijk(ghost_nodes_inner[i].ijk);
    ghost_nodes_inner_color[i] = color[ijk[2]][ijk[1]][ijk[0]];
  }

  Color.RestoreDataPointerAndInsert();

  // If there is a single processor, "unionization" is not needed --> we are done.
  if(mpi_size==1)
    return mycolor;

  //Color.WriteToVTRFile("Sign0.vtr", "color0");

  //---------------------------------
  // Part II. Unionize colors (for multi-processor runs)
  //---------------------------------
  color = Color.GetDataPointer();
  // ----------------------
  // II.1. Each subdomain checks for equivalent colors with neighbors
  // ----------------------
  set<Int3> equiv; //[0]: local color;  [1]: owner id;  [2]: color from owner
  vector<bool> unique_color(mycolor+1, true); //whether color is uniquely inside subdomain

  for(int i=0; i<(int)ghost_nodes_inner.size(); i++) {
    Int3 &ijk(ghost_nodes_inner[i].ijk);
    if(ghost_nodes_inner_color[i] == 0) {//occluded
      assert(color[ijk[2]][ijk[1]][ijk[0]]==0); //should be occluded from owner
      continue;
    }
    if(color[ijk[2]][ijk[1]][ijk[0]] == 0) {//occluded
      assert(ghost_nodes_inner_color[i] == 0);
      continue;
    }
  
    equiv.insert(Int3(ghost_nodes_inner_color[i], ghost_nodes_inner[i].owner_proc, color[ijk[2]][ijk[1]][ijk[0]]));
    unique_color[ghost_nodes_inner_color[i]] = false; 
  }

  // now, add colors that are unique to this subdomain (e.g., a closure that is entirely inside the subdomain)
  for(int c = 1; c <= mycolor; c++)
    if(unique_color[c])
      equiv.insert(Int3(c, mpi_rank, c));

/*
  //Debug (check for inconsistencies)
  auto last_it = equiv.end();
  last_it--;
  if(equiv.size()>1) { 
    for(auto it = equiv.begin(); it != last_it; it++) {
      auto it_next = it;
      it_next++;
      for(auto it2 = it_next; it2 != equiv.end(); it2++) {
        if((*it)[0] == (*it2)[0] && (*it)[1] == (*it2)[1] && (*it)[2] != (*it2)[2]) {
          fprintf(stdout,"\033[0;31m*** Error: Detected inconsistent colors. Proc.%d (me): %d = Proc.%d, %d or %d.\n\033[0m",
                         mpi_rank, (*it)[0], (*it)[1], (*it)[2], (*it2)[2]);
          exit(-1);
        }
        if((*it)[0] != (*it2)[0] && (*it)[1] == (*it2)[1] && (*it)[2] == (*it2)[2]) {
          fprintf(stdout,"\033[0;31m*** Error: Detected inconsistent colors. Proc.%d (me): %d or %d = Proc.%d, %d.\n\033[0m",
                         mpi_rank, (*it)[0], (*it2)[0], (*it)[1], (*it)[2]);
          exit(-1);
        }
      }
    } 
  }
*/

  // ----------------------
  // II.2. Proc. 0 collects data from everyone
  // ----------------------
  int worker = 0; //ID of the processor that does this job
  int my_data_count = 4*equiv.size();
  vector<int> my_data(my_data_count,0);
  int counter = 0;
  for(auto it = equiv.begin(); it != equiv.end(); it++) {
    my_data[counter++] = mpi_rank;
    my_data[counter++] = (*it)[0];
    my_data[counter++] = (*it)[1];
    my_data[counter++] = (*it)[2];
  }

  // talk
  vector<int> counts;
  if(mpi_rank == worker) { //I am the receiver
    counts.resize(mpi_size,-1);
    MPI_Gather(&my_data_count, 1, MPI_INT, counts.data(), 1, MPI_INT, worker, comm);
  } else //I am a sender
    MPI_Gather(&my_data_count, 1, MPI_INT, NULL, 1, MPI_INT, worker, comm);

  // work
  vector<int> displacements;
  counter = 0;
  if(mpi_rank == worker) {
    displacements.resize(mpi_size,-1); 
    for(int i=0; i<(int)displacements.size(); i++) {
      displacements[i] = counter;
      counter += counts[i];
    }
  }

  // talk again
  vector<int> all_data;
  if(mpi_rank == worker) { //receiver
    all_data.resize(counter);
    MPI_Gatherv(my_data.data(), my_data_count, MPI_INT, all_data.data(), counts.data(), displacements.data(), MPI_INT, worker, comm);
/*
    for(int i=0; i<counter/4; i++)
      fprintf(stdout,"[%d] %d = [%d] %d.\n", all_data[4*i], all_data[4*i+1], all_data[4*i+2], all_data[4*i+3]);
*/
  } else {//sender
    MPI_Gatherv(my_data.data(), my_data_count, MPI_INT, NULL, NULL, NULL, MPI_INT, worker, comm);
  }


  // ----------------------
  // II.3. Proc. 0 unionizes colors
  // ----------------------
  vector<set<pair<int, int> > > loc2glob;  //loc2glob[i]: i-th proc.;  loc2glob[i][j]: j-th color-pair (old->new)
  int current_color(-1);
  if(mpi_rank == worker) {
    
    loc2glob.resize(mpi_size);

    int num_pairs = counter/4;
    vector<bool> unionized(num_pairs, false);

    current_color = 0;

    set<pair<int,int> > glob2loc;

    for(int i=0; i<num_pairs; i++) {

      if(unionized[i])
        continue;

      current_color++; //color starts at 1, not 0, which is for occluded

      glob2loc.clear(); //a new color

      glob2loc.insert(std::make_pair(all_data[4*i], all_data[4*i+1]));
      glob2loc.insert(std::make_pair(all_data[4*i+2], all_data[4*i+3]));

      loc2glob[all_data[4*i]].insert(make_pair(all_data[4*i+1], current_color));
      loc2glob[all_data[4*i+2]].insert(make_pair(all_data[4*i+3], current_color));
      unionized[i] = true;
      bool new_member = true; //has a new pair that has "current_color"

      while(new_member==true) {

        new_member = false;

        for(int j=0; j<num_pairs; j++) {
        
          if(unionized[j])
            continue;

          pair<int,int> p1(all_data[4*j],   all_data[4*j+1]); //p1 and p2 should have the same final color
          pair<int,int> p2(all_data[4*j+2], all_data[4*j+3]);

          if(glob2loc.find(p1) != glob2loc.end()) {
            glob2loc.insert(p2);
            loc2glob[p2.first].insert(make_pair(p2.second, current_color));
            unionized[j] = true;
            new_member = true;
          } else if(glob2loc.find(p2) != glob2loc.end()) {
            glob2loc.insert(p1);
            loc2glob[p1.first].insert(make_pair(p1.second, current_color));
            unionized[j] = true;
            new_member = true;
          }
        }
      }
    } 
/*
    for(int i=0; i<mpi_size; i++)
      for(auto it = loc2glob[i].begin(); it != loc2glob[i].end(); it++)
        fprintf(stdout,"[%d] loc %d --> glob %d.\n", i, it->first, it->second);
*/
  }

  // ----------------------
  // II.4. Everyone gets the final colors
  // ----------------------
  // talk again to figure out package sizes
  if(mpi_rank == worker) { 

    for(int i=0; i<mpi_size; i++)
      counts[i] = 2*loc2glob[i].size();

    MPI_Scatter(counts.data(), 1, MPI_INT, &my_data_count, 1, MPI_INT, worker, comm);
  } else 
    MPI_Scatter(NULL, 1, MPI_INT, &my_data_count, 1, MPI_INT, worker, comm);

 
  assert(my_data_count/2 == mycolor);

  // worker prepares packages
  if(mpi_rank == worker) {
    counter = 0;
    for(int i=0; i<mpi_size; i++) {
      displacements[i] = counter;

      for(auto it = loc2glob[i].begin(); it != loc2glob[i].end(); it++) {
        all_data[counter++] = it->first;
        all_data[counter++] = it->second;
      }
    }
  }

  // talk again --> each subdomain gets the final colors
  if(mpi_rank == worker) //sender
    MPI_Scatterv(all_data.data(), counts.data(), displacements.data(), MPI_INT, my_data.data(), my_data_count, MPI_INT, worker, comm);
  else //receiver
    MPI_Scatterv(NULL, NULL, NULL, MPI_INT, my_data.data(), my_data_count, MPI_INT, worker, comm);


  // ----------------------
  // II.5. Everyone applies the final colors
  // ----------------------
  vector<int> old2new(mycolor+1, -1);
  old2new[0] = 0; //occluded --> 0
  int old_color(-1), new_color(-1);
  for(int i=0; i<my_data_count/2; i++) {
    old_color = my_data[2*i];
    new_color = my_data[2*i+1];
    assert(old_color>0 && old_color<=mycolor);
    assert(old2new[old_color]==-1); //no conflicts.
    old2new[old_color] = new_color;
    //fprintf(stdout,"[%d] old color %d --> new color %d.\n", mpi_rank, old_color, new_color);
  } 
    
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++)
        color[k][j][i] = old2new[color[k][j][i]];

  Color.RestoreDataPointerAndInsert();


  // ----------------------
  // III. Finalization
  // ----------------------
  // final talk --> everyone should know (and return) the total number of non-zero colors
  MPI_Bcast(&current_color, 1, MPI_INT, worker, comm);

  return current_color;

}


//-----------------------------------------------------------------------------------

