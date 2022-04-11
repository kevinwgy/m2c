#include<FloodFill.h>
#include<Utils.h>

using std::vector;
using std::set;

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

  Vec3D***     ob = (Vec3D***) Obs.GetDataPointer();
  double*** color = Color.GetDataPointer();

/*
  if(mpi_rank==1) {
  fprintf(stderr,"ii0 = %d, jj0= %d, kk0 = %d, iimax_in = %d, jjmax_in = %d, kkmax_in = %d.\n", ii0, jj0, kk0, iimax_in, jjmax_in, kkmax_in);
  for(int k=kk0; k<kkmax_in; k++)
    for(int j=jj0; j<jjmax_in; j++)
      for(int i=ii0; i<iimax_in; i++) {
        if(ob[k][j][i][0] != non_obstruction_flag)
          fprintf(stderr,"[%d][%d][%d] --/-- [%d][%d][%d].\n", k,j,i-1, k,j,i);
        if(ob[k][j][i][1] != non_obstruction_flag)
          fprintf(stderr,"[%d][%d][%d] --/-- [%d][%d][%d].\n", k,j-1,i, k,j,i);
        if(ob[k][j][i][2] != non_obstruction_flag)
          fprintf(stderr,"[%d][%d][%d] --/-- [%d][%d][%d].\n", k-1,j,i, k,j,i);
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
    assert(seed != Int3(iimax_in-1,jjmax_in-1,kkmax_in-1)); //everyone else are decided. this is the last one...
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
  Color.RestoreDataPointerAndInsert();

  Color.WriteToVTRFile("Sign0.vtr", "color0");
  //---------------------------------
  // Part II. Unionize colors
  //---------------------------------
  int iter;
  int numColors(-1);
  for(iter=0; iter<mpi_size; iter++) {

    color = Color.GetDataPointer();

    //loop through internal ghost nodes who just got a new value from its owner. Find colors
    //that need to changed
    int i,j,k;
    std::map<int,int> old2new; //old color to new color
    for(auto it = ghost_nodes_inner.begin(); it != ghost_nodes_inner.end(); it++) {
      i = it->ijk[0];
      j = it->ijk[1];
      k = it->ijk[2];
      if(Color.IsHere(i+1,j,k,false) && //left boundary
         ob[k][j][i+1][0] == non_obstruction_flag && //not obstructed
         color[k][j][i+1] > color[k][j][i]) { //colors do not match, and mine needs to be reset
        if(old2new.find(color[k][j][i+1]) == old2new.end() || old2new[color[k][j][i+1]] > color[k][j][i])
          old2new[color[k][j][i+1]] = color[k][j][i];
      }
      if(Color.IsHere(i-1,j,k,false) && //right boundary
         ob[k][j][i][0] == non_obstruction_flag && //not obstructed
         color[k][j][i-1] > color[k][j][i]) { //colors dos not match, and mine needs to be reset
        if(old2new.find(color[k][j][i-1]) == old2new.end() || old2new[color[k][j][i-1]] > color[k][j][i])
          old2new[color[k][j][i-1]] = color[k][j][i]; 
      }
      if(Color.IsHere(i,j+1,k,false) && //bottom boundary
         ob[k][j+1][i][1] == non_obstruction_flag && //not obstructed
         color[k][j+1][i] > color[k][j][i]) { //colors dos not match, and mine needs to be reset
        if(old2new.find(color[k][j+1][i]) == old2new.end() || old2new[color[k][j+1][i]] > color[k][j][i])
          old2new[color[k][j+1][i]] = color[k][j][i]; 
      }
      if(Color.IsHere(i,j-1,k,false) && //top boundary
         ob[k][j][i][1] == non_obstruction_flag && //not obstructed
         color[k][j-1][i] > color[k][j][i]) { //colors dos not match, and mine needs to be reset
        if(old2new.find(color[k][j-1][i]) == old2new.end() || old2new[color[k][j-1][i]] > color[k][j][i])
          old2new[color[k][j-1][i]] = color[k][j][i]; 
      }
      if(Color.IsHere(i,j,k+1,false) && //back boundary
         ob[k+1][j][i][2] == non_obstruction_flag && //not obstructed
         color[k+1][j][i] > color[k][j][i]) { //colors dos not match, and mine needs to be reset
        if(old2new.find(color[k+1][j][i]) == old2new.end() || old2new[color[k+1][j][i]] > color[k][j][i])
          old2new[color[k+1][j][i]] = color[k][j][i]; 
      }
      if(Color.IsHere(i,j,k-1,false) && //front boundary
         ob[k][j][i][2] == non_obstruction_flag && //not obstructed
         color[k-1][j][i] > color[k][j][i]) { //colors dos not match, and mine needs to be reset
        if(old2new.find(color[k-1][j][i]) == old2new.end() || old2new[color[k-1][j][i]] > color[k][j][i])
          old2new[color[k-1][j][i]] = color[k][j][i]; 
      }
    }

    // Check if we are done...
    int total_color_change = old2new.size();
    MPI_Allreduce(MPI_IN_PLACE, &total_color_change, 1, MPI_INT, MPI_SUM, comm);

    if(total_color_change == 0) { //yeah!
      for(int k=kk0_in; k<kkmax_in; k++)
        for(int j=jj0_in; j<jjmax_in; j++)
          for(int i=ii0_in; i<iimax_in; i++)
            numColors = std::max(numColors, (int)color[k][j][i]); //find numColors

      MPI_Allreduce(MPI_IN_PLACE, &numColors, 1, MPI_INT, MPI_MAX, comm);

      Color.RestoreDataPointerToLocalVector();
      break;
    } 

    // Go over the nodes and update color
    if(!old2new.empty()) {
      for(int k=kk0_in; k<kkmax_in; k++)
        for(int j=jj0_in; j<jjmax_in; j++)
          for(int i=ii0_in; i<iimax_in; i++) {
            if(old2new.find(color[k][j][i]) != old2new.end())
              color[k][j][i] = old2new[color[k][j][i]];
          }

    }
    Color.RestoreDataPointerAndInsert();

  } 

  if(iter==mpi_size) {
    print_error("*** Error: Flood-fill failed after %d iterations.\n", iter);
    exit_mpi();
  }

  Obs.RestoreDataPointerToLocalVector();

  return numColors;
}


//-----------------------------------------------------------------------------------

