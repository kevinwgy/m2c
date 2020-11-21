#ifndef _SPACEVARIABLE_
#define _SPACEVARIABLE_
#include <petscmdma.h>

// -----------------------------------------------
// This class stores all the DM's
// -----------------------------------------------
class DataManagers2D {

public:
  // Note: To avoid ambiguities, ghosted means DM_BOUNDARY_GHOSTED. We do NOT use DM_BOUNDARY_MIRROR, 
  //       _PERIODIC, _TWISTED, ... Also, we always use DMDA_STENCIL_BOX for "stencil type"
  //
  static DM ghosted1_1dof; //ghosted"1" --> stencil width is 1 
  static DM ghosted1_2dof;  
  static DM ghosted1_3dof;  
  static DM ghosted1_4dof;  
  static DM ghosted1_5dof;  

public:
  DataManagers2D();
  DataManagers2D(MPI_Comm comm, int NX, int NY);
  ~DataManagers2D();

  void CreateAllDataManagers(MPI_Comm comm, int NX, int NY);
  void DestroyAllDataManagers(); //need to call this before "PetscFinalize()".

};

// -----------------------------------------------
// Defines a space variable
// -----------------------------------------------
class SpaceVariable2D {
  MPI_Comm   comm; 
  DM*        dm;   
  Vec        globalVec;  //each process only stores a local portion, without ghost
  Vec        localVec; //local portion of globalVec (separate memory allocation), with ghost layer
  void*      array; //user should only edit "array", which is a pointer to localVec. this class handles 
                    //array <-> localVec <-> globalVec
  int        dof;

  int        NX, NY; //global number of grid points in x and y directions
  int        nProcX, nProcY; //number of processors in X and Y directions

  int        i0, j0; //lower-left corner of the actual subdomain
  int        imax, jmax; //upper-right corner of the actual subdomain
  int        nx, ny; //width of the actual subdomain in x and y directions

  bool       ghosted; //whether da is ghosted.
  int        ghost_width; //width of the ghost layer (usually 1)
  int        ghost_i0, ghost_j0; //lower-left corner of the ghost layer
  int        ghost_imax, ghost_jmax; //upper-right corner of the ghost layer
  int        ghost_nx, ghost_ny; //width of the ghosted subdomain in x and y directions

public:
  SpaceVariable2D(MPI_Comm &comm_, DM *dm_);
  ~SpaceVariable2D();

  void *GetDataPointer(); 
  void RestoreDataPointer(); 
  void AssembleInsert();
  void AssembleAdd();
  void Destroy(); //should be called before PetscFinalize!

  void GetCornerIndices(int *i0_, int *j0_, int *imax_=0, int *jmax_=0) {
    *i0_ = i0; *j0_ = j0; 
    if(imax_) *imax_0 = imax; 
    if(jmax_) *jmax_jmax;
  }

  void GetSize(int *nx_, int *ny_) {*nx_ = nx;  *ny_ = ny;}

  void GetGhostedCornerIndices(int *i0_, int *j0_, int *imax_=0, int *jmax_=0) {
    *i0_ = ghost_i0; *j0_ = ghost_j0;
    if(imax_) *imax_ = ghost_imax;
    if(jmax_) *jmax_ = ghost_jmax;
  }

  void GetGhostedSize(int *ghost_nx_, int *ghost_ny_) {*ghost_nx_ = ghost_nx; *ghost_ny_ = ghost_ny;}

  int  NumGhostLayers() {return ghost_width;}
  int  NumDOF() {return dof;}

  int  NumProcs(int *nProcX_, int *nProcY_) {*nProcX_ = nProcX; *nProcY_ = nProcY;}
  int  GetGlobalSize(int *NX_, int *NY_) {*NX_ = NX; *NY_ = NY;}
 
};

#endif
