#ifndef _SPACEVARIABLE_
#define _SPACEVARIABLE_
#include <petscdmda.h>

/*******************************************
 * This class stores all the DM's
 *******************************************
*/
class DataManagers3D {

public:
  //! Note: To avoid ambiguities, ghosted means DM_BOUNDARY_GHOSTED. We do NOT use DM_BOUNDARY_MIRROR, 
  //!       _PERIODIC, _TWISTED, ... Also, we always use DMDA_STENCIL_BOX for "stencil type"
  //
  static DM ghosted1_1dof; //!< ghosted"1" --> stencil width is 1 
  static DM ghosted1_2dof;  
  static DM ghosted1_3dof;  
  static DM ghosted1_5dof;  

public:
  DataManagers3D();
  DataManagers3D(MPI_Comm comm, int NX, int NY, int NZ);
  ~DataManagers3D();

  int CreateAllDataManagers(MPI_Comm comm, int NX, int NY, int NZ);
  void DestroyAllDataManagers(); //!< need to call this before "PetscFinalize()".

};

/*******************************************
 * Defines a space variable
 *******************************************
 */
class SpaceVariable3D {
  MPI_Comm&  comm; 
  DM*        dm;   
  Vec        globalVec;  //!< each process only stores a local portion, without ghost
  Vec        localVec; //!< local portion of globalVec (separate memory allocation), with ghost layer
  double***  array; /**< user should only edit "array", which is a pointer to localVec. this class handles 
                      *  array <-> localVec <-> globalVec */
  int        dof;

  int        NX, NY, NZ; //!< global number of grid points in x, y and z directions
  int        nProcX, nProcY, nProcZ; //!< number of processors in x, y and z directions

  int        i0, j0, k0; //!< lower-left corner of the actual subdomain
  int        imax, jmax, kmax; //!< upper-right corner of the actual subdomain
  int        nx, ny, nz; //!< width of the actual subdomain in x, y, and z directions

  bool       ghosted; //!< whether da is ghosted.
  int        ghost_width; //!< width of the ghost layer (usually 1)
  int        ghost_i0, ghost_j0, ghost_k0; //!< lower-left corner of the ghost layer
  int        ghost_imax, ghost_jmax, ghost_kmax; //!< upper-right corner of the ghost layer
  int        ghost_nx, ghost_ny, ghost_nz; //!< width of the ghosted subdomain in x, y and z directions

public:
  SpaceVariable3D(MPI_Comm &comm_, DM *dm_);
  ~SpaceVariable3D();

  double*** GetDataPointer(); 

  /** The following two functions involve MPI communications
   *  Note that only the data in the real domain gets "communicated" (i.e. inserted or added)
   *  Data in the ghost boundary (i.e. outside the physical domain) do not participate 
   *  in any communications
   */
  void RestoreDataPointerAndInsert();
  void RestoreDataPointerAndAdd();

  void RestoreDataPointerToLocalVector(); //!< caution: does not update globalVec
  void Destroy(); //!< should be called before PetscFinalize!

  void StoreMeshCoordinates(SpaceVariable3D &coordinates);
  void WriteToVTRFile(const char *filename); //!< write vector to file (w/o mesh info)

  inline void GetCornerIndices(int *i0_, int *j0_, int *k0_, int *imax_=0, int *jmax_=0, int *kmax_=0) {
    *i0_ = i0; *j0_ = j0; *k0_ = k0; 
    if(imax_) *imax_ = imax; 
    if(jmax_) *jmax_ = jmax;
    if(kmax_) *kmax_ = kmax;
  }

  inline void GetSize(int *nx_, int *ny_, int *nz_) {*nx_ = nx;  *ny_ = ny; *nz_ = nz;}

  inline void GetGhostedCornerIndices(int *i0_, int *j0_, int *k0_, int *imax_=0, int *jmax_=0, int *kmax_=0) {
    *i0_ = ghost_i0; *j0_ = ghost_j0; *k0_ = ghost_k0;
    if(imax_) *imax_ = ghost_imax;
    if(jmax_) *jmax_ = ghost_jmax;
    if(kmax_) *kmax_ = ghost_kmax;
  }

  inline void GetGhostedSize(int *ghost_nx_, int *ghost_ny_, int *ghost_nz_) {
    *ghost_nx_ = ghost_nx; *ghost_ny_ = ghost_ny; *ghost_nz_ = ghost_nz;}

  inline int  NumGhostLayers() {return ghost_width;}
  inline int  NumDOF() {return dof;}

  inline void NumProcs(int *nProcX_, int *nProcY_, int *nProcZ_) {*nProcX_ = nProcX; *nProcY_ = nProcY; *nProcZ_ = nProcZ;}
  inline void GetGlobalSize(int *NX_, int *NY_, int *NZ_) {*NX_ = NX; *NY_ = NY; *NZ_ = NZ;}
 
  inline Vec& GetRefToGlobalVec() {return globalVec;}

  inline bool OutsidePhysicalDomain(int i, int j, int k) {return (i<0 || i>=NX || j<0 || j>=NY || k<0 || k>=NZ);}

  //! operators
  void AXPlusB(double a, double b, bool workOnGhost = false); //!< self = a*self + b;
  void AXPlusBY(double a, double b, SpaceVariable3D &y, bool workOnGhost = false); //!< self = a*self + b*vector_y
  void SetConstantValue(double a, bool workOnGhost = false); //!< set value to a
  
};

#endif
