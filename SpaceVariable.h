/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _SPACEVARIABLE_
#define _SPACEVARIABLE_
#include <petscdmda.h>
#include <vector>

/*******************************************
 * This class stores all the DM's
 *******************************************
*/
class DataManagers3D {

public:
  //! Note: To avoid ambiguities, ghosted means DM_BOUNDARY_GHOSTED. We do NOT use DM_BOUNDARY_MIRROR, 
  //!       _PERIODIC, _TWISTED, ... Also, we always use DMDA_STENCIL_BOX for "stencil type"
  //
  DM ghosted1_1dof; //!< ghosted"1" --> stencil width is 1 
  DM ghosted1_2dof;  
  DM ghosted1_3dof;  
  DM ghosted1_4dof;  
  DM ghosted1_5dof;  
  DM ghosted1_6dof;  
  DM ghosted1_9dof;  

  DM ghosted2_1dof;
  DM ghosted2_3dof;

public:
  DataManagers3D();
  DataManagers3D(MPI_Comm comm, int NX, int NY, int NZ);
  ~DataManagers3D();

  int CreateAllDataManagers(MPI_Comm comm, int NX, int NY, int NZ);
  void DestroyAllDataManagers(); //!< need to call this before "PetscFinalize()".

};

/*******************************************
 * Defines a space variable
 * Note: Upon initialization, SpaceVariable3D
 *       is filled with 0.
 *******************************************
 */
class SpaceVariable3D {

  MPI_Comm*  comm; 
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
  int        internal_ghost_i0, internal_ghost_j0, internal_ghost_k0; //!< include internal ghosts (i.e. inside physical domain)
  int        internal_ghost_imax, internal_ghost_jmax, internal_ghost_kmax; //!< include internal ghosts
  int        ghost_nx, ghost_ny, ghost_nz; //!< width of the ghosted subdomain in x, y and z directions

  int        numNodes0; //number of interior nodes
  int        numNodes1; //number of interior nodes + internal ghost nodes
  int        numNodes2; //number of interior nodes + internal & external ghost nodes

public:
  SpaceVariable3D(MPI_Comm &comm_, DM *dm_);
  SpaceVariable3D(); //must be followed by a call to function Setup(...)
  ~SpaceVariable3D();

  void Setup(MPI_Comm &comm_, DM *dm_);

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
  void WriteToVTRFile(const char *filename, const char *varname = NULL); //!< write vector to file (w/o mesh info)
  void SetOutputVariableName(const char *name); //!< give it a name, which will show up in output files (VTR/VTK)

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

  inline void GetInternalGhostedCornerIndices(int *i0_, int *j0_, int *k0_, int *imax_=0, int *jmax_=0, int *kmax_=0) {
    *i0_ = internal_ghost_i0; *j0_ = internal_ghost_j0; *k0_ = internal_ghost_k0;
    if(imax_) *imax_ = internal_ghost_imax;
    if(jmax_) *jmax_ = internal_ghost_jmax;
    if(kmax_) *kmax_ = internal_ghost_kmax;
  }

  inline void GetGhostedSize(int *ghost_nx_, int *ghost_ny_, int *ghost_nz_) {
    *ghost_nx_ = ghost_nx; *ghost_ny_ = ghost_ny; *ghost_nz_ = ghost_nz;}

  inline void GetInternalGhostedSize(int *ghost_nx_, int *ghost_ny_, int *ghost_nz_) {
    *ghost_nx_ = internal_ghost_imax - internal_ghost_i0; 
    *ghost_ny_ = internal_ghost_jmax - internal_ghost_j0; 
    *ghost_nz_ = internal_ghost_kmax - internal_ghost_k0;}

  inline int  NumGhostLayers() {return ghost_width;}
  inline int  NumDOF() {return dof;}

  inline void NumProcs(int *nProcX_, int *nProcY_, int *nProcZ_) {*nProcX_ = nProcX; *nProcY_ = nProcY; *nProcZ_ = nProcZ;}
  inline void GetGlobalSize(int *NX_, int *NY_, int *NZ_) {*NX_ = NX; *NY_ = NY; *NZ_ = NZ;}
 
  inline Vec& GetRefToGlobalVec() {return globalVec;}

  inline bool OutsidePhysicalDomain(int i, int j, int k) {return (i<0 || i>=NX || j<0 || j>=NY || k<0 || k>=NZ);}

  double CalculateGlobalMin(int mydof = 0, bool workOnGhost = false);

  double CalculateGlobalMax(int mydof = 0, bool workOnGhost = false);

  inline int NumInternalNodes() {return numNodes0;}
  inline int NumNodesIncludingInternalGhosts() {return numNodes1;}
  inline int NumNodesIncludingGhosts() {return numNodes2;}

  inline bool OutsidePhysicalDomainAndUnpopulated(int i, int j, int k)
  {
    int count = 0;
    if(i<0 || i>=NX) count++;
    if(j<0 || j>=NY) count++;
    if(k<0 || k>=NZ) count++;
    if(count>1)
      return true;
    return false;
  } 

  inline int BoundaryType(int i, int j, int k) //0~interior, 1~boundary-face, 2~boundary-edge, 3~boundary-corner
  {return int(i<0) + int(i>=NX) + int(j<0) + int(j>=NY) + int(k<0) + int(k>=NZ);}

  inline bool IsHere(int i, int j, int k, bool include_ghost = false)
  {
    if(include_ghost) 
      return i>=ghost_i0 && i<ghost_imax && j>=ghost_j0 && j<ghost_jmax && k>=ghost_k0 && k<ghost_kmax;
    else
      return i>=i0 && i<imax && j>=j0 && j<jmax && k>=k0 && k<kmax;
    return false;
  }

  inline bool IsHereOrInternalGhost(int i, int j, int k)
  {
    return i>=internal_ghost_i0 && i<internal_ghost_imax && 
           j>=internal_ghost_j0 && j<internal_ghost_jmax && 
           k>=internal_ghost_k0 && k<internal_ghost_kmax;
  }

  inline bool IsHereOrConnectedGhost(int i, int j, int k)
  {return int(i<i0 || i>=imax) + int(j<j0 || j>=jmax) + int(k<k0 || k>=kmax) <= 1;}

  //! operators
  void AXPlusB(double a, double b, bool workOnGhost = false); //!< self = a*self + b;
  void AXPlusBY(double a, double b, SpaceVariable3D &y, bool workOnGhost = false); //!< self = a*self + b*vector_y
  void AXPlusBY(double a, double b, SpaceVariable3D &y, std::vector<int>& Xindices,
                std::vector<int>& Yindices, bool workOnGhost = false); //!< customized version
  void SetConstantValue(double a, bool workOnGhost = false); //!< set value to a
  
};

#endif
