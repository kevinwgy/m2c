/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _IONIZATION_OPERATOR_H_
#define _IONIZATION_OPERATOR_H_

#include<SpaceVariable.h>
#include<VarFcnBase.h>
#include<IoData.h>
#include<SahaEquationSolver.h>

/****************************************************************
 * Class IonizationOperator is responsible for computing and
 * outputing ionization variables. Currently, only the Saha equation
 * is implemented. Can be extended in future.
 * *************************************************************/

struct IonizationOperator {

  IoData& iod;

  std::vector<SahaEquationSolver*> saha; //one "saha" for each materialid

  //! Variables for output
  SpaceVariable3D Zav; //!< average charge per heavy particle (dimensionless), 
  SpaceVariable3D Nh; //!< number density of heavy particles or nuclei
  SpaceVariable3D Ne; //!< number density of electrons (=Zav*Nh) 
  std::map<int, SpaceVariable3D*> AlphaRJ; //!< molar fraction of species/element j, charge r 
  int max_charge_in_output; //!< max charge number in output files (i.e. dim of each SpaceVariable3D in AlphaRJ)

public:

  IonizationOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_, std::vector<VarFcnBase*> &varFcn_);
  ~IonizationOperator();

  void Destroy();

  Vec3D ComputeIonizationAtOnePoint(int id, double rho, double p);

  void ComputeIonization(SpaceVariable3D &V, SpaceVariable3D &ID);

  // Get reference to results (e.g., for outputing)
  SpaceVariable3D& GetReferenceToZav() {return Zav;}
  SpaceVariable3D& GetReferenceToNh() {return Nh;}
  SpaceVariable3D& GetReferenceToNe() {return Ne;}
  std::map<int, SpaceVariable3D*>& GetReferenceToAlphaRJ() {return AlphaRJ;}

};

#endif
