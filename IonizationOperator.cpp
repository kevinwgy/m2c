/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<IonizationOperator.h>
#include<NonIdealSahaEquationSolver.h>
#include<Vector5D.h>

extern int verbose;

//-----------------------------------------------------------------------

IonizationOperator::IonizationOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_, 
                                      std::vector<VarFcnBase*> &varFcn_)
                  : iod(iod_), Zav(comm_, &(dm_all_.ghosted1_1dof)),
                    Nh(comm_, &(dm_all_.ghosted1_1dof)), Ne(comm_, &(dm_all_.ghosted1_1dof))
{
  // Create Saha solvers
  if(iod.ion.materialMap.dataMap.size() == 0) {
    print_error("*** Error: Unable to create ionization operator. No models specified in the input file.\n");
    exit_mpi();
  }
  int nMat = varFcn_.size();
  saha.resize(nMat, NULL);
  for(auto it = iod.ion.materialMap.dataMap.begin(); it != iod.ion.materialMap.dataMap.end(); it++) {
    if(it->first<0 || it->first>=nMat) {
      print_error("*** Error: Ionization model specified for an unknown material id (%d).\n", it->first);
      exit_mpi();
    }
    print("- Initializing Saha Equation solver for material %d.\n", it->first);
    saha[it->first] = new SahaEquationSolver(*(it->second), iod, varFcn_[it->first], &comm_);
    switch (it->second->type) {
      case MaterialIonizationModel::SAHA_IDEAL :
        print("- Initializing ideal Saha Equation solver for material %d.\n", it->first);
        saha[it->first] = new SahaEquationSolver(*(it->second), iod, varFcn_[it->first], &comm_);
        break;
      case MaterialIonizationModel::SAHA_NONIDEAL :
        print("- Initializing non-ideal Saha Equation solver for material %d.\n", it->first);
        saha[it->first] = new NonIdealSahaEquationSolver(*(it->second), iod, varFcn_[it->first], &comm_);
        break;
      default :
        print_error("*** Error: Cannot initialize ionization solver for material %d. Unknown type.\n",
                    it->first);
        exit_mpi();
    }
  }
  for(int i=0; i<(int)saha.size(); i++) { //create dummy solvers for materials w/o ionization model
    if(saha[i] == NULL)
      saha[i] = new SahaEquationSolver(iod, varFcn_[i]); 
  }

  // Create AlphaRJ depending on user-requested outputs
  for(int iSpecies=0; iSpecies<OutputData::MAXSPECIES; iSpecies++) {
    if(iod.output.molar_fractions[iSpecies] == OutputData::ON) {
      switch (iod.output.max_charge_number) {
        case 0:
          AlphaRJ[iSpecies] = new SpaceVariable3D(comm_, &(dm_all_.ghosted1_2dof)); break;
        case 1:
          AlphaRJ[iSpecies] = new SpaceVariable3D(comm_, &(dm_all_.ghosted1_3dof)); break;
        case 2:
          AlphaRJ[iSpecies] = new SpaceVariable3D(comm_, &(dm_all_.ghosted1_4dof)); break;
        case 3:
          AlphaRJ[iSpecies] = new SpaceVariable3D(comm_, &(dm_all_.ghosted1_5dof)); break;
        case 4:
          AlphaRJ[iSpecies] = new SpaceVariable3D(comm_, &(dm_all_.ghosted1_6dof)); break;
        default:
          print_error("*** Error: Max. charge number (for output) must be between 0 and 4. Found %d.\n",
                      iod.output.max_charge_number); 
          exit_mpi();
      }
    }
  }
  max_charge_in_output = AlphaRJ.empty() ? -1 : iod.output.max_charge_number;

}

//-----------------------------------------------------------------------

IonizationOperator::~IonizationOperator()
{
  for(auto it = saha.begin(); it != saha.end(); it++)
    if(*it)
      delete *it;
  for(auto it = AlphaRJ.begin(); it != AlphaRJ.end(); it++)
    if(it->second)
      delete it->second;
}
 

//-----------------------------------------------------------------------

void
IonizationOperator::Destroy()
{
  Zav.Destroy();
  Nh.Destroy();
  Ne.Destroy();
  for(auto it = AlphaRJ.begin(); it != AlphaRJ.end(); it++)
    it->second->Destroy();
}

//-----------------------------------------------------------------------

Vec3D
IonizationOperator::ComputeIonizationAtOnePoint(int id, double rho, double p)
{
  double v[5] = {rho, 0.0, 0.0, 0.0, p}; //velocities are not needed
  std::map<int, vector<double> > nodal_alphas; //empty map
  Vec3D result; //(zav, nh, he)

  saha[id]->Solve(v, result[0], result[1], result[2], nodal_alphas);

  return result;
}

//-----------------------------------------------------------------------

void
IonizationOperator::ComputeIonization(SpaceVariable3D &V, SpaceVariable3D &ID)
{

  Vec5D***   v  = (Vec5D***)V.GetDataPointer();
  double*** id  = ID.GetDataPointer();
  double*** zav = Zav.GetDataPointer();
  double*** nh  = Nh.GetDataPointer();
  double*** ne  = Ne.GetDataPointer();

  std::map<int, double***> alphas;
  std::map<int, vector<double> > nodal_alphas;

  for(auto it = AlphaRJ.begin(); it != AlphaRJ.end(); it++) {
    alphas[it->first] = it->second->GetDataPointer();
    nodal_alphas[it->first] = vector<double>(max_charge_in_output+2);
  }

  // Main loop
  int k0, kmax, j0, jmax, i0, imax;
  ID.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        saha[id[k][j][i]]->Solve(v[k][j][i], zav[k][j][i], nh[k][j][i], ne[k][j][i], nodal_alphas);

        auto it0 = nodal_alphas.begin();
        for(auto it = alphas.begin(); it != alphas.end(); it++) {
          for(int p=0; p<max_charge_in_output+2; p++) 
            (it->second)[k][j][i*(max_charge_in_output+2)+p] = (it0->second)[p];
          it0++;
        }

      }


  V.RestoreDataPointerToLocalVector();
  ID.RestoreDataPointerToLocalVector();
  Zav.RestoreDataPointerAndInsert();
  Nh.RestoreDataPointerAndInsert();
  Ne.RestoreDataPointerAndInsert();
  for(auto it = AlphaRJ.begin(); it != AlphaRJ.end(); it++)
    it->second->RestoreDataPointerAndInsert();

  if(verbose>=2)
    print("- Computed ionizations.\n");
}

//-----------------------------------------------------------------------




