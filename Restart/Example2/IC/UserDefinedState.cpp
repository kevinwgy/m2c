/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include "UserDefinedState.h"
#include "Vector3D.h"
#include "Vector5D.h"
#include "MathTools/trilinear_interpolation.h"

#include <vtkSmartPointer.h>
#include <vtkXMLRectilinearGridReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <iostream>
#include <cassert>
#include <string>
#include <numeric>
#include <set>

using std::pair;
using std::array;
using std::vector;
using std::string;
using std::set;

//------------------------------------------------------------
// This is an example of initializing the state variables (V, ID, Phi)
// using an M2C solution snapshot, possibly obtained on a different mesh.
// To compile this example, run "cmake ." followed by "make".
//------------------------------------------------------------

class MyStateCalculator : public UserDefinedState{

public:
  void GetUserDefinedState(int i0, int j0, int k0, int imax, int jmax, int kmax,
                           Vec3D*** coords, Vec5D*** v, double*** id, vector<double***> phi);
};

//------------------------------------------------------------

//-------------------------------------------------------------------------

template<typename T>
T Interpolate(vector<T> val, vector<bool>& valid, Vec3D &coeff) //!< value (instead of ref) of val is passed in
{
  int nValid = 8;
  for(int i=0; i<8; i++)
    if(!valid[i]) {
      val[i] = 0.0;
      nValid--;
    }
  assert(nValid>0);

  if(nValid<8) {
      T avg = (1.0/(double)nValid)*std::accumulate(val.begin(), val.end(), T(0.0));
      for(int i=0; i<8; i++)
        if(!valid[i])
          val[i] = avg;
  }

  return MathTools::trilinear_interpolation(coeff, val[0], val[1], val[2], val[3], val[4], val[5], val[6], val[7]);
}

//------------------------------------------------------------

void
MyStateCalculator::GetUserDefinedState(int i0, int j0, int k0, int imax, int jmax, int kmax,
                                       Vec3D*** coords, Vec5D*** v, double*** id, vector<double***> phi)
{

  // Note:
  // 1. Each processor core calls this functio with different inputs.
  // 2. Override the original values of v, id, and possibly phi. 
  // 3. If phi is updated, make sure it is consistent with id.
  // 4. The state variables (v) at each node are: rho,u,v,w,p.
  // 3. Do NOT change coords.

  //-------------------------------------------------------------
  double tol = 0.0025; //!< problem-specific (see usage below)
  int VAPOR_ID  = 1;
  int LIQUID_ID = 0;
  int INACTIVE_ID = 2; //!< problem-specific (inactive material id)
  //-------------------------------------------------------------

  //-------------------------------------------------------------
  // Read a binary solution file in VTR (VTK Rectilinear) format
  //-------------------------------------------------------------
  std::string filename = "IC/solution.vtr"; //Relative path from the place simulation is run.

  // Create a VTK reader for VTR files
  vtkSmartPointer<vtkXMLRectilinearGridReader> reader = vtkSmartPointer<vtkXMLRectilinearGridReader>::New();
  reader->SetFileName(filename.c_str());
  reader->Update();  // Read the file

  // Get the rectilinear grid from the reader
  vtkSmartPointer<vtkRectilinearGrid> grid = reader->GetOutput();
  if(!grid) {
    std::cerr << "Error reading VTR file." << std::endl;
    return;
  }

  // -----------------------
  // Extract coordinates
  // -----------------------
  vtkDataArray* xCoords = grid->GetXCoordinates();
  vtkDataArray* yCoords = grid->GetYCoordinates();
  vtkDataArray* zCoords = grid->GetZCoordinates();
  if (!xCoords || !yCoords || !zCoords) {
    std::cerr << "Error retrieving coordinates." << std::endl;
    return;
  }
  vector<double> xyz[3];
  for(int i=0; i<xCoords->GetNumberOfTuples(); i++)
    xyz[0].push_back(xCoords->GetComponent(i,0));
  for(int i=0; i<yCoords->GetNumberOfTuples(); i++)
    xyz[1].push_back(yCoords->GetComponent(i,0));
  for(int i=0; i<zCoords->GetNumberOfTuples(); i++)
    xyz[2].push_back(zCoords->GetComponent(i,0));


  // -----------------------
  // Extract state variables
  // -----------------------
  int nVar = grid->GetPointData()->GetNumberOfArrays();
  assert(nVar == 7); //'density', 'velocity', 'pressure', 'materialid', 'levelset0', 'temperature', 'laser_radiance'
  auto density = grid->GetPointData()->GetArray(0);
  assert((string)density->GetName() == "density");
  auto velocity = grid->GetPointData()->GetArray(1);
  assert((string)velocity->GetName() == "velocity");
  auto pressure = grid->GetPointData()->GetArray(2);
  assert((string)pressure->GetName() == "pressure");
  auto matid = grid->GetPointData()->GetArray(3);
  assert((string)matid->GetName() == "materialid");
  auto levelset = grid->GetPointData()->GetArray(4);
  assert((string)levelset->GetName() == "levelset0");
  

  //-------------------------------------------------------------
  // Specify v, id, and phi through interpolation
  //-------------------------------------------------------------
  Int3 ijk(0);
  Vec3D coeff;
  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        if(id[k][j][i] == INACTIVE_ID)
          continue; //this node is covered by solid body; nothing to be done.

        //-------------------------------
        // Find relative location & coordinates
        //-------------------------------
        bool found[3] = {false,false,false}; 
        for(int dim=0; dim<3; dim++) {
          double x = coords[k][j][i][dim];
          if(x<xyz[dim].front()) {
            if(x>=xyz[dim].front()-tol) {
              found[dim] = true;
              ijk[dim] = 0;
              coeff[dim] = 0.0;
            }
          }
          else if(x>=xyz[dim].back()) {
            if(x<xyz[dim].back()+tol) {
              found[dim] = true;
              ijk[dim] = xyz[dim].size()-1;
              coeff[dim] = 0.0;
            }
          }
          else { //x is in [min, max)
            for(int p=1; p<(int)xyz[dim].size(); p++) {
              if(x<xyz[dim][p]) {
                found[dim] = true;
                ijk[dim] = p-1;
                coeff[dim] = (x - xyz[dim][p-1]) / (xyz[dim][p] - xyz[dim][p-1]);
                break;
              }
            }
          }
        }
        
        if(!found[0] || !found[1] || !found[2])
          continue; //skip this node

        
        //-------------------------------
        // Get neighbors (some of them may be invalid / out of domain)
        //-------------------------------
        vector<bool> valid(8, false);
        valid[0] = true;                                //(i,  j,  k)
        valid[1] = ijk[0]+1 < xyz[0].size();            //(i+1,j,  k)
        valid[2] = ijk[1]+1 < xyz[1].size();            //(i,  j+1,k)
        valid[3] = valid[1] && valid[2];                //(i+1,j+1,k)
        valid[4] = ijk[2]+1 < xyz[2].size();              //(i,  j,  k+1)
        valid[5] = valid[4] && valid[1];                //(i+1,j,  k+1) 
        valid[6] = valid[4] && valid[2];                //(i,  j+1,k+1)
        valid[7] = valid[4] && valid[3];                //(i+1,j+1,k+1)

        vector<Vec5D> v_neigh(8, 0.0);
        vector<double> id_neigh(8, 0); 
        vector<double> phi_neigh(8, 0.0);

        int counter = 0;
        set<int> ids;
        for(int kk=0; kk<2; kk++)
          for(int jj=0; jj<2; jj++)
            for(int ii=0; ii<2; ii++) {
              if(!valid[counter])
                continue;

              int index = (ijk[2]+kk)*(xyz[0].size()*xyz[1].size())
                        + (ijk[1]+jj)*xyz[0].size() + ijk[0]+ii;

              int myid = matid->GetTuple(index)[0];

              if(myid == INACTIVE_ID){
                valid[counter] = false;
                continue;
              }

              ids.insert(myid);

              v_neigh[counter][0] = density->GetTuple(index)[0];
              v_neigh[counter][1] = velocity->GetTuple(index)[0];
              v_neigh[counter][2] = velocity->GetTuple(index)[1];
              v_neigh[counter][3] = velocity->GetTuple(index)[2];
              v_neigh[counter][4] = pressure->GetTuple(index)[0];
              id_neigh[counter]   = myid;
              phi_neigh[counter]  = levelset->GetTuple(index)[0];

              counter++;
            }

        if(ids.empty())
          continue; //no valid neighbor. nothing can be done.

        if(ids.size()==1) {//far from interface, straightforward interpolation
          id[k][j][i] = *ids.begin(); //the only valid id in the neighborhood
          v[k][j][i]  = Interpolate(v_neigh, valid, coeff);
          (phi[0])[k][j][i] = Interpolate(phi_neigh, valid, coeff);
        }
        else {
          (phi[0])[k][j][i] = Interpolate(phi_neigh, valid, coeff);
          if((phi[0])[k][j][i]<0) //inside bubble
            id[k][j][i] = VAPOR_ID;
          else
            id[k][j][i] = LIQUID_ID; //outside bubble (and not covered by solid)

          bool has_same_id_neighbor = false;
          vector<bool> same_id(8, true);
          for(int p=0; p<8; p++) {
            if(id_neigh[p] != id[k][j][i])
              same_id[p] = false;
            else
              has_same_id_neighbor = true;
          }
          assert(has_same_id_neighbor);
   
          v[k][j][i]  = Interpolate(v_neigh, same_id, coeff); //avoid crossing interface

          //interpolate pressure across interface
          vector<double> p_neigh(8, 0.0);
          for(int p=0; p<8; p++)
            p_neigh[p] = v_neigh[p][4];

          v[k][j][i][4]  = Interpolate(p_neigh, valid, coeff); //avoid crossing interface
        }
      }
}


//------------------------------------------------------------
// The class factory (Note: Do NOT change these functions except the word "MyStateCalculator".)
extern "C" UserDefinedState* Create() {
  return new MyStateCalculator; //TODO: If you've changed the name of the derived class, this needs to be updated.
}

extern "C" void Destroy(UserDefinedState* uds) {
  delete uds;
}

//------------------------------------------------------------

