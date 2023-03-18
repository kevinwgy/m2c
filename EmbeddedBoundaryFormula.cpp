/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include<EmbeddedBoundaryFormula.h>
#include<cassert>

//if the verification script is activated, turn following line on
//#include<trilinear_interpolation.h>

//-----------------------------------------------------------------------------

EmbeddedBoundaryFormula::EmbeddedBoundaryFormula(Operation operation_, double eps_)
                       : operation(operation_), scenario(NODE), constant(0.0),
                         eps(eps_)
{ }

//-----------------------------------------------------------------------------

EmbeddedBoundaryFormula::~EmbeddedBoundaryFormula()
{ }

//-----------------------------------------------------------------------------

void
EmbeddedBoundaryFormula::SpecifyFormula(Operation op, ImageScenario sc, std::vector<Int3> &node_, 
                                        std::vector<double> &coeff_, double constant_)
{
  assert(node_.size() == coeff_.size());

  operation = op;
  scenario  = sc;
  node      = node_;
  coeff     = coeff_;
  constant  = constant_;

}

//-----------------------------------------------------------------------------

void
EmbeddedBoundaryFormula::BuildMirroringFormula(Int3& ghost, Int3& image, Vec3D& xi0, double constant_)
{
  operation = MIRRORING;
  constant = constant_;

  node.clear();
  coeff.clear();

  //---------------------------------------------------------------------
  // Step 1: Check xi0 to identify degenerate cases (node, edge, face)
  //---------------------------------------------------------------------
  int involved[2][2][2];
  for(int k=0; k<2; k++)
    for(int j=0; j<2; j++)
      for(int i=0; i<2; i++)
        involved[k][j][i] = 1; //by default, all the 8 nodes are involved

  std::vector<double> xi;

  if(fabs(xi0[0])<=eps || fabs(xi0[0]-1.0)<eps) {
    int i = fabs(xi0[0])<=eps ? 1 : 0; 
    for(int k=0; k<2; k++)
      for(int j=0; j<2; j++)
        involved[k][j][i] = 0;
  } else 
    xi.push_back(xi0[0]);
    
  if(fabs(xi0[1])<=eps || fabs(xi0[1]-1.0)<eps) {
    int j = fabs(xi0[1])<=eps ? 1 : 0; 
    for(int k=0; k<2; k++)
      for(int i=0; i<2; i++)
        involved[k][j][i] = 0;
  } else
    xi.push_back(xi0[1]);
 
  if(fabs(xi0[2])<=eps || fabs(xi0[2]-1.0)<eps) {
    int k = fabs(xi0[2])<=eps ? 1 : 0; 
    for(int j=0; j<2; j++)
      for(int i=0; i<2; i++)
        involved[k][j][i] = 0;
  } else
    xi.push_back(xi0[2]);


  //---------------------------------------------------------------------
  // Step 2: Deal with difference cases separately
  //---------------------------------------------------------------------
  std::vector<Int3> node_tmp;
  for(int k=0; k<2; k++)
    for(int j=0; j<2; j++)
      for(int i=0; i<2; i++)
        if(involved[k][j][i])
          node_tmp.push_back(image + Int3(i,j,k));

  switch(node_tmp.size()) {
    case 1: //node
      assert(xi.size() == 0);
      node = node_tmp;  
      coeff.push_back(1.0);   
      scenario = NODE;
      break;
    case 2: //edge
      assert(xi.size() == 1);
      FinalizeMirroringFormulaEdge(ghost, node_tmp, xi);
      break; 
    case 4: //face
      assert(xi.size() == 2);
      FinalizeMirroringFormulaFace(ghost, node_tmp, xi);
      break;
    case 8: //element (interior)
      assert(xi.size() == 3);
      FinalizeMirroringFormulaElem(ghost, node_tmp, xi);
      break;
    default:
      fprintf(stdout,"*** Error: Detected error in the construction of mirroring formulas. Ghost node: %d %d %d.\n",
                     ghost[0], ghost[1], ghost[2]);
      break;
  }
}

//-----------------------------------------------------------------------------

void
EmbeddedBoundaryFormula::BuildLinearExtrapolationFormula(Int3& ghost, Int3& image, Vec3D& xi0,
                                                         double constant_)
{
  //everything is the same as BuildMirroringFormula, except the "operation" type;
  BuildMirroringFormula(ghost, image, xi0, constant_); 
  operation = LINEAR_EXTRAPOLATION;
}

//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------

void
EmbeddedBoundaryFormula::FinalizeMirroringFormulaEdge(Int3 &ghost, std::vector<Int3> &node_tmp, 
                                                      std::vector<double>& xi)
{
  if(ghost != node_tmp[0] && ghost != node_tmp[1]) {
    scenario = EDGE_REMOTE;
    node = node_tmp;
    coeff.push_back(1-xi[0]);
    coeff.push_back(xi[0]);
    //Formula: v(ghost) = coeff[0]*v(node[0]) + coeff[1]*v(node[1])
  }  
  else {
    scenario = EDGE_SHARED;
    if(ghost == node_tmp[0])
      node.push_back(node_tmp[1]);        
    else
      node.push_back(node_tmp[0]);
    coeff.push_back(1.0);
    //Formula: v(ghost) = coeff[0]*v(node[0])   (where coeff[0] = 1.0)
  }
}

//-----------------------------------------------------------------------------

void
EmbeddedBoundaryFormula::FinalizeMirroringFormulaFace(Int3 &ghost, std::vector<Int3> &node_tmp, 
                                                      std::vector<double>& xi)
{
  std::vector<double> coeff_tmp;
  coeff_tmp.push_back((1-xi[1])*(1-xi[0]));
  coeff_tmp.push_back((1-xi[1])*xi[0]);
  coeff_tmp.push_back(xi[1]*(1-xi[0]));
  coeff_tmp.push_back(xi[1]*xi[0]);

  scenario = FACE_REMOTE;
  double denom = 0;
  for(int i=0; i<(int)node_tmp.size(); i++) {
    if(node_tmp[i] != ghost) {
      node.push_back(node_tmp[i]);
      coeff.push_back(coeff_tmp[i]);
    } else {
      denom = 1.0 - coeff_tmp[i];
      scenario = FACE_SHARED;
    }
  }
  if(scenario == FACE_SHARED) {
    assert(denom != 0.0);
    for(int i=0; i<(int)coeff.size(); i++)
      coeff[i] /= denom; 
  }
  //Formula: v(ghost) = sum ( coeff[i]*v(node[i]) )
}

//-----------------------------------------------------------------------------

void
EmbeddedBoundaryFormula::FinalizeMirroringFormulaElem(Int3 &ghost, std::vector<Int3> &node_tmp, 
                                                      std::vector<double>& xi)
{
  std::vector<double> coeff_tmp(8, 0.0);
  coeff_tmp[0] = (1-xi[0])*(1-xi[1])*(1-xi[2]);
  coeff_tmp[1] = xi[0]*(1-xi[1])*(1-xi[2]);
  coeff_tmp[2] = (1-xi[0])*xi[1]*(1-xi[2]);
  coeff_tmp[3] = xi[0]*xi[1]*(1-xi[2]);
  coeff_tmp[4] = (1-xi[0])*(1-xi[1])*xi[2];
  coeff_tmp[5] = xi[0]*(1-xi[1])*xi[2];
  coeff_tmp[6] = (1-xi[0])*xi[1]*xi[2];
  coeff_tmp[7] = xi[0]*xi[1]*xi[2];

  scenario = ELEMENT_REMOTE;
  double denom = 0;
  for(int i=0; i<(int)node_tmp.size(); i++) {
    if(node_tmp[i] != ghost) {
      node.push_back(node_tmp[i]);
      coeff.push_back(coeff_tmp[i]);
    } else {
      denom = 1.0 - coeff_tmp[i];
      scenario = ELEMENT_SHARED;
    }
  }
  if(scenario == ELEMENT_SHARED) {
    assert(denom != 0.0);
    for(int i=0; i<(int)coeff.size(); i++)
      coeff[i] /= denom; 
  }
  //Formula: v(ghost) = sum ( coeff[i]*v(node[i]) )

/*
  //verification
  std::vector<double> coeff2(8, 0.0);
  coeff2[0] = MathTools::trilinear_interpolation<double>(1,0,0,0,0,0,0,0,xi.data()); 
  coeff2[1] = MathTools::trilinear_interpolation<double>(0,1,0,0,0,0,0,0,xi.data()); 
  coeff2[2] = MathTools::trilinear_interpolation<double>(0,0,1,0,0,0,0,0,xi.data()); 
  coeff2[3] = MathTools::trilinear_interpolation<double>(0,0,0,1,0,0,0,0,xi.data()); 
  coeff2[4] = MathTools::trilinear_interpolation<double>(0,0,0,0,1,0,0,0,xi.data()); 
  coeff2[5] = MathTools::trilinear_interpolation<double>(0,0,0,0,0,1,0,0,xi.data()); 
  coeff2[6] = MathTools::trilinear_interpolation<double>(0,0,0,0,0,0,1,0,xi.data()); 
  coeff2[7] = MathTools::trilinear_interpolation<double>(0,0,0,0,0,0,0,1,xi.data());
  for(int i=0; i<8; i++)
    assert(fabs(coeff2[i]-coeff_tmp[i])<1.0e-14);
*/
}

//-----------------------------------------------------------------------------

double
EmbeddedBoundaryFormula::Evaluate(double*** v, double vin)
{

  double res = 0.0;

  if(operation==MIRRORING) {
    assert(vin==0.0); //should be specified as "constant"
    res = constant;
    for(int i=0; i<(int)node.size(); i++)   
      res += coeff[i]*v[node[i][2]/*k*/][node[i][1]/*j*/][node[i][0]/*i*/];
  }
  else if(operation==LINEAR_EXTRAPOLATION) {
    for(int i=0; i<(int)node.size(); i++)   
      res += coeff[i]*v[node[i][2]/*k*/][node[i][1]/*j*/][node[i][0]/*i*/];
    res = (vin - constant*res)/(1.0-constant); //constant should be at most 0.5.
  }

  return res;
}

//-----------------------------------------------------------------------------

Vec3D
EmbeddedBoundaryFormula::Evaluate3D(Vec3D*** v, Vec3D vin)
{

  Vec3D res(0.0);

  if(operation==MIRRORING) {
    assert(vin[0]==vin[1] && vin[1]==vin[2] && vin[2]==0.0); //should be specified as "constant"
    res = constant;
    for(int i=0; i<(int)node.size(); i++)   
      res += coeff[i]*v[node[i][2]/*k*/][node[i][1]/*j*/][node[i][0]/*i*/];
  }
  else if(operation==LINEAR_EXTRAPOLATION) {
    for(int i=0; i<(int)node.size(); i++)   
      res += coeff[i]*v[node[i][2]/*k*/][node[i][1]/*j*/][node[i][0]/*i*/];
    res = (vin - constant*res)/(1.0-constant); //constant should be at most 0.5.
  }

  return res;
}

//-----------------------------------------------------------------------------

double
EmbeddedBoundaryFormula::Evaluate(std::vector<double>& v, double vin)
{
  assert(v.size() == node.size());
  double res = 0.0; 

  if(operation==MIRRORING) {
    assert(vin==0.0); //should be specified as "constant"
    res = constant;
    for(int i=0; i<(int)v.size(); i++)   
      res += coeff[i]*v[i];
  }
  else if(operation==LINEAR_EXTRAPOLATION) {
    for(int i=0; i<(int)v.size(); i++)   
      res += coeff[i]*v[i];
    res = (vin - constant*res)/(1.0-constant); //constant should be at most 0.5.
  }

  return res;
}

//-----------------------------------------------------------------------------

Vec3D
EmbeddedBoundaryFormula::Evaluate3D(std::vector<Vec3D>& v, Vec3D vin)
{
  assert(v.size() == node.size());
  Vec3D res(0.0);

  if(operation==MIRRORING) {
    assert(vin[0]==vin[1] && vin[1]==vin[2] && vin[2]==0.0); //should be specified as "constant"
    res = constant;
    for(int i=0; i<(int)v.size(); i++)   
      res += coeff[i]*v[i];
  }
  else if(operation==LINEAR_EXTRAPOLATION) {
    for(int i=0; i<(int)v.size(); i++)   
      res += coeff[i]*v[i];
    res = (vin - constant*res)/(1.0-constant); //constant should be at most 0.5.
  }

  return res;
}

//-----------------------------------------------------------------------------




//-----------------------------------------------------------------------------






