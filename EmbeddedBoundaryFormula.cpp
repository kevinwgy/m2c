#include<EmbeddedBoundaryFormula.h>
#include<trilinear_interpolation.h>
#include<cassert>

//-----------------------------------------------------------------------------

EmbeddedBoundaryFormula::EmbeddedBoundaryFormula(double eps_)
                       : operation(MIRRORING), scenario(NODE), constant(0.0),
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
  assert(nodes_.size() == coeff_.size());

  operation = op;
  scenario  = sc;
  node      = node_;
  coeff     = coeff_;
  constant  = constant_;

}

//-----------------------------------------------------------------------------

void
EmbeddedBoundaryFormula::BuildMirroringFormula(Int3& ghost, Int3& image, Vec3D& xi0)
{
  operation = MIRRORING;
  constant = 0.0;

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
      fprintf(stderr,"*** Error: Detected error in the construction of mirroring formulas. Ghost node: %d %d %d.\n",
                     ghost[0], ghost[1], ghost[2]);
      break;
  }
}

//-----------------------------------------------------------------------------

void
EmbeddedBoundaryFormula::FinalizeMirroringFormulaEdge(Int3 &ghost, std::vector<Int3> &node_tmp, std::vector<double>& xi)
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
    coeff.push_back(xi[0]);
    //Formula: v(ghost) = coeff[0]*
  }
}

//-----------------------------------------------------------------------------






