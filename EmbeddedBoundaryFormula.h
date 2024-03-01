/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _EMBEDDED_BOUNDARY_FORMULA_H_
#define _EMBEDDED_BOUNDARY_FORMULA_H_

#include<Vector3D.h>
#include<vector>

//------------------------------------------------
// class EmbeddedBoundaryFormula is responsible for
// building and applying an affine function (usually
// trilinear linear interpolation or extrapolation) 
// to update a ghost node.
//
// Note: We treat the degnerate cases separately, 
// where the image of the ghost node is a node, 
// along an edge, or on a face. This is NOT ncessary.
// One could use the trilinear interpolation formulas
// associated with "ELEMENT_SHARED" and "ELEMENT_REMOTE"
// to handle all the cases.
//------------------------------------------------

class EmbeddedBoundaryFormula {

public:

  enum Operation {MIRRORING = 0, LINEAR_EXTRAPOLATION = 1};

  enum ImageScenario {NODE = 0, EDGE_SHARED = 1, EDGE_REMOTE = 2,
                      FACE_SHARED = 3, FACE_REMOTE = 4, 
                      ELEMENT_SHARED = 5, ELEMENT_REMOTE = 6};

private: 

  double eps; //geometry tolerance (length)

  Operation operation;
  ImageScenario scenario;

  std::vector<Int3> node;
  std::vector<double> coeff; //!< Note: sizes of node and coeff are the same
  double constant; //!< MIRRORING: constant in affine function; EXTRAP: abs((ghost-interface)/(ghost-image))
  
public:

  EmbeddedBoundaryFormula(Operation operation_ = MIRRORING, double eps_ = 1.0e-12);
  ~EmbeddedBoundaryFormula();

  //! Explicitly specify the formula
  void SpecifyFormula(Operation op, ImageScenario sc, std::vector<Int3> &node_, std::vector<double> &coeff_,
                      double constant_ = 0.0);

  //! Build the mirroring formula
  void BuildMirroringFormula(Int3& ghost, Int3& image, Vec3D& xi, double constant_ = 0.0);

  //! Build the linear extrapolation formula
  void BuildLinearExtrapolationFormula(Int3& ghost, Int3& image, Vec3D& xi, double constant_); 

  //! Get Info
  Operation            GetOperationType() {return operation;}
  ImageScenario        GetImageScenario() {return scenario;}
  int                  GetSize()          {return node.size();}
  std::vector<Int3>&   GetNodes()         {return node;}
  std::vector<double>& GetCoefficients()  {return coeff;}
  double               GetConstant()      {return constant;}

  //! Apply the formula
  double Evaluate(double*** v, double vin=0.0); //!< the 3D data structure (stored in SpaceVariable3D)
  Vec3D Evaluate3D(Vec3D*** v, Vec3D vin=0.0);

  //! Apply the formula
  double Evaluate(std::vector<double>& v, double vin=0.0); //!< v must have the same size as "node", and the correct order.
  Vec3D Evaluate3D(std::vector<Vec3D>& v, Vec3D vin=0.0);

private:

  //! Functions used in BuildMirroringFormula
  void FinalizeMirroringFormulaEdge(Int3 &ghost, std::vector<Int3> &node_tmp, std::vector<double>& xi);
  void FinalizeMirroringFormulaFace(Int3 &ghost, std::vector<Int3> &node_tmp, std::vector<double>& xi);
  void FinalizeMirroringFormulaElem(Int3 &ghost, std::vector<Int3> &node_tmp, std::vector<double>& xi);

};

#endif
