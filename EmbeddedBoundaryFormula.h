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

  double eps; //geometry tolerance (length)

  enum Operation {MIRRORING = 0, SIZE = 1} operation;

  enum ImageScenario {NODE = 0, EDGE_SHARED = 1, EDGE_REMOTE = 2,
                      FACE_SHARED = 3, FACE_REMOTE = 4, 
                      ELEMENT_SHARED = 5, ELEMENT_REMOTE = 6, SIZE = 7};
  ImageScenario scenario;

  std::vector<Int3> node;
  std::vector<double> coeff; //!< Note: sizes of node and coeff are usually NOT the same
  double constant; //!< constant term in the affine function (often 0).
  
public:

  EmbeddedBoundaryFormula(double eps_ = 1.0e-12);
  ~EmbeddedBoundaryFormula();

  //! Explicitly specify the formula
  void SpecifyFormula(Operation op, ImageScenario sc, std::vector<Int3> &node_, std::vector<double> &coeff_,
                      double constant_ = 0.0);

  //! Build the formula
  void BuildMirroringFormula(Int3& ghost, Int3& image, Vec3D& xi);

  //! Get Info
  Operation            GetOperationType() {return operation;}
  ImageScenario        GetImageScenario() {return scenario;}
  int                  GetSize()          {return N;}
  std::vector<Int3>&   GetNodes()         {return node;}
  std::vector<double>& GetCoefficients()  {return coeff;}
  double               GetConstant()      {return constant;}

  //! Apply the formula
  double Evaluate(double*** v); //!< the 3D data structure (stored in SpaceVariable3D)
  void Evaluate(double*** vin, double* vout, int dof);

  //! Apply the formula
  template <typename T>
  T Evaluate(std::vector<T>& v); //!< v must have size N, and the correct order.

private:

  //! Functions used in BuildMirroringFormula
  void FinalizeMirroringFormulaEdge(Int3 &ghost, std::vector<Int3> &node_tmp, std::vector<double>& xi);
  void FinalizeMirroringFormulaFace(Int3 &ghost, std::vector<Int3> &node_tmp, std::vector<double>& xi);
  void FinalizeMirroringFormulaElem(Int3 &ghost, std::vector<Int3> &node_tmp, std::vector<double>& xi);

};

#endif
