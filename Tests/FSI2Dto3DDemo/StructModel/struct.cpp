#include <iostream>
#include <gmsh.h>
#include <string>
#include <math.h>
#include <vector>
#include <map>

static double AngleOfRevolution = M_PI/4;

int main(void) {
	/* 		Create a planar BSpline surface and 
 *  			revolve it by 30 degrees...........		*/

	std::string modelName("CakeSlice");
	gmsh::initialize();
	gmsh::model::add(modelName);

	/* 
 *		
 *		Only inner surface control points and thickness are taken as
 *		design variables for now. Point 2 has 2 degrees of freedom, i.e.
 *		it can move in x and y. Point 3 can move only in x with x>100.
 *
 *		Outer surface will have the same shape as inner surface but at an
 *		offset equal to the thickness of the chamber.
 *
 *		Capsule rotated such that x-y plane lies in the middle... required
 *		for axisymmetric analysis with m2c.
 *
 * 	*/

	// fixed quantities
	double lx = 330;
	double ly = 230;
	double dx = 10.0;
	
	// variable quantities
	//double var1 = var[0];  // x location of Point 2
	//double var2 = var[1];  // y location of Point 2
	//double var3 = var[2];  // x location of Point 3 -- y fixed
	//double var4 = var[3];   // chamber wall thickness
	double var1 = 330;  // x location of Point 2
	double var2 = 230;  // y location of Point 2
	double var3 = 100;  // x location of Point 3 -- y fixed
	double var4 = 30;   // chamber wall thickness

	/*		Inner Surface		 */
	double sina = std::sin(AngleOfRevolution/2);
	double cosa = std::cos(AngleOfRevolution/2);
	gmsh::model::occ::addPoint(    lx,         0,          0, dx,  1);
	gmsh::model::occ::addPoint(  var1, var2*cosa, -var2*sina, dx,  2);
	gmsh::model::occ::addPoint(  var3,   ly*cosa,   -ly*sina, dx,  3);
	gmsh::model::occ::addPoint(var3/2,   ly*cosa,   -ly*sina, dx,  4); // additional point placed to make the top side flatter
	gmsh::model::occ::addPoint(     0,   ly*cosa,   -ly*sina, dx,  5);

	/*		Outer Surface		 */
	gmsh::model::occ::addPoint(  lx+var4,                0,                 0, dx, 6);
	gmsh::model::occ::addPoint(var1+var4, (var2+var4)*cosa, -(var2+var4)*sina, dx, 7);
	gmsh::model::occ::addPoint(     var3,   (ly+var4)*cosa,   -(ly+var4)*sina, dx, 8);
	gmsh::model::occ::addPoint(   var3/2,   (ly+var4)*cosa,   -(ly+var4)*sina, dx, 9); // additional point
	gmsh::model::occ::addPoint(        0,   (ly+var4)*cosa,   -(ly+var4)*sina, dx,10);

	int degreeU = 3; // parameter along the mold-line of capsule
	int degreeV = 1; // parameter along the thickness direction

	std::vector<int> pointTags = {1,2,3,4,5,6,7,8,9,10};

	// weights 
	std::vector<double> weight(10, 1.0);
	weight[1] = weight[6] = sqrt(2)/2.;

	std::vector<double> knotsU = {0,0.5,1}; // obtained from http://nurbscalculator.in/
	std::vector<double> knotsV = {0,1};
	std::vector<int> multU = {4,1,4};
	std::vector<int> multV = {2,2};

	// create base surface
	gmsh::model::occ::addBSplineSurface(pointTags, 5, 1, degreeU, degreeV, weight, knotsU, knotsV, multU, multV);

	gmsh::vectorpair capsule;
	gmsh::model::occ::revolve({{2,1}}, 0, 0, 0, 1, 0, 0, AngleOfRevolution, capsule);
	gmsh::model::occ::synchronize();

	//// copy surfaces for contact
	//gmsh::vectorpair contactInner, contactOuter;
	//gmsh::model::occ::copy({{2,2}}, contactInner);
	//gmsh::model::occ::copy({{2,4}}, contactOuter);
	//gmsh::model::occ::synchronize();

	gmsh::model::addPhysicalGroup(2, {3}, 1); // x=0
	gmsh::model::addPhysicalGroup(2, {1}, 2); // lmpc
	gmsh::model::addPhysicalGroup(2, {5}, 3); // lmpc
	gmsh::model::addPhysicalGroup(2, {2}, 4); // inner surface
	gmsh::model::addPhysicalGroup(2, {4}, 5); // outer surface
	//gmsh::model::addPhysicalGroup(2, {contactInner[0].second}, 6); // inner face sheet
	//gmsh::model::addPhysicalGroup(2, {contactOuter[0].second}, 7); // outer face sheet
	gmsh::model::addPhysicalGroup(3, {1}, 1);
	gmsh::option::setNumber("Mesh.Smoothing", 20);
	gmsh::option::setNumber("Mesh.SmoothRatio", 3);

	gmsh::model::occ::synchronize();

	gmsh::model::mesh::generate(3);

	gmsh::write(modelName+".msh");
	
	return 0;
}
