// INSTRUCTIONS
// 1. Run in command prompt with the following:
// gmsh -3 -optimize -format msh4 -o xxx.mesh this_file.geo
// 2. Run in command prompt with the following:
// gmsh2aerof ...

// Unit: mm
// INPUTS
xmin = -9.1;
xmid = -9;
xmax = -5.0;
ymin = 0.0;
ymax = 0.6;
zmin = -0.02;
zmax = 0.02;
dx = 0.2;

Point(1) = {xmin, ymin, zmin, dx};
Point(2) = {xmax, ymin, zmin, dx};
Point(3) = {xmid, ymax, zmin, dx};
Point(4) = {xmin, ymax, zmin, dx};
Point(5) = {xmin, ymin, zmax, dx};
Point(6) = {xmax, ymin, zmax, dx};
Point(7) = {xmid, ymax, zmax, dx};
Point(8) = {xmin, ymax, zmax, dx};
Line(1) = {1, 4};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {1, 5};
Line(5) = {4, 8};
Line(6) = {3, 7};
Line(7) = {2, 6};
Line(8) = {5, 8};
Line(9) = {8, 7};
Line(10)= {7, 6};
Line Loop(1) = {1,5,-8,-4};
Line Loop(2) = {2,6,-9,-5};
Line Loop(3) = {3,7,-10,-6};
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};

// Define boundary conditions with GMSH tags
//  1 = OutletFixedSurface_1
//  2 = SymmetrySurface_2
//  3 = StickFixedSurface_3
//  4 = InletFixedSurface_4

// Symmetry
Physical Surface(1) = {1,2,3};

// Saves mesh in .MESH format for conversion to AERO-F format
//Mesh.SaveElementTagType=2;
