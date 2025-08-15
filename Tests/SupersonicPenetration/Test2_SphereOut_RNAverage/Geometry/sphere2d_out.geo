// INSTRUCTIONS
// 1. Run in command prompt with the following:
// gmsh -3 -optimize -format msh4 -o xxx.mesh this_file.geo
// 2. Run in command prompt with the following:
// gmsh2aerof ...

SetFactory("OpenCASCADE");

// Unit: mm
// INPUTS
zmin = -0.02;
zmax = 0.02;
x0 = -3.5;
r = 1.0;
dx = 0.05;
Point(1) = {x0,   0, zmin, dx};
Point(2) = {x0+r, 0, zmin, dx};
Point(3) = {x0-r, 0, zmin, dx};

Circle(1) = {3,1,2};

out[] = Extrude {0, 0, zmax-zmin} {Line{1};};

Physical Surface(1) = {out[1]};

//Mesh.SaveElementTagType=2;
