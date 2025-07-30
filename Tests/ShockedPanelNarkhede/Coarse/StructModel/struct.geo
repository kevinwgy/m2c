// Unit: mm 

// INPUTS
dx = 0.5;
L = 50;
w = 1; // through plane thickness
t = 1; // in plane thickness

// beam
Point(1) = { 0.0, 0.0, -w, dx};
Point(2) = { 0.0,   L, -w, dx}; 
Point(3) = { 0.0,   L,  w, dx}; 
Point(4) = { 0.0, 0.0,  w, dx};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {1, 4};

Line Loop(1) = {1, 2, 3, -4};
Plane Surface(1) = {1};
Transfinite Curve{1} = 201;
Transfinite Curve{2} = 2;
Transfinite Curve{3} = 201;
Transfinite Curve{4} = 2;
Transfinite Surface{1};
Recombine Surface{1};

// extrude
Ex[] = Extrude {t, 0, 0} {
  Surface{1};
  Layers{2};
  Recombine;
};

// symmetry surfaces (top and bottom)
Physical Surface(1) = {Ex[2], Ex[4]};

// fixed surface
Physical Surface(2) = {Ex[5]};

// wetted surfaces (FS interface)
Physical Surface(3) = {1, Ex[3], Ex[0]};

// volume
Physical Volume(1) = {Ex[1]}; 

// Smoothing the final mesh
Mesh.Smoothing = 20;
Mesh.SmoothRatio = 3;
