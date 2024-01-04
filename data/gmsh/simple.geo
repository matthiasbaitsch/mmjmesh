SetFactory("OpenCASCADE");

// Points
Point(1) = {0, 0, 0, 0};
Point(2) = {1, 0, 0, 0};
Point(3) = {0.75, 1, 0, 0};
Point(4) = {0.25, 1, 0, 0};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Connect lines
Curve Loop(1) = {1, 2, 3, 4};

// Surface
Plane Surface(1) = {1};

// Mesh
MeshSize{:} = 0.75;
Mesh 2;

// Save
Save "simple.msh";




