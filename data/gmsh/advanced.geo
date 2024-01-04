SetFactory("OpenCASCADE");

// Rectangles
Rectangle(1) = {0, 0, 0, 1, 1, 0};
Rectangle(2) = {1, 0, 0, 1, 1, 0};

// Glue rectangles together
Coherence;

// Physical names
Physical Curve("c1") = {3, 4, 7};
Physical Curve("c2") = {1, 5};
Physical Curve("c3") = {6};
Physical Surface("s1") = {1};
Physical Surface("s2") = {2};

// Mesh
MeshSize{:} = 0.5;
Mesh 2;

// Save
Save "advanced.msh";



