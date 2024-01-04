SetFactory("OpenCASCADE");

// Rectangles
Rectangle(1) = {0, 0, 0, 0.5, 1, 0};
Rectangle(2) = {0.5, 0.3, 0, 0.25, 0.6, 0};
Rectangle(3) = {0.75, 0, 0, 0.5, 1, 0};
Rectangle(4) = {1.25, 0.5, 0, 0.5, 0.5, 0};

// Glue rectangles together
Coherence;

// Physical names
Physical Curve("c1") = {7, 17, 18};
Physical Curve("c2") = {5, 19, 14, 15};
Physical Surface("s1") = {1};
Physical Surface("s3") = {4};

// Create mesh
MeshSize{:} = 0.05;
Mesh 2;

Save "complex-g2.msh";



