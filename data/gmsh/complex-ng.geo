SetFactory("OpenCASCADE");

// rectangles
Rectangle(1) = {0, 0, 0, 0.5, 1, 0};
Rectangle(2) = {0.5, 0.3, 0, 0.25, 0.6, 0};
Rectangle(3) = {0.75, 0, 0, 0.5, 1, 0};
Rectangle(4) = {1.25, 0.5, 0, 0.5, 0.5, 0};

// glue rectangles together
Coherence;

// mesh
MeshSize{:} = 0.05;
Mesh 2;

Save "complex-ng.msh";



