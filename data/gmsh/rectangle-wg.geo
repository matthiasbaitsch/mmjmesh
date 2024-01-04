SetFactory("OpenCASCADE");

Rectangle(1) = {0, 0, 0, 2, 1, 0};
Physical Curve("C1") = {1, 2, 4};
Physical Curve("C2") = {1, 3};
Physical Curve("C3") = {2, 3};
Physical Surface("S1") = {1};
MeshSize{:} = 0.5;
Mesh 2;
Save "rectangle-wg.msh";
