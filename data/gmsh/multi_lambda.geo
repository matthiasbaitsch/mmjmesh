SetFactory("OpenCASCADE");

Rectangle(1) = {0, 0, 0, 1, 1, 0};
Rectangle(2) = {1, 0, 0, 1, 1, 0};
Rectangle(3) = {2, 0, 0, 1, 1, 0};
Rectangle(4) = {0, 1, 0, 1, 1, 0};
Rectangle(5) = {1, 1, 0, 1, 1, 0};
Rectangle(6) = {2, 1, 0, 1, 1, 0};
Rectangle(7) = {0, 2, 0, 1, 1, 0};
Rectangle(8) = {1, 2, 0, 1, 1, 0};
Rectangle(9) = {2, 2, 0, 1, 1, 0};
Coherence;

Physical Curve("ΓD") = {4, 13, 20, 9, 16, 23};
Physical Surface("s1") = {1};
Physical Surface("s2") = {2};
Physical Surface("s3") = {3};
Physical Surface("s4") = {4};
Physical Surface("s5") = {5};
Physical Surface("s6") = {6};
Physical Surface("s7") = {7};
Physical Surface("s8") = {8};
Physical Surface("s9") = {9};

MeshSize{:} = 0.1;
Mesh 2;

Save "multi_lambda.msh";

