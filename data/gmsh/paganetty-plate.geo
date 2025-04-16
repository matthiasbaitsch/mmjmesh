// Gmsh project created on Sun Apr 13 17:25:25 2025
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
Point(2) = {31, 0, 0, 1.0};
Point(3) = {31, -7, 0, 1.0};
Point(4) = {16.5, -7, 0, 1.0};
Point(5) = {0, -9, 0, 1.0};
Point(6) = {11, -9, 0, 1.0};
Point(7) = {27, -17, 0, 1.0};
Point(8) = {23.385, -20.795, 0, 1.0};
Point(9) = {26, -7, 0, 1.0};
Point(10) = {22, -7, 0, 1.0};
Point(11) = {23.37963268, -13.55138567, 0, 1.0};
Point(12) = {20.48294524, -10.79290749, 0, 1.0};
Point(13) = {0, -4, 0, 1.0};
Point(14) = {5, -9, 0, 1.0};
//+
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {7, 8};
Line(4) = {8, 6};
Line(5) = {6, 14};
Line(6) = {14, 5};
Line(7) = {5, 13};
Line(8) = {13, 1};
Line(9) = {11, 12};
Line(10) = {12, 10};
Line(11) = {10, 9};
Circle(12) = {3, 4, 7};
Circle(13) = {9, 4, 11};
//+
Curve Loop(1) = {1, 2, 12, 3, 4, 5, 6, 7, 8};
Curve Loop(2) = {10, 11, 13, 9};
Plane Surface(1) = {1, 2};
//+
Physical Curve("free", 14) = {6, 7, 10, 9, 11};
Physical Curve("fixed", 15) = {5, 4, 3, 12, 2, 1, 8, 13};
Physical Point("column", 17) = {5, 10, 12};
Physical Surface("plate", 16) = {1};
Recombine Surface {1};
//+
Mesh 2;
//+

