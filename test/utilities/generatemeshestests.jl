using Test
using DomainSets
using DomainSets: ×
using LinearAlgebra

using MMJMesh
using MMJMesh.Meshes
using MMJMesh.Utilities
using MMJMesh.Geometries
using MMJMesh.Mathematics

m = makemeshoninterval(0.0, 1.2, 10)
@test nedges(m) == 10

a = 3
m = makemeshonrectangle(9.0, 4.5, 2a, a)
@test nfaces(m) == 2a * a
@test coordinates(m, 1) == [0, 0]

m = makemeshonrectangle(9.0, 4.5, 2a, a, TRIANGLE)
@test nfaces(m) == 2 * 2a * a

m = makemeshonrectangle((1 .. 3) × (4 .. 5), 2)
@test nfaces(m) == 2
@test coordinates(m, 1) == [1, 4]
@test coordinates(m, 6) == [3, 5]

m = makemeshonrectangle((1 .. 3) × (4 .. 5), 2, 2)
@test nfaces(m) == 4
@test coordinates(m, 1) == [1, 4]
@test coordinates(m, 9) == [3, 5]

cartesianfrompolar(x) = x[1] * [cos(x[2]), sin(x[2])]
m = makemeshonrectangle((1 .. 2) × (0 .. (3 / 2)π), 10, 70, gmap=cartesianfrompolar)
@test coordinates(m, 11 * 71) ≈ [0, -2]
@test geometrytype(face(m, 1)) == GeometricObjectI

m = makemeshonrectangle(QHat, 2, 3, gmap=AffineMapping(Diagonal([2, 3]), [2, 1]));
@test coordinates(m, 12) == [4, 4]
@test geometrytype(face(m, 1)) == Box
