using Test
using LinearAlgebra

using MMJMesh
using MMJMesh.Meshes
using MMJMesh.Utilities
using MMJMesh.Geometries
using MMJMesh.Mathematics

m = Mesh(0.0 .. 1.2, 10)
@test nedges(m) == 10

a = 3
m = Mesh((0 .. 9.0) × (0 .. 4.5), 2a, a)
@test nfaces(m) == 2a * a
@test coordinates(m, 1) == [0, 0]

m = Mesh((0 .. 9.0) × (0 .. 4.5), 2a, a, meshtype=TRIANGLE)
@test nfaces(m) == 2 * 2a * a

m = Mesh((1 .. 3) × (4 .. 5), 2)
@test nfaces(m) == 2
@test coordinates(m, 1) == [1, 4]
@test coordinates(m, 6) == [3, 5]

m = Mesh((1 .. 3) × (4 .. 5), 2, 2)
@test nfaces(m) == 4
@test coordinates(m, 1) == [1, 4]
@test coordinates(m, 9) == [3, 5]

m = Mesh(4.0, 2.0, 4, meshtype=CRISSCROSS)
@test nfaces(m) == 4 * 8

cartesianfrompolar(x) = x[1] * [cos(x[2]), sin(x[2])]
m = Mesh((1 .. 2) × (0 .. (3 / 2)π), 10, 70, gmap=cartesianfrompolar)
@test coordinates(m, 11 * 71) ≈ [0, -2]
@test geometrytype(face(m, 1)) == GeometricObjectI

m = Mesh(QHat, 2, 3, gmap=AffineMapping(Diagonal([2, 3]), [2, 1]));
@test coordinates(m, 12) == [4, 4]
@test geometrytype(face(m, 1)) == Box
