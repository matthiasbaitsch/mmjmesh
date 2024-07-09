using Test
using DomainSets
using DomainSets: ×

using MMJMesh
using MMJMesh.Meshes
using MMJMesh.Utilities


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
