using Test

using MMJMesh
using MMJMesh.Meshes
using MMJMesh.Utilities

m = makemeshoninterval(0.0, 1.2, 10)
@test nedges(m) == 10

a = 3
m = makemeshonrectangle(9.0, 4.5, 2a, a)
@test nfaces(m) == 2a * a
