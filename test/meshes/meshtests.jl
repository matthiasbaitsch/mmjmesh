module MeshTests

using Test
using MMJMesh
using MMJMesh.Meshes
using MMJMesh.Utilities

a = 3
m = makemeshonrectangle(9.0, 4.5, 2a, a)

@test nentities(m, 0) == (a + 1) * (2a + 1)
@test nentities(m, 1) == (a + 1) * 2a + a * (2a + 1)
@test nentities(m, 2) == 2a^2

end