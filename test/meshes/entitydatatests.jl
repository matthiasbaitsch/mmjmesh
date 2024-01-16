using Test

using MMJMesh
using MMJMesh.Meshes
using MMJMesh.Groups
using MMJMesh.Utilities


# -------------------------------------------------------------------------------------------------
# Associate data
# -------------------------------------------------------------------------------------------------
a = 5
m = makemeshonrectangle(9.0, 4.5, 2a, a)

m.groups[:g1] = FaceGroup([1, 2, 3, 13, 14, 7, 15])
m.groups[:g2] = FaceGroup([2, 13, 7, 5])
m.groups[:g3] = EdgeGroup([1, 2, 3, 13, 14, 7, 15])

e13 = element(m, 13)
e42 = element(m, 42)

# Associate with everything
m.data[:foo] = 42
@test m.data[:foo] == 42
# @test e42.data[:foo] == 42

# Associate with group
# m.data[:bar, :g1] = 13
# @test m.data[:bar, :g1] == 13

#@test e13.data[:bar] == 42
#@test e42.data[:bar] == nothing