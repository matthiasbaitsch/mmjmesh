using Test

using MMJMesh
using MMJMesh.Meshes
using MMJMesh.Utilities


# -------------------------------------------------------------------------------------------------
# Associate data
# -------------------------------------------------------------------------------------------------
a = 5
m = makemeshonrectangle(9.0, 4.5, 2a, a)

# Groups
m.groups[:g1] = FaceGroup([1, 2, 3, 13, 14, 7, 15])
m.groups[:g2] = FaceGroup([45, 56])

# Elements
e13 = element(m, 13)
e42 = element(m, 42)
e45 = element(m, 45)

# Associate with everything
m.data[:foo] = 42
@test m.data[:foo] == 42
@test e42.data[:foo] == 42

# Associate with group
m.data[:bar, :g1] = 13
m.data[:bar, :g2] = 81
@test e13.data[:bar]  == 13
@test isnothing(e42.data[:bar])
@test e45.data[:bar]  == 81