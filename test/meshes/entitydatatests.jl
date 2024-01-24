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

# Entities
e1 = edge(m, 1)
f13 = element(m, 13)
f42 = element(m, 42)
f45 = element(m, 45)

# Associate with everything
m.data[:foo] = 42
@test m.data[:foo] == 42
@test f42.data[:foo] == 42

# Associate with group
m.data[:bar, :g1] = 13
m.data[:bar, :g2] = 81
@test f13.data[:bar]  == 13
@test isnothing(f42.data[:bar])
@test f45.data[:bar]  == 81

# Associate with faces
m.data[:baz, :faces] = 91

@test f13.data[:baz] == 91
@test isnothing(e1.data[:baz])