using Test

using MMJMesh.Meshes
using MMJMesh.Groups
using MMJMesh.Utilities
import MMJMesh.Meshes: collectgroups


# -------------------------------------------------------------------------------------------------
# Entities and groups
# -------------------------------------------------------------------------------------------------
a = 5
m = makemeshonrectangle(9.0, 4.5, 2a, a)
m.groups[:g1] = FaceGroup([1, 2, 3, 13, 14, 7, 15])
m.groups[:g2] = FaceGroup([2, 13, 7, 5])
m.groups[:g3] = EdgeGroup([1, 2, 3, 13, 14, 7, 15])

e13 = element(m, 13)
e42 = element(m, 42)
@test e13 ∈ m.groups[:g1]
@test e13 ∉ m.groups[:g3]
@test e42 ∉ m.groups[:g1]


# -------------------------------------------------------------------------------------------------
# Boundary groups
# -------------------------------------------------------------------------------------------------
a = 5
m = makemeshonrectangle(9.0, 4.5, 2a, a)

@test m.groups.entries[:boundarynodes] == nothing
@test m.groups.entries[:boundaryedges] == nothing
@test length(m.groups[:boundarynodes]) == 30
@test length(m.groups[:boundaryedges]) == 30


# -------------------------------------------------------------------------------------------------
# Collectgroups
# -------------------------------------------------------------------------------------------------
a = 5
m = makemeshonrectangle(9.0, 4.5, 2a, a)

m.groups[:g1] = FaceGroup([1, 2, 3, 13, 14, 7, 15])
m.groups[:g2] = FaceGroup([2, 13, 7, 5])

fg = collectgroups(m, d=2)
@test fg[1] == [:g1]
@test fg[2] == [:g2, :g1]

ng = collectgroups(m, d=0, predefined=true)
@test ng[1] == [:boundarynodes]
