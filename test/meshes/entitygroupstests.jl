using Test

using MMJMesh
using MMJMesh.Meshes
using MMJMesh.Utilities


# -------------------------------------------------------------------------------------------------
# EntityGroupCollection
# -------------------------------------------------------------------------------------------------

# Test predefined only
a = 5
m = makemeshonrectangle(9.0, 4.5, 2a, a)
@test groupnames(m.groups) == []
@test groupnames(m.groups, d=1) == []
@test groupnames(m.groups, d=1, predefined=true) == [:boundaryedges]
@test groupnames(m.groups, predefined=true) |> sort == 
                                                    [:boundaryedges, :boundaryfaces, :boundarynodes]
@test !hasgroups(m.groups, 0)
@test !hasgroups(m.groups, 1)
@test !hasgroups(m.groups, 2)
@test hasgroups(m.groups, 0, predefined=true)
@test hasgroups(m.groups, 1, predefined=true)

# Test with more groups
m.groups[:g1] = NodeGroup([1, 2, 3, 13, 14, 7, 15])
m.groups[:g2] = FaceGroup([2, 13, 7, 5])
m.groups[:g3] = EdgeGroup([4, 5, 9])
@test groupnames(m.groups) |> sort == [:g1, :g2, :g3]
@test groupnames(m.groups, predefined=true) |> sort == 
                                    [:boundaryedges, :boundaryfaces, :boundarynodes, :g1, :g2, :g3]
@test groupnames(m.groups, d=1, predefined=true) |> sort == [:boundaryedges, :g3]

# Collect groups
a = 5
m = makemeshonrectangle(9.0, 4.5, 2a, a)
m.groups[:g1] = FaceGroup([1, 2, 3, 13, 14, 7, 15])
m.groups[:g2] = FaceGroup([2, 13, 7, 5])
fg = collectgroups(m, d=2)
@test fg[1] == [:g1]
@test fg[2] == [:g2, :g1]
ng = collectgroups(m, d=0, predefined=true)
@test ng[1] == [:boundarynodes]


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

@test isnothing(m.groups.entries[:boundarynodes])
@test isnothing(m.groups.entries[:boundaryedges])
@test length(m.groups[:boundarynodes]) == 30
@test length(m.groups[:boundaryedges]) == 30

