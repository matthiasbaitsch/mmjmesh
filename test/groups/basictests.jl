using Test

using MMJMesh.Meshes
using MMJMesh.Groups
using MMJMesh.Utilities

import MMJMesh.Groups: dim, groupnames, hasgroups
import MMJMesh.Meshes: collectgroups

# Test predefined only
a = 5
m = makemeshonrectangle(9.0, 4.5, 2a, a)
@test groupnames(m.groups) == []
@test groupnames(m.groups, d=1) == []
@test groupnames(m.groups, d=1, predefined=true) == [:boundaryedges]
@test groupnames(m.groups, predefined=true) |> sort == [:boundaryedges, :boundaryfaces, :boundarynodes]
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
@test groupnames(m.groups, predefined=true) |> sort == [:boundaryedges, :boundaryfaces, :boundarynodes, :g1, :g2, :g3]
@test groupnames(m.groups, d=1, predefined=true) |> sort == [:boundaryedges, :g3]
