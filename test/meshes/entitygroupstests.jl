using Test

using MMJMesh
using MMJMesh.Meshes
using MMJMesh.Meshes: _collectgroups, _idvector
using MMJMesh.Utilities


# -------------------------------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------------------------------

# _idvector
function validate(values, ids)
      @test length(values) == length(ids)
      for i = eachindex(values), j = eachindex(values)
            @test cmp(values[i], values[j]) == cmp(ids[i], ids[j])
      end
      return true
end

# IDs for numbers
values = [1.3, 2.1, 5.3, 1.3, 9.1, 2.1]
v2 = copy(values)
ids = _idvector(values)
@test v2 == values
@test validate(values, ids)

# IDs for vectors of symbols
values = [[], [:a], [:a], [:b], [:a, :b]]
ids = _idvector(values)
@test validate(values, ids)

# Collect groups
a = 5
m = makemeshonrectangle(9.0, 4.5, 2a, a)
m.groups[:g1] = FaceGroup([1, 2, 3, 13, 14, 7, 15])
m.groups[:g2] = FaceGroup([2, 13, 7, 5])
fg = _collectgroups(m, d=2)
@test fg[1] == [:g1]
@test fg[2] == [:g2, :g1]
ng = _collectgroups(m, d=0, predefined=true)
@test ng[1] |> sort == [:boundarynodes, :nodes]

# ID vector
@test length(_idvector(fg)) == nfaces(m)


# -------------------------------------------------------------------------------------------------
# Entities and groups
# -------------------------------------------------------------------------------------------------

a = 5
m = makemeshonrectangle(9.0, 4.5, 2a, a)
m.groups[:g1] = FaceGroup([1, 2, 3, 13, 14, 7, 15])
m.groups[:g2] = FaceGroup([2, 13, 7, 5])
m.groups[:g3] = EdgeGroup([1, 2, 3, 13, 14, 7, 15])

@test edim(m.groups[:g1]) == 2
@test edim(m.groups[:g2]) == 2
@test edim(m.groups[:g3]) == 1

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


# -------------------------------------------------------------------------------------------------
# GroupCollection
# -------------------------------------------------------------------------------------------------

# Test predefined only
a = 5
m = makemeshonrectangle(9.0, 4.5, 2a, a)
@test groupnames(m) == []
@test groupnames(m, d=1) == []
@test groupnames(m, d=1, predefined=true) == [:edges, :boundaryedges]
@test groupnames(m, predefined=true) |> sort ==
      [:boundaryedges, :boundaryfaces, :boundarynodes, :edges, :faces, :nodes, :solids]
@test !hasgroups(m, d=0)
@test !hasgroups(m, d=1)
@test !hasgroups(m, d=2)
@test hasgroups(m, d=0, predefined=true)
@test hasgroups(m, d=1, predefined=true)

# Test with more groups
m.groups[:g1] = NodeGroup([1, 2, 3, 13, 14, 7, 15])
m.groups[:g2] = FaceGroup([2, 13, 7, 5])
m.groups[:g3] = EdgeGroup([4, 5, 9])
@test groupnames(m) |> sort == [:g1, :g2, :g3]
@test groupnames(m, predefined=true) |> sort ==
      [:boundaryedges, :boundaryfaces, :boundarynodes, :edges, :faces, :g1, :g2, :g3, :nodes, :solids]
@test groupnames(m, d=1, predefined=true) |> sort == [:boundaryedges, :edges, :g3]

# Group by entity
f5 = face(m, 5)
f7 = face(m, 7)
f33 = face(m, 33)
@test groups(node(m, 3)) == [:g1, :boundarynodes, :nodes]
@test groups(f5) == [:g2, :faces]
@test groups(f7) == [:g2, :faces]
@test groups(f33) == [:faces]
@test group(f5) == :g2
@test group(f7) == :g2
@test group(f33) == :faces


# -------------------------------------------------------------------------------------------------
# Dimension groups
# -------------------------------------------------------------------------------------------------

g = m.groups[:edges]
for f ∈ edges(m)
      @test f ∈ g
end

g = m.groups[:faces]
for f ∈ faces(m)
      @test f ∈ g
end
