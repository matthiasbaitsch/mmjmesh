using Test

using MMJMesh
using MMJMesh.Meshes
using MMJMesh.Meshes: _entitygroupnames, _idvector
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

values = [1.3, 2.1, 5.3, 1.3, 9.1, 2.1]
v2 = copy(values)
ids = _idvector(values)
@test v2 == values
@test validate(values, ids)

values = [[], [:a], [:a], [:b], [:a, :b]]
ids = _idvector(values)
@test validate(values, ids)

# _entitygroupnames
a = 5
m = Mesh((0 .. 9.0) × (0 .. 4.5), 2a, a)
definegroup!(m, 2, :g1, [1, 2, 3, 13, 14, 7, 15])
definegroup!(m, 2, :g2, [2, 13, 7, 5])
fg = _entitygroupnames(m, d=2)
@test fg[1] == [:g1]
@test fg[2] == [:g1, :g2]
ng = _entitygroupnames(m, d=0, predefined=true)
@test ng[1] |> sort == [:b1, :b4, :boundarynodes, :nodes]
@test length(_idvector(fg)) == nfaces(m)


# -------------------------------------------------------------------------------------------------
# Entities and groups
# -------------------------------------------------------------------------------------------------

a = 5
m = Mesh((0 .. 9.0) × (0 .. 4.5), 2a, a)
definegroup!(m, 2, :g1, [1, 2, 3, 13, 14, 7, 15])
definegroup!(m, 2, :g2, [2, 13, 7, 5])
definegroup!(m, 1, :g3, [1, 2, 3, 13, 14, 7, 15])

@test mesh(group(m, :g1)) === m
@test mesh(group(m, :g2)) === m
@test mesh(group(m, :g3)) === m
@test mesh(group(m, :elements)) === m

for gn = groupnames(m, predefined=true)
      @test name(group(m, gn)) == gn
end

@test edim(group(m, :g1)) == 2
@test edim(group(m, :g2)) == 2
@test edim(group(m, :g3)) == 1

e13 = element(m, 13)
e42 = element(m, 42)
@test e13 ∈ group(m, :g1)
@test e13 ∉ group(m, :g3)
@test e42 ∉ group(m, :g1)

@test nodeindices(group(m, :g2)) ==
      [2, 3, 5, 6, 7, 8, 13, 14, 15, 16, 17, 18, 19, 25, 26]


# -------------------------------------------------------------------------------------------------
# Boundary groups
# -------------------------------------------------------------------------------------------------

a = 5
m = Mesh((0 .. 9.0) × (0 .. 4.5), 2a, a)
@test isnothing(m.groups.entries[:boundarynodes])
@test isnothing(m.groups.entries[:boundaryedges])
@test length(group(m, :boundarynodes)) == 30
@test length(group(m, :boundaryedges)) == 30


# -------------------------------------------------------------------------------------------------
# GroupCollection
# -------------------------------------------------------------------------------------------------

# Test predefined only
a = 5
m = Mesh((0 .. 9.0) × (0 .. 4.5), 2a, a)
@test groupnames(m) == [:b1, :b2, :b3, :b4]
@test groupnames(m, d=1) == []
@test groupnames(m, d=1, predefined=true) == [:boundaryedges, :edges]
@test groupnames(m, predefined=true) |> sort ==
      [:b1, :b2, :b3, :b4, :boundaryedges, :boundaryfaces, :boundarynodes, :edges, :elements, :faces, :nodes, :solids]
@test hasgroups(m, d=0)
@test !hasgroups(m, d=1)
@test !hasgroups(m, d=2)
@test hasgroups(m, d=0, predefined=true)
@test hasgroups(m, d=1, predefined=true)

# Test with more groups
definegroup!(m, 0, :g1, [1, 2, 3, 13, 14, 7, 15])
definegroup!(m, 2, :g2, [2, 13, 7, 5])
definegroup!(m, 1, :g3, [4, 5, 9])
@test groupnames(m) |> sort == [:b1, :b2, :b3, :b4, :g1, :g2, :g3]
@test groupnames(m, predefined=true) ==
      [
      :b1, :b2, :b3, :b4,
      :boundaryedges, :boundaryfaces, :boundarynodes, :edges, :elements,
      :faces, :g1, :g2, :g3, :nodes, :solids
]
@test groupnames(m, d=1, predefined=true) == [:boundaryedges, :edges, :g3]

# Group by entity
f5 = face(m, 5)
f7 = face(m, 7)
f33 = face(m, 33)
@test groupnames(node(m, 3)) |> sort == [:b1, :boundarynodes, :g1, :nodes]
@test groupnames(f5) |> sort == [:elements, :faces, :g2]
@test groupnames(f7) |> sort == [:elements, :faces, :g2]
@test groupnames(f33) == [:elements, :faces]
@test groupname(f5) == :g2
@test groupname(f7) == :g2
@test groupname(f33) == :elements


# -------------------------------------------------------------------------------------------------
# Dimension groups
# -------------------------------------------------------------------------------------------------

g = group(m, :edges)
for f ∈ edges(m)
      @test f ∈ g
end

g = group(m, :faces)
for f ∈ faces(m)
      @test f ∈ g
end


# -------------------------------------------------------------------------------------------------
# Dimension groups
# -------------------------------------------------------------------------------------------------

@test indices(group(m, :g2), 0) == [2, 3, 5, 6, 7, 8, 13, 14, 15, 16, 17, 18, 19, 25, 26]
@test indices(group(m, :g2), 1) == [2, 5, 6, 7, 10, 12, 14, 15, 16, 18, 20, 21, 22, 35, 37, 38]
