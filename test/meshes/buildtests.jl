using Test
using IntervalSets

using MMJMesh
using MMJMesh.Meshes
using MMJMesh.Mathematics


# -------------------------------------------------------------------------------------------------
# Example 1
# -------------------------------------------------------------------------------------------------

nn = 20
m = Mesh(1, 2)

i1 = addnode!(m, [0, 7])
i2 = addnode!(m, [5, 11])
is1 = addnodes!(m, [1, 2], [12, 12])
is2 = addnodes!(m, 4.0 * MappingFromComponents(Cos(0 .. π), Sin(0 .. π)), nn)
is3 = addnodes!(m, 4.5 * MappingFromComponents(Cos(0 .. π), Sin(0 .. π)), nn)
is4 = addnodes!(m, [3, -1], [-3, -1], nn)

es1 = addelement!(m, 1, 2, group=:g1)
es2 = addelement!(m, 1, 3, group=:g1)
es3 = addelements!(m, is2, is3, group=:conn)
es4 = addelements!(m, is2[1:end-1], is3 .+ 1, group=:diag)
es5 = addelements!(m, is3[1:end-1], is2 .+ 1, group=:diag)
es6 = addelements!(m, is2[1:end-1], is2 .+ 1, group=:chord)
es7 = addelements!(m, is3[1:end-1], is3 .+ 1, group=:chord)
es8 = addelements!(m, is2, is4, group=:g1)

@test nnodes(m) == 64
@test i1 == 1
@test i2 == 2
@test id(node(m, i1)) == 1
@test id(node(m, i2)) == 2

@test es7 == 80:98
@test group(edge(m, 64)) == :chord


# -------------------------------------------------------------------------------------------------
# Example 2
# -------------------------------------------------------------------------------------------------

nx = 8
ny = 5
m = Mesh(2, 2)

nds = addnodes!(m, [[x, y] for y = range(0 .. 2(ny - 1), ny) for x = range(0 .. 2(nx - 1), nx)])

for i = 1:ny-1
    s = (i - 1) * nx + 1
    addelements!(m, nds[s:s+nx-2], nds[s+1:end], nds[s+nx+1:end], nds[s+nx:end])
end

@test nnodes(m) == 40
@test nfaces(m) == 28