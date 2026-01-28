using Test

using MMJMesh
using MMJMesh.Meshes


# Set up
coords = [0.0 1.0 2.0 0.1 0.9 1.9; 0.0 0.1 0.0 0.9 1.0 0.9]
elts = [[1, 2, 5, 4], [2, 3, 6], [2, 6, 5]]
m = Mesh(coords, elts, 2)
definegroup!(:n1, nodes(m, [2, 4, 5]))
definegroup!(:g1, edges(m, [2, 3]))

# Number of entities
@test nnodes(m) == 6
@test nedges(m) == 8
@test nfaces(m) == 3

# Iterate over entities
@test sum(index(e) for e = nodes(m)) == sum(1:6)
@test sum(index(e) for e = edges(m)) == sum(1:8)
@test sum(index(e) for e = faces(m)) == sum(1:3)
@test sum(index(e) for e = entities(m)) == sum(1:6) + sum(1:8) + sum(1:3)

# Groups
@test edgeindices(m, :n1) == group(m, :g1)
@test edgeindices(m, :n1, select=all) == group(m, :g1)
@test nodeindices(m, :g1) == group(m, :n1)
@test nodeindices(m, :g1, select=any) == group(m, :n1)

# Coordinates
@test coordinates(m) == coords
@test coordinates(m) == coords
@test coordinates(m, 2) == coords[:, 2]
@test coordinates(m, 2, 1) == coords[1, 2]
@test coordinates(m, :g1) == coords[:, group(m, :n1)]
