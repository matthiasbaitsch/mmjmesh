using Test

using MMJMesh.Meshes


# Set up
coords = [0.0 1.0 2.0 0.1 0.9 1.9; 0.0 0.1 0.0 0.9 1.0 0.9]
elts = [[1, 2, 5, 4], [2, 3, 6], [2, 6, 5]]
m = Mesh(coords, elts, 2)


# Test
@test nnodes(m) == 6
@test nedges(m) == 8
@test nfaces(m) == 3

