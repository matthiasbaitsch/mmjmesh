using Test
using MMJMesh.Meshes

coords = [0.0 1.0 2.0 0.1 0.9 1.9; 0.0 0.1 0.0 0.9 1.0 0.9]
elts = [[1, 2, 5, 4], [2, 3, 6], [2, 6, 5]]
m = Mesh(coords, elts, 2)

@test nnodes(m) == 6
@test nedges(m) == 8
@test nfaces(m) == 3

e2 = edge(m, 2)
f1 = face(m, 1)
@test nodeIdxs(f1) == [1, 2, 5, 4]
@test edgeIdxs(f1) == [1, 2, 3, 4]
@test faceIdxs(f1) == [3]
@test length(e2) == hypot(0.1, 0.9)


elements(m)