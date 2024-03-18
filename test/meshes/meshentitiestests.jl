using Test

using MMJMesh.Meshes
using MMJMesh.Geometries


# Set up
coords = [0.0 1.0 2.0 0.1 0.9 1.9; 0.0 0.1 0.0 0.9 1.0 0.9]
elts = [[1, 2, 5, 4], [2, 3, 6], [2, 6, 5]]
m = Mesh(coords, elts, 2)


# -------------------------------------------------------------------------------------------------
# Node
# -------------------------------------------------------------------------------------------------

n5 = node(m, 5)
@test nedges(n5) == 3
@test coordinates(n5) == [0.9, 1.0]
@test coordinate(n5, 1) == 0.9
@test coordinate(n5, 2) == 1.0


# -------------------------------------------------------------------------------------------------
# Edge
# -------------------------------------------------------------------------------------------------

e2 = edge(m, 2)
@test length(e2) == hypot(0.1, 0.9)

g = geometry(e2)
p = parametrization(g)
@test p(-1) == coordinates(e2, 1)
@test p(1) == coordinates(e2, 2)

# -------------------------------------------------------------------------------------------------
# Face
# -------------------------------------------------------------------------------------------------

f1 = face(m, 1)
@test nnodes(f1) == 4
@test nodeindices(f1) == [1, 2, 5, 4]
@test edgeindices(f1) == [1, 2, 3, 4]
@test faceindices(f1) == [3]
@test coordinates(f1) == [0.0 1.0 0.9 0.1; 0.0 0.1 1.0 0.9]
@test coordinates(f1, 1) == [0.0; 0.0]
@test coordinates(f1, 2) == [1.0; 0.1]
@test coordinates(f1, 4) == [0.1; 0.9]

f3 = face(m, 3)
@test nnodes(f3) == 3
@test nodeindices(f3) == [2, 6, 5]
@test coordinates(f3) == [1.0  1.9  0.9; 0.1 0.9 1.0]
@test coordinates(f3, 1) == [1.0, 0.1]
@test coordinates(f3, 2) == [1.9, 0.9]
@test coordinates(f3, 3) == [0.9, 1.0]
