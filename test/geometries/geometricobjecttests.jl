using Test

using MMJMesh
using MMJMesh.Meshes
using MMJMesh.Geometries

coords = [0.0 1.0 2.0 0.1 0.9 1.9; 0.0 0.1 0.0 0.9 1.0 0.9]
elts = [[1, 2, 5, 4], [2, 3, 6], [2, 6, 5]]
m = Mesh(coords, elts, 2)

e2 = edge(m, 2)
F = parametrization(geometry(e2))
@test F(-1) == coordinates(e2, 1)
@test F(1) == coordinates(e2, 2)

f1 = face(m, 1)
F = parametrization(geometry(f1))
@test F(-1, -1) == coordinates(f1, 1)
@test F(1, -1) == coordinates(f1, 2)
@test F(1, 1) == coordinates(f1, 3)
@test F(-1, 1) == coordinates(f1, 4)

# TODO test triangles

