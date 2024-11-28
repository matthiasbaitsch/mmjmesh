using Test
using LinearAlgebra

using MMJMesh
using MMJMesh.Meshes
using MMJMesh.Geometries

coords = [0.0 1.0 2.0 0.1 0.9 1.9; 0.0 0.1 0.0 0.9 1.0 0.9]
elts = [[1, 2, 5, 4], [2, 3, 6], [2, 6, 5]]
m = Mesh(coords, elts, 2)

# Edge 3
e2 = edge(m, 2)
g2 = geometry(e2)
F = parametrization(g2)

@test measure(g2) ≈ sqrt(0.82)
@test F(-1) == coordinates(e2, 1)
@test F(1) == coordinates(e2, 2)

# Face 1 (quadrangle)
f1 = face(m, 1)
gf1 = geometry(f1)
F = parametrization(gf1)

@test measure(gf1) == 0.81
@test F(-1, -1) == coordinates(f1, 1)
@test F(1, -1) == coordinates(f1, 2)
@test F(1, 1) == coordinates(f1, 3)
@test F(-1, 1) == coordinates(f1, 4)

# Face 2 (triangle) TODO: Shape functions on triangle
f2 = face(m, 2)
gf2 = geometry(f2)
# F = parametrization(gf2)

@test measure(gf2) ≈ 0.445
# @test F(-1, -1) == coordinates(f2, 1)
# @test F(1, -1) == coordinates(f2, 2)
# @test F(1, 1) == coordinates(f2, 3)
