using Test

using MMJMesh
using MMJMesh.Geometries

# Basic functionality
p1 = Point(0, 0)
p2 = Point(2, 3)
b = Box(p1, p2)
@test p1 ∈ b
@test p2 ∈ b

# Parametrization
phi = parametrization(Box([1, 2], [5, 7]))
@test phi(-1, -1) == [1, 2]
@test phi(1, -1) == [5, 2]
@test phi(1, 1) == [5, 7]
@test phi(-1, 1) == [1, 7]