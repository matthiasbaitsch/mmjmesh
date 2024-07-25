using Test
using IntervalSets

import DomainSets
using DomainSets: ×

using MMJMesh
using MMJMesh.Geometries

# Basic functionality
p1 = Point(0, 0)
p2 = Point(2, 3)
b = Box(p1, p2)
@test p1 ∈ b
@test p2 ∈ b

# Parametrization 2D
phi = parametrization(Box([1, 2], [5, 7]))
@test phi(-1, -1) == [1, 2]
@test phi(1, -1) == [5, 2]
@test phi(1, 1) == [5, 7]
@test phi(-1, 1) == [1, 7]

# Parametrization 3D
phi = parametrization(Box([1, 2, 3], [5, 7, 9]))
@test phi(-1, -1, -1) == [1, 2, 3]
@test phi(1, 1, 1) == [5, 7, 9]

# From DomainSets.Rectangle
b = Box((1 .. 2) × (2 .. 7))
@test minimum(b) == Point(1, 2)
@test maximum(b) == Point(2, 7)