using Test

using MMJMesh
using MMJMesh.Geometries


p1 = Point(0, 0)
p2 = Point(2, 3)
b = Box(p1, p2)

@test p1 ∈ b
@test p2 ∈ b

