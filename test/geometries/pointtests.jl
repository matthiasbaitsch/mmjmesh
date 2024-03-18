using Test

using MMJMesh
using MMJMesh.Geometries

p1 = Point(1, 2, 3)
p2 = Point(3, 2, 1)
p3 = Point(10, 20, 30)

@test p2 - p1 == [2, 0, -2]
@test p2 ≤ p2
@test p2 ≥ p2
@test p1 < p3
@test p3 > p1
@test !(p1 < p2)
@test !(p1 > p2)

