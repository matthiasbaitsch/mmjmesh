using Test

using MMJMesh
using MMJMesh.Mathematics


u = linesegment([1, 2, 3], [3, 2, 1])
@test domain(u) == IHat
@test u(-1) == [1, 2, 3]
@test u(0) == [2, 2, 2]
@test u(1) == [3, 2, 1]
