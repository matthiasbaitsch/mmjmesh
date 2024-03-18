using Test

using MMJMesh
using MMJMesh.Geometries


dt = 1
pts = [0 1; 1 0]

o = GeometricObjectI{dt}(pts)
p = parametrization(o)
@test p(-1) == pts[:, 1]
@test p(1) == pts[:, 2]


