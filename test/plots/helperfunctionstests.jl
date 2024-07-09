using Test
using DomainSets
using LinearAlgebra

using MMJMesh.Mathematics
using MMJMesh.Plots: _collectlines, _mergemeshes, _getnintervals, sample2d

x1, x2 = _collectlines([[[0, 0], [1, 0], [2, 0]], [[0, 1], [2, 1]]])
@test isequal(x1, [0, 1, 2, NaN, 0, 2, NaN])
@test isequal(x2, [0, 0, 0, NaN, 1, 1, NaN])


m1 = sample2d(x -> 1, domain=QHat, npoints=1)
m2 = sample2d(x -> 2, domain=QHat, npoints=1, gmap=AffineMapping(Diagonal([1, 1]), [2, 0]))
m = _mergemeshes([m1, m2])
@test m[1] == [
    -1.0 1.0 -1.0 1.0 1.0 3.0 1.0 3.0
    -1.0 -1.0 1.0 1.0 -1.0 -1.0 1.0 1.0
    1.0 1.0 1.0 1.0 2.0 2.0 2.0 2.0
]
@test m[2] == [1 2 4; 1 4 3; 5 6 8; 5 8 7]



n1, n2 = _getnintervals((1 .. 2) × (3 .. 5), 10)
@test n1 == 5
@test n2 == 10