@testitem "Helperfunctions" begin

    using Test
    using LinearAlgebra: Diagonal

    using MMJMesh
    using MMJMesh.Mathematics
    using MMJMesh.Plots: _collectlines, _mergemeshes, _getnintervals, sample2d, sample2dlines
    using DomainSets: ×

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

    # Triangle
    x, t = sample2d(x -> x[1] + 2x[2], domain=THat, npoints=2)
    @test x ≈ [
        0.0 0.5 1.0 0.0 0.5 0.0
        0.0 0.0 0.0 0.5 0.5 1.0
        0.0 0.5 1.0 1.0 1.5 2.0
    ]
    @test t == [1 2 4; 2 5 4; 2 3 5; 4 5 6]

    lines = sample2dlines(x -> 0.0, domain=THat, npoints=2, nmeshlines=1)
    @test length(lines) == 6
    @test lines[1][1] ≈ [0, 0, 0] && lines[1][end] ≈ [1, 0, 0]       # y = 0
    @test lines[end][1] ≈ [1, 0, 0] && lines[end][end] ≈ [0, 1, 0]   # x + y = 1


    n1, n2 = _getnintervals((1 .. 2) × (3 .. 5), 10)
    @test n1 == 5
    @test n2 == 10

end
