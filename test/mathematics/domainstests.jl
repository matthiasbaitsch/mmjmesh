using Test
using Symbolics
using IntervalSets
using DomainSets: ×

using MMJMesh.Mathematics


# -------------------------------------------------------------------------------------------------
# Basic functionality
# -------------------------------------------------------------------------------------------------

@test dimension(R) == 1
@test (-1 .. 1) ∩ (-10:10) == -1:1
@test (-10:10) ∩ (-1 .. 1) == -1:1
@test isfinite(IHat)
@test !isfinite(R)
@test !isfinite(R⁺)
@test 0 ∈ R
@test 0 ∉ R⁺
@test 0 ∈ R⁺₀

@test dimension(R2) == 2
@test dimension(R3) == 3
@test isfinite(QHat)
@test !isfinite(R2)

@test [1, 2] ∈ R2
@test [1, 2] ∉ R3
@test [1, 2, 3] ∈ R3
@test [1, 2, 3] ∉ R2


# -------------------------------------------------------------------------------------------------
# Points on domain
# -------------------------------------------------------------------------------------------------

# Interval
@test points(IHat, :corners) == [-1, 1]
@test points(IHat, :sides, 3) == []
@test points(IHat, :interior, 3) == [-1 / 2, 0, 1 / 2]

# Rectangle
@assert points(QHat, :corners) == [[-1, -1], [1, -1], [1, 1], [-1, 1]]
@assert points(QHat, :sides, 2) ≈ [
    [-1 / 3, -1], [1 / 3, -1],
    [1, -1 / 3], [1, 1 / 3],
    [1 / 3, 1], [-1 / 3, 1],
    [-1, 1 / 3], [-1, -1 / 3]
]
@test points(QHat, :interior, 2) ≈ [
    [-1 / 3, -1 / 3], [1 / 3, -1 / 3],
    [-1 / 3, 1 / 3], [1 / 3, 1 / 3]
]

K = (1 .. 2) × (3 .. 5)
@assert points(K, :corners) == [[1, 3], [2, 3], [2, 5], [1, 5]]
@assert points(K, :sides, 2) ≈ [
    [4 / 3, 3], [5 / 3, 3],
    [2, 11 / 3], [2, 13 / 3],
    [5 / 3, 5], [4 / 3, 5],
    [1, 13 / 3], [1, 11 / 3]
]
@test points(K, :interior, 2) ≈ [
    [4 / 3, 11 / 3], [5 / 3, 11 / 3],
    [4 / 3, 13 / 3], [5 / 3, 13 / 3]
]

@variables a, b
K = (0 .. a) × (0 .. b)
@test isequal(points(K, :corners), [[0, 0], [a, 0], [a, b], [0, b]])
@test isequal(points(K, :sides, 1), [[a / 2, 0], [a, b / 2], [a / 2, b], [0, b / 2]])
@test isequal(
    points(K, :interior, 2),
    [[a / 3, b / 3], [2a / 3, b / 3], [a / 3, 2b / 3], [2a / 3, 2b / 3]]
)
