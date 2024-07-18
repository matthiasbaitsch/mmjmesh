using Test
using IntervalSets

using MMJMesh.Mathematics

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

