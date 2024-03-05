using Test
using LinearAlgebra

using MMJMesh
using MMJMesh.Mathematics
using MMJMesh.Plots
import MMJMesh.Plots: X1, X2, W1, W2, sampleadaptive


# Test numerical integration formulae
function nint(f, a, b)
    w = b - a
    y1 = f.(a .+ w * X1)
    y2 = f.(a .+ w * X2)
    w * dot(W1, y1), w * dot(W2, y2)
end
check(a, b) = a[1] ≈ b[1] && a[2] ≈ b[2]

@test check(nint(x -> x, 0, 2), (2.0, 2.0))
@test check(nint(x -> x^2, 0, 2), (8.0 / 3.0, 8.0 / 3.0))


# Test helper functions
# @test isnan(safeeval(sin, 1.0 / 0.0))
# @test valuerange(sin, 0, 1, 21) - 2 < 1e-4
# @test abs(valuerange(x -> sin(1 / x), 0, 1, 21) - 2) < 0.5


# Adaptive sampling
x, y = sampleadaptive(Cos(), 0, 2π, ir=true)
# @test length(x) == 1
