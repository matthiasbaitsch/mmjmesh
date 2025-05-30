using Test

using MMJMesh
using MMJMesh.Mathematics
using MMJMesh.Plots
import MMJMesh.Plots: sample1d, SafeEval


# -------------------------------------------------------------------------------------------------
# Basic functionality with f(x) = -1 + x^2
# -------------------------------------------------------------------------------------------------

f = Polynomial(-1, 0, 1)
xref = [
    0.0 0.11756712243832271 0.23513424487664542 0.36309740085764786 0.4910605568386503 0.6292282154250157 0.7673958740113812 0.8836979370056905 1.0
]'

# Number of points
for d ∈ 0:3, np ∈ 3:6
    p = sample1d(f, 0, 1, npoints=np, maxrecursion=d, maxangle=0)
    @test size(p, 2) == (np - 1) * 2^d + 1
end

# Compare to reference values
x, p = sample1d(f, 0, 1, maxrecursion=1, maxangle=0, rp=true)
@test p[1, :] ≈ xref
@test f.(x) ≈ p[2, :]


# Test scaling factor
p1 = sample1d(f, 0, 1, maxrecursion=20)
p2 = sample1d(41 * f, 0, 1, maxrecursion=20)
p3 = sample1d(41 * f, 0, 1, maxrecursion=20, yscale=1.0 / 41.0)

@test size(p1, 2) < size(p2, 2)
@test p1[1, :] == p3[1, :]

# Insert root
p1 = sample1d(f, -2, 2, maxrecursion=1)
p2 = sample1d(f, -2, 2, maxrecursion=1, ir=true)
x3, p3 = sample1d(f, -2, 2, maxrecursion=1, ir=true, rp=true)

@test size(p2, 2) == size(p1, 2) + 2
@test 0 ∈ p2[2, :]
@test p1[:, 1] == p2[:, 1]
@test p1[:, end] == p2[:, end]
@test p2 == p3


# -------------------------------------------------------------------------------------------------
# Special cases
# -------------------------------------------------------------------------------------------------

f = 1 / Polynomial(0, 1)
sample1d(f, 0, 1)


# -------------------------------------------------------------------------------------------------
# Parametric curve
# -------------------------------------------------------------------------------------------------

u = ParametricCurve(Sin(), Cos())
x, pts = sample1d(u, 0, 2π, rp=true)
@test length(x) == 257
@test size(pts, 2) == 257


# -------------------------------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------------------------------

se = SafeEval(Sin())
@test isnan(se(Inf))

