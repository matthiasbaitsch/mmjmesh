using Test

using IntervalSets

using MMJMesh
using MMJMesh.Tests
using MMJMesh.Mathematics
using MMJMesh.Mathematics: _makebreakpoints, indexofpieceat

# Helpers
@test _makebreakpoints(1:2, 0 .. 3) == 0:3
@test _makebreakpoints(0:3, 0 .. 3) == 0:3
@test _makebreakpoints(-1:4, 0 .. 3) == 0:3

# Testees
pw1 = PiecewiseFunction([0, 5, 12, 15], [Sin(), Cos(), Sin()])
pw2 = PiecewiseFunction([1.1, 6.5, 11.1], [Polynomial(3, 2, 1), Polynomial(1, 2, 3)])

# Basic
@test domain(pw1) == 0 .. 15
@test npieces(pw1) == 3
@test pois(pw1) == [5, 12]

# Index of x
@test indexofpieceat(pw1, -1) == 1
@test indexofpieceat(pw1, 0) == 1
@test indexofpieceat(pw1, 0.5) == 1
@test indexofpieceat(pw1, 5) == 1
@test indexofpieceat(pw1, 5.5) == 2
@test indexofpieceat(pw1, 12) == 2
@test indexofpieceat(pw1, 12.5) == 3
@test indexofpieceat(pw1, 15) == 3
@test indexofpieceat(pw1, 15.5) == 3

# Function value and derivatives
@test pw1(1) == sin(1)
@test derivativeat(pw1, 1) ≈ cos(1)
@test derivativeat(pw1, 1, 2) ≈ -sin(1)
@test pw1(6) == cos(6)
@test derivativeat(pw1, 6) ≈ -sin(6)
@test derivativeat(pw1, 6, 2) ≈ -cos(6)

# Integration
@test integrate(pw1, 1 .. 2) == integrate(Sin(), 1 .. 2)
@test integrate(pw1, 1 .. 6) == integrate(Sin(), 1 .. 5) + integrate(Cos(), 5 .. 6)
@test integrate(pw1, 1 .. 14) ==
      integrate(Sin(), 1 .. 5) + integrate(Cos(), 5 .. 12) + integrate(Sin(), 12 .. 14)
@test integrate(pw1, 0 .. 15) ==
      integrate(Sin(), 0 .. 5) + integrate(Cos(), 5 .. 12) + integrate(Sin(), 12 .. 15)

# Operations
@test (pw1 + Sin()) isa PiecewiseFunction
@test (Sin() + pw1) isa PiecewiseFunction
@test (pw1 * Sin()) isa PiecewiseFunction
@test (Sin() * pw1) isa PiecewiseFunction
@test (π * pw1) isa PiecewiseFunction
@test (pw1 + 3Sin()) isa PiecewiseFunction
@test (3Sin() + pw1) isa PiecewiseFunction

@test validateoperations(pw1, Sin())
@test validateoperations(pw1, pw1)
@test validateoperations(pw1, pw2)
@test validateoperations(pw2, pw1)

# Linear interpolation
x = [1, 2, 5]
y = [2, 1, 6]
f = interpolate(x, y)

for i ∈ eachindex(x)
    @test f(x[i]) ≈ y[i]
end