using Test
using Symbolics

using MMJMesh
using MMJMesh.Mathematics
using MMJMesh.Mathematics: _simplify, _isintegervalue, _integerize, _integerize!

include("validatemappings.jl")


# -------------------------------------------------------------------------------------------------
# Real coefficients
# -------------------------------------------------------------------------------------------------

x = [1, 2, 3]
f = MPolynomial([3 1; 1 1; 0 2], [-2.0, 3.0])
g = MPolynomial([3 1 1; 1 1 2; 2 2 3], [1.0, 2.0, 4.0])

@test f == f
@test f(x) == 50
@test (3.1 * f)(x) == 3.1 * 50
@test (f + g)(x) == f(x) + g(x)
@test domaintype(f) == SVector{3,<:Real}
@test domain(MPolynomial([1 2; 2 1], [1, 2])) == R2

# 2x₁²x₂⁵+3x₁³x₂²
f = MPolynomial([2 3; 5 2], [2, 3])
x = SVector(1.0, 2.0)
@test derivative(f, [1, 0]) == MPolynomial([1 2; 5 2], [4, 9])
@test derivative(f, [2, 0]) == MPolynomial([0 1; 5 2], [4, 18])
@test derivative(f, [1, 1]) == MPolynomial([1 2; 4 1], [20, 18])
@test derivativeat(f, x, [1, 1]) == derivative(f, [1, 1])(x)

gradf = derivative(f, 1)
@test gradf[1] == derivative(f, [1, 0])
@test gradf[2] == derivative(f, [0, 1])

hf = derivative(f, 2)
@test hf[1][1] == derivative(f, [2, 0])
@test hf[1][2] == derivative(f, [1, 1])
@test hf[2][1] == derivative(f, [1, 1])
@test hf[2][2] == derivative(f, [0, 2])

# Operations
x = SVector(1.0, 2.0)
f = MPolynomial([2 2; 3 4], [4, 1])
g = MPolynomial([1 4 1; 6 3 4], [6, 5, 4])

h = f * g
@test h(x) == f(x) * g(x)
@test typeof(h) == typeof(g)

h = f + g
@test h(x) == f(x) + g(x)
@test typeof(h) == typeof(g)

# Check simplify
f = MPolynomial([1 2; 1 2], [4, 1])
g = MPolynomial([1 1; 1 1], [6, 5])

@test f * g == MPolynomial([3 2; 3 2], [11, 44])
@test f + g == MPolynomial([2 1; 2 1], [1, 15])

# Antiderivative
ns = [1, 2, 3]
f = MPolynomial([1 2 3; 6 5 4; 1 2 3], [5, 4, 3])
ff = derivative(antiderivative(f, ns), ns)

@test f.p.exponents == ff.p.exponents
@test f.p.coefficients ≈ ff.p.coefficients


# -------------------------------------------------------------------------------------------------
# Multivariate monomials
# -------------------------------------------------------------------------------------------------

ps = mmonomials(2, 1, type=Int)
@test ps[1] == MPolynomial([0; 0;;], [1])
@test ps[2] == MPolynomial([1; 0;;], [1])
@test ps[3] == MPolynomial([0; 1;;], [1])
@test ps[4] == MPolynomial([1; 1;;], [1])
@test ps[1](0, 0) isa Integer

ps = mmonomials(2, 1, type=BigInt)
@test ps[1].p.coefficients isa Vector{BigInt}


# -------------------------------------------------------------------------------------------------
# Symbolic coefficients
# -------------------------------------------------------------------------------------------------

# Declarations
@variables a, b

# Integration
f = MPolynomial([1 3; 1 2], [1, 3])
@test integrate(f, 1 .. 2, 5 .. 9) == 2307

# Simplification
@test ([0; 0;;], [0]) == _simplify([1 2; 2 1], [0, 0])
@test isequal(([3; 1;;], [a]), _simplify([3 2; 1 1], [a, 0]))
@test isequal(([0; 0;;], [0]), _simplify([3 2; 1 1], [0 * a, 0]))

# Make integers
u = [2.0 + 3.0a, b]

@test _isintegervalue(1)
@test _isintegervalue(2 // 2)
@test _isintegervalue(1.0)
@test !_isintegervalue(1.1)
@test !_isintegervalue(3 // 2)

@test string(_integerize(1.0a + a)) == "2a"
@test string(_integerize(2.0a)) == "2a"
@test string(_integerize(12 // 6 * a)) == "2a"
@test string(_integerize(2.0a + 6 // 3)) == "2 + 2a"
@test string(_integerize(0.0a)) == "0"

@test string(_integerize!(u)) == string([2 + 3a, b])
@test string(u) == string([2 + 3a, b])

# Symbolic coefficients and parameters
f = MPolynomial([1 2; 2 1], [a, b]);
g = MPolynomial([3 2; 2 3], [7, 2]);

@test isequal(f(2, 3), 18a + 12b)
@test isequal((a * f)(2, 3), 18a^2 + 12a * b)
@test isequal((f * g)(2, 3), (g * f)(2, 3))
@test isequal((f * f)(1, 2), 16 * a^2 + 16 * a * b + 4 * b^2)
@test isequal(g(a, b), 7(a^3) * (b^2) + 2(a^2) * (b^3))
@test isequal(ValueAtLF([a, b])(f), 2(a^2) * (b^2))

# Integerization of coefficients
f = MPolynomial([1 2; 2 1], [2.0a, b + 6 // 3]);
@test string(f.p.coefficients) == string([2 + b, 2a])

# Integration with rationals
f = MPolynomial([1; 0;;], [0 * a + 1])
F = antiderivative(f, [1, 1])
@test string(F.p.coefficients[1]) == "1//2"

f = MPolynomial([1; 0;;], [1])
F = antiderivative(f, [1, 1])
@test F.p.coefficients[1] == 0.5

# TODO generic test