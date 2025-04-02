using Test
using Symbolics
using StaticArrays
using LinearAlgebra

using MMJMesh
using MMJMesh.MMJBase
using MMJMesh.Mathematics
using MMJMesh.Mathematics: MPolynomial2,
      _monomialsat, _monomialsderivativeat, _monomialsderivative,
      _lt,
      _subscript, _superscript, _prettymonomial, _factorialpower

# Uncomment in order to work with this file
# include("Validate.jl")
using .Validate

# -------------------------------------------------------------------------------------------------
# Helper functions tests
# -------------------------------------------------------------------------------------------------

# Factorial power
@test _factorialpower(3, 1) == 3
@test _factorialpower(3, 2) == 6
@test _factorialpower(3, 3) == 6
@test _factorialpower(3, 3) == 6
@test _factorialpower(3, 4) == 0

# Monomial evaluation with [x₁x₂⁵, x₁³x₂⁴] at (3, 2)
@test _monomialsat(SA[1 3; 5 4], SA[3, 2]) == [3^1 * 2^5, 3^3 * 2^4]

# Monomial derivative evaluation with [x₁x₂⁵, x₁³x₂⁴] at (3, 2)
@test _monomialsderivativeat(SA[1 3; 5 4], SA[3, 2], [0, 0]) == [3 * 2^5, 3^3 * 2^4]
@test _monomialsderivativeat(SA[1 3; 5 4], SA[3, 2], -[0, 0]) == [3 * 2^5, 3^3 * 2^4]
@test _monomialsderivativeat(SA[1 3; 5 4], SA[3, 2], [1, 1]) == [5 * 2^4, 3 * 3^2 * 4 * 2^3]
@test _monomialsderivativeat(SA[1 3; 5 4], SA[3, 2], [2, 1]) == [0, 6 * 3 * 4 * 2^3]
@test _monomialsderivativeat(SA[1 3; 5 4], SA[3, 2], [1, 2]) == [20 * 2^3, 3 * 3^2 * 12 * 2^2]
@test _monomialsderivativeat(SA[1 3; 5 4], SA[3, 2], -[1, 1]) ==
      [1 / 2 * 3^2 * 1 / 6 * 2^6, 1 / 4 * 3^4 * 1 / 5 * 2^5]
@test _monomialsderivativeat(SA[1 3; 5 4], SA[3, 2], -[2, 1]) ==
      [1 / 6 * 3^3 * 1 / 6 * 2^6, 1 / 20 * 3^5 * 1 / 5 * 2^5]
@test _monomialsderivativeat(SA[1 3; 5 4], SA[3, 2], -[1, 2]) ==
      [1 / 2 * 3^2 * 1 / 42 * 2^7, 1 / 4 * 3^4 * 1 / 30 * 2^6]

# Monomial derivatives with [x₁x₂⁵, x₁³x₂⁴]
@test _monomialsderivative(SA[1 3; 5 4], [0, 0]) == ([1.0, 1.0], [1 3; 5 4])
@test _monomialsderivative(SA[1 3; 5 4], [1, 0]) == ([1, 3], [0 2; 5 4])
@test _monomialsderivative(SA[1 3; 5 4], [0, 1]) == ([5, 4], [1 3; 4 3])
@test _monomialsderivative(SA[1 3; 5 4], [1, 2]) == ([1 * 5 * 4, 3 * 4 * 3], [0 2; 3 2])
@test _monomialsderivative(SA[1 3; 5 4], -[0, 0]) == ([1.0, 1.0], [1 3; 5 4])
@test _monomialsderivative(SA[1 3; 5 4], -[1, 0]) == ([1 / 2, 1 / 4], [2 4; 5 4])
@test _monomialsderivative(SA[1 3; 5 4], -[0, 1]) == ([1 / 6, 1 / 5], [1 3; 6 5])
@test _monomialsderivative(SA[1 3; 5 4], [-1, 2]) == ([1 / 2 * 5 * 4, 1 / 4 * 4 * 3], [2 4; 3 2])
@test _monomialsderivative(SA[1 3; 5 4], -[1, 2]) ==
      ([1 / 2 * 1 / 6 * 1 / 7, 1 / 4 * 1 / 5 * 1 / 6], [2 4; 7 6])

# Exponent comparison
@test _lt([1, 2], [2, 2])
@test !_lt([2, 2], [2, 2])
@test !_lt([2, 2], [1, 2])
@test _lt([1, 4], [4, 1])
@test !_lt([1, 4], [1, 4])
@test !_lt([4, 1], [1, 4])

# Super- and subscripts
@test _subscript(3) == "₃"
@test _subscript(43) == "₄₃"
@test _superscript(3) == "³"
@test _superscript(43) == "⁴³"

# Pretty print monomial
@test _prettymonomial([0, 2, 1]) == "x₁⁰x₂²x₃¹"


# -------------------------------------------------------------------------------------------------
# Basic tests
# -------------------------------------------------------------------------------------------------

# Polynomials
# 3x₁³x₂⁴ + 2x₁¹x₂⁵
u = MPolynomial2([1 3; 5 4], [2, 3])
# x₁³x₂⁷⋅[3, 4] + x₁²x₂⁶⋅[2, 5] + x₁¹x₂⁵⋅[1, 6]
v = MPolynomial2([1 2 3; 5 6 7], [1 2 3; 6 5 4])
# x₁³x₂⁷ ⋅ [9 8 7; 5 8 1; 8 2 4; 6 5 4] + x₁¹x₂⁵ ⋅ [1 2 3; 6 5 4; 4 1 9; 4 3 2]
w = MPolynomial2([1 3; 5 7], [1 2 3; 6 5 4; 4 1 9; 4 3 2;;; 9 8 7; 5 8 1; 8 2 4; 6 5 4])

# Test parameters
x = [3, 2]
@variables a, b
ns1 = SA[1 0; 0 1; 3 3]
ns2 = SArray{Tuple{2,3,2},Int}([2 1 3; 1 0 2;;; 0 1 2; 1 2 3])

# Comparison
@test u === u
@test u === MPolynomial2([1 3; 5 4], [2, 3])

# Simplify in constructor: Duplicates
@test u === MPolynomial2([1 3 3 1; 5 4 4 5], [1, 2, 1, 1])
@test v === MPolynomial2([1 2 3 3 2 1; 5 6 7 7 6 5], [0 1 2 1 1 1; 5 4 3 1 1 1])
@test w === MPolynomial2(
      [1 3 1; 5 7 5],
      [0 1 2; 5 4 3; 3 0 8; 3 2 1;;; 9 8 7; 5 8 1; 8 2 4; 6 5 4;;; 1 1 1; 1 1 1; 1 1 1; 1 1 1]
)

# Simplify in constructor: Zeros
@test nterms(MPolynomial2([1 3 3 1; 5 4 4 5], [1, 1, 1, -1])) == 1
@test nterms(MPolynomial2([1 2 3 3 2 1; 5 6 7 7 6 5], [0 1 2 -2 1 1; 5 4 3 -3 1 1])) == 2
@test nterms(MPolynomial2(
      [1 3 4; 5 7 6],
      [0 1 2; 5 4 3; 3 0 8; 3 2 1;;; 9 8 7; 5 8 1; 8 2 4; 6 5 4;;; 0 0 0; 0 0 0; 0 0 0; 0 0 0]
)) == 2

# Properties
@test exponents(u) == [3 1; 4 5]
@test coefficients(u) == [3, 2]
@test exponents(v) == [3 2 1; 7 6 5]
@test coefficients(v) == [3 2 1; 4 5 6]
@test nterms(u) == 2
@test nterms(v) == 3

# Components
@test [v[1](x), v[2](x)] == v(x)
@test MappingFromComponents(components(v)...)(x) == v(x)
@test domain(MPolynomial2([1 3; 5 4], [2 3; 3 2], QHat)[1]) == QHat

# Values
@test valueat(u, x) == 2 * 3 * 2^5 + 3 * 3^3 * 2^4
@test valueat(v, x) == 3 * 2^5 * [1, 6] + 3^2 * 2^6 * [2, 5] + 3^3 * 2^7 * [3, 4]
@test valueat(w, x) ==
      3^3 * 2^7 * [9 8 7; 5 8 1; 8 2 4; 6 5 4] + 3 * 2^5 * [1 2 3; 6 5 4; 4 1 9; 4 3 2]

# Mathematica code for reference values
# u[x1_, x2_] := 3 x1^3 x2^4 + 2 x1 x2^5
# Derivative[0, 1][u][3, 2];
# Derivative[1, 0][u][3, 2];
# Derivative[2, 3][u][3, 2];
# Derivative[4, 1][u][3, 2];
# Derivative[0, -1][u][3, 2] // N;
# Derivative[-1, 0][u][3, 2] // N;
# N[Derivative[-2, -3][u][3, 2], 16];
# N[Derivative[-4, -1][u][3, 2], 16];
# v[x1_, x2_] := x1^3 x2^7 {3, 4} + x1^2 x2^6 {2, 5} + x1 x2^5 {1, 6}
# Derivative[0, 1][v][3, 2];
# Derivative[1, 0][v][3, 2];
# Derivative[2, 3][v][3, 2];
# Derivative[4, 1][v][3, 2];
# N[Derivative[0, -1][v][3, 2], 16];
# Derivative[-1, 0][v][3, 2];
# N[Derivative[-2, -3][v][3, 2], 16];
# Derivative[-4, -1][v][3, 2] // N

# Evaluation of partial derivatives
@test derivativeat(u, x, [0, 1]) == 3072
@test derivativeat(u, x, [1, 0]) == 1360
@test derivativeat(u, x, [2, 3]) == 2592
@test derivativeat(u, x, [4, 1]) == 0
@test derivativeat(v, x, [0, 1]) == [39984, 58464]
@test derivativeat(v, x, [1, 0]) == [11168, 15936]
@test derivativeat(v, x, [2, 3]) == [185280, 251520]
@test derivativeat(v, x, [4, 1]) == [0, 0]

# Evaluation of antiderivatives
@test derivativeat(u, x, -[0, 1]) ≈ 582.4
@test derivativeat(u, x, -[1, 0]) == 1260
@test derivativeat(u, x, -[2, 3]) ≈ 29.0742857
@test derivativeat(u, x, -[4, 1]) == 93.18857142857143
@test derivativeat(v, x, -[0, 1]) == [2953.142857142857, 4470.857142857143]
@test derivativeat(v, x, -[1, 0]) == [9072, 14112]
@test derivativeat(v, x, -[2, 3]) ≈ [68.98285714285714, 123.9771428571429]
@test derivativeat(v, x, -[4, 1]) == [345.6, 648]

# Evaluation of composed derivatives
@test derivativeat(u, x, ns1) ==
      [derivativeat(u, x, [1, 0]), derivativeat(u, x, [0, 1]), derivativeat(u, x, [3, 3])]
@test derivativeat(u, x, ns2) ==
      [derivativeat(u, x, [2, 0]) derivativeat(u, x, [1, 1]) derivativeat(u, x, [3, 2])
      derivativeat(u, x, [1, 1]) derivativeat(u, x, [0, 2]) derivativeat(u, x, [2, 3])]
@test derivativeat(v, x, ns1) ==
      hcat(derivativeat(v, x, [1, 0]), derivativeat(v, x, [0, 1]), derivativeat(v, x, [3, 3]))

# Evaluation of specific derivatives
@test derivativeat(u, x, 1) == [derivativeat(u, x, [1, 0]), derivativeat(u, x, [0, 1])]
@test derivativeat(v, x, 1) ==
      hcat(derivativeat(v, x, [1, 0]), derivativeat(v, x, [0, 1]))

# Partial derivative functions
@test derivative(u, [1, 0])(x) == derivativeat(u, x, [1, 0])
@test derivative(u, [0, 1])(x) == derivativeat(u, x, [0, 1])
@test derivative(u, [1, 3])(x) == derivativeat(u, x, [1, 3])
@test derivative(v, [1, 0])(x) == derivativeat(v, x, [1, 0])
@test derivative(v, [0, 1])(x) == derivativeat(v, x, [0, 1])
@test derivative(v, [1, 3])(x) == derivativeat(v, x, [1, 3])

# Composed derivative functions
@test derivative(u, ns1)(x) == derivativeat(u, x, ns1)
@test derivative(u, ns2)(x) == derivativeat(u, x, ns2)
@test derivative(v, ns1)(x) == derivativeat(v, x, ns1)

# Specific order derivative functions
@test derivative(u, 1)(x) == [derivativeat(u, x, [1, 0]), derivativeat(u, x, [0, 1])]
@test derivative(u, 1)(x) == derivativeat(u, x, 1)
@test derivative(v, 1)(x) == derivativeat(v, x, 1)

# Operations
u2 = MPolynomial2([3 1; 2 3], [7, 1])
v2 = MPolynomial2([3 1 2; 4 3 2], [4 3 1; 1 4 1])
@test nterms(u + u) == 2
@test nterms(v + v) == 3
@test (u + u)(x) == 2u(x)
@test (u - u)(x) == 0
@test (v + v)(x) == 2v(x)
@test (v - v)(x) == [0, 0]
@test (u + u2)(x) == u(x) + u2(x)
@test (v + v2)(x) == v(x) + v2(x)
@test (3u)(x) == 3u(x)
@test (u * 3)(x) == 3u(x)
@test (3v)(x) == 3v(x)
@test (v * 3)(x) == 3v(x)
@test (u * u2)(x) == u(x) * u2(x)

# Dot product and multiplication by matrix
@test ([4, 1] ⋅ v)(x) == [4, 1] ⋅ v(x)
@test ([4 1; 5 2; 9 1] * v)(x) == [4 1; 5 2; 9 1] * v(x)

# Monomials
@test domain(mmonomials2(2, 2)) == R^2
@test domain(mmonomials2(2, 2, QHat)) == QHat
@test domain(mmonomials2(2, 2, QHat)[1]) == QHat


# -------------------------------------------------------------------------------------------------
# Old tests: Real coefficients
# -------------------------------------------------------------------------------------------------

x = [1, 2, 3]
f = MPolynomial2([3 1; 1 1; 0 2], [-2.0, 3.0], (1 .. 3) × (0 .. 2π) × (2 .. 8))
g = MPolynomial2([3 1 1; 1 1 2; 2 2 3], [1.0, 2.0, 4.0], (1.1 .. 5.2) × (0.1 .. 0.4) × (1.8 .. 1.9))

@test validate(f)
@test validate(g)

@test f == f
@test f === f
@test f(x) == 50
@test (3.1 * f)(x) == 3.1 * 50
@test (f + g)(x) == f(x) + g(x)
@test domaintype(f) == InR3
@test domain(MPolynomial2([1 2; 2 1], [1, 2])) == R2

# 2x₁²x₂⁵+3x₁³x₂²
f = MPolynomial2([2 3; 5 2], [2, 3])
x = SVector(1.0, 2.0)
@test derivative(f, [1, 0]) == MPolynomial2([1 2; 5 2], [4, 9])
@test derivative(f, [2, 0]) == MPolynomial2([0 1; 5 2], [4, 18])
@test derivative(f, [1, 1]) == MPolynomial2([1 2; 4 1], [20, 18])
@test derivativeat(f, x, [1, 1]) == derivative(f, [1, 1])(x)
@test derivativeat(f, [1, 2], [1, 1]) == derivative(f, [1, 1])(x)

gradf = derivative(f, 1)
@test gradf[1] == derivative(f, [1, 0])
@test gradf[2] == derivative(f, [0, 1])

hf = derivative(f, 2)
@test hf[1, 1] == derivative(f, [2, 0])
@test hf[1, 2] == derivative(f, [1, 1])
@test hf[2, 1] == derivative(f, [1, 1])
@test hf[2, 2] == derivative(f, [0, 2])

# Operations
x = SVector(1.0, 2.0)
f = MPolynomial2([2 2; 3 4], [4, 1])
g = MPolynomial2([1 4 1; 6 3 4], [6, 5, 4])

h = f * g
@test h(x) == f(x) * g(x)
@test typeof(h) <: PolynomialRnToR

h = f + g
@test h(x) == f(x) + g(x)
@test typeof(h) <: PolynomialRnToR

# Check simplify
f = MPolynomial2([1 2; 1 2], [4, 1])
g = MPolynomial2([1 1; 1 1], [6, 5])

@test f * g == MPolynomial2([3 2; 3 2], [11, 44])
@test f + g == MPolynomial2([2 1; 2 1], [1, 15])

# # Antiderivative
ns = [1, 2, 3]
f = MPolynomial2([1 2 3; 6 5 4; 1 2 3], [5, 4, 3])
@test f ≈ derivative(antiderivative(f, ns), ns)


f = MPolynomial2([3 2; 1 9], [1, 2])
@test degree(f) == 11
@test degree(f, 1) == 3
@test degree(f, 2) == 9
@test degrees(f) == [degree(f, 1), degree(f, 2)]


# -------------------------------------------------------------------------------------------------
# Old tests: Multivariate monomials
# -------------------------------------------------------------------------------------------------

ps = mmonomials2(2, 1, type=Int)
@test ps[4] == MPolynomial2([0; 0;;], [1])
@test ps[3] == MPolynomial2([0; 1;;], [1])
@test ps[2] == MPolynomial2([1; 0;;], [1])
@test ps[1] == MPolynomial2([1; 1;;], [1])

@test coefficients(ps[1]) isa AbstractVector{Int}
@test ps[1](0, 0) isa Integer

ps = mmonomials2(2, 1, type=BigInt)
@test coefficients(ps[1]) isa AbstractVector{BigInt}


# -------------------------------------------------------------------------------------------------
# Integration
# -------------------------------------------------------------------------------------------------

f = MPolynomial2([1 3; 1 2], [1, 3])
@test integrate(f, 1 .. 2, 5 .. 9) == 2307


# -------------------------------------------------------------------------------------------------
# Symbolic coefficients
# -------------------------------------------------------------------------------------------------

# Multiplication with symbol
@test isequal(coefficients(a * MPolynomial2([1 3; 1 2], [1, 3])), [3a, a])
@test isequal(coefficients(a * MPolynomial2([1 3 4; 1 2 5], [1 0 0; 2 3 0])), [0 a; 3a 2a])

@test isequal(
      coefficients(
            a * MPolynomial2([1 3 4; 1 2 5], [1 0 0; 2 3 0;;; 3 2 1; 5 4 9;;; 0 0 0; 0 0 0])
      ),
      [3a 2a 1a; 5a 4a 9a;;; a 0 0; 2a 3a 0]
)

# Symbolic coefficients and parameters
f = MPolynomial2([1 2; 2 1], [a, b])
g = MPolynomial2([3 2; 2 3], [7, 2])

@test isequal(f(2, 3), 18a + 12b)
@test isequal((a * f)(2, 3), 18a^2 + 12a * b)
@test isequal((f * g)(2, 3), (g * f)(2, 3))
@test isequal((f * f)(1, 2), 16 * a^2 + 16 * a * b + 4 * b^2)
@test isequal(g(a, b), 7(a^3) * (b^2) + 2(a^2) * (b^3))
@test isequal(ValueAtLF([a, b])(f), 2(a^2) * (b^2))

# Integerization
f = MPolynomial2([1 2; 2 1], [2.0a, b + 6 // 3]) |> integerize
@test string(coefficients(f)) == "Num[2 + b, 2a]"

# Rationalization
f = MPolynomial2([1; 0;;], [1])
F = antiderivative(f, [1, 1]) |> rationalize
@test coefficient(F, 1) == 0.5
@test string(coefficient(F, 1)) == "1//2"

