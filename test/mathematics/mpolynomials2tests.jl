using Test
using Symbolics
using StaticArrays

using MMJMesh
using MMJMesh.Mathematics
using MMJMesh.Mathematics: MPolynomial2, _monomialsat, _monomialsderivativeat, _lt,
      _subscript, _superscript, _prettymonomial, _factorialpower
using MMJMesh.Mathematics: _simplify, _isintegervalue, _integerize, _integerize!


# -------------------------------------------------------------------------------------------------
# Helper functions tests
# -------------------------------------------------------------------------------------------------

# Factorial power
@test _factorialpower(3, 1) == 3
@test _factorialpower(3, 2) == 6
@test _factorialpower(3, 3) == 6
@test _factorialpower(3, 3) == 6
@test _factorialpower(3, 4) == 0

# Monomial evaluation with [x₁x₃⁵, x₁³x₂⁴] at (3, 2)
@test _monomialsat(SA[1 3; 5 4], SA[3, 2]) == [3^1 * 2^5, 3^3 * 2^4]

# Monomial derivative evaluation with [x₁x₃⁵, x₁³x₂⁴] at (3, 2)
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
u = MPolynomial2([1 3; 5 4], [2, 3])              # 3x₁³x₂⁴ + 2x₁¹x₂⁵
v = MPolynomial2([1 2 3; 5 6 7], [1 2 3; 6 5 4])  # x₁³x₂⁷⋅[3, 4] + x₁²x₂⁶⋅[2, 5] + x₁¹x₂⁵⋅[1, 6]

# Properties
@test exponents(u) == [3 1; 4 5]
@test coefficients(u) == [3, 2]
@test exponents(v) == [3 2 1; 7 6 5]
@test coefficients(v) == [3 2 1; 4 5 6]
@test nterms(u) == 2
@test nterms(v) == 3

# Test point
x = [3, 2]

# Values
@test valueat(u, x) == 2 * 3 * 2^5 + 3 * 3^3 * 2^4
@test valueat(v, x) == 3 * 2^5 * [1, 6] + 3^2 * 2^6 * [2, 5] + 3^3 * 2^7 * [3, 4]

# Derivatives
@test derivativeat(u, x, [0, 1]) == 3072
@test derivativeat(u, x, [1, 0]) == 1360
@test derivativeat(u, x, [2, 3]) == 2592
@test derivativeat(u, x, [4, 1]) == 0
@test derivativeat(v, x, [0, 1]) == [39984, 58464]
@test derivativeat(v, x, [1, 0]) == [11168, 15936]
@test derivativeat(v, x, [2, 3]) == [185280, 251520]
@test derivativeat(v, x, [4, 1]) == [0, 0]

# Antiderivatives
@test derivativeat(u, x, -[0, 1]) ≈ 582.4
@test derivativeat(u, x, -[1, 0]) == 1260
@test derivativeat(u, x, -[2, 3]) ≈ 29.0742857
@test derivativeat(u, x, -[4, 1]) == 93.18857142857143
@test derivativeat(v, x, -[0, 1]) == [2953.142857142857, 4470.857142857143]
@test derivativeat(v, x, -[1, 0]) == [9072, 14112]
@test derivativeat(v, x, -[2, 3]) ≈ [68.98285714285714, 123.9771428571429]
@test derivativeat(v, x, -[4, 1]) == [345.6, 648]

# Mathematica code for above reference values
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

# Composed derivatives
@test derivativeat(u, x, SA[1 0; 0 1; 3 3]) ==
      [derivativeat(u, x, [1, 0]), derivativeat(u, x, [0, 1]), derivativeat(u, x, [3, 3])]
@test derivativeat(u, x, SArray{Tuple{2,3,2},Int}([2 1 3; 1 0 2;;; 0 1 2; 1 2 3])) ==
      [derivativeat(u, x, [2, 0]) derivativeat(u, x, [1, 1]) derivativeat(u, x, [3, 2])
      derivativeat(u, x, [1, 1]) derivativeat(u, x, [0, 2]) derivativeat(u, x, [2, 3])]
@test derivativeat(v, x, SA[1 0; 0 1; 3 3]) ==
      hcat(derivativeat(v, x, [1, 0]), derivativeat(v, x, [0, 1]), derivativeat(v, x, [3, 3]))


# using BenchmarkTools
# @btime derivativeat(u, x, SA[1 0; 0 1; 3 3]);
# @btime derivativeat(v, x, SA[1 0; 0 1; 3 3]);
