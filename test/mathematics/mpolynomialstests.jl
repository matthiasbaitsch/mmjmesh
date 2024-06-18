using Test
using Symbolics

using MMJMesh
using MMJMesh.Mathematics
using MMJMesh.Mathematics: _simplify, _isintegervalue, _integerize, _integerize!

include("validatemappings.jl")


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