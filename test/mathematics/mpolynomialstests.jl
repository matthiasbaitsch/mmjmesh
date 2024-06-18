using Test
using Symbolics

using MMJMesh
using MMJMesh.Mathematics


include("validatemappings.jl")


# Integration
f = MPolynomial([1 3; 1 2], [1, 3])
@test integrate(f, 1 .. 2, 5 .. 9) == 2307

# Simplification
f = MPolynomial([1 2; 2 1], [1.0, 0.0])
@test length(f.p.coefficients) == 1

# Symbolic variables
@variables a, b
f = MPolynomial([1 2; 2 1], [a, b]);
g = MPolynomial([3 2; 2 3], [7, 2]);
@test simplify(f(2, 3) == 18a + 12b) == true
@test simplify((a * f)(2, 3) == 18a^2 + 12a * b) == true
@test simplify((f * g)(2, 3) == (g * f)(2, 3)) == true
@test simplify((f * f)(1, 2) == 16 * a^2 + 16 * a * b + 4 * b^2) == true
@test simplify(g(a, b) == 7(a^3) * (b^2) + 2(a^2) * (b^3)) == true
@test simplify(ValueAtLF([a, b])(f) == 2(a^2) * (b^2)) == true

