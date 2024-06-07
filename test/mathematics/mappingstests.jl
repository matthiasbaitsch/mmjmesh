using Test
using Symbolics
using StaticArrays
using IntervalSets
using LinearAlgebra

using MMJMesh
using MMJMesh.Mathematics

include("validatemappings.jl")


# -------------------------------------------------------------------------------------------------
# Similar
# -------------------------------------------------------------------------------------------------

rand(5, 5) * monomials(0:4)
Base.similar(monomials(0:4), Polynomial, Base.OneTo(5))
@test length(similar([Sin(), Cos()], Sin)) == 2


# -------------------------------------------------------------------------------------------------
# Matrix vector product
# -------------------------------------------------------------------------------------------------

m1 = Sin()
m2 = Cos()

# Needs union
[0 0; 0 0] * [m1, m1]
[1.1 2.2; 3.1 5.8] * [m1, m1]
[1.1 2.2; 3.1 5.8] * [m1, m2]


# -------------------------------------------------------------------------------------------------
# Interval domains
# -------------------------------------------------------------------------------------------------

@test 3 ∈ R
@test 0 ∉ RPlus
@test -3 ∉ RPlus
@test 0 ∈ R0Plus

@test intersect(R, [1, 2, 3]) == [1, 2, 3]
@test intersect(RPlus, [-1, 1, 2, 3]) == [1, 2, 3]
@test R ∩ [1, 2, 3] == [1, 2, 3]
@test RPlus ∩ [-1, 1, 2, 3] == [1, 2, 3]


# -------------------------------------------------------------------------------------------------
# Mapping from components
# -------------------------------------------------------------------------------------------------

# Parametric curve
g = MappingFromComponents(Sin(), Cos())

@test g[1] == Sin()
@test valueat(g, 1.0) == [sin(1.0), cos(1.0)]
@test derivativeat(g, 1.0, 1) == [cos(1), -sin(1)]
@test derivativeat(g, 1.0, 2) ≈ [-sin(1), -cos(1)]
@test derivative(g, 1)(1.0) ≈ derivativeat(g, 1.0, 1)
@test derivative(g, 2)(1.0) ≈ derivativeat(g, 1.0, 2)
# TODO generic test

# Vector field R2 -> R2
x = SVector(1.0, 2.0)
f1 = MPolynomial([1 2; 2 1], [1, 2])
f2 = ProductFunction(Sin(), Cos())
g = MappingFromComponents(f1, f2)

@test valueat(g, x) == [f1(x), f2(x)]
@test g(x) == [f1(x), f2(x)]
# TODO generic test


# -------------------------------------------------------------------------------------------------
# Ad hoc mapping
# -------------------------------------------------------------------------------------------------

f = makefunction(x -> sin(x[1]) * sin(x[2]), 0 .. 5π, 0 .. 5π)
@test valueat(f, [1.0, 2.0]) == sin(1) * sin(2)


# -------------------------------------------------------------------------------------------------
# Functions R → R
# -------------------------------------------------------------------------------------------------

# Sine and Cosine
s = Sin(0 .. 2π)
c = Cos(0 .. 2π)
validate(s)
validate(antiderivative(s))
validate(c)
validate(antiderivative(c))
@test antiderivative(Sin(), 2) == -Sin()

# Polynomial
p = Polynomial(4, 6, 1, 9, 2, -1)
@test degree(p) == 5
validate(p)
validate(antiderivative(p))

Polynomial([1, 2, 3]) == Polynomial(1, 2, 3)

# Roots and domain
@test roots(Polynomial([-1, 0, 1])) == [-1, 1]
@test roots(Polynomial([-1, 0, 1], RPlus)) == [1]
@test roots(Polynomial([-1, 0, 1], -5 .. 0)) == [-1]
@test roots(Polynomial([-1, 0, 1]), -5 .. 0) == [-1]
@test roots(Polynomial([-1, 0, 1]), 4 .. 5) == []

# Lagrange
p = [0, 1, 3, 4, 5.5]
l4 = lagrangepolynomials(p, 0 .. 5.5)
validate.(l4, atol=1e-4)
for i ∈ 1:5, j ∈ 1:5
    @test isapprox(l4[i](p[j]), i == j ? 1 : 0, atol=1e-14)
end

# From roots
c = [1, 2, 3]
p = fromroots(c)
@test roots(p) ≈ c

# Monomials
@test Polynomial([1, 2, 3, 4, 5]) == [1, 2, 3, 4, 5]' * monomials(0:4)


# -------------------------------------------------------------------------------------------------
# Functions Rn -> R
# -------------------------------------------------------------------------------------------------

# Product function, two parameters
f = ProductFunction(Sin(0 .. 1), Cos(0.5 .. 5))
x = SVector(1.0, 2.0)

@test domain(f) == (0 .. 1) × (0.5 .. 5)
@test valueat(f, x) == sin(1) * cos(2)
@test derivativeat(f, x, [1, 0]) ≈ cos(1.0) * cos(2.0)
@test derivativeat(f, x, [2, 1]) ≈ sin(1.0) * sin(2.0)
@test derivative(f, [2, 1])(x) ≈ derivativeat(f, x, [2, 1])
@test derivativeat(f, x, 1) == [cos(1) * cos(2), -sin(1) * sin(2)]
@test derivative(f, 1)(x) == derivativeat(f, x, 1)
@test derivativeat(f, x, 2) ≈ [
    -sin(1)*cos(2) -cos(1)*sin(2)
    -cos(1)*sin(2) -sin(1)*cos(2)
]
@test derivative(f, 2)(x) ≈ derivativeat(f, x, 2)
@test derivativeat(f, x, 1) isa SVector{2,Float64}
@test derivativeat(f, x, 2) isa SMatrix{2,2,Float64}

@test derivativeat(f, Vector(x)) == derivativeat(f, x)

@test gradientat(f, x) == derivativeat(f, x, 1)
@test gradientat(f, Vector(x)) == derivativeat(f, x, 1)
@test gradientat(f, x...) == derivativeat(f, x, 1)

@test hessianat(f, x) == derivativeat(f, x, 2)
@test hessianat(f, Vector(x)) == derivativeat(f, x, 2)
@test hessianat(f, x...) == derivativeat(f, x, 2)

@test laplacian(f)(x) ≈ -2 * f(x)
@test laplacianat(f, x) ≈ -2 * f(x)
@test laplacianat(f, Vector(x)) == laplacianat(f, x)
@test laplacianat(f, x...) == laplacianat(f, x)

@test ∇(f) == derivative(f, 1)
@test H(f) == derivative(f, 2)
@test Δ(f)(x) ≈ laplacianat(f, x)

# TODO Generic test


# Product function three variables
f = ProductFunction(Sin(1 .. 2), Cos(3 .. 4), Polynomial([2, 3], 6 .. 7))
x = SVector(1.0, 2.0, 3.0)

@test domain(f) == (1 .. 2) × (3 .. 4) × (6 .. 7)
@test valueat(f, x) == sin(1) * cos(2) * 11
@test derivativeat(f, x, 1) ≈ [cos(1) * cos(2) * 11, -sin(1) * sin(2) * 11, sin(1) * cos(2) * 3]
@test derivativeat(f, x, 2) ≈ [
    -sin(1)*cos(2)*11 -cos(1)*sin(2)*11 cos(1)*cos(2)*3
    -cos(1)*sin(2)*11 -sin(1)*cos(2)*11 -sin(1)*sin(2)*3
    cos(1)*cos(2)*3 -sin(1)*sin(2)*3 0
]
@test derivative(f, 2)(x) ≈ derivativeat(f, x, 2)
@test derivativeat(f, x, 1) isa SVector{3,Float64}
@test derivativeat(f, x, 2) isa SMatrix{3,3,Float64}
# TODO Generic test


# Antiderivative
ns = [3, 2]
f = ProductFunction(Sin(), Cos())
@test derivative(antiderivative(f, ns), ns) == f


# Integrals
f = ProductFunction(Sin(), Cos())
@test integrate(f, 1 .. 2, 5 .. 9) ≈ 1.31133267192572

f = MPolynomial([1 3; 1 2], [1, 3])
@test integrate(f, 1 .. 2, 5 .. 9) == 2307


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


# -------------------------------------------------------------------------------------------------
# Parametric curves R → Rn
# -------------------------------------------------------------------------------------------------

c = ParametricCurve(Sin(), Cos())
@test c(0.0) == [0, 1]
@test derivativeat(c, 0) ≈ [1, 0]
@test derivativeat(c, 0, 2) ≈ [0, -1]
@test c'(0) ≈ [1, 0]
@test c''(0) ≈ [0, -1]
validate(c)

n = UnitNormal(c)
@test n(0.0) ≈ [0, 1]
@test norm(n(0.0)) ≈ 1
@test abs(n(0.0) ⋅ c'(0.0)) < 1e-16

f1 = Polynomial(1, 2, 3) * c
f2 = c * Polynomial(1, 2, 3)
@test f1 == f2
@test f1(π) == f2(π)
@test codomaintype(f1) <: SVector
@test codomaintype(f2) <: SVector


# -------------------------------------------------------------------------------------------------
# Operations
# -------------------------------------------------------------------------------------------------

m1 = Sin()
m2 = Sin()
m3 = Cos()
m4 = Polynomial([0, 0, 1], 0 .. 4)

# Multiply by scalar
@test m1 * 2 == 2 * m1
@test 2 * m1 == 2 * m1
@test 2 * m1 == 2 * m2
@test 2 * m1 + 3 * m2 == 5 * m1
@test (-1) * ((-1) * m1) == m1
@test (-1) * (-m1) == m1
validate(2 * m1, rtol=1e-4)

# Add and subtract
@test m1 + zero(m1) == m1
@test m1 + m1 == 2.0 * m1
@test 2 * m1 + m1 == 3 * m1
@test m1 + 2 * m1 == 3 * m1
@test derivative(m1 + m3) == Cos() - Sin()
@test antiderivative(m1 + m3) == -Cos() + Sin()
@test m1 - m1 == zero(m1)
@test m3 - m1 == Cos() - Sin()
validate(m1 + m3, rtol=1e-4)

# Composition m5 = Sin(x^2)
m5 = m1 ∘ m4
@test m5' == (Cos() ∘ Polynomial([0, 0, 1], 0 .. 4)) * Polynomial([0, 2], 0 .. 4)
@test isapprox(m5(√π), 0, atol=1e-15)
validate(m5, rtol=1e-4)

# Product
p = m1 * m2
validate(p)

# Quotient
p = m1 / m3
@test p(0.2) ≈ tan(0.2)
validate(p, rtol=1e-4)
f = 2 / Polynomial(0, 1)
@test f(0.2) ≈ 10
@test pois(f) == [0]


# -------------------------------------------------------------------------------------------------
# Mappings to Rn
# -------------------------------------------------------------------------------------------------

f1 = MPolynomial([1 2 3; 3 2 1], [6, 5, 4])
f2 = ProductFunction(Sin(), Cos())

u1 = MappingFromComponents(f1, f2)
u2 = MappingFromComponents(f2, f1)


u = u1 ⋅ u2
@test u(2, 3) == 2 * f1(2, 3) * f2(2, 3)

u = ∇(f1) ⋅ ∇(f2)
@test u(2, 3) ≈ derivative(f1, [1, 0])(2, 3) * cos(2) * cos(3) - derivative(f1, [0, 1])(2, 3) * sin(2) * sin(3)


# -------------------------------------------------------------------------------------------------
# Mappings Rn → R
# -------------------------------------------------------------------------------------------------

# Multivariate polynomial
x = [1, 2, 3]
f = MPolynomial([3 1; 1 1; 0 2], [-2.0, 3.0])
g = MPolynomial([3 1 1; 1 1 2; 2 2 3], [1.0, 2.0, 4.0])

@test f == f
@test f(x) == 50
@test (3.1 * f)(x) == 3.1 * 50
@test (f + g)(x) == f(x) + g(x)
@test domaintype(f) == SVector{3,<:Real}
@test domain(MPolynomial([1 2; 2 1], [1, 2])) == R2

ps = mmonomials(2, 1)
@test ps[1] == MPolynomial([0; 0;;], [1])
@test ps[2] == MPolynomial([1; 0;;], [1])
@test ps[3] == MPolynomial([0; 1;;], [1])
@test ps[4] == MPolynomial([1; 1;;], [1])

# 2x₁²x₂⁵+3x₁³x₂²
f = MPolynomial([2 3; 5 2], [2, 3])
x = SVector(1.0, 2.0)
@test derivative(f, [1, 0]) == MPolynomial([1 2; 5 2], [4, 9])
@test derivative(f, [2, 0]) == MPolynomial([0 1; 5 2], [4, 18])
@test derivative(f, [1, 1]) == MPolynomial([1 2; 4 1], [20, 18])
@test derivativeat(f, x, [1, 1]) == derivative(f, [1, 1])(x)

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
# Vector fields
# -------------------------------------------------------------------------------------------------

x = SVector(1.0, 2.0)
v = MappingFromComponents(ProductFunction(Sin(), Cos()), ProductFunction(Cos(), Sin()))

@test divergenceat(v, x) ≈ cos(1) * cos(2) + cos(1) * cos(2)
@test divergenceat(v, 1.0, 2.0) == divergenceat(v, x)
@test divergenceat(v, [1.0, 2.0]) == divergenceat(v, x)
@test divergence(v)(x) ≈ divergenceat(v, x)
@test div(v) == divergence(v)
