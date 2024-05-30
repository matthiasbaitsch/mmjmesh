using Test
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
# Functions R → R
# -------------------------------------------------------------------------------------------------

# Sine and Cosine
s = Sin(0 .. 2π)
c = Cos(0 .. 2π)
validate(s)
validate(antiderivative(s))
validate(c)
validate(antiderivative(c))

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
# Mappings R → Rn
# -------------------------------------------------------------------------------------------------

c = ParametricCurve(Sin(), Cos())
@test c(0.0) == [0, 1]
@test derivativeat(c, 0) ≈ [1, 0]
@test derivativeat(c, 0, 2) ≈ [0, -1]
@test c'(0) ≈ [1, 0]
@test c''(0) ≈ [0, -1]
validate(c)

n = UnitNormal(c)
@test norm(n(0)) ≈ 1
@test abs(n(0) ⋅ c'(0)) < 1e-16

f1 = Polynomial(1, 2, 3) * c
f2 = c * Polynomial(1, 2, 3)
@test f1 == f2
@test f1(π) == f2(π)
@test codomaintype(f1) <: SVector
@test codomaintype(f2) <: SVector


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
@test domaintype(f) == SVector{3,Float64}

# Product function, two parameters
f = ProductFunction(Sin(0 .. 1), Cos(0.5 .. 5))
x = SVector(1.0, 2.0)
@test valueat(f, x) == sin(1) * cos(2)
@test gradientat(f, x) == [cos(1) * cos(2), -sin(1) * sin(2)]
# TODO @test ∇(f)(x) == gradientat(f, x)
@test hessianat(f, x) ≈ [
    -sin(1)*cos(2) -cos(1)*sin(2)
    -cos(1)*sin(2) -sin(1)*cos(2)
]
# TODO @test H(f)(x) == hessian(f, x)
@test gradientat(f, x) isa SVector{2,Float64}
@test hessianat(f, x) isa MMatrix{2,2,Float64}
@test domain(f) == (0 .. 1) × (0.5 .. 5)
# TODO Generic test

# Product function three variables
f = ProductFunction(Sin(1 .. 2), Cos(3 .. 4), Polynomial([2, 3], 6 .. 7))
x = SVector(1.0, 2.0, 3.0)
@test valueat(f, x) == sin(1) * cos(2) * 11
@test gradientat(f, x) ≈ [cos(1) * cos(2) * 11, -sin(1) * sin(2) * 11, sin(1) * cos(2) * 3]
@test hessianat(f, x) ≈ [
    -sin(1)*cos(2)*11 -cos(1)*sin(2)*11 cos(1)*cos(2)*3
    -cos(1)*sin(2)*11 -sin(1)*cos(2)*11 -sin(1)*sin(2)*3
    cos(1)*cos(2)*3 -sin(1)*sin(2)*3 0
]
# TODO @test H(f)(x) == hessian(f, x)
@test gradientat(f, x) isa SVector{3,Float64}
@test hessianat(f, x) isa MMatrix{3,3,Float64}
@test domain(f) == (1 .. 2) × (3 .. 4) × (6 .. 7)
# TODO Generic test
