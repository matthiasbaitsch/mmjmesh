using Test
using StaticArrays
using LinearAlgebra
using DomainSets: ×, Rectangle, ProductDomain

using MMJMesh
using MMJMesh.Mathematics
using MMJMesh.Mathematics: derivativetype, dimension

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
# Derivative types
# -------------------------------------------------------------------------------------------------

T1 = InR
T2 = InR2
T3 = InR3

@test derivativetype(T1, T2) === InR2
@test derivativetype(T2, T3, 1) == InRⁿˣᵐ{3,2}
@test derivativetype(T2, T3, 2) == SArray{Tuple{3,2,2},<:Real}
@test derivativetype(T3, T2, 1) == InRⁿˣᵐ{2,3}
@test derivativetype(T3, T2, 2) == SArray{Tuple{2,3,3},<:Real}

@test derivativetype(Sin()) == InR
@test derivativetype(Sin(), 2) == InR
@test derivativetype(MappingFromComponents(Sin(), Cos())) == InR2
@test derivativetype(MappingFromComponents(Sin(), Cos()), 2) == InR2
@test derivativetype(ProductFunction(Sin(), Cos())) == InR2
@test derivativetype(ProductFunction(Sin(), Cos()), 2) == InRⁿˣᵐ{2,2}
@test derivativetype(
    MappingFromComponents(
        ProductFunction(Sin(), Cos(), Sin()),
        ProductFunction(Cos(), -Sin(), -Sin())
    )
) == SMatrix{2,3,<:Real}


# -------------------------------------------------------------------------------------------------
# Zero function
# -------------------------------------------------------------------------------------------------

z = Zero{InR2,InR,QHat}()
@test z(1, 1) == 0
@test derivative(z) == Zero{InR2,InR2,QHat}()
@test derivative(z, 2) == Zero{InR2,InRⁿˣᵐ{2,2},QHat}()
@test derivativeat(z, InR2([1, 1])) == [0, 0]
@test derivativeat(z, [1, 1]) == [0, 0]
@test derivativeat(z, [1, 1], 2) == [0 0; 0 0]


# -------------------------------------------------------------------------------------------------
# One function
# -------------------------------------------------------------------------------------------------

o1 = One{InR,R}()
@test One(Sin()) === o1
@test o1(1) == 1
@test derivativeat(o1, 1) == 0
@test validate(o1, atol=1e-8)

o3 = One{InRⁿ{3},R3}()
@test valueat(o3, [1, 2, 3]) == 1
@test valueat(o3, SVector{3,Float64}([1, 2, 3])) == 1
@test o3(1, 2, 3) == 1
@test derivativeat(o3, [1, 2, 3]) == [0, 0, 0]
@test o3'(1, 2, 3) == [0, 0, 0]
@test o3''(1, 2, 3) == zeros(3, 3)
@test (3 * o3)(1, 2, 3) == 3


# -------------------------------------------------------------------------------------------------
# Mapping from components
# -------------------------------------------------------------------------------------------------

# Parametric curve
g = MappingFromComponents(Sin(), Cos())

@test validate(g)
@test g[1] == Sin()
@test valueat(g, 1.0) == [sin(1.0), cos(1.0)]
@test derivativeat(g, 1.0, 1) == [cos(1), -sin(1)]
@test derivativeat(g, 1.0, 2) ≈ [-sin(1), -cos(1)]
@test derivative(g, 1)(1.0) ≈ derivativeat(g, 1.0, 1)
@test derivative(g, 2)(1.0) ≈ derivativeat(g, 1.0, 2)

# To matrix
m1 = MappingFromComponents(Sin(), 2 * Cos())
m2 = MappingFromComponents(3 * Cos(), 4 * Sin())
m3 = MappingFromComponents(m1, m2)
@test typeof(m3) <: AbstractMapping{InR,InRⁿˣᵐ{2,2}}
@test m3(1) == stack([m1(1), m2(1)])'

# Vector field R2 -> R2
x = SVector(1.0, 2.0)
d = (0 .. 3) × (1 .. 5)
f1 = MPolynomial([1 2; 2 1], [1, 2], d)
f2 = ProductFunction(Sin(0 .. 3), Cos(1 .. 5))
g = MappingFromComponents(f1, f2)

@test validate(g)
@test valueat(g, x) == [f1(x), f2(x)]
@test g(x) == valueat(g, x)
@test divergenceat(g, x) == derivativeat(f1, x, [1, 0]) + derivativeat(f2, x, [0, 1])
@test divergence(g)(x) == divergenceat(g, x)
@test div(g)(x) == divergenceat(g, x)
@test derivativeat(g, x) == [
    derivativeat(f1, x, [1, 0]) derivativeat(f1, x, [0, 1])
    derivativeat(f2, x, [1, 0]) derivativeat(f2, x, [0, 1])
]
@test jacobian(g)(x) == derivativeat(g, x)
@test jacobianat(g, x) == derivativeat(g, x)


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
@test validate(s)
@test validate(antiderivative(s))
@test validate(c)
@test validate(antiderivative(c))
@test antiderivative(Sin(), 2) == -Sin()

# Exp
f = Exp(0 .. 2π)
@test validate(f)
@test validate(antiderivative(f))

# Identity
id = Identity(0 .. 5)
@test id(3) == 3
@test validate(id, atol=1e-6)

# Polynomials
p = Polynomial(4, 6, 1, 9, 2, -1)
@test degree(p) == 5
@test validate(p)
# XXX @test validate(antiderivative(p))

# Constructor
@test Polynomial([1, 2, 3]) == Polynomial(1, 2, 3)

# Roots and domain
@test roots(Polynomial([-1, 0, 1])) == [-1, 1]
@test roots(Polynomial([-1, 0, 1], RPlus)) == [1]
@test roots(Polynomial([-1, 0, 1], -5 .. 0)) == [-1]
@test roots(Polynomial([-1, 0, 1]), -5 .. 0) == [-1]
@test roots(Polynomial([-1, 0, 1]), 4 .. 5) == []

# Lagrange
p = [0, 1, 3, 4, 5.5]
l4 = lagrangepolynomials(p, 0 .. 5.5)
for l ∈ l4
    @test validate(l, atol=1e-4)
end
for i ∈ 1:5, j ∈ 1:5
    @test isapprox(l4[i](p[j]), i == j ? 1 : 0, atol=1e-14)
end

# From roots
c = [1, 2, 3]
p = fromroots(c)
@test roots(p) ≈ c

# Monomials
@test Polynomial([1, 2, 3, 4, 5]) == [1, 2, 3, 4, 5]' * monomials(0:4)

# Affine function from one interval into another
f = affinefunction(1 .. 2, 5 .. 1)
@test degree(f) == 1
@test f(1) == 5
@test f(2) == 1
@test integrate(f, 1 .. 2) == 3

# Affine map
f = AffineMapping(2, 3)
@test f(2) == 7
@test derivativeat(f, 2) == 2
@test derivativeat(f, 2, 2) == 0
@test derivativeat(f, 2, 3) == 0
@test f'(2) == 2

u = AffineMapping([-1 1; 2 3], [1, 2], QHat)
@test u(1, 2) == [2, 10]
@test derivativeat(u, [1, 2]) == [-1 1; 2 3]
@test derivativeat(u, [1, 2], 2) == zeros(2, 2, 2)
# TODO Derivative when needed. Requires product of something with mapping to R

f = ProductFunction(Sin(), Sin())
g = AffineMapping(Diagonal([π, π]), zeros(2))
h = f ∘ g
@test isapprox(h([1, 1]), 0.0, atol=1e-14)


# -------------------------------------------------------------------------------------------------
# Functions Rn -> R
# -------------------------------------------------------------------------------------------------

## _nn function
using MMJMesh.Mathematics: _nn
@test _nn(2, 1, 1) == [2, 0]
@test _nn(2, 1, 2, 1, 2) == [2, 2]

# Product function, two parameters
f = ProductFunction(Sin(0 .. 1), Cos(0.5 .. 5))
x = SVector(1.0, 2.0)

@test validate(f)
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

gradf = derivative(f, 1)
@test gradf[1] == derivative(f, [1, 0])
@test gradf[2] == derivative(f, [0, 1])

hf = derivative(f, 2)
@test hf[1][1] == derivative(f, [2, 0])
@test hf[1][2] == derivative(f, [1, 1])
@test hf[2][1] == derivative(f, [1, 1])
@test hf[2][2] == derivative(f, [0, 2])

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


# Product function three variables
f = ProductFunction(Sin(1 .. 2), Cos(3 .. 4), Polynomial([2, 3], 6 .. 7))
x = SVector(1.0, 2.0, 3.0)

@test validate(f)
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


# Antiderivative
ns = [3, 2]
f = ProductFunction(Sin(), Cos())
@test derivative(antiderivative(f, ns), ns) == f


# Integrals
f = ProductFunction(Sin(), Cos())
@test integrate(f, 1 .. 2, 5 .. 9) ≈ 1.31133267192572
@test integrate(f, (1 .. 2) × (5 .. 9)) ≈ 1.31133267192572


# -------------------------------------------------------------------------------------------------
# Parametric curves
# -------------------------------------------------------------------------------------------------

c = ParametricCurve(Sin(), Cos())
@test c(0.0) == [0, 1]
@test derivativeat(c, 0) ≈ [1, 0]
@test derivativeat(c, 0, 2) ≈ [0, -1]
@test c'(0) ≈ [1, 0]
@test c''(0) ≈ [0, -1]
@test validate(c)

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
# Operations and simplification rules
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
@test validate(2 * m1, rtol=1e-4)

# Add and subtract
@test m1 + zero(m1) == m1
@test m1 + m1 == 2.0 * m1
@test 2 * m1 + m1 == 3 * m1
@test m1 + 2 * m1 == 3 * m1
@test derivative(m1 + m3) == Cos() - Sin()
@test antiderivative(m1 + m3) == -Cos() + Sin()
@test m1 - m1 == zero(m1)
@test m3 - m1 == Cos() - Sin()
@test validate(m1 + m3, rtol=1e-4)


# Add or subtract constant
f = Exp()

@test (2 + f)(1) == 2 + ℯ
@test (f + 2)(1) == ℯ + 2
@test (2 - f)(1) == 2 - ℯ
@test (f - 2)(1) == ℯ - 2

@test validate(f + 2)
@test validate(2 + f)
@test validate(f - 2)
@test validate(2 - f)


# Composition m5 = Sin(x^2)
m5 = m1 ∘ m4
@test m5' == (Cos() ∘ Polynomial([0, 0, 1], 0 .. 4)) * Polynomial([0, 2], 0 .. 4)
@test isapprox(m5(√π), 0, atol=1e-15)
@test validate(m5, rtol=1e-4)

# Product
p = m1 * m2
@test validate(p)


# Quotient
p = m1 / m3
@test p(0.2) ≈ tan(0.2)
@test validate(p, rtol=1e-4)
f = 2 / Polynomial(0, 1)
@test f(0.2) ≈ 10
@test pois(f) == [0]

# Multiply by one function
f1 = Sin()
u1 = MappingFromComponents(Sin(), Cos(), Sin())
f3 = MPolynomial([4 3 2; 3 2 1; 2 1 0], [1, 2, 3])
@test o1 * f1 === f1
@test f1 * o1 === f1
@test One(f1) === o1
@test One(Sin{R}) === o1
@test o1 * u1 === u1
@test u1 * o1 === u1
@test o3 * f3 === f3
@test f3 * o3 === f3


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
@test u(2, 3) ≈ derivative(f1, [1, 0])(2, 3) * cos(2) * cos(3) -
                derivative(f1, [0, 1])(2, 3) * sin(2) * sin(3)

@test isequal(∂x * f1, ∂x(f1))


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
