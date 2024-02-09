using Test
using StaticArrays
using IntervalSets
using FiniteDifferences

using MMJMesh
using MMJMesh.Mathematics

include("validatemappings.jl")


# -------------------------------------------------------------------------------------------------
# AllOf
# -------------------------------------------------------------------------------------------------

RR = AllOf{Real}()
R3 = AllOf{SVector{3,Real}}()
Z3 = AllOf{SVector{3,Int}}()
a = [1, 2, 6]
b = [1, 2, 7]

@test 3 ∈ RR
@test π ∈ RR
x = SA[1, 2, 8]
@test x ∈ R3
@test x ∈ Z3
x = SA[1.1, 2, 8]
@test x ∈ R3
@test x ∉ Z3


@test intersect(RR, RR) === RR
@test intersect(RR, RR, RR) === RR
@test intersect(RR, a) == a
@test intersect(RR, a, b) == [1, 2]
@test intersect(RR, (a,)...) == [1, 2, 6]
@test intersect(a, RR) == a
@test_throws MethodError intersect(a, b, RR) == a
@test intersect(RR, a, b) == [1, 2]
@test union(RR, RR) === RR
@test union(RR, a) === RR
@test union(RR, a, b) === RR
@test union(a, RR) === RR
@test_throws MethodError union(a, b, RR) === RR


# -------------------------------------------------------------------------------------------------
# Matrix vector product
# -------------------------------------------------------------------------------------------------

m1 = Sin()
m2 = Cos()

# Needs union
[0 0; 0 0] * [m1, m1]
[1.1 2.2; 3.1 5.8] * [m1, m1]
[1.1 2.2; 3.1 5.8] * [m1, m2]

# Needs similar
rand(5, 5) * monomials(0:4)


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
p = Polynomial([4, 6, 1, 9, 2, -1])
@test degree(p) == 5
validate(p)
validate(antiderivative(p))

# Roots in domain XXX
# roots(Polynomial([-1, 0, 1]))

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
# XXX @test roots(p) ≈ c

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
# XXX validate(p, rtol=1e-4)

