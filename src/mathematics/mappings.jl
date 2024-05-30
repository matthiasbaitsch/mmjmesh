import Polynomials as P
import FixedPolynomials as FP

using IntervalSets
using StaticArrays
using DomainSets: ×


"""
Implementation of the concept of mappings as elements of a vector space.

There are many things missing

- Nearly all elementary functions

- Image of function

- Roots of composed functions

- Inverse functions

- Printing is very basic

- ...

"""


# -------------------------------------------------------------------------------------------------
# General concept of mapping
# -------------------------------------------------------------------------------------------------

"""
    AbstractMapping{DT,CT,D}

Abstract type for mappings from a domain `X` to a codomain `Y`. In this implementation `X` is
represented by the type `DT` and `Y` by the type `CT`. The third parameter `D` is a subset of
`DT` and specifies the actual domain, use `Any` if all elements of `DT` are in the domain.
"""
abstract type AbstractMapping{DT,CT,D} <: Function end


# Fundamental operations

"""
    valueat(m::AbstractMapping{DT, CT}, x::DT) -> CT

Evaluate the mapping `m` at the point `x`.
"""
function valueat end

"""
    derivativeat(m::AbstractMapping{DT, CT}, x::DT, n::Integer = 1) -> dt(DT, CT, n)

Evaluate the (generalized) `n`th derivative of the mapping `m` at the point `x`.
"""
function derivativeat end

"""
    derivative(m::AbstractMapping{DT, CT}, n::Integer = 1) -> AbstractMapping{DT, dt(DT, CT, n)}

The (generalized) `n`th derivative of the mapping `m`.
"""
function derivative end


# Extract parameters
domaintype(::AbstractMapping{DT,CT,D}) where {DT,CT,D} = DT
codomaintype(::AbstractMapping{DT,CT,D}) where {DT,CT,D} = CT
domain(::AbstractMapping{DT,CT,D}) where {DT,CT,D} = D


# Shorthand notations for derivative and evaluation
Base.adjoint(m::AbstractMapping) = derivative(m)
(m::AbstractMapping{DT})(x::T) where {DT,T<:DT} = valueat(m, x)


# Make similar array 
function Base.similar(a::Vector, ::Type{T}, dims::AbstractUnitRange...) where {T<:AbstractMapping}
    if !isempty(dims)
        return similar(a, AbstractMapping, Base.to_shape(dims))
    else
        return Vector{AbstractMapping}(undef, length(a))
    end
end


# Neutral element w.r.t. addition
struct Zero{DT,CT,D} <: AbstractMapping{DT,CT,D} end
valueat(::Zero{DT,CT,D}, x::DT) where {DT,CT,D} = zero(DT)
derivativeat(::Zero{DT,CT,D}, x::DT) where {DT,CT,D} = DT(0)
derivative(z::Zero) = z
Base.zero(::AbstractMapping{DT,CT,D}) where {DT,CT,D} = Zero{DT,CT,D}()
Base.show(io::IO, ::Zero) = print(io, "0")


# TODO: One

# TODO: Add value

# TODO: Transform a f(b(x + c)) + d ???


# -------------------------------------------------------------------------------------------------
# Construction of mappings from mappings
# -------------------------------------------------------------------------------------------------

# Sum
struct Sum{DT,CT,D} <: AbstractMapping{DT,CT,D}
    m1::AbstractMapping{DT,CT}
    m2::AbstractMapping{DT,CT}
    Sum(
        m1::AbstractMapping{DT,CT,D1}, m2::AbstractMapping{DT,CT,D2}
    ) where {DT,CT,D1,D2} = new{DT,CT,D1 ∩ D2}(m1, m2)
end

valueat(s::Sum, x) = valueat(s.m1, x) + valueat(s.m2, x)
derivative(f::Sum, n::Integer=1) = Sum(derivative(f.m1, n), derivative(f.m2, n))
derivativeat(f::Sum, x, n::Integer=1) = derivativeat(f.m1, x, n) + derivativeat(f.m2, x, n)
Base.show(io::IO, s::Sum) = print(io, s.m1, " + ", s.m2)
Base.:(==)(m1::Sum{DT,CT,D}, m2::Sum{DT,CT,D}) where {DT,CT,D} =
    (m1.m1 == m2.m1 && m1.m2 == m2.m2)

# Number times mapping
struct ScaledMapping{DT,CT,D} <: AbstractMapping{DT,CT,D}
    a::Real
    m::AbstractMapping
    ScaledMapping(a::Real, m::AbstractMapping{DT,CT,D}) where {DT,CT,D} = new{DT,CT,D}(a, m)
end

valueat(s::ScaledMapping, x) = s.a * valueat(s.m, x)
derivative(f::ScaledMapping, n::Integer=1) = f.a * derivative(f.m, n)
derivativeat(f::ScaledMapping, x, n::Integer=1) = f.a * derivativeat(f.m, x, n)
Base.show(io::IO, s::ScaledMapping) = print(io, s.a, " * ", s.m)
Base.:(==)(m1::ScaledMapping{DT,CT,D}, m2::ScaledMapping{DT,CT,D}) where {DT,CT,D} =
    (m1.a == m2.a && m1.m == m2.m)

function ct(t1, t2)
    if t1 <: Real && t2 <: Real
        return t1
    elseif t1 <: Real && t2 <: StaticVector
        return t2
    elseif t1 <: StaticVector && t2 <: Real
        return t1
    end
    error()
end

# Mapping times mapping, only useful if `*` is defined for type `CT`
struct ProductMapping{DT,CT,D} <: AbstractMapping{DT,CT,D}
    m1::AbstractMapping
    m2::AbstractMapping
    ProductMapping(
        m1::AbstractMapping{DT,CT1,D1}, m2::AbstractMapping{DT,CT2,D2}
    ) where {DT,CT1,CT2,D1,D2} = new{DT,ct(CT1, CT2),D1 ∩ D2}(m1, m2)
end

valueat(p::ProductMapping, x) = valueat(p.m1, x) * valueat(p.m2, x)
_derivative(f::ProductMapping) = f.m1' * f.m2 + f.m1 * f.m2'
derivativeat(f::ProductMapping, x, n::Integer=1) = sum(
    [binomial(n, k) * derivativeat(f.m1, x, k) * derivativeat(f.m2, x, n - k) for k in 0:n]
)
Base.show(io::IO, m::ProductMapping) = print(io, m.m1, " * ", m.m2)
Base.:(==)(m1::ProductMapping{DT,CT,D}, m2::ProductMapping{DT,CT,D}) where {DT,CT,D} =
    (m1.m1 == m2.m1 && m1.m2 == m2.m2)

# Mapping divided by mapping, only useful if `/` is defined for type `CT`
struct QuotientMapping{DT,CT,D} <: AbstractMapping{DT,CT,D}
    m1::AbstractMapping{DT,CT}
    m2::AbstractMapping{DT,CT}
    QuotientMapping(
        m1::AbstractMapping{DT,CT,D1}, m2::AbstractMapping{DT,CT,D2}
    ) where {DT,CT,D1,D2} = new{DT,CT,D1 ∩ D2}(m1, m2)
end

valueat(q::QuotientMapping, x) = valueat(q.m1, x) / valueat(q.m2, x)
_derivative(f::QuotientMapping) = (f.m1' * f.m2 - f.m1 * f.m2') / (f.m2 * f.m2)
_derivativeat(f::QuotientMapping, x) =
    (derivativeat(f.m1, x) * f.m2(x) - f.m1(x) * derivativeat(f.m2, x)) / f.m2(x)^2
Base.show(io::IO, m::QuotientMapping) = print(io, "(", m.m1, ") / (", m.m2, ")")
Base.:(==)(m1::QuotientMapping{DT,CT,D}, m2::QuotientMapping{DT,CT,D}) where {DT,CT,D} =
    (m1.m1 == m2.m1 && m1.m2 == m2.m2)


# Composition of mappings m1 ∘ m2
struct ComposedMapping{DT,CT,D} <: AbstractMapping{DT,CT,D}
    m1::AbstractMapping
    m2::AbstractMapping
    ComposedMapping(
        m1::AbstractMapping{DT2,CT,D1}, m2::AbstractMapping{DT,DT2,D}
    ) where {DT,DT2,CT,D,D1} = new{DT,CT,D}(m1, m2)
end
Base.:∘(m1::AbstractMapping, m2::AbstractMapping) = ComposedMapping(m1, m2)

valueat(c::ComposedMapping, x) = valueat(c.m1, valueat(c.m2, x))
_derivative(f::ComposedMapping) = (f.m1' ∘ f.m2) * f.m2'
_derivativeat(f::ComposedMapping, x) = derivativeat(f.m1, valueat(f.m2, x)) * derivativeat(f.m2, x)
Base.show(io::IO, m::ComposedMapping) = print(io, m.m1, " ∘ ", m.m2)
Base.:(==)(m1::ComposedMapping{DT,CT,D}, m2::ComposedMapping{DT,CT,D}) where {DT,CT,D} =
    (m1.m1 == m2.m1 && m1.m2 == m2.m2)


# Generic implementation of higher order derivatives
function derivative(f::AbstractMapping, n::Integer=1)
    n == 0 && return f
    n == 1 && return _derivative(f)
    return derivative(derivative(f, n - 1))
end

function derivativeat(f::AbstractMapping, x, n::Integer=1)
    n == 0 && return f(x)
    n == 1 && return _derivativeat(f, x)
    return derivative(f, n)(x)
end


# -------------------------------------------------------------------------------------------------
# Interval domains
# -------------------------------------------------------------------------------------------------

# Common domains
const R = -Inf .. Inf
const RPlus = Interval{:open,:open}(0, Inf)
const R0Plus = Interval{:closed,:open}(0, Inf)
const IHat = -1.0 .. 1.0
const ReferenceInterval = IHat

# Select elements in interval
Base.intersect(s::AbstractInterval, a::AbstractVector{T}) where {T} = T[x for x in a if x ∈ s]


# -------------------------------------------------------------------------------------------------
# Mappings from R
# -------------------------------------------------------------------------------------------------

# Base type for functions from R
const MappingFromR{CT,D} = AbstractMapping{Real,CT,D}

"""
    pois(m::MappingFromR) -> Real[]

Return the points of interest of the mapping `m`.
"""
pois(::MappingFromR) = Real[]

"""
    roots(m::MappingFromR, I::Union{Nothing,Interval}=nothing) -> Real[]

Return the roots of the mapping `m` either in `I` or the domain of `m`.
"""
function roots(m::MappingFromR, I::Union{Nothing,Interval}=nothing)
    isnothing(I) && return roots(m, domain(m))
    return I ∩ _roots(m)
end

# Introduce ConstructedMapping interface?
pois(p::Sum) = pois(p.m1) ∪ pois(p.m2)
pois(s::ScaledMapping) = pois(s.m)
pois(p::ProductMapping) = pois(p.m1) ∪ pois(p.m2)
pois(q::QuotientMapping) = pois(q.m1) ∪ pois(q.m2) ∪ roots(q.m2)


# -------------------------------------------------------------------------------------------------
# Domains in Rn
# -------------------------------------------------------------------------------------------------

const ReferenceQuadrilateral = ReferenceInterval × ReferenceInterval


# -------------------------------------------------------------------------------------------------
# Mappings from Rn
# -------------------------------------------------------------------------------------------------

const MappingFromRn{N,CT,D} = AbstractMapping{SVector{N,Float64},CT,D}

valueat(m::MappingFromRn{N}, x::T...) where {N,T<:Real} = valueat(m, SVector{N,Float64}(x))
(m::MappingFromRn{N})(x::T...) where {N,T<:Real} = valueat(m, SVector{N,Float64}(x))


# -------------------------------------------------------------------------------------------------
# Functions into R
# -------------------------------------------------------------------------------------------------

const FunctionToR{DT,D} = AbstractMapping{DT,Real,D}


# -------------------------------------------------------------------------------------------------
# Functions R → R
# -------------------------------------------------------------------------------------------------

const FunctionRToR{D} = AbstractMapping{Real,Real,D}


# Fundamental operations

"""
    antiderivative(f::FunctionRToR) -> FunctionRToR

An antiderivative of the function `f`.
"""
function antiderivative end


"""
    integrate(f::FunctionRToR, I::Interval) -> Real

The integral of the function `f` over the interval `I`.
"""
function integrate(f::FunctionRToR, I::Interval)
    F = antiderivative(f)
    return F(rightendpoint(I)) - F(leftendpoint(I))
end


# Rules for antiderivative
antiderivative(f::ScaledMapping{Real,Real}) = f.a * antiderivative(f.m)
antiderivative(f::Sum{Real,Real}) = antiderivative(f.m1) + antiderivative(f.m2)


# -------------------------------------------------------------------------------------------------
# Mappings R → Rn
# -------------------------------------------------------------------------------------------------

const FunctionRToRn{N,D} = AbstractMapping{Real,SVector{N,Float64},D}

# TODO Generalize to mapping from components
# Parametric curve from components
struct ParametricCurve{N,D} <: AbstractMapping{Real,SVector{N,Float64},D}
    components::Vector{FunctionRToR{D}}
end

function ParametricCurve(components::FunctionRToR{D}...) where {D}
    n = length(components)
    cc = Vector{FunctionRToR{D}}(undef, n)
    for i ∈ eachindex(components)
        cc[i] = components[i]
    end
    return ParametricCurve{n,D}(cc)
end

valueat(u::ParametricCurve{N,D}, x::Real) where {N,D} =
    SVector{N}([c(x) for c ∈ u.components])
derivativeat(u::ParametricCurve{N,D}, x::Real, n::Integer=1) where {N,D} =
    SVector{N}([derivativeat(c, x, n) for c ∈ u.components])
_derivative(u::ParametricCurve{N,D}, n::Integer=1) where {N,D} =
    ParametricCurve{N,D}([derivative(c, n) for c ∈ u.components])
Base.show(io::IO, u::ParametricCurve) =
    print(io, "ParametricCurve[$(u.components)]")
Base.:(==)(c1::ParametricCurve{N,D}, c2::ParametricCurve{N,D}) where {N,D} =
    (c1.components == c2.components)


# Unit normal
struct UnitNormal{D} <: AbstractMapping{Real,SVector{2,Float64},D}
    u::AbstractMapping{Real,SVector{2,Float64},D}
end

function valueat(u::UnitNormal, x)
    t = derivativeat(u.u, x)
    n = SA[-t[2], t[1]]
    return n / norm(n)
end


# -------------------------------------------------------------------------------------------------
# Functions Rn → R
# -------------------------------------------------------------------------------------------------

const FunctionRnToR{N,D} = AbstractMapping{SVector{N,Float64},Real,D}
(f::FunctionRnToR{N})(x::AbstractVector) where {N} = valueat(f, SVector{N}(x))

gradientat(f::FunctionRnToR, x::AbstractArray) = derivativeat(f, x, 1)
hessianat(f::FunctionRnToR, x::AbstractArray) = derivativeat(f, x, 2)

∇(f::FunctionRnToR) = derivative(f, 1)
H(f::FunctionRnToR) = derivative(f, 2)

# Multivariate polynomial
struct MPolynomial{N,D} <: FunctionRnToR{N,D}
    p::FP.Polynomial
    MPolynomial(exponents::Matrix{Int}, coefficients::AbstractVector{T}, d=Any) where {T} =
        new{size(exponents, 1),d}(FP.Polynomial(exponents, coefficients))
end

valueat(p::MPolynomial, x::AbstractVector) = p.p(x)
# derivative(p::Polynomial{D}, n::Integer=1) where {D} = Polynomial(Polynomials.derivative(p.p, n), D)
# derivativeat(p::Polynomial, x::Real, n::Integer=1) = derivative(p, n)(x) # TODO: Efficient implementation
# antiderivative(p::Polynomial{D}, n::Integer=1) where {D} = Polynomial(Polynomials.integrate(p.p, n), D)
Base.show(io::IO, p::MPolynomial) = print(io, p.p)
Base.:(==)(p1::MPolynomial{N,D}, p2::MPolynomial{N,D}) where {N,D} = p1.p == p2.p

function Base.:(+)(p1::MPolynomial, p2::MPolynomial)
    e1 = FP.exponents(p1.p)
    c1 = FP.coefficients(p1.p)
    e2 = FP.exponents(p2.p)
    c2 = FP.coefficients(p2.p)
    m = Dict{AbstractArray{Int},Float64}()

    for (e, c) ∈ zip(eachcol(e1), c1)
        m[e] = c
    end

    for (e, c) ∈ zip(eachcol(e2), c2)
        if e ∈ keys(m)
            m[e] += c
        else
            m[e] = c
        end
    end

    return MPolynomial(stack(keys(m)), collect(values(m)), domain(p1))
end

Base.:(*)(a::Real, p::MPolynomial{N,D}) where {N,D} =
    MPolynomial(FP.exponents(p.p), a * FP.coefficients(p.p), D)


# Multivariate monomials
function mmonomials(n::Integer, p::Integer, dom; predicate=nothing)
    if n == 2
        if isnothing(predicate)
            predicate = (p1, p2) -> p1 <= p && p2 <= p
        end
        return reshape(
            [
                MPolynomial([p1; p2;;], [1.0], dom)
                for p1 in 0:p, p2 in 0:p
                if predicate(p1, p2)
            ],
            :
        )
    else
        error("Not implemented yet")
    end
end


# Product of functions, one for each parameter
struct ProductFunction{N,D} <: FunctionRnToR{N,D}
    factors
    ProductFunction(fs::FunctionRToR...) = new{length(fs),ProductDomain(domain.(fs)...)}(fs)
end

valueat(f::ProductFunction{N}, x::AbstractArray) where {N} =
    prod([f.factors[i](x[i]) for i ∈ 1:N])

function derivativeat(
    f::ProductFunction{N}, x::AbstractArray, n::Integer=1
) where {N}
    if n == 1
        return SVector{N}([
            derivativeat(f.factors[i], x[i]) *
            prod([valueat(f.factors[j], x[j]) for j ∈ 1:N if i ≠ j])
            for i ∈ 1:N
        ])
    elseif n == 2
        H = ones(MMatrix{N,N})
        d = [derivativeat(f.factors[i], x[i], n) for n ∈ 0:2, i ∈ 1:N]
        for i = 1:N, j = 1:N, k = 1:N
            if i == j == k
                H[i, j] *= d[3, k]
            elseif i == k || j == k
                H[i, j] *= d[2, k]
            else
                H[i, j] *= d[1, k]
            end
        end
        return H
    else
        throw(DomainError("Derivatives of order > 1 not implemented yet"))
    end
end


# -------------------------------------------------------------------------------------------------
# Special functions
# -------------------------------------------------------------------------------------------------


# Sine function
struct Sin{D} <: FunctionRToR{D}
    Sin(d=R) = new{d}()
end

valueat(::Sin, x::Real) = sin(x)
derivativeat(::Sin, x, n::Integer=1) = sin(x + n * π / 2)
derivative(::Sin{D}, n::Integer=1) where {D} = (-1.0)^(n ÷ 2) * (mod(n, 2) == 0 ? Sin(D) : Cos(D))
antiderivative(::Sin{D}) where {D} = -Cos(D)
Base.show(io::IO, ::Sin) = print(io, "sin(x)")


# Cosine function
struct Cos{D} <: FunctionRToR{D}
    Cos(d=R) = new{d}()
end

valueat(::Cos, x) = cos(x)
derivativeat(::Cos, x, n::Integer=1) = cos(x + n * π / 2)
derivative(::Cos{D}, n::Integer=1) where {D} =
    (-1.0)^((n + 1) ÷ 2) * (mod(n, 2) == 0 ? Cos(D) : Sin(D))
antiderivative(::Cos{D}) where {D} = Sin(D)
Base.show(io::IO, ::Cos) = print(io, "cos(x)")


# Polynomials
struct Polynomial{D} <: FunctionRToR{D}
    p::P.Polynomial
end

Polynomial(p::P.Polynomial, d=R) = Polynomial{d}(p)
Polynomial(c::AbstractArray, d=R) = Polynomial{d}(P.Polynomial(c))
Polynomial(c::T...; d=R) where {T<:Real} = Polynomial{d}(P.Polynomial(c))

_roots(p::Polynomial) = P.roots(p.p)
degree(p::Polynomial) = P.degree(p.p)
fromroots(r::AbstractArray{<:Real}, d=R) = Polynomial(P.fromroots(r), d)
valueat(p::Polynomial, x::Real) = p.p(x)
derivative(p::Polynomial{D}, n::Integer=1) where {D} = Polynomial(P.derivative(p.p, n), D)
derivativeat(p::Polynomial, x::Real, n::Integer=1) = derivative(p, n)(x) # TODO: Efficient implementation
antiderivative(p::Polynomial{D}, n::Integer=1) where {D} = Polynomial(P.integrate(p.p, n), D)
Base.show(io::IO, p::Polynomial) = print(io, p.p)
Base.:(==)(p1::Polynomial{D1}, p2::Polynomial{D2}) where {D1,D2} = (D1 == D2 && p1.p == p2.p)

Base.:+(p1::Polynomial{D1}, p2::Polynomial{D2}) where {D1,D2} = Polynomial(p1.p + p2.p, D1 ∩ D2)
Base.:*(p1::Polynomial{D1}, p2::Polynomial{D2}) where {D1,D2} = Polynomial(p1.p * p2.p, D1 ∩ D2)
Base.:*(a::Real, p::Polynomial{D}) where {D} = Polynomial(a * p.p, D)


# Lagrange polynomials
function lagrangepolynomials(c::AbstractArray{<:Real}, d=R)
    indices = 1:length(c)
    normalize(f, x) = 1 / f(x) * f
    return [normalize(fromroots(c[filter(j -> j != i, indices)], d), c[i]) for i in indices]
end


# Monomials
function monomials(p::AbstractArray{<:Integer}, d=R)
    coeffs(pp) = [i == pp + 1 ? 1 : 0 for i in 1:pp+1]
    return [Polynomial(coeffs(n), d) for n in p]
end


# -------------------------------------------------------------------------------------------------
# Operators and simplification rules
# -------------------------------------------------------------------------------------------------

struct AdHocMapping{DT,CT,D} <: AbstractMapping{DT,CT,D}
    m::Function
end
AdHocMapping(m::Function, dt=Real, ct=Real, d=Any) = AdHocMapping{dt,ct,d}(m)

valueat(m::AdHocMapping{DT}, x::T) where {DT,T<:DT} = m.m(x)


# -------------------------------------------------------------------------------------------------
# Operators and simplification rules
# -------------------------------------------------------------------------------------------------


# +
Base.:+(m1::AbstractMapping, m2::AbstractMapping) = Sum(m1, m2)
Base.:+(z::Zero, ::Zero) = z
Base.:+(::Zero, m::AbstractMapping) = m
Base.:+(::Zero, m::ScaledMapping) = m
Base.:+(m::AbstractMapping, ::Zero) = m
Base.:+(m1::ScaledMapping, m2::AbstractMapping) = m2 + m1
Base.:+(m1::AbstractMapping, m2::ScaledMapping) =
    m1 == m2.m ? (1 + m2.a) * m2.m : Sum(m1, m2)
Base.:+(m1::ScaledMapping, m2::ScaledMapping) =
    m1.m === m2.m ? (m1.a + m2.a) * m1.m : Sum(m1, m2)
Base.:+(m1::T, m2::T) where {T<:AbstractMapping} = m1 == m2 ? 2.0 * m1 : Sum(m1, m2)


# -
Base.:-(m::AbstractMapping) = -1.0 * m
Base.:-(m1::AbstractMapping, m2::AbstractMapping) = m1 + (-m2)


# *
Base.:*(m1::AbstractMapping{DT,CT1}, m2::AbstractMapping{DT,CT2}) where {DT,CT1,CT2} = ProductMapping(m1, m2)
Base.:*(m1::AbstractMapping{DT,CT1}, m2::AbstractMapping{DT,CT2}) where {DT,CT1,CT2<:Real} = ProductMapping(m2, m1)

Base.:*(a::Real, m::AbstractMapping) = a == 1 ? m : a == 0 ? zero(m) : ScaledMapping(a, m)
Base.:*(m::AbstractMapping, a::Real) = a * m
Base.:*(a::Real, m::ScaledMapping) = (a * m.a) * m.m


# /
Base.:/(m1::AbstractMapping, m2::AbstractMapping) = QuotientMapping(m1, m2)
Base.:/(a::Real, m2::AbstractMapping) = QuotientMapping(Polynomial(a), m2)


# -------------------------------------------------------------------------------------------------
# Convenience functions
# -------------------------------------------------------------------------------------------------

makefunction(f::Function, d::Rectangle) = AdHocMapping(f, SVector{2,Float64}, Real, d)
makefunction(f::Function, xrange::Interval, yrange::Interval) = makefunction(f, xrange × yrange)