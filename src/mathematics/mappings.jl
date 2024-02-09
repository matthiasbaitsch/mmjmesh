using IntervalSets
using StaticArrays

import Polynomials
import CairoMakie as cm


# XXX ???
using MMJMesh.MMJBase


# XXX move
export AllOf, AbstractMapping, AbstractMappingFromR, AbstractFunctionToR
export domaintype, codomaintype, domain, valueat, derivativeat, derivative, R, IHat
export antiderivative, integrate, sample, plot, pois, roots
export Sin, Cos, Polynomial, fromroots, lagrangepolynomials, monomials, degree


# -------------------------------------------------------------------------------------------------
# Set of all elements of type T
# -------------------------------------------------------------------------------------------------

struct AllOf{T} end
Base.eltype(::AllOf{T}) where {T} = T

Base.in(x, ::AllOf{T}) where {T} = typeof(x) <: T
Base.in(x::AbstractArray, ::AllOf{T}) where {T<:SVector} =
    return eltype(x) <: eltype(T) && length(x) == length(T)

function Base.intersect(a::AllOf, itrs...)  
    isempty(itrs) && return a
    length(itrs) == 1 && return itrs[1]
    length(itrs) == 2 && return intersect(itrs[1], itrs[2])
    return intersect(itrs[1], itrs[2:end]...)
end


Base.intersect(s, ::AllOf) = s
Base.intersect(a::AllOf{T}, ::AllOf{T}) where {T} = a

Base.intersect(a::Tuple{AllOf{T}}) where {T} = a[1]


Base.union(a::AllOf, itrs...) = a
Base.union(_, a::AllOf) = a
Base.union(a::AllOf{T}, ::AllOf{T}) where {T} = a


# -------------------------------------------------------------------------------------------------
# General concept of mapping
# -------------------------------------------------------------------------------------------------

"""
    AbstractMapping{DT,CT,D}

Abstract type for mappings from a domain `X` to a codomain `Y`. In this implementation `X` is
represented by the type `DT` and `Y` by the type `CT`. The third parameter `D` is a set of elements 
of `DT` and can be used to incorporate additional information about the domain of interest.
"""
abstract type AbstractMapping{DT,CT,D} end


# Fundamental operations

"""
    valueat(m::AbstractMapping{DT, CT}, x::DT) -> CT

Evaluate the mapping `m` at the point `x`.
"""
function valueat end

"""
    derivativeat(m::AbstractMapping{DT, CT}, x::DT, n::Int = 1) -> dt(DT, CT, n)

Evaluate the (generalized) `n`th derivative of the mapping `m` at the point `x`.
"""
function derivativeat end

"""
    derivative(m::AbstractMapping{DT, CT}, n::Int = 1) -> AbstractMapping{DT, dt(DT, CT, n)}

The (generalized) `n`th derivative of the mapping `m`.
"""
function derivative end


# Extract parameters
domaintype(::AbstractMapping{DT,CT,D}) where {DT,CT,D} = DT
codomaintype(::AbstractMapping{DT,CT,D}) where {DT,CT,D} = CT
domain(::AbstractMapping{DT,CT,D}) where {DT,CT,D} = D


# Shorthand notations for derivative and evaluation
Base.adjoint(m::AbstractMapping) = derivative(m)
(m::AbstractMapping{DT})(x::DT) where {DT} = valueat(m, x)


# Make similar array 
Base.similar(a::AbstractArray, ::Type{T}, dims::Base.OneTo...) where {T<:AbstractMapping} =
    similar(a, AbstractMapping, Base.to_shape(dims))


# Zero element w.r.t. addition
struct Zero{DT,CT,D} <: AbstractMapping{DT,CT,D} end
valueat(::Zero{DT,CT,D}, x::DT) where {DT,CT,D} = zero(DT)
derivativeat(::Zero{DT,CT,D}, x::DT) where {DT,CT,D} = DT(0)
derivative(z::Zero) = z
Base.zero(::AbstractMapping{DT,CT,D}) where {DT,CT,D} = Zero{DT,CT,D}()
Base.show(io::IO, ::Zero) = print(io, "0")


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
derivative(f::Sum, n::Int=1) = Sum(derivative(f.m1, n), derivative(f.m2, n))
derivativeat(f::Sum, x, n::Int=1) = derivativeat(f.m1, x, n) + derivativeat(f.m2, x, n)
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
derivative(f::ScaledMapping, n::Int=1) = f.a * derivative(f.m, n)
derivativeat(f::ScaledMapping, x, n::Int=1) = f.a * derivativeat(f.m, x, n)
Base.show(io::IO, s::ScaledMapping) = print(io, s.a, " * ", s.m)
Base.:(==)(m1::ScaledMapping{DT,CT,D}, m2::ScaledMapping{DT,CT,D}) where {DT,CT,D} =
    (m1.a == m2.a && m1.m == m2.m)

# Mapping times mapping, only useful if `*` is defined for type `CT`
struct ProductMapping{DT,CT,D} <: AbstractMapping{DT,CT,D}
    m1::AbstractMapping{DT,CT}
    m2::AbstractMapping{DT,CT}
    ProductMapping(
        m1::AbstractMapping{DT,CT,D1}, m2::AbstractMapping{DT,CT,D2}
    ) where {DT,CT,D1,D2} = new{DT,CT,D1 ∩ D2}(m1, m2)
end

valueat(p::ProductMapping, x) = valueat(p.m1, x) * valueat(p.m2, x)
_derivative(f::ProductMapping) = f.m1' * f.m2 + f.m1 * f.m2'
derivativeat(f::ProductMapping, x, n::Int=1) = sum(
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
function derivative(f::AbstractMapping, n::Int=1)
    n == 0 && return f
    n == 1 && return _derivative(f)
    return derivative(derivative(f, n - 1))
end

function derivativeat(f::AbstractMapping, x, n::Int=1)
    n == 0 && return f(x)
    n == 1 && return _derivativeat(f, x)
    return derivative(f, n)(x)
end


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
Base.:*(m1::AbstractMapping, m2::AbstractMapping) = ProductMapping(m1, m2)
Base.:*(a::Real, m::AbstractMapping) = a == 1 ? m : a == 0 ? zero(m) : ScaledMapping(a, m)
Base.:*(m::AbstractMapping, a::Real) = a * m
Base.:*(a::Real, m::ScaledMapping) = (a * m.a) * m.m

# /
Base.:/(m1::AbstractMapping, m2::AbstractMapping) = QuotientMapping(m1, m2)


# -------------------------------------------------------------------------------------------------
# Mappings from R
# -------------------------------------------------------------------------------------------------

# Common domains
const R = -Inf .. Inf
const RPlus = Interval{:open,:open}(0, Inf)
const R0Plus = 0 .. Inf
const IHat = -1.0 .. 1.0

# Base type for functions from R
const AbstractMappingFromR{Real,CT,D} = AbstractMapping{Real,CT,D}

"""
    pois(m::AbstractMappingFromR) -> Real[]

Return the points of interest of the mapping `m`.
"""
pois(::AbstractMappingFromR) = Real[]

"""
    roots(m::AbstractMappingFromR, I::Union{Nothing,Interval}=nothing) -> Real[]

Return the roots of the mapping `m` either in `I` or the domain of `m`.
"""
function roots(m::AbstractMappingFromR, I::Union{Nothing,Interval}=nothing)
    isnothing(I) && return roots(m, domain(m))
    return I ∩ _roots(m)
end

# Introduce ConstructedMapping interface?
pois(p::Sum) = pois(p.m1) ∪ pois(p.m2)
pois(s::ScaledMapping) = pois(s.m)
pois(p::ProductMapping) = pois(p.m1) ∪ pois(p.m2)
pois(q::QuotientMapping) = pois(q.m1) ∪ pois(q.m2) ∪ roots(q.m2)


# -------------------------------------------------------------------------------------------------
# Functions into R
# -------------------------------------------------------------------------------------------------

const AbstractFunctionToR{DT,Real,D} = AbstractMapping{DT,Real,D}


# -------------------------------------------------------------------------------------------------
# Functions R → R
# -------------------------------------------------------------------------------------------------

const AbstractFunctionRToR{D} = AbstractMapping{Real,Real,D}


# Fundamental operations

"""
    antiderivative(f::AbstractFunctionRToR) -> AbstractFunctionRToR

An antiderivative of the function `f`.
"""
function antiderivative end


"""
    integrate(f::AbstractFunctionRToR, I::Interval) -> Real

The integral of the function `f` over the interval `I`.
"""
function integrate(f::AbstractFunctionRToR, I::Interval)
    F = antiderivative(f)
    return F(rightendpoint(I)) - F(leftendpoint(I))
end


# Rules for antiderivative
antiderivative(f::ScaledMapping{Real,Real}) = f.a * antiderivative(f.m)
antiderivative(f::Sum{Real,Real}) = antiderivative(f.m1) + antiderivative(f.m2)


# TODO: Adaptive sampling
function sample(f::AbstractFunctionRToR, n=100)
    dom = domain(f)
    a, b = width(dom) != Inf ? endpoints(dom) : (0, 5)
    x = LinRange(a, b, n)
    y = f.(x)
    x, y
end

# TODO: Recipe
function plot(f::AbstractFunctionRToR, fs::AbstractFunctionRToR...)
    n = 100
    f = cm.lines(sample(f, n)...)
    for g ∈ fs
        cm.lines!(sample(g, n)...)
    end
    f
end


# -------------------------------------------------------------------------------------------------
# Special functions
# -------------------------------------------------------------------------------------------------


# Sine function
struct Sin{D} <: AbstractFunctionRToR{D}
    Sin(d=R) = new{d}()
end

valueat(::Sin, x::Real) = sin(x)
derivativeat(::Sin, x, n::Int=1) = sin(x + n * π / 2)
derivative(::Sin{D}, n::Int=1) where {D} = (-1.0)^(n ÷ 2) * (mod(n, 2) == 0 ? Sin(D) : Cos(D))
antiderivative(::Sin{D}) where {D} = -Cos(D)
Base.show(io::IO, ::Sin) = print(io, "sin(x)")


# Cosine function
struct Cos{D} <: AbstractFunctionRToR{D}
    Cos(d=R) = new{d}()
end

valueat(::Cos, x) = cos(x)
derivativeat(::Cos, x, n::Int=1) = cos(x + n * π / 2)
derivative(::Cos{D}, n::Int=1) where {D} =
    (-1.0)^((n + 1) ÷ 2) * (mod(n, 2) == 0 ? Cos(D) : Sin(D))
antiderivative(::Cos{D}) where {D} = Sin(D)
Base.show(io::IO, ::Cos) = print(io, "cos(x)")


# Polynomials
struct Polynomial{D} <: AbstractFunctionRToR{D}
    p::Polynomials.Polynomial
end

Polynomial(p::Polynomials.Polynomial, d=R) = Polynomial{d}(p)
Polynomial(c::AbstractArray, d=R) = Polynomial{d}(Polynomials.Polynomial(c))


# Basic operations
_roots(p::Polynomial{D}) where {D} = Polynomials.roots(p.p)
degree(p::Polynomial{D}) where {D} = Polynomials.degree(p.p)
fromroots(r::AbstractArray{<:Real}, d=R) = Polynomial(Polynomials.fromroots(r), d)
valueat(p::Polynomial, x::Real) = p.p(x)
derivative(p::Polynomial{D}, n::Int=1) where {D} = Polynomial(Polynomials.derivative(p.p, n), D)
derivativeat(p::Polynomial, x::Real, n::Int=1) = derivative(p, n)(x)
antiderivative(p::Polynomial{D}, n::Int=1) where {D} = Polynomial(Polynomials.integrate(p.p, n), D)
Base.show(io::IO, p::Polynomial) = print(io, p.p)
Base.:(==)(p1::Polynomial{D1}, p2::Polynomial{D2}) where {D1,D2} = (D1 == D2 && p1.p == p2.p)


# Rules
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
function monomials(p::AbstractArray{Int}, d=R)
    coeffs(pp) = [i == pp + 1 ? 1 : 0 for i in 1:pp+1]
    return [Polynomial(coeffs(n), d) for n in p]
end
