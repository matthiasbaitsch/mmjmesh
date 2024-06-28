using IntervalSets
using StaticArrays
using LinearAlgebra
using DomainSets: ×, Rectangle, ProductDomain

import Polynomials as P

using MMJMesh.MMJBase


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
abstract type AbstractMapping{DT,CT,D} end


# Fundamental operations

"""
    valueat(m::AbstractMapping{DT, CT}, x::DT) -> CT

Evaluate the mapping `m` at the point `x`.
"""
function valueat end

"""
    derivativeat(m::AbstractMapping{DT, CT}, x::DT, n::Integer = 1) -> dt(DT, CT, n)
    derivativeat(m::AbstractMapping{DT, CT}, x::DT, ns::AbstractArray{Integer}) -> CT

- Evaluates the (generalized) `n`-th derivative at `x`.

- Evaluates the (generalized) `ns`-th partial derivative at `x`. For instance `derivativeat(f, x, [2, 1])` computes ``f_{xxy}(x, y)``.
"""
function derivativeat end

"""
    derivative(m::AbstractMapping{DT, CT}, n::Integer = 1) -> AbstractMapping{DT, dt(DT, CT, n)}
    derivative(m::AbstractMapping{DT, CT}, ns::AbstractArray{Integer}) -> AbstractMapping{DT, CT}

- The (generalized) `n`-th derivative of `m`.

- The (generalized) `ns`-th partial derivative of `m`.
"""
function derivative end


domaintype(::AbstractMapping{DT,CT,D}) where {DT,CT,D} = DT
codomaintype(::AbstractMapping{DT,CT,D}) where {DT,CT,D} = CT
domain(::AbstractMapping{DT,CT,D}) where {DT,CT,D} = D


# Shorthand notations for derivative and evaluation
Base.adjoint(m::AbstractMapping) = derivative(m)
(m::AbstractMapping{DT})(x::T) where {DT,T<:DT} = valueat(m, x)


# Make similar array 
function Base.similar(a::Vector, ::Type{T}, dims::AbstractUnitRange...) where {T<:AbstractMapping}
    TT = AbstractMapping{domaintype(a[1]),codomaintype(a[1])}
    if !isempty(dims)
        return similar(a, TT, Base.to_shape(dims))
    else
        return Vector{TT}(undef, length(a))
    end
end


# Neutral element w.r.t. addition
struct Zero{DT,CT,D} <: AbstractMapping{DT,CT,D} end
valueat(::Zero{DT,CT,D}, x::DT) where {DT,CT,D} = zero(CT)
derivative(z::Zero, _...) = z
antiderivative(z::Zero, _...) = z
derivativeat(::Zero{DT,CT,D}, x::DT, _...) where {DT,CT,D} = CT(0)
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

valueat(s::Sum{DT}, x::DT) where {DT} = valueat(s.m1, x) + valueat(s.m2, x)
derivative(f::Sum, n::Integer=1) = Sum(derivative(f.m1, n), derivative(f.m2, n))
derivativeat(f::Sum{DT}, x::DT, n::Integer=1) where {DT} =
    derivativeat(f.m1, x, n) + derivativeat(f.m2, x, n)
Base.show(io::IO, s::Sum) = print(io, s.m1, " + ", s.m2)
Base.:(==)(m1::Sum{DT,CT,D}, m2::Sum{DT,CT,D}) where {DT,CT,D} = (m1.m1 == m2.m1 && m1.m2 == m2.m2)

# Number times mapping
struct ScaledMapping{DT,CT,D} <: AbstractMapping{DT,CT,D}
    a::Float64
    m::AbstractMapping
    ScaledMapping(a::Real, m::AbstractMapping{DT,CT,D}) where {DT,CT,D} = new{DT,CT,D}(a, m)
end

valueat(s::ScaledMapping{DT}, x::DT) where {DT} = s.a * valueat(s.m, x)
derivative(f::ScaledMapping, n::Integer=1) = f.a * derivative(f.m, n)
derivativeat(f::ScaledMapping{DT}, x::DT, n::Integer=1) where {DT} = f.a * derivativeat(f.m, x, n)
Base.show(io::IO, s::ScaledMapping) = print(io, s.a, " * ", s.m)
Base.:(==)(m1::ScaledMapping{DT,CT,D}, m2::ScaledMapping{DT,CT,D}) where {DT,CT,D} =
    (m1.a == m2.a && m1.m == m2.m)

# Mapping times mapping, only useful if `*` is defined for type `CT`
struct ProductMapping{DT,CT,D} <: AbstractMapping{DT,CT,D}
    m1::AbstractMapping
    m2::AbstractMapping
    ProductMapping(
        m1::AbstractMapping{DT,CT1,D1}, m2::AbstractMapping{DT,CT2,D2}
    ) where {DT,CT1,CT2,D1,D2} = new{DT,Base.promote_op(*, CT1, CT2),D1 ∩ D2}(m1, m2)
end

valueat(p::ProductMapping{DT}, x::DT) where {DT} = valueat(p.m1, x) * valueat(p.m2, x)
_derivative(f::ProductMapping) = f.m1' * f.m2 + f.m1 * f.m2'
derivativeat(f::ProductMapping{DT}, x::DT, n::Integer=1) where {DT} = sum(
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

valueat(q::QuotientMapping{DT}, x::DT) where {DT} = valueat(q.m1, x) / valueat(q.m2, x)
_derivative(f::QuotientMapping) = (f.m1' * f.m2 - f.m1 * f.m2') / (f.m2 * f.m2)
_derivativeat(f::QuotientMapping{DT}, x::DT) where {DT} =
    (derivativeat(f.m1, x) * f.m2(x) - f.m1(x) * derivativeat(f.m2, x)) / f.m2(x)^2
Base.show(io::IO, m::QuotientMapping) = print(io, "(", m.m1, ") / (", m.m2, ")")
Base.:(==)(m1::QuotientMapping{DT,CT,D}, m2::QuotientMapping{DT,CT,D}) where {DT,CT,D} =
    (m1.m1 == m2.m1 && m1.m2 == m2.m2)


# Composition of mappings m1 ∘ m2
struct ComposedMapping{DT,CT,D} <: AbstractMapping{DT,CT,D}
    m1::AbstractMapping
    m2::AbstractMapping
    ComposedMapping(
        m1::AbstractMapping{DT1,CT1,D1}, m2::AbstractMapping{DT2,CT2,D2}
    ) where {DT1,CT1,D1,DT2,CT2<:DT1,D2} = new{DT2,CT1,D2}(m1, m2)
end
Base.:∘(m1::AbstractMapping, m2::AbstractMapping) = ComposedMapping(m1, m2)

valueat(c::ComposedMapping{DT}, x::DT) where {DT} = valueat(c.m1, valueat(c.m2, x))
_derivative(f::ComposedMapping) = (f.m1' ∘ f.m2) * f.m2'
_derivativeat(f::ComposedMapping{DT}, x::DT) where {DT} =
    derivativeat(f.m1, valueat(f.m2, x)) * derivativeat(f.m2, x)
Base.show(io::IO, m::ComposedMapping) = print(io, m.m1, " ∘ ", m.m2)
Base.:(==)(m1::ComposedMapping{DT,CT,D}, m2::ComposedMapping{DT,CT,D}) where {DT,CT,D} =
    (m1.m1 == m2.m1 && m1.m2 == m2.m2)


# Mapping from components
struct MappingFromComponents{DT,CT,D} <: AbstractMapping{DT,CT,D}
    components

    MappingFromComponents(components::AbstractMapping{DT,CT,D}...) where {DT,CT,D} =
        new{DT,SVector{length(components),CT},D}(collect(components))
end

valueat(m::MappingFromComponents{DT,CT,D}, x::DT) where {DT,CT,D} =
    SVector{size(CT, 1)}([valueat(m.components[i], x) for i ∈ eachindex(m.components)])
derivativeat(m::MappingFromComponents{DT,CT,D}, x::DT, n::Integer=1) where {DT,CT,D} =
    SVector{size(CT, 1)}([derivativeat(m.components[i], x, n) for i ∈ eachindex(m.components)])
derivative(m::MappingFromComponents, n::Integer=1) =
    MappingFromComponents([derivative(c, n) for c ∈ m.components]...)

function valueat(m::MappingFromComponents{DT,SVector{N,SVector{M,Float64}},D}, x::DT) where {DT,N,M,D}
    v = MMatrix{N,M,Float64}(undef)
    for i = 1:N
        v[i, :] = valueat(m.components[i], x)
    end
    v
end

Base.getindex(m::MappingFromComponents, i::Integer) = m.components[i]
Base.show(io::IO, m::MappingFromComponents) = print(io, "MappingFromComponents[$(m.components)]")
Base.:(==)(m1::MappingFromComponents, m2::MappingFromComponents) = m1.components == m2.components


# Generic implementation of higher-order derivatives

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
# Domains
# -------------------------------------------------------------------------------------------------

# Common domains
const R = -Inf .. Inf
const R2 = R^2
const R3 = R^3
const RPlus = Interval{:open,:open}(0, Inf)
const R0Plus = Interval{:closed,:open}(0, Inf)

const IHat = -1.0 .. 1.0
const ReferenceInterval = IHat

const QHat = IHat × IHat
const ReferenceQuadrilateral = QHat

# Select elements in interval
Base.intersect(s::AbstractInterval, a::AbstractVector{T}) where {T} = T[x for x in a if x ∈ s]



# -------------------------------------------------------------------------------------------------
# Mappings from R
# -------------------------------------------------------------------------------------------------

# Base type for functions from R
const MappingFromR{CT,D} = AbstractMapping{Real,CT,D}

"""
    pois(m::MappingFromR) -> Float64[]

Return the points of interest of the mapping `m`.
"""
pois(::MappingFromR) = Float64[]

# Introduce ConstructedMapping interface?
pois(p::Sum) = pois(p.m1) ∪ pois(p.m2)
pois(s::ScaledMapping) = pois(s.m)
pois(p::ProductMapping) = pois(p.m1) ∪ pois(p.m2)
pois(q::QuotientMapping) = pois(q.m1) ∪ pois(q.m2) ∪ roots(q.m2)

"""
    roots(m::MappingFromR, I::Union{Nothing,Interval}=nothing) -> Float64[]

Return the roots of the mapping `m` either in `I` or the domain of `m`.
"""
function roots(m::MappingFromR, I::Union{Nothing,Interval}=nothing)
    isnothing(I) && return roots(m, domain(m))
    return I ∩ _roots(m)
end


# -------------------------------------------------------------------------------------------------
# Mappings from Rn
# -------------------------------------------------------------------------------------------------

const MappingFromRn{N,CT,D} = AbstractMapping{SVector{N,<:Real},CT,D}

# Parameter type conversion rules
(f::MappingFromRn{N})(x::Vector{<:Real}) where {N} = valueat(f, SVector{N}(x))
(m::MappingFromRn{N})(x::Real...) where {N} = valueat(m, SVector{N}(x))

valueat(f::MappingFromRn{N}, x::Vector{<:Real}) where {N} = valueat(f, SVector{N}(x))
valueat(m::MappingFromRn{N}, x::Real...) where {N} = valueat(m, SVector{N}(x))

derivativeat(m::MappingFromRn{N}, x::Vector{<:Real}, n::Integer=1) where {N} =
    derivativeat(m, SVector{N}(x), n)

derivativeat(
    f::MappingFromRn{N}, x::Vector{<:Real}, ns::AbstractArray{<:Integer}
) where {N} = derivativeat(f, SVector{N}(x), ns)


# -------------------------------------------------------------------------------------------------
# Functions to R
# -------------------------------------------------------------------------------------------------

const FunctionToR{DT,D} = AbstractMapping{DT,Float64,D}

LinearAlgebra.dot(f1::FunctionToR{DT}, f2::FunctionToR{DT}) where {DT} = f1 * f2

# -------------------------------------------------------------------------------------------------
# Mappings to Rn
# -------------------------------------------------------------------------------------------------

const MappingToRn{DT,N,D} = AbstractMapping{DT,SVector{N,Float64},D}

LinearAlgebra.dot(u1::MappingToRn{DT,N,D}, u2::MappingToRn{DT,N,D}) where {DT,N,D} =
    sum([u1[i] * u2[i] for i = 1:N])


# -------------------------------------------------------------------------------------------------
# Functions R → R
# -------------------------------------------------------------------------------------------------

const FunctionRToR{D} = AbstractMapping{Real,Float64,D}


# Fundamental operations

"""
    antiderivative(f::FunctionRToR, n=1) -> FunctionRToR

An `n`-th antiderivative of the function `f`.
"""
function antiderivative(f::FunctionRToR, n::Integer=1)
    n == 0 && return f
    n == 1 && return _antiderivative(f)
    return antiderivative(antiderivative(f, n - 1))
end


"""
    _antiderivative(::FunctionRToR) -> FunctionRToR

Fallback for function types which do not provide higher order antiderivatives.
"""
function _antiderivative(::FunctionRToR)
    @notimplemented
end


"""
    integrate(f::FunctionRToR, I::Interval) -> Float64

The integral of the function `f` over the interval `I`.
"""
function integrate(f::FunctionRToR, I::Interval)
    F = antiderivative(f)
    return F(rightendpoint(I)) - F(leftendpoint(I))
end


# Rules for antiderivative
antiderivative(f::ScaledMapping{Real,Float64}, n::Integer) = f.a * antiderivative(f.m, n)
antiderivative(f::Sum{Real,Float64}, n::Integer) = antiderivative(f.m1, n) + antiderivative(f.m2, n)


# -------------------------------------------------------------------------------------------------
# Parametric curves R → Rn
# -------------------------------------------------------------------------------------------------

const ParametricCurve{N,D} = AbstractMapping{Real,SVector{N,Float64},D}

ParametricCurve(components::FunctionRToR...) = MappingFromComponents(components...)

# Unit normal
struct UnitNormal{D} <: ParametricCurve{2,D}
    u::ParametricCurve
    UnitNormal(u::ParametricCurve{2,D}) where {D} = new{D}(u)
end

function valueat(u::UnitNormal, x::Real)
    t = derivativeat(u.u, x)
    l = norm(t)
    return SA[-t[2]/l, t[1]/l]
end


# -------------------------------------------------------------------------------------------------
# Functions Rn → R
# -------------------------------------------------------------------------------------------------

const FunctionRnToR{N,D} = AbstractMapping{SVector{N,<:Real},Float64,D}


# Operator
struct Operator
    op
end
(op::Operator)(x) = op.op(x)
Base.:(*)(op::Operator, x) = op.op(x)


# Partial derivatives of functions R2 -> R
const ∂x = Operator(f -> derivative(f, [1, 0]))
const ∂y = Operator(f -> derivative(f, [0, 1]))
const ∂xx = Operator(f -> derivative(f, [2, 0]))
const ∂yy = Operator(f -> derivative(f, [0, 2]))
const ∂xy = Operator(f -> derivative(f, [1, 1]))


# Helper function 
function _nn(::MappingFromRn{N}, is::Integer...) where {N}
    n = MVector{N,Int}(zeros(N))
    for i ∈ is
        n[i] += 1
    end
    return n
end


# Derivatives
function derivative(f::FunctionRnToR{N}, n::Integer=1) where {N}
    if n == 1
        return MappingFromComponents([derivative(f, _nn(f, i)) for i = 1:N]...)
    elseif n == 2
        return MappingFromComponents([derivative(derivative(f, _nn(f, i))) for i = 1:N]...)
    else
        throw(DomainError("Derivatives of order > 1 not implemented yet"))
    end
end

function derivativeat(f::FunctionRnToR{N}, x::SVector{N}, n::Integer=1) where {N}
    if n == 1
        return SVector{N}([derivativeat(f, x, _nn(f, i)) for i = 1:N])
    elseif n == 2
        return SMatrix{N,N}([derivativeat(f, x, _nn(f, i, j)) for i = 1:N, j = 1:N])
    else
        throw(DomainError("Derivatives of order > 1 not implemented yet"))
    end
end

derivativeat(
    f::FunctionRnToR{N,D}, x::SVector{N}, ns::AbstractArray{<:Integer}
) where {N,D} = valueat(derivative(f, ns), x)


# Differential operators
const gradient(f::FunctionRnToR) = derivative(f, 1)
const hessian(f::FunctionRnToR) = derivative(f, 2)
const laplacian(f::FunctionRnToR{N}) where {N} = sum([derivative(f, _nn(f, i, i)) for i = 1:N])
const ∇(f::FunctionRnToR) = gradient(f)
const H(f::FunctionRnToR) = hessian(f)
const Δ(f::FunctionRnToR) = laplacian(f)

gradientat(f::FunctionRnToR, x::SVector{N}) where {N} = derivativeat(f, x, 1)
hessianat(f::FunctionRnToR, x::SVector{N}) where {N} = derivativeat(f, x, 2)
laplacianat(f::FunctionRnToR{N}, x::SVector{N}) where {N} =
    sum([derivativeat(f, x, _nn(f, i, i)) for i = 1:N])

gradientat(f::FunctionRnToR{N}, x::Vector{<:Real}) where {N} = gradientat(f, SVector{N,Float64}(x))
gradientat(m::FunctionRnToR{N}, x::Real...) where {N} = gradientat(m, SVector{N,Float64}(x))
hessianat(f::FunctionRnToR{N}, x::Vector{<:Real}) where {N} = hessianat(f, SVector{N,Float64}(x))
hessianat(m::FunctionRnToR{N}, x::Real...) where {N} = hessianat(m, SVector{N,Float64}(x))
laplacianat(f::FunctionRnToR{N}, x::Vector{<:Real}) where {N} = laplacianat(f, SVector{N,Float64}(x))
laplacianat(m::FunctionRnToR{N}, x::Real...) where {N} = laplacianat(m, SVector{N,Float64}(x))


# Integration
function integrate(f::FunctionRnToR{2}, I1::Interval, I2::Interval)
    a = leftendpoint(I1)
    b = rightendpoint(I1)
    c = leftendpoint(I2)
    d = rightendpoint(I2)
    F = antiderivative(f, [1, 1])
    return F(a, c) + F(b, d) - F(a, d) - F(b, c)
end

integrate(f::FunctionRnToR{2}, d::Rectangle) = integrate(f, component(d, 1), component(d, 2))


"""
    ProductFunction(fs::FunctionRToR...)

Construct product function such as ``h(x, y) = f(x) g(y)``.
"""
struct ProductFunction{N,D} <: FunctionRnToR{N,D}
    factors::Vector{FunctionRToR}

    ProductFunction(fs::FunctionRToR...) =
        new{length(fs),ProductDomain(domain.(fs)...)}(collect(fs))
end

valueat(f::ProductFunction{N}, x::SVector{N}) where {N} =
    prod([f.factors[i](x[i]) for i ∈ 1:N])
derivativeat(f::ProductFunction{N}, x::SVector{N}, n::AbstractArray{<:Integer}) where {N} =
    prod([derivativeat(f.factors[i], x[i], n[i]) for i in 1:N])
derivative(f::ProductFunction{N}, n::AbstractArray{<:Integer}) where {N} =
    ProductFunction([derivative(f.factors[i], n[i]) for i in 1:N]...)
antiderivative(f::ProductFunction{N}, ns::AbstractArray{<:Integer}) where {N} =
    ProductFunction([antiderivative(f.factors[i], ns[i]) for i = 1:N]...)
Base.show(io::IO, f::ProductFunction) = print(io, "ProductFunction with factors $(f.factors)")
Base.:(==)(f1::ProductFunction, f2::ProductFunction) = f1.factors == f2.factors


# -------------------------------------------------------------------------------------------------
# Vector fields
# -------------------------------------------------------------------------------------------------

const VectorField{N,D} = AbstractMapping{SVector{N,<:Real},SVector{N,Float64},D}

divergence(v::VectorField{N}) where {N} = sum([derivative(v[i], _nn(v, i)) for i = 1:N])
divergenceat(v::VectorField{N}, x::SVector{N,<:Real}) where {N} =
    sum([derivativeat(v[i], x, _nn(v, i)) for i = 1:N])
divergenceat(v::VectorField{N}, x::Real...) where {N} = divergenceat(v, SVector{N}(x))
divergenceat(v::VectorField{N}, x::Vector{<:Real}) where {N} = divergenceat(v, SVector{N}(x))

Base.div(v::VectorField) = divergence(v)


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
_antiderivative(::Sin{D}) where {D} = -Cos(D)
Base.show(io::IO, ::Sin) = print(io, "sin(x)")


# Cosine function
struct Cos{D} <: FunctionRToR{D}
    Cos(d=R) = new{d}()
end

valueat(::Cos, x::Real) = cos(x)
derivativeat(::Cos, x, n::Integer=1) = cos(x + n * π / 2)
derivative(::Cos{D}, n::Integer=1) where {D} =
    (-1.0)^((n + 1) ÷ 2) * (mod(n, 2) == 0 ? Cos(D) : Sin(D))
_antiderivative(::Cos{D}) where {D} = Sin(D)
Base.show(io::IO, ::Cos) = print(io, "cos(x)")


# Polynomials
struct Polynomial{D} <: FunctionRToR{D}
    p::P.Polynomial
end

Polynomial(p::P.Polynomial, d=R) = Polynomial{d}(p)
Polynomial(c::AbstractArray, d=R) = Polynomial{d}(P.Polynomial(c))
Polynomial(c::Real...; d=R) = Polynomial{d}(P.Polynomial(c))

_roots(p::Polynomial) = P.roots(p.p)
degree(p::Polynomial) = P.degree(p.p)
fromroots(r::AbstractArray{<:Real}, d=R) = Polynomial(P.fromroots(r), d)
valueat(p::Polynomial, x::Real) = p.p(x)
derivative(p::Polynomial{D}, n::Integer=1) where {D} = Polynomial(P.derivative(p.p, n), D)
derivativeat(p::Polynomial, x::Real, n::Integer=1) = derivative(p, n)(x) # TODO: Efficient implementation
antiderivative(p::Polynomial{D}, n::Integer=1) where {D} = Polynomial(P.integrate(p.p, n), D)
Base.show(io::IO, p::Polynomial) = print(io, p.p)
Base.:(==)(p1::Polynomial{D1}, p2::Polynomial{D2}) where {D1,D2} = (D1 == D2 && p1.p == p2.p)

Base.:(+)(p1::Polynomial{D1}, p2::Polynomial{D2}) where {D1,D2} = Polynomial(p1.p + p2.p, D1 ∩ D2)
Base.:(*)(p1::Polynomial{D1}, p2::Polynomial{D2}) where {D1,D2} = Polynomial(p1.p * p2.p, D1 ∩ D2)
Base.:(*)(a::Real, p::Polynomial{D}) where {D} = Polynomial(a * p.p, D)


# Lagrange polynomials
function lagrangepolynomials(c::AbstractArray{<:Float64}, d=R)
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
# Ad hoc mapping
# -------------------------------------------------------------------------------------------------

struct AdHocMapping{DT,CT,D} <: AbstractMapping{DT,CT,D}
    m::Function
end
AdHocMapping(m::Function, dt=Real, ct=Float64, d=R) = AdHocMapping{dt,ct,d}(m)

valueat(m::AdHocMapping, x::Real) = m.m(x)
valueat(m::AdHocMapping, x::SVector{N,<:Number}) where {N} = m.m(x)

makefunction(f::Function, d::Rectangle) = AdHocMapping(f, SVector{2,<:Real}, Float64, d)
makefunction(f::Function, xrange::Interval, yrange::Interval) = makefunction(f, xrange × yrange)


# -------------------------------------------------------------------------------------------------
# Operators and simplification rules
# -------------------------------------------------------------------------------------------------

# +
Base.:(+)(m1::AbstractMapping, m2::AbstractMapping) = Sum(m1, m2)
Base.:(+)(z::Zero, ::Zero) = z
Base.:(+)(::Zero, m::AbstractMapping) = m
Base.:(+)(::Zero, m::ScaledMapping) = m
Base.:(+)(m::AbstractMapping, ::Zero) = m
Base.:(+)(m1::ScaledMapping, m2::AbstractMapping) = m2 + m1
Base.:(+)(m1::AbstractMapping, m2::ScaledMapping) =
    m1 == m2.m ? (1 + m2.a) * m2.m : Sum(m1, m2)
Base.:(+)(m1::ScaledMapping, m2::ScaledMapping) =
    m1.m === m2.m ? (m1.a + m2.a) * m1.m : Sum(m1, m2)
Base.:(+)(m1::T, m2::T) where {T<:AbstractMapping} = m1 == m2 ? 2.0 * m1 : Sum(m1, m2)


# -
Base.:-(m::AbstractMapping) = -1.0 * m
Base.:-(m1::AbstractMapping, m2::AbstractMapping) = m1 + (-m2)


# *
Base.:(*)(m1::AbstractMapping{DT,CT1}, m2::AbstractMapping{DT,CT2}) where {DT,CT1,CT2} =
    ProductMapping(m1, m2)
Base.:(*)(m1::AbstractMapping{DT,CT1}, m2::AbstractMapping{DT,CT2}) where {DT,CT1,CT2<:Real} =
    ProductMapping(m2, m1)

Base.:(*)(a::Real, m::AbstractMapping) = a == 1 ? m : a == 0 ? zero(m) : ScaledMapping(a, m)
Base.:(*)(m::AbstractMapping, a::Real) = a * m
Base.:(*)(a::Real, m::ScaledMapping) = (a * m.a) * m.m


# /
Base.:(/)(m1::AbstractMapping, m2::AbstractMapping) = QuotientMapping(m1, m2)
Base.:(/)(a::Real, m2::AbstractMapping) = QuotientMapping(Polynomial(a), m2)

