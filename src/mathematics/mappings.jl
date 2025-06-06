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


## Type declaration 

"""
    AbstractMapping{DT,CT,D}

Abstract type for mappings from a domain type `DT` to a codomain type `CT`. The third parameter `D` 
is a set of `DT`s and specifies the actual domain, use `Any` if all elements of `DT` are in the 
domain. TODO: Use DomainSets.Full
"""
abstract type AbstractMapping{DT,CT,D} end


## Fundamental operations

"""
    valueat(m::AbstractMapping{DT, CT}, x::DT) -> CT

Evaluate the mapping `m` at `x`.
"""
function valueat end

"""
    derivativeat(m::AbstractMapping{DT, CT}, x::DT, n::Integer = 1) -> dt(DT, CT, n)
    derivativeat(m::AbstractMapping{DT, CT}, x::DT, ns::AbstractArray{Integer}) -> CT

- Evaluates the (generalized) `n`-th derivative at `x` where `dt(DT, CT, n)` is the derivative
  type.

- Evaluates the (generalized) `ns`-th partial derivative at `x`. For instance 
  `derivativeat(f, x, [2, 1])` computes ``f_{xxy}(x, y)``.
"""
function derivativeat end

"""
    derivative(m::AbstractMapping{DT, CT}, n::Integer = 1) -> AbstractMapping{DT, dt(DT, CT, n)}
    derivative(m::AbstractMapping{DT, CT}, ns::AbstractArray{Integer}) -> AbstractMapping{DT, CT}

- The (generalized) `n`-th derivative of `m`.

- The (generalized) `ns`-th partial derivative of `m`.
"""
function derivative end


## Fallback implementation of higher-order derivatives

function derivative(f::AbstractMapping, n::Integer=1)
    n == 0 && return f
    n == 1 && return _derivative(f)
    return derivative(derivative(f, n - 1))
end

function derivativeat(f::AbstractMapping{DT}, x::DT, n::Integer=1) where {DT}
    n == 0 && return f(x)
    n == 1 && return _derivativeat(f, x)
    return derivative(f, n)(x)
end


## Shorthand notations for derivative and evaluation

Base.adjoint(m::AbstractMapping) = derivative(m, 1)
(m::AbstractMapping{DT})(x::DT) where {DT} = valueat(m, x)


## Matrix times vector of functions

Base.promote_op(matprod, ::Type{<:Real}, t::Type{<:AbstractMapping}) = t

function Base.:(*)(A::AbstractMatrix{T}, x::AbstractVector{<:AbstractMapping{DT,CT,D}}) where {T,DT,CT,D}
    TS = AbstractMapping{DT,CT,D}
    mul!(similar(x, TS, axes(A, 1)), A, x)
end


# -------------------------------------------------------------------------------------------------
# Types of domain, codomain and derivatives
# -------------------------------------------------------------------------------------------------

""" Type of elements in the domain of the mapping. """
domaintype(::Type{<:AbstractMapping{DT}}) where {DT} = DT
domaintype(m::AbstractMapping) = domaintype(typeof(m))

""" Type of elements in the codomain of the mapping. """
codomaintype(::Type{<:AbstractMapping{DT,CT}}) where {DT,CT} = CT
codomaintype(m::AbstractMapping) = codomaintype(typeof(m))

""" Domain of the mapping. """
domain(::Type{<:AbstractMapping{DT,CT,D}}) where {DT,CT,D} = D
domain(m::AbstractMapping) = domain(typeof(m))

"""
    derivativetype(m, n)
    derivativetype(DT, CT, n)

Type of the `n`-th order derivative for the mapping m or domain type `DT` and codomain
type `CT`. For instance, the second derivative of a function R2 to R is a 2x2 matrix.
"""
derivativetype(f::T, n::Integer=1) where {T<:AbstractMapping} = derivativetype(typeof(f), n)
derivativetype(::Type{<:AbstractMapping{DT,CT}}, n::Integer) where {DT,CT} =
    derivativetype(DT, CT, n)
derivativetype(::Type{InR}, t, ::Integer=1) = t

function derivativetype(::Type{InRⁿ{N}}, ::Type{InR}, n::Int) where {N}
    n == 1 && return InRⁿ{N}
    n == 2 && return InRᵐˣⁿ{N,N}
    error("Not implemented yet")
end

function derivativetype(::Type{InRⁿ{N}}, ::Type{InRⁿ{M}}, n::Int) where {N,M}
    n == 1 && return InRᵐˣⁿ{M,N}
    n == 2 && return SArray{Tuple{M,N,N},<:Real}
    error("Not implemented yet")
end

"""
	degree(f)
	degree(f, i)

Polynomial degree of f. Returns the

- number -1 if `f`` is the zero function

- largest exponent if `f` is a polynomial

- largest sum of exponents if `f` is a multivariate polynomial

- `Inf` in all other cases

For a multivariate mapping `f`, the invocation `degree(f, i)` returns the polynomial degree 
in the i-th parameter.
"""
degree(::AbstractMapping, i=1) = Inf


# -------------------------------------------------------------------------------------------------
# Neutral elements
# -------------------------------------------------------------------------------------------------

""" Neutral element w.r.t. addition. """
struct Zero{DT,CT,D} <: AbstractMapping{DT,CT,D} end

degree(::Zero) = -1
valueat(::Zero{DT,CT}, x::DT) where {DT,CT} = zero(CT)
derivative(::Zero{DT,CT,D}, n::Int=1) where {DT,CT,D} = Zero{DT,derivativetype(DT, CT, n),D}()
derivativeat(::Zero{DT,CT}, x::DT, n::Int=1) where {DT,CT} = zero(derivativetype(DT, CT, n))

Base.zero(::Type{<:AbstractMapping{DT,CT,D}}) where {DT,CT,D} = Zero{DT,CT,D}()
Base.zero(m::AbstractMapping) = zero(typeof(m))
Base.show(io::IO, ::Zero) = print(io, "0(x)")
Base.isequal(::Zero{DT,CT,D}, ::Zero{DT,CT,D}) where {DT,CT,D} = true


""" Neutral element w.r.t. multiplication. """
struct One{DT,D} <: AbstractMapping{DT,InR,D} end

One(::Type{<:AbstractMapping{DT,CT,D}}) where {DT,CT,D} = One{DT,D}()
One(m::AbstractMapping) = One(typeof(m))

degree(::One) = 0
valueat(::One{DT}, x::DT) where {DT} = 1.0
derivative(::One{DT,D}, n::Integer=1) where {DT,D} = Zero{DT,derivativetype(DT, InR, n),D}()
derivativeat(::One{DT}, x::DT, n::Integer=1) where {DT} = zero(derivativetype(DT, InR, n))

Base.one(::Type{<:AbstractMapping{DT,CT,D}}) where {DT,CT,D} = One{DT,D}()
Base.one(m::AbstractMapping) = one(typeof(m))
Base.show(io::IO, ::One) = print(io, "1(x)")
Base.isequal(::One{DT,D}, ::One{DT,D}) where {DT,D} = true


# -------------------------------------------------------------------------------------------------
# Identity and constant mappings
# -------------------------------------------------------------------------------------------------

""" Identity function. """
struct Identity{DT,D} <: AbstractMapping{DT,DT,D} end

Identity(d::T) where {T<:AbstractInterval} = Identity{InR,d}()

degree(::Identity) = 1
valueat(::Identity{DT}, x::DT) where {DT} = x
_derivative(::Identity{DT,D}) where {DT,D} = One{DT,D}()
_derivativeat(::Identity{DT}, x::DT, n::Int=1) where {DT} = one(derivativetype(DT, DT, 1))
antiderivative(::Identity{InR,D}, n::Int=1) where {D} =
    1 / factorial(n + 1) * Polynomial(_monomial_coeffs(n + 1), D)
Base.show(io::IO, ::Identity) = print(io, "identity(x)")


""" An independent variable of a function. """
parameter(d::AbstractInterval=R) = Identity(d)


""" Constant mapping. """
struct ConstantMapping{DT,CT,D} <: AbstractMapping{DT,CT,D}
    c::Any
    ConstantMapping(c, dt, d) = new{dt,typeof(c),d}(c)
end

valueat(m::ConstantMapping{DT}, ::DT) where {DT} = m.c


# -------------------------------------------------------------------------------------------------
# Mappings from R
# -------------------------------------------------------------------------------------------------

""" Mappings from the real numbers. """
const MappingFromR{CT,D} = AbstractMapping{InR,CT,D}


"""
    pois(m::MappingFromR) -> Float64[]

Return the points of interest of the mapping `m`.
"""
pois(::MappingFromR) = Float64[]


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

""" Mappings from ``R^n``. """
const MappingFromRn{N,CT,D} = AbstractMapping{InRⁿ{N},CT,D}


# -------------------------------------------------------------------------------------------------
# Functions to R (scalar valued)
# -------------------------------------------------------------------------------------------------

""" Function to ``R``. """
const FunctionToR{DT,D} = AbstractMapping{DT,InR,D}

LinearAlgebra.dot(f1::FunctionToR{DT}, f2::FunctionToR{DT}) where {DT} = f1 * f2


# -------------------------------------------------------------------------------------------------
# Mappings to Rn (vector valued)
# -------------------------------------------------------------------------------------------------

const MappingToRn{DT,N,D} = AbstractMapping{DT,InRⁿ{N},D}

components(u::MappingToRn{DT,N}) where {DT,N} = [u[i] for i = 1:N]

Base.length(::MappingToRn{DT,N}) where {DT,N} = N
Base.eltype(::MappingToRn{DT,N,D}) where {DT,N,D} = FunctionToR{DT,D}
Base.getindex(::MappingToRn, i::Integer) = @notimplemented
Base.iterate(u::MappingToRn, idx=1) = idx > length(u) ? nothing : (u[idx], idx + 1)


# -------------------------------------------------------------------------------------------------
# Mappings to Rmxn (matrix valued)
# -------------------------------------------------------------------------------------------------

const MappingToRmxn{DT,M,N,D} = AbstractMapping{DT,InRᵐˣⁿ{M,N},D}

Base.transpose(::MappingToRmxn) = @notimplemented
Base.:(*)(::MappingToRmxn, x::RealVecOrMat) = @notimplemented


# -------------------------------------------------------------------------------------------------
# Functions R → R
# -------------------------------------------------------------------------------------------------

""" Function ``R \\to R``. """
const FunctionRToR{D} = AbstractMapping{InR,InR,D}


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
    if isempty(I)
        return 0.0
    end
    try
        F = antiderivative(f)
        return F(rightendpoint(I)) - F(leftendpoint(I))
    catch _
        return quadgk(f, leftendpoint(I), rightendpoint(I))[1]
    end
end


# -------------------------------------------------------------------------------------------------
# Parametric curves R → Rn
# -------------------------------------------------------------------------------------------------

""" Parametric curve``R \\to R^n``. """
const ParametricCurve{N,D} = AbstractMapping{InR,InRⁿ{N},D}

ParametricCurve(components::FunctionRToR...) = MappingFromComponents(components...)


## Unit normal

struct UnitNormal{D} <: ParametricCurve{2,D}
    u::ParametricCurve{2,D}
    UnitNormal(u::ParametricCurve{2,D}) where {D} = new{D}(u)
end

function valueat(u::UnitNormal, x::InR)
    t = derivativeat(u.u, x)
    l = norm(t)
    return SA[-t[2]/l, t[1]/l]
end


# -------------------------------------------------------------------------------------------------
# Functions Rn → R
# -------------------------------------------------------------------------------------------------

""" Function ``R^n \\to R`` of several real variables. """
const FunctionRnToR{N,D} = AbstractMapping{InRⁿ{N},InR,D}


## Helper functions

"""
    _nn(n, indices...)

Converts sequence of variable indices to vector of length `n` indicating orders
of partial derivatives. For instance, calling `_nn(2, 1, 2, 1, 1)` representing 
``w_{xyxx}`` returns `[3, 1]` which corresponds to ``w_{xxxy}``.
"""
function _nn(n::Integer, indices::Integer...)
    n = zeros(Int, n)
    for i ∈ indices
        n[i] += 1
    end
    return n
end

"""
    _derivative(f::FunctionRnToR{N}, n::Integer)

Default implementation of `n`-th order derivative of multivariate function `f`. Currently
implemented for 

- `n = 1`: Returns the gradient function

- `n = 2`: Returns the Hessian function
"""
function _derivative(f::FunctionRnToR{N}, n::Integer) where {N}
    if n == 1
        return MappingFromComponents([derivative(f, _nn(N, i)) for i = 1:N]...)
    elseif n == 2
        return MappingFromComponents([derivative(derivative(f, _nn(N, i))) for i = 1:N]...)
    else
        throw(DomainError("Derivatives of order > 1 not implemented yet"))
    end
end

"""
    _derivativeat(f::FunctionRnToR{N}, x::InRⁿ{N}, n::Integer)

Default implementation for the evaluation `n`-th order derivative of multivariate function `f`
at position `x`. Currently implemented for

- `n = 1`: Returns the gradient of `f` at `x`

- `n = 2`: Returns the Hessian of `f` at `x`
"""
function _derivativeat(f::FunctionRnToR{N}, x::InRⁿ{N}, n::Integer) where {N}
    if n == 1
        return SVector{N}([derivativeat(f, x, _nn(N, i)) for i = 1:N])
    elseif n == 2
        return SMatrix{N,N}([derivativeat(f, x, _nn(N, i, j)) for i = 1:N, j = 1:N])
    else
        throw(DomainError("Derivatives of order > 1 not implemented yet"))
    end
end

"""
    _derivativeat(f::FunctionRnToR{N,D}, x::InRⁿ{N}, ns::AbstractArray{<:Integer})

Default implementation for the evaluation of the `ns` partial derivative of `f` at `x`. Note
that this is an inefficient implementation which computes the derivative function first.
"""
_derivativeat(
    f::FunctionRnToR{N,D}, x::InRⁿ{N}, ns::AbstractArray{<:Integer}
) where {N,D} = valueat(derivative(f, ns), x)


## Differential operators

# Operator struct for operator matrices
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

# Gradient, Hessian, Laplacian
""" Gradient of function `f` """
const gradient(f::FunctionRnToR) = derivative(f, 1)

""" Hessian of function `f` """
const hessian(f::FunctionRnToR) = derivative(f, 2)

""" Laplacian of function `f` """
const laplacian(f::FunctionRnToR{N}) where {N} = sum([derivative(f, _nn(N, i, i)) for i = 1:N])

""" `∇(f)` returns the gradient of function `f` """
const ∇(f::FunctionRnToR) = gradient(f)

""" `H(f)` returns the Hessian of function `f` """
const H(f::FunctionRnToR) = hessian(f)

""" `Δ(f)` returns the Laplacian of function `f` """
const Δ(f::FunctionRnToR) = laplacian(f)

gradientat(f::FunctionRnToR, x::InRⁿ{N}) where {N} = derivativeat(f, x, 1)
hessianat(f::FunctionRnToR, x::InRⁿ{N}) where {N} = derivativeat(f, x, 2)
laplacianat(f::FunctionRnToR{N}, x::InRⁿ{N}) where {N} =
    sum([derivativeat(f, x, _nn(N, i, i)) for i = 1:N])

# Integrate functions R2 -> R over rectangle
function integrate(f::FunctionRnToR{2}, I1::Interval, I2::Interval)
    a = leftendpoint(I1)
    b = rightendpoint(I1)
    c = leftendpoint(I2)
    d = rightendpoint(I2)
    F = antiderivative(f, [1, 1])
    return F(a, c) + F(b, d) - F(a, d) - F(b, c)
end
integrate(f::FunctionRnToR{2}, d::DomainSets.Rectangle) =
    integrate(f, DomainSets.component(d, 1), DomainSets.component(d, 2))


# -------------------------------------------------------------------------------------------------
# Mappings Rn to Rm
# -------------------------------------------------------------------------------------------------

""" Vector valued mapping ''R^n \\to R^m'' with ``n,m > 1`` of several real variables. """
const MappingRnToRm{N,M,D} = AbstractMapping{InRⁿ{N},InRⁿ{M},D}

"""
    jacobian(u::MappingRnToRm)

Jacobian of the mapping `u`.
"""
jacobian(u::MappingRnToRm) = derivative(u)

"""
    jacobianat(u::MappingRnToRm, x::InRⁿ)

Evaluates the jacobian of the mapping `u` at point `x`.
"""
jacobianat(u::MappingRnToRm, x::InRⁿ{N}) where {N} = derivativeat(u, x)


# -------------------------------------------------------------------------------------------------
# Vector fields
# -------------------------------------------------------------------------------------------------

""" Vector field ''R^n \\to R^n'' with ``n > `1`. """
const VectorField{N,D} = MappingRnToRm{N,N,D}

"""
    divergence(v::VectorField)

Divergence of the vector field `v`.
"""
divergence(v::VectorField{N}) where {N} = sum([derivative(v[i], _nn(N, i)) for i = 1:N])

"""
    divergenceat(v::VectorField, x::InRⁿ)

Evaluates the divergence of the vector field `v` at point `x`.
"""
divergenceat(v::VectorField{N}, x::InRⁿ{N}) where {N} =
    sum([derivativeat(v[i], x, _nn(N, i)) for i = 1:N])

"""
    div(v)

Shorthand notation for `divergence(v)`
"""
Base.div(v::VectorField) = divergence(v)


# -------------------------------------------------------------------------------------------------
# Construction of mappings from mappings
# -------------------------------------------------------------------------------------------------


## Sum of mappings

struct Sum{DT,CT,D} <: AbstractMapping{DT,CT,D}
    m1::AbstractMapping{DT,CT}
    m2::AbstractMapping{DT,CT}
    Sum(
        m1::AbstractMapping{DT,CT,D1}, m2::AbstractMapping{DT,CT,D2}
    ) where {DT,CT,D1,D2} = new{DT,CT,D1 ∩ D2}(m1, m2)
end

pois(p::Sum) = pois(p.m1) ∪ pois(p.m2)
valueat(s::Sum{DT}, x::DT) where {DT} = valueat(s.m1, x) + valueat(s.m2, x)
derivative(f::Sum, n::Integer=1) = Sum(derivative(f.m1, n), derivative(f.m2, n))
derivativeat(f::Sum{DT}, x::DT, n::Integer=1) where {DT} =
    derivativeat(f.m1, x, n) + derivativeat(f.m2, x, n)
antiderivative(f::Sum{InR,InR}, n::Integer=1) =
    antiderivative(f.m1, n) + antiderivative(f.m2, n)

Base.show(io::IO, s::Sum) = print(io, s.m1, " + ", s.m2)
Base.:(==)(m1::Sum{DT,CT,D}, m2::Sum{DT,CT,D}) where {DT,CT,D} = (m1.m1 == m2.m1 && m1.m2 == m2.m2)


## Number times mapping

struct ScaledMapping{DT,CT,D} <: AbstractMapping{DT,CT,D}
    a::Real
    m::AbstractMapping
    ScaledMapping(a::Real, m::AbstractMapping{DT,CT,D}) where {DT,CT,D} = new{DT,CT,D}(a, m)
end

pois(s::ScaledMapping) = pois(s.m)
degree(s::ScaledMapping) = degree(s.m)
valueat(s::ScaledMapping{DT}, x::DT) where {DT} = s.a * valueat(s.m, x)
derivative(f::ScaledMapping, n::Integer=1) = f.a * derivative(f.m, n)
derivativeat(f::ScaledMapping{DT}, x::DT, n::Integer=1) where {DT} = f.a * derivativeat(f.m, x, n)
antiderivative(f::ScaledMapping{InR,InR}, n::Integer=1) =
    f.a * antiderivative(f.m, n)


Base.show(io::IO, s::ScaledMapping) = print(io, s.a, " * ", s.m)
Base.:(==)(m1::ScaledMapping{DT,CT,D}, m2::ScaledMapping{DT,CT,D}) where {DT,CT,D} =
    (isequal(m1.a, m2.a) && m1.m == m2.m)


## Mapping times mapping, only useful if `*` is defined for type `CT`

struct ProductMapping{DT,CT,D} <: AbstractMapping{DT,CT,D}
    m1::AbstractMapping
    m2::AbstractMapping
end

ProductMapping(
    m1::AbstractMapping{DT,CT1,D1}, m2::AbstractMapping{DT,CT2,D2}
) where {DT,CT1,CT2,D1,D2} = ProductMapping{DT,Base.promote_op(*, CT1, CT2),D1 ∩ D2}(m1, m2)

ProductMapping(m::AbstractMapping{DT,CT,D1}, f::AbstractMapping{DT,InR,D2}) where {DT,CT,D1,D2} =
    ProductMapping{DT,CT,D1 ∩ D2}(f, m)

pois(p::ProductMapping) = pois(p.m1) ∪ pois(p.m2)
valueat(p::ProductMapping{DT}, x::DT) where {DT} = valueat(p.m1, x) * valueat(p.m2, x)
_derivative(f::ProductMapping) = f.m1' * f.m2 + f.m1 * f.m2'
derivativeat(f::ProductMapping{DT}, x::DT, n::Integer=1) where {DT} = sum(
    [binomial(n, k) * derivativeat(f.m1, x, k) * derivativeat(f.m2, x, n - k) for k in 0:n]
)

Base.show(io::IO, m::ProductMapping) = print(io, m.m1, " * ", m.m2)
Base.:(==)(m1::ProductMapping{DT,CT,D}, m2::ProductMapping{DT,CT,D}) where {DT,CT,D} =
    (m1.m1 == m2.m1 && m1.m2 == m2.m2)
Base.promote_op(*, ::Type{InR}, ::Type{InRⁿ{N}}) where {N} = InRⁿ{N}
Base.promote_op(*, ::Type{InR}, ::Type{InR}) = InR


## Mapping divided by mapping, only useful if `/` is defined for type `CT`

struct QuotientMapping{DT,CT,D} <: AbstractMapping{DT,CT,D}
    m1::AbstractMapping{DT,CT}
    m2::AbstractMapping{DT,CT}
    QuotientMapping(
        m1::AbstractMapping{DT,CT,D1}, m2::AbstractMapping{DT,CT,D2}
    ) where {DT,CT,D1,D2} = new{DT,CT,D1 ∩ D2}(m1, m2)
end

pois(q::QuotientMapping) = pois(q.m1) ∪ pois(q.m2) ∪ roots(q.m2)
valueat(q::QuotientMapping{DT}, x::DT) where {DT} = valueat(q.m1, x) / valueat(q.m2, x)
_derivative(f::QuotientMapping) = (f.m1' * f.m2 - f.m1 * f.m2') / (f.m2 * f.m2)
_derivativeat(f::QuotientMapping{DT}, x::DT) where {DT} =
    (derivativeat(f.m1, x) * f.m2(x) - f.m1(x) * derivativeat(f.m2, x)) / f.m2(x)^2

Base.show(io::IO, m::QuotientMapping) = print(io, "(", m.m1, ") / (", m.m2, ")")
Base.:(==)(m1::QuotientMapping{DT,CT,D}, m2::QuotientMapping{DT,CT,D}) where {DT,CT,D} =
    (m1.m1 == m2.m1 && m1.m2 == m2.m2)


## Composition of mappings m1 ∘ m2

struct ComposedMapping{DT,CT,D} <: AbstractMapping{DT,CT,D}
    m1::AbstractMapping
    m2::AbstractMapping
    ComposedMapping(
        m1::AbstractMapping{DT1,CT1,D1}, m2::AbstractMapping{DT2,CT2,D2}
    ) where {DT1,CT1,D1,DT2,CT2<:DT1,D2} = new{DT2,CT1,D2}(m1, m2)
end
Base.:(∘)(m1::AbstractMapping, m2::AbstractMapping) = ComposedMapping(m1, m2)

valueat(c::ComposedMapping{DT}, x::DT) where {DT} = valueat(c.m1, valueat(c.m2, x))
_derivative(f::ComposedMapping) = (f.m1' ∘ f.m2) * f.m2'
_derivativeat(f::ComposedMapping{DT}, x::DT) where {DT} =
    derivativeat(f.m1, valueat(f.m2, x)) * derivativeat(f.m2, x)

function antiderivative(f::ComposedMapping{InR,InR}, n::Integer=1)
    g = f.m1
    h = f.m2

    if degree(h) == 1
        return (1 / derivativeat(h, 0))^n * (antiderivative(g, n) ∘ h)
    else
        @notimplemented
    end
end

Base.show(io::IO, m::ComposedMapping) = print(io, m.m1, " ∘ ", m.m2)
Base.:(==)(m1::ComposedMapping{DT,CT,D}, m2::ComposedMapping{DT,CT,D}) where {DT,CT,D} =
    (m1.m1 == m2.m1 && m1.m2 == m2.m2)


## Product of functions R to R

struct ProductFunction{N,D} <: FunctionRnToR{N,D}
    factors::Vector{FunctionRToR}

    ProductFunction(fs::FunctionRToR...) =
        new{length(fs),ProductDomain(domain.(fs)...)}(collect(fs))
end

valueat(f::ProductFunction{N}, x::InRⁿ{N}) where {N} =
    prod([f.factors[i](x[i]) for i ∈ 1:N])

derivativeat(f::ProductFunction{N}, x::InRⁿ{N}, n::AbstractArray{<:Integer}) where {N} =
    prod([derivativeat(f.factors[i], x[i], n[i]) for i in 1:N])
derivativeat(f::ProductFunction{N}, x::InRⁿ{N}, n::Integer=1) where {N} =
    _derivativeat(f, x, n)

derivative(f::ProductFunction{N}, n::AbstractArray{<:Integer}) where {N} =
    ProductFunction([derivative(f.factors[i], n[i]) for i in 1:N]...)
derivative(f::ProductFunction{N}, n::Integer=1) where {N} = _derivative(f, n)

antiderivative(f::ProductFunction{N}, ns::AbstractArray{<:Integer}) where {N} =
    ProductFunction([antiderivative(f.factors[i], ns[i]) for i = 1:N]...)

Base.show(io::IO, f::ProductFunction) = print(io, "ProductFunction with factors $(f.factors)")
Base.:(==)(f1::ProductFunction, f2::ProductFunction) = f1.factors == f2.factors


## Mappings D -> R^n or D -> R^nxm from components

"""
    MappingFromComponents(c1, c2, ..., cn)

Construct a mapping ``X \\to Z`` from Mappings ``c_k: X \\to Y, k = 1, \\dots, n`` 
such that ``Z = Y^n``.
"""
struct MappingFromComponents{DT,CT,D,N} <: AbstractMapping{DT,CT,D}
    components::SVector{N}
    MappingFromComponents(components::SVector{N,<:FunctionToR{DT,D}}) where {N,DT,D} =
        new{DT,InRⁿ{N},D,N}(components)
    MappingFromComponents(components::SVector{M,<:MappingToRn{DT,N,D}}) where {M,DT,N,D} =
        new{DT,InRᵐˣⁿ{M,N},D,M}(components)
end

MappingFromComponents(components::AbstractMapping...) = MappingFromComponents(collect(components))
MappingFromComponents(components::Vector{<:AbstractMapping}) =
    MappingFromComponents(SVector{length(components)}(components))

# Value at
valueat(m::MappingFromComponents{DT,CT,D}, x::DT) where {DT,CT,D} =
    SVector{size(CT, 1)}([valueat(m.components[i], x) for i ∈ eachindex(m.components)])

function valueat(m::MappingFromComponents{DT,InRᵐˣⁿ{N,M}}, x::DT) where {DT,N,M}
    v = MMatrix{N,M,Float64}(undef)
    for i = 1:N
        v[i, :] = valueat(m.components[i], x)
    end
    return v
end

# Derivatives
derivative(m::MappingFromComponents, n::Integer=1) =
    MappingFromComponents([derivative(c, n) for c ∈ m.components]...)
derivativeat(m::MappingFromComponents{DT,CT,D}, x::DT, n::Integer=1) where {DT,CT,D} =
    SVector{size(CT, 1)}([derivativeat(m.components[i], x, n) for i ∈ eachindex(m.components)])

function derivativeat(
    m::MappingFromComponents{InRⁿ{N},InRⁿ{M}}, x::InRⁿ{N}, n::Integer=1
) where {N,M}
    @assert n == 1
    J = MMatrix{M,N,Float64}(undef)
    for i = 1:M
        J[i, :] = derivativeat(m.components[i], x)
    end
    return J
end

# Dot and matrix products
LinearAlgebra.dot(x::RealVec, m::MappingFromComponents{DT,InRⁿ{N}}) where {DT,N} =
    dot(x, m.components)
Base.:(*)(A::RealMat, m::MappingFromComponents{DT,InRⁿ{N}}) where {DT,N} =
    MappingFromComponents(A * m.components)

# From Base
Base.eltype(m::MappingFromComponents) = eltype(m.components)
Base.length(m::MappingFromComponents) = length(m.components)
Base.getindex(m::MappingFromComponents, i::Integer) = m.components[i]
Base.show(io::IO, m::MappingFromComponents) = print(io, "MappingFromComponents[$(m.components)]")
Base.:(==)(m1::MappingFromComponents, m2::MappingFromComponents) = m1.components == m2.components


# -------------------------------------------------------------------------------------------------
# Special functions
# -------------------------------------------------------------------------------------------------


## Sine function

struct Sin{D} <: FunctionRToR{D}
    Sin(d=R) = new{d}()
end

Sin(::Identity{DT,D}) where {DT,D} = Sin(D)
Sin(f::T) where {T<:FunctionToR} = Sin() ∘ f
Base.sin(f::T) where {T<:FunctionToR} = Sin(f)

valueat(::Sin, x::InR) = sin(x)
derivativeat(::Sin, x::InR, n::Integer=1) = sin(x + n * π / 2)
derivative(::Sin{D}, n::Integer=1) where {D} = (-1.0)^(n ÷ 2) * (mod(n, 2) == 0 ? Sin(D) : Cos(D))
_antiderivative(::Sin{D}) where {D} = -Cos(D)
Base.show(io::IO, ::Sin) = print(io, "sin(x)")


## Cosine function

struct Cos{D} <: FunctionRToR{D}
    Cos(d=R) = new{d}()
end

Cos(::Identity{DT,D}) where {DT,D} = Cos(D)
Cos(f::T) where {T<:FunctionToR} = Cos() ∘ f
Base.cos(f::T) where {T<:FunctionToR} = Cos(f)

valueat(::Cos, x::InR) = cos(x)
derivativeat(::Cos, x::InR, n::Integer=1) = cos(x + n * π / 2)
derivative(::Cos{D}, n::Integer=1) where {D} =
    (-1.0)^((n + 1) ÷ 2) * (mod(n, 2) == 0 ? Cos(D) : Sin(D))
_antiderivative(::Cos{D}) where {D} = Sin(D)
Base.show(io::IO, ::Cos) = print(io, "cos(x)")


## Exp function

struct Exp{D} <: FunctionRToR{D}
    Exp(d=R) = new{d}()
end

Exp(::Identity{DT,D}) where {DT,D} = Exp(D)
Exp(f::T) where {T<:FunctionRToR} = Exp() ∘ f
Base.exp(f::T) where {T<:FunctionRToR} = Exp(f)

valueat(::Exp, x::InR) = exp(x)
derivativeat(::Exp, x::InR, ::Integer=1) = exp(x)
derivative(f::Exp, ::Integer=1) = f
antiderivative(f::Exp, ::Integer=1) = f
Base.show(io::IO, ::Exp) = print(io, "exp(x)")


"""
    Polynomial([a0, a1, a2, ...], d=R)

Construct polynomial `a0 + a1 x + a2 x^2 + ...` defined on domain `d`.
"""
struct Polynomial{D} <: FunctionRToR{D}
    p::Polynomials.Polynomial
end

Polynomial(p::Polynomials.Polynomial, d=R) = Polynomial{d}(p)
Polynomial(c::AbstractArray, d=R) = Polynomial{d}(Polynomials.Polynomial(c))
Polynomial(c::Real...; d=R) = Polynomial{d}(Polynomials.Polynomial(c))

coefficients(p::Polynomial) = Polynomials.coeffs(p.p)

_roots(p::Polynomial) = Polynomials.roots(p.p)
degree(p::Polynomial) = Polynomials.degree(p.p)
fromroots(roots::AbstractArray{<:Real}, d=R) = Polynomial(Polynomials.fromroots(roots), d)
valueat(p::Polynomial, x::InR) = p.p(x)
derivative(p::Polynomial{D}, n::Integer=1) where {D} =
    Polynomial(Polynomials.derivative(p.p, n), D)
derivativeat(p::Polynomial, x::InR, n::Integer=1) = derivative(p, n)(x)
antiderivative(p::Polynomial{D}, n::Integer=1) where {D} =
    Polynomial(Polynomials.integrate(p.p, n), D)

Base.show(io::IO, p::Polynomial) = print(io, p.p)
Base.:(+)(p1::Polynomial{D1}, p2::Polynomial{D2}) where {D1,D2} = Polynomial(p1.p + p2.p, D1 ∩ D2)
Base.:(*)(a::Num, p::Polynomial{D}) where {D} = Polynomial(a * p.p, D)
Base.:(*)(a::Real, p::Polynomial{D}) where {D} = Polynomial(a * p.p, D)
Base.:(*)(p1::Polynomial{D1}, p2::Polynomial{D2}) where {D1,D2} = Polynomial(p1.p * p2.p, D1 ∩ D2)
Base.:(==)(p1::Polynomial{D1}, p2::Polynomial{D2}) where {D1,D2} = (D1 == D2 && p1.p == p2.p)


"""
    lagrangepolynomials(x, D=R)

Lagrange polynomials for nodes ``x_0, \\dots, x_k`` defined on the domain `D`.
"""
function lagrangepolynomials(c::AbstractArray{<:Float64}, d=R)
    indices = 1:length(c)
    normalize(f, x) = 1 / f(x) * f
    return [normalize(fromroots(c[filter(j -> j != i, indices)], d), c[i]) for i in indices]
end

_monomial_coeffs(n) = [i == n + 1 for i = 1:n+1]

"""
    monomials(p, D=R)

Monomials of degree ``p_1, \\dots, p_n`` defined on the domain `D`.
"""
function monomials(p::AbstractArray{<:Integer}, d=R)
    return MappingFromComponents([Polynomial(_monomial_coeffs(n), d) for n in p])
end


"""
    affinefunction(I₁::Interval, I₂::Interval)

Affine function ``f: I_1 \\to I_2`` with ``f(I_1) = I_2``.
"""
function affinefunction(I1::Interval, I2::Interval)
    a, b = endpoints(I1)
    c, d = endpoints(I2)
    s = (d - c) / (b - a)
    return Polynomial(c - a * s, s)
end


"""
    AffineMapping(A, b)

Creates the affine map `A*x + b`.
"""
struct AffineMapping{DT,CT,D} <: AbstractMapping{DT,CT,D}

    A
    b

    function AffineMapping(a::Real, b::Real, d=R)
        @assert dimension(d) == 1
        return new{InR,InR,d}(a, b)
    end

    function AffineMapping(A::AbstractMatrix, b::AbstractVector, d=nothing)
        n = size(A, 2)
        m = size(A, 1)
        dd = !isnothing(d) ? d : R^n
        @assert length(b) == m
        @assert dimension(dd) == m
        return new{InRⁿ{n},InRⁿ{m},d}(A, b)
    end
end

degree(::AffineMapping) = 1

valueat(m::AffineMapping{DT,CT}, x::DT) where {DT,CT} = CT(m.A * x + m.b)

function derivativeat(m::AffineMapping{DT}, ::DT, n::Integer=1) where {DT}
    n == 1 && return m.A
    return zero(derivativetype(m, n))
end

_derivative(m::AffineMapping{DT,CT,D}) where {DT,CT,D} = ConstantMapping(m.A, DT, D)


# -------------------------------------------------------------------------------------------------
# Ad hoc mapping
# -------------------------------------------------------------------------------------------------

struct AdHocMapping{DT,CT,D} <: AbstractMapping{DT,CT,D}
    m::Function
end
AdHocMapping(m::Function, dt=Real, ct=Float64, d=R) = AdHocMapping{dt,ct,d}(m)

valueat(m::AdHocMapping{DT}, x::DT) where {DT} = m.m(x)

makefunction(f::Function, d::DomainSets.Rectangle) = AdHocMapping(f, InR2, Float64, d)
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
Base.:(+)(m1::AbstractMapping{DT}, m2::AbstractMapping{DT}) where {DT} =
    m1 == m2 ? 2.0 * m1 : Sum(m1, m2)

# + polynomial special
_addscaled(f::Polynomial, a::Real, ::One{InR,D}) where {D} = Polynomial(a, d=D) + f
_addscaled(f::Polynomial, a::Real, ::Identity{InR,D}) where {D} = Polynomial(0, a, d=D) + f
Base.:(+)(f::Polynomial, ::Identity{InR,D}) where {D} = f + Polynomial(0, 1, d=D)
Base.:(+)(::Identity{InR,D}, f::Polynomial) where {D} = f + Polynomial(0, 1, d=D)

# + one special
Base.:(+)(f::MappingFromR{InR}, ::One{InR,D}) where {D} = f + Polynomial(1, d=D)
Base.:(+)(::One{InR,D}, f::MappingFromR{InR}) where {D} = f + Polynomial(1, d=D)

# + one/identity special
_addscaled(a1::Real, ::One{InR,D1}, a2::Real, ::Identity{InR,D2}) where {D1,D2} =
    Polynomial([a1, a2], D1 ∩ D2)
_addscaled(a1::Real, id::Identity, a2::Real, one::One) = _addscaled(a2, one, a1, id)
_addscaled(one::One, a2::Real, id::Identity) = _addscaled(1, one, a2, id)

# -
Base.:-(m::AbstractMapping) = -1.0 * m
Base.:-(m1::AbstractMapping{DT}, m2::AbstractMapping{DT}) where {DT} = m1 + (-m2)

# *
function Base.:(*)(a::Real, m::AbstractMapping)
    if a == 1
        return m
    elseif a == 0
        return zero(m)
    else
        return ScaledMapping(a, m)
    end
end

Base.:(*)(a::Num, m::AbstractMapping) = ScaledMapping(a, m)
Base.:(*)(m::AbstractMapping, a::Real) = a * m
Base.:(*)(::One{DT}, m::AbstractMapping{DT}) where {DT} = m
Base.:(*)(m::AbstractMapping{DT}, ::One{DT}) where {DT} = m
Base.:(*)(a::Real, m::ScaledMapping) = (a * m.a) * m.m
Base.:(*)(::One{DT}, m::ScaledMapping{DT}) where {DT} = m
Base.:(*)(m::ScaledMapping{DT}, ::One{DT}) where {DT} = m
Base.:(*)(m1::ScaledMapping{DT}, m2::ScaledMapping{DT}) where {DT} = m1.a * m2.a * (m1.m * m2.m)
Base.:(*)(m1::ScaledMapping{DT}, m2::AbstractMapping{DT}) where {DT} = m1.a * (m1.m * m2)
Base.:(*)(m1::AbstractMapping{DT}, m2::ScaledMapping{DT}) where {DT} = m2.a * (m1 * m2.m)
Base.:(*)(m1::AbstractMapping{DT}, m2::AbstractMapping{DT}) where {DT} =
    ProductMapping(m1, m2)

# /
Base.:(/)(m1::AbstractMapping, m2::AbstractMapping) = QuotientMapping(m1, m2)
Base.:(/)(a::Real, m2::AbstractMapping) = QuotientMapping(Polynomial(a), m2)
Base.:(/)(m2::AbstractMapping, a::Real) = 1 / a * m2

# ∘
Base.:(∘)(m1::ScaledMapping, m2::AbstractMapping) = m1.a * (m1.m ∘ m2)

# ^
Base.:(^)(::Identity{InR,D}, n::Int) where {D} = Polynomial(_monomial_coeffs(n), D)

# Add or subtract a Number
Base.:(+)(f::Identity{Real,D}, a::Real) where {D} = a == 0 ? f : Polynomial([a, 1], D)
Base.:(+)(f::FunctionToR, a::Real) = f + a * One(f)
Base.:(+)(a::Real, f::FunctionToR) = f + a
Base.:(-)(a::Real, f::FunctionToR) = a + (-f)
Base.:(-)(f::FunctionToR, a::Real) = f + (-a)

# Scaled mapping special
Base.:(+)(m1::ScaledMapping{DT,CT}, ::Zero{DT,CT}) where {DT,CT} = m1

# Two scaled functions
_addscaled(a1::Real, m1::AbstractMapping, a2::Real, m2::AbstractMapping) =
    m1 == m2 ? (a1 + a2) * m1 : Sum(a1 * m1, a2 * m2)

Base.:(+)(m1::ScaledMapping{DT,CT}, m2::ScaledMapping{DT,CT}) where {DT,CT} =
    _addscaled(m1.a, m1.m, m2.a, m2.m)

# Scaled and other function
_addscaled(m1::AbstractMapping, a::Real, m2::AbstractMapping) =
    m1 == m2 ? (1 + a) * m1 : Sum(m1, a * m2)

Base.:(+)(f1::AbstractMapping{DT,CT}, f2::ScaledMapping{DT,CT}) where {DT,CT} =
    _addscaled(f1, f2.a, f2.m)
Base.:(+)(f1::ScaledMapping{InR,InR}, f2::One{InR}) = _addscaled(f2, f1.a, f1.m)
Base.:(+)(m1::ScaledMapping{DT,CT}, m2::AbstractMapping{DT,CT}) where {DT,CT} = m2 + m1

# Enable dot product
LinearAlgebra.dot(a::Real, m::AbstractMapping) = a * m
LinearAlgebra.dot(m::AbstractMapping, a::Real) = a * m
LinearAlgebra.dot(a::Num, m::AbstractMapping) = a * m
LinearAlgebra.dot(m::AbstractMapping, a::Num) = a * m
