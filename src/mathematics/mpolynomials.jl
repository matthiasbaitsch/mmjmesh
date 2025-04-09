# -------------------------------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------------------------------

"""
    _factorialpower(m, n)

Return the product ``m (m - 1) ... (m - n)
"""
_factorialpower(m::Integer, n::Integer) = prod([m - i for i = 0:n-1])

function _monomialsat(exponents::SMatrix{N,NT}, x::InRⁿ{N}) where {N,NT}
    c = Vector{eltype(x)}(undef, NT)
    for i = 1:NT
        v = 1
        for j = 1:N
            v *= Base.FastMath.pow_fast(x[j], exponents[j, i])
        end
        c[i] = v
    end
    return SizedVector{NT}(c)
end

function _monomialsderivativeat(
    exponents::StaticMatrix{N,NT}, x::InRⁿ{N}, ns::IntegerVec
) where {N,NT}
    @assert length(ns) == N
    c = MVector{NT,Float64}(undef)
    for i = 1:NT
        v = 1
        for j = 1:N
            e = exponents[j, i]
            d = ns[j]
            if d < 0
                v *= 1 / _factorialpower(e - d, -d) * Base.FastMath.pow_fast(x[j], e - d)
            elseif d <= e
                v *= _factorialpower(e, d) * Base.FastMath.pow_fast(x[j], e - d)
            else
                v = 0
            end
        end
        c[i] = v
    end
    return c
end

function _monomialsderivative(
    exponents::StaticMatrix{N,NT}, ns::IntegerVec
) where {N,NT}
    factors = MVector{NT,Float64}(undef)
    nexponents = MMatrix{N,NT,Int}(undef)

    for i = 1:NT
        v = 1.0
        for j = 1:N
            e = exponents[j, i]
            d = ns[j]
            if d < 0
                v *= 1 / _factorialpower(e - d, -d)
                nexponents[j, i] = e - d
            elseif d <= e
                v *= _factorialpower(e, d)
                nexponents[j, i] = e - d
            else
                v = 0
                nexponents[j, i] = 0
            end
        end
        factors[i] = v
    end

    return factors, nexponents
end


# Construction helpers

_dd(exponents) = R^size(exponents, 1)

function _lt(a::AbstractArray, b::AbstractArray)
    sa = sum(a)
    sb = sum(b)
    sa != sb && return sa < sb
    for (ea, eb) = zip(a, b)
        ea != eb && return ea < eb
    end
    return false
end

# Sort handling
_arrange(c::AbstractVector, p::AbstractVector) = c[p]
_arrange(c::AbstractMatrix, p::AbstractVector) = c[:, p]
_arrange(c::AbstractArray{T,3}, p::AbstractVector) where T = c[:, :, p]

# Duplicates handling
function _squeeze(duplicates, coefficients::AbstractVector)
    idx = 0
    nd = sum(duplicates)
    nt = length(coefficients)
    sc = similar(coefficients, eltype(coefficients), nt - nd)
    for i = eachindex(coefficients)
        if duplicates[i]
            sc[idx] += coefficients[i]
        else
            idx += 1
            sc[idx] = coefficients[i]
        end
    end
    return sc
end

function _squeeze(duplicates, coefficients::AbstractMatrix)
    idx = 0
    nd = sum(duplicates)
    nt = size(coefficients, 2)
    sc = similar(coefficients, eltype(coefficients), (size(coefficients, 1), nt - nd))
    for i = 1:nt
        if duplicates[i]
            sc[:, idx] += coefficients[:, i]
        else
            idx += 1
            sc[:, idx] = coefficients[:, i]
        end
    end
    return sc
end

function _squeeze(duplicates, coefficients::AbstractArray)
    idx = 0
    nd = sum(duplicates)
    nt = size(coefficients, 3)
    sc = similar(
        coefficients,
        eltype(coefficients),
        (size(coefficients, 1), size(coefficients, 2), nt - nd)
    )
    for i = 1:nt
        if duplicates[i]
            sc[:, :, idx] += coefficients[:, :, i]
        else
            idx += 1
            sc[:, :, idx] = coefficients[:, :, i]
        end
    end
    return sc
end

# Zero handling
_findnonzeros(coefficients::AbstractVector) = .!isequal.(coefficients, 0)
_findnonzeros(coefficients::AbstractMatrix) = [!isequal(norm(col), 0) for col = eachcol(coefficients)]
_findnonzeros(coefficients::AbstractArray{T,3}) where {T} =
    [!isequal(norm(coefficients[:, :, i]), 0) for i = axes(coefficients, 3)]
_dropzeros(coefficients::AbstractVector, nonzeros) = coefficients[nonzeros]
_dropzeros(coefficients::AbstractMatrix, nonzeros) = coefficients[:, nonzeros]
_dropzeros(coefficients::AbstractArray, nonzeros) = coefficients[:, :, nonzeros]

# Normalize
function _normalize(exponents, coefficients)
    @assert size(exponents, 2) == size(coefficients, ndims(coefficients))
    nt = size(exponents, 2)
    perm = sort!(collect(1:nt), lt=(i, j) -> _lt(exponents[:, i], exponents[:, j]), rev=true)
    exponents = exponents[:, perm]
    coefficients = _arrange(coefficients, perm)

    # Duplicates
    duplicates = [false; [exponents[:, i] == exponents[:, i+1] for i = 1:nt-1]]
    if any(duplicates)
        exponents = exponents[:, .!duplicates]
        coefficients = _squeeze(duplicates, coefficients)
    end

    # Zeros
    nonzeros = _findnonzeros(coefficients)
    if !all(nonzeros)
        exponents = exponents[:, nonzeros]
        coefficients = _dropzeros(coefficients, nonzeros)
    end

    return exponents, coefficients
end

# Operations helpers
_appendcoefficients(a::AbstractVector, b::AbstractVector) = [a; b]
_appendcoefficients(a::AbstractMatrix, b::AbstractMatrix) = [a b]

# Evaluation helpers
_mulcm(c::StaticVector, m::StaticVector) = c ⋅ m
_mulcm(c::StaticMatrix, m::StaticVector) = c * m
_mulcm(c::StaticArray, m::StaticVector) = sum(m[i] * c[:, :, i] for i = eachindex(m))

# Misc helpers
function _id(n)
    id = MMatrix{n,n,Int}(undef)
    id .= 0
    id[diagind(id)] .= 1
    return id
end

# Helper functions to print

_abs(a) = abs(a)
_abs(a::Num) = a
_signchar(c::Real) = c isa Num || c >= 0 ? " + " : " - "

function _subscript(i::Integer)
    @assert i >= 0
    i > 9 && return join([_subscript(d) for d = reverse(digits(i))])
    return string('₀' + i)
end

function _superscript(i::Integer)
    @assert i >= 0
    i > 9 && return join([_superscript(d) for d = reverse(digits(i))])
    i == 1 && return "¹"
    i == 2 && return "²"
    i == 3 && return "³"
    return string('⁰' + i)
end

function _prettymonomial(exponents)
    a = ""
    for (i, e) = enumerate(exponents)
        a *= "x$(_subscript(i))$(_superscript(e))"
    end
    return a
end


# -------------------------------------------------------------------------------------------------
# Implementation
# -------------------------------------------------------------------------------------------------

struct MPolynomial{N,CT,D,NT,S} <: AbstractMapping{InRⁿ{N},CT,D}
    exponents::SMatrix{N,NT,Int}
    coefficients::SArray{S}
end

function MPolynomial(exponents::IntegerMat, coefficients::AbstractVector, D=_dd(exponents))
    exponents, coefficients = _normalize(exponents, coefficients)
    N = size(exponents, 1)
    NT = length(coefficients)
    S = Tuple{NT}
    return MPolynomial{N,InR,D,NT,S}(exponents, coefficients)
end

function MPolynomial(exponents::IntegerMat, coefficients::AbstractMatrix, D=_dd(exponents))
    exponents, coefficients = _normalize(exponents, coefficients)
    N = size(exponents, 1)
    M = size(coefficients, 1)
    NT = size(coefficients, 2)
    S = Tuple{M,NT}
    return MPolynomial{N,InRⁿ{M},D,NT,S}(exponents, coefficients)
end

function MPolynomial(exponents::IntegerMat, coefficients::AbstractArray{T,3}, D=_dd(exponents)) where {T}
    exponents, coefficients = _normalize(exponents, coefficients)
    N = size(exponents, 1)
    P = size(coefficients, 1)
    Q = size(coefficients, 2)
    NT = size(coefficients, 3)
    S = Tuple{P,Q,NT}
    return MPolynomial{N,InRᵐˣⁿ{P,Q},D,NT,S}(exponents, coefficients)
end


## Typedefs

const PolynomialRnToR{N,D,NT} = MPolynomial{N,InR,D,NT,Tuple{NT}}
const PolynomialRnToRm{N,M,D,NT} = MPolynomial{N,InRⁿ{M},D,NT,Tuple{M,NT}}
const PolynomialRnToRpxq{N,P,Q,D,NT} = MPolynomial{N,InRᵐˣⁿ{P,Q},D,NT,Tuple{P,Q,NT}}


## Properties

exponents(p::MPolynomial) = p.exponents
exponents(p::MPolynomial, idx::Integer) = p.exponents[:, idx]
coefficients(p::MPolynomial) = p.coefficients
coefficient(p::PolynomialRnToR, idx::Integer) = p.coefficients[idx]
coefficient(p::PolynomialRnToRm, idx::Integer) = p.coefficients[:, idx]
coefficient(p::PolynomialRnToRpxq, idx::Integer) = p.coefficients[:, :, idx]
nterms(::Type{<:MPolynomial{N,CT,D,NT}}) where {N,CT,D,NT} = NT
nterms(p::MPolynomial) = nterms(typeof(p))
degree(p::MPolynomial) = sum(exponents(p, 1))
degree(p::MPolynomial, c::Integer) = maximum(exponents(p)[c, :])
degrees(p::MPolynomial) = vec(maximum(exponents(p), dims=2))
Base.getindex(p::PolynomialRnToRm{N,M,D}, i::Integer) where {N,M,D} =
    MPolynomial(exponents(p), coefficients(p)[i, :], D)
Base.getindex(p::PolynomialRnToRpxq{N,P,Q,D}, i::Integer, j::Integer) where {N,P,Q,D} =
    MPolynomial(exponents(p), coefficients(p)[i, j, :], D)


## Evaluation

# Value
valueat(p::MPolynomial{N}, x::InRⁿ{N}) where {N} =
    _mulcm(coefficients(p), _monomialsat(exponents(p), x))

# Partial derivative
derivativeat(p::MPolynomial{N}, x::InRⁿ{N}, ns::IntegerVec) where {N} =
    _mulcm(coefficients(p), _monomialsderivativeat(exponents(p), x, ns))

# To R: Vector of partial derivatives (e.g. gradient)
function derivativeat(
    p::PolynomialRnToR{N}, x::InRⁿ{N}, ns::StaticMatrix{ND,N,<:Integer}
) where {N,ND}
    d = Vector{Base.promote_op(*, eltype(coefficients(p)), Float64)}(undef, ND)
    for i = 1:ND
        d[i] = derivativeat(p, x, ns[i, :])
    end
    return d
end

# To R: Matrix of partial derivatives (e.g. Hessian)
function derivativeat(
    p::PolynomialRnToR{N}, x::InRⁿ{N}, ns::StaticArray{Tuple{ND,MD,N},<:Integer}
) where {N,ND,MD}
    d = MArray{Tuple{ND,MD},Float64}(undef)
    for i = 1:ND, j = 1:MD
        d[i, j] = derivativeat(p, x, ns[i, j, :])
    end
    return d
end

# To R: Derivative of specific order (shorthand notation)
function derivativeat(
    p::MPolynomial{N}, x::InRⁿ{N}, n::Integer
) where {N}
    n == 1 && return derivativeat(p, x, _id(N))
    @notimplemented
end

# To Rm: Matrix of partial derivatives (e.g. Jacobian)  
function derivativeat(
    p::PolynomialRnToRm{N,M}, x::InRⁿ{N}, ns::StaticMatrix{ND,N,<:Integer}
) where {N,M,ND}
    d = MMatrix{M,ND,Float64}(undef)
    for i = 1:ND
        d[:, i] = derivativeat(p, x, ns[i, :])
    end
    return d
end


## Derivatives

# To R: Partial derivative
function derivative(p::PolynomialRnToR, ns::IntegerVec)
    a, e = _monomialsderivative(exponents(p), ns)
    return MPolynomial(e, a .* coefficients(p), domain(p))
end

# To R: Vector of partial derivatives (e.g. gradient)
function derivative(
    p::PolynomialRnToR{N,D,NT}, ns::StaticMatrix{ND,N,<:Integer}
) where {N,D,NT,ND}
    idx = 1
    e = MMatrix{N,NT * ND,Int}(undef)
    c = Matrix{eltype(coefficients(p))}(undef, ND, NT * ND)
    c .= 0
    for i = 1:ND
        ai, ei = _monomialsderivative(exponents(p), ns[i, :])
        e[:, idx:(idx+NT-1)] = ei
        c[i, idx:(idx+NT-1)] = ai .* coefficients(p)
        idx += NT
    end
    return MPolynomial(e, c, domain(p))
end

# To R: Matrix of partial derivatives (e.g. Hessian)
function derivative(
    p::PolynomialRnToR{N,D,NT}, ns::StaticArray{Tuple{ND,MD,N},<:Integer}
) where {N,D,NT,ND,MD}
    idx = 1
    e = MMatrix{N,NT * ND * MD,Int}(undef)
    c = MArray{Tuple{ND,MD,NT * ND * MD},eltype(coefficients(p))}(undef)
    c .= 0
    for i = 1:ND, j = 1:MD
        ai, ei = _monomialsderivative(exponents(p), ns[i, j, :])
        e[:, idx:(idx+NT-1)] = ei
        c[i, j, idx:(idx+NT-1)] = ai .* coefficients(p)
        idx += NT
    end
    return MPolynomial(e, c, domain(p))
end

# Derivative of specific order
# TODO: Move to mappings.jl
function derivative(
    p::MPolynomial{N}, n::Integer
) where {N}
    n == 1 && return derivative(p, _id(N))

    if n == 2
        ns = MArray{Tuple{N,N,N},Int}(zeros(N, N, N))
        for i = 1:N
            ns[i, :, i] .+= 1
        end
        for i = 1:N
            ns[:, i, i] .+= 1
        end
        return derivative(p, ns)
    end
    @notimplemented
end

# To Rm: Partial derivative as function Rn -> Rm
function derivative(p::PolynomialRnToRm, ns::IntegerVec)
    a, e = _monomialsderivative(exponents(p), ns)
    return MPolynomial(e, (a .* coefficients(p)')', domain(p))
end

# To Rm: Matrix of partial derivatives (e.g. Jacobian)
function derivative(
    p::PolynomialRnToRm{N,M,D,NT}, ns::StaticMatrix{ND,N,<:Integer}
) where {N,M,D,NT,ND}
    idx = 1
    e = MMatrix{N,NT * ND,Int}(undef)
    c = MArray{Tuple{M,ND,NT * ND},eltype(coefficients(p))}(undef)
    c .= 0
    for i = 1:ND
        ai, ei = _monomialsderivative(exponents(p), ns[i, :])
        e[:, idx:(idx+NT-1)] = ei
        c[:, i, idx:(idx+NT-1)] = (ai .* coefficients(p)')'
        idx += NT
    end
    return MPolynomial(e, c, domain(p))
end

antiderivative(p::MPolynomial, ns::IntegerVec) = derivative(p, -ns)


## Operations

Base.transpose(p::PolynomialRnToRpxq) =
    MPolynomial(exponents(p), permutedims(coefficients(p), (2, 1, 3)), domain(p))

function Base.:(+)(p1::MPolynomial{N,CT,D}, p2::MPolynomial{N,CT,D}) where {N,CT,D}
    e = [exponents(p1) exponents(p2)]
    c = _appendcoefficients(coefficients(p1), coefficients(p2))
    return MPolynomial(e, c, D)
end

function Base.:(*)(p1::PolynomialRnToR{N,D}, p2::PolynomialRnToR{N,D}) where {N,D}
    idx = 1
    nt1 = nterms(p1)
    nt2 = nterms(p2)
    t1 = eltype(coefficients(p1))
    t2 = eltype(coefficients(p2))
    tt = promote_type(t1, t2)
    e = MMatrix{N,nt1 * nt2,Int}(undef)
    c = Vector{tt}(undef, nt1 * nt2)
    for i = 1:nt1, j = 1:nt2
        e[:, idx] = exponents(p1, i) + exponents(p2, j)
        c[idx] = coefficient(p1, i) * coefficient(p2, j)
        idx += 1
    end
    return MPolynomial(e, c, D)
end

Base.:(*)(a::Num, p::MPolynomial) = MPolynomial(exponents(p), a * coefficients(p), domain(p))
Base.:(*)(a::Real, p::MPolynomial) = MPolynomial(exponents(p), a * coefficients(p), domain(p))


## Algebraic operations

LinearAlgebra.dot(a::AbstractVector, p::PolynomialRnToRm) =
    MPolynomial(exponents(p), vec(coefficients(p)' * a))

Base.:(*)(a::AbstractMatrix, p::PolynomialRnToRm) =
    MPolynomial(exponents(p), a * coefficients(p), domain(p))

function Base.:(*)(p::PolynomialRnToRpxq{N,P,Q,D,NT}, x::RealVec) where {N,P,Q,D,NT}
    c = coefficients(p)
    @assert length(x) == Q
    return MPolynomial(exponents(p), [c[i, :, j] ⋅ x for i = 1:P, j = 1:NT], D)
end

function Base.:(*)(p::PolynomialRnToRpxq{N,P,Q,D,NT}, A::RealMat) where {N,P,Q,D,NT}
    c = coefficients(p)
    @assert size(A, 1) == Q
    return MPolynomial(exponents(p), stack([c[:, :, j] * A for j = 1:NT], dims=3), D)
end

function Base.:(*)(A::RealMat, p::PolynomialRnToRpxq{N,P,Q,D,NT}) where {N,P,Q,D,NT}
    c = coefficients(p)
    @assert size(A, 2) == P
    return MPolynomial(exponents(p), stack([A * c[:, :, j] for j = 1:NT], dims=3), D)
end


## Manipulation

Base.rationalize(p::MPolynomial) =
    MPolynomial(exponents(p), rationalize.(coefficients(p)), domain(p))
MMJMesh.MMJBase.integerize(p::MPolynomial) =
    MPolynomial(exponents(p), integerize.(coefficients(p)), domain(p))


## Show

function Base.show(io::IO, p::PolynomialRnToR)
    for i = 1:nterms(p)
        c = coefficient(p, i)
        i > 1 && print(io, _signchar(c))
        print(io, _abs(c), _prettymonomial(exponents(p, i)))
    end
end

function Base.show(io::IO, p::PolynomialRnToRm)
    for i = 1:nterms(p)
        i > 1 && print(io, " + ")
        print(io, _prettymonomial(exponents(p, i)), " ⋅ ", coefficient(p, i))
    end
end

function Base.show(io::IO, p::PolynomialRnToRpxq)
    for i = 1:nterms(p)
        i > 1 && print(io, " + ")
        print(io, _prettymonomial(exponents(p, i)), " ⋅ ", coefficient(p, i))
    end
end


## Compare

Base.:(==)(p1::MPolynomial, p2::MPolynomial) =
    (p1.exponents == p2.exponents && p1.coefficients == p2.coefficients)

Base.isapprox(
    p1::MPolynomial, p2::MPolynomial; atol::Real=0, rtol::Real=atol > 0 ? 0 : √eps(Float64)
) =
    (
        p1.exponents == p2.exponents &&
        isapprox(coefficients(p1), coefficients(p2); atol=atol, rtol=rtol)
    )


## Monomials

"""
	mmonomials(n::Integer, p::Integer, dom=R^n, predicate=(...) -> true; type=Float64)

Generate multivariate momonials of `n` components up to degree `p`. Optionally, a domain, a 
predicate and a coefficient type can be specified
"""
function mmonomials(n::Integer, p::Integer, dom=R^n, predicate=(ps...) -> true; type=Float64)
    if n == 2
        e = tomatrix(
            sort!(
                [[p1, p2] for p1 in 0:p, p2 in 0:p if predicate(p1, p2)],
                lt=_lt, rev=true
            )
        )
        c = zeros(type, size(e, 2), size(e, 2))
        c[diagind(c)] .= 1
        return MPolynomial(e, c, dom)
    else
        error("Not implemented yet")
    end
end
