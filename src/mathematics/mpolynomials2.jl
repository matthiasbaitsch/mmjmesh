# -------------------------------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------------------------------

"""
    _factorialpower(m, n)

Return the product ``m (m - 1) ... (m - n)
"""
_factorialpower(m::Integer, n::Integer) = prod([m - i for i = 0:n-1])

function _monomialsat(exponents::SMatrix{N,NT}, x::InRⁿ{N}) where {N,NT}
    c = MVector{NT,Float64}(undef)
    for i = 1:NT
        v = 1
        for j = 1:N
            v *= Base.FastMath.pow_fast(x[j], exponents[j, i])
        end
        c[i] = v
    end
    return c
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

# Helper functions to print

_signchar(c::Real) = c >= 0 ? " + " : " - "

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


# Construction helpers

_dd(exponents) = InRⁿ{size(exponents, 1)}

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

# Zero handling
_findnonzeros(coefficients::AbstractVector) = coefficients .!= 0
_findnonzeros(coefficients::AbstractMatrix) = [norm(col) != 0 for col = eachcol(coefficients)]
_dropzeros(coefficients::AbstractVector, nonzeros) = coefficients[nonzeros]
_dropzeros(coefficients::AbstractMatrix, nonzeros) = coefficients[:, nonzeros]

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


# -------------------------------------------------------------------------------------------------
# Implementation
# -------------------------------------------------------------------------------------------------

struct MPolynomial2{N,CT,D,NT,S} <: AbstractMapping{InRⁿ{N},CT,D}
    exponents::SMatrix{N,NT,Int}
    coefficients::SArray{S}
end

function MPolynomial2(exponents::IntegerMat, coefficients::AbstractVector, D=_dd(exponents))
    exponents, coefficients = _normalize(exponents, coefficients)
    N = size(exponents, 1)
    NT = length(coefficients)
    S = Tuple{NT}
    return MPolynomial2{N,InR,D,NT,S}(exponents, coefficients)
end

function MPolynomial2(exponents::IntegerMat, coefficients::AbstractMatrix, D=_dd(exponents))
    exponents, coefficients = _normalize(exponents, coefficients)
    N = size(exponents, 1)
    M = size(coefficients, 1)
    NT = size(coefficients, 2)
    S = Tuple{M,NT}
    return MPolynomial2{N,InRⁿ{M},D,NT,S}(exponents, coefficients)
end


## Typedefs

const PolynomialRnToR{N,D,NT} = MPolynomial2{N,InR,D,NT,Tuple{NT}}
const PolynomialRnToRm{N,M,D,NT} = MPolynomial2{N,InRⁿ{M},D,NT,Tuple{M,NT}}


## Properties

exponents(p::MPolynomial2) = p.exponents
exponents(p::MPolynomial2, idx::Integer) = p.exponents[:, idx]
coefficients(p::MPolynomial2) = p.coefficients
coefficient(p::PolynomialRnToR, idx::Integer) = p.coefficients[idx]
coefficient(p::PolynomialRnToRm, idx::Integer) = p.coefficients[:, idx]
nterms(::Type{<:MPolynomial2{N,CT,D,NT}}) where {N,CT,D,NT} = NT
nterms(p::MPolynomial2) = nterms(typeof(p))
Base.getindex(p::PolynomialRnToRm, i::Integer) = MPolynomial2(exponents(p), coefficients(p)[i, :])


## Evaluation

_mulop(::PolynomialRnToR) = dot
_mulop(::PolynomialRnToRm) = *

valueat(p::MPolynomial2{N}, x::InRⁿ{N}) where {N} =
    _mulop(p)(coefficients(p), _monomialsat(exponents(p), x))

derivativeat(p::MPolynomial2{N}, x::InRⁿ{N}, ns::IntegerVec) where {N} =
    _mulop(p)(coefficients(p), _monomialsderivativeat(exponents(p), x, ns))

function derivativeat(
    p::PolynomialRnToR{N}, x::InRⁿ{N}, ns::StaticMatrix{ND,N,<:Integer}
) where {N,ND}
    d = MVector{ND,Float64}(undef)
    for i = 1:ND
        d[i] = derivativeat(p, x, ns[i, :])
    end
    return d
end

function derivativeat(
    p::PolynomialRnToR{N}, x::InRⁿ{N}, ns::StaticArray{Tuple{ND,MD,N},<:Integer}
) where {N,ND,MD}
    d = MArray{Tuple{ND,MD},Float64}(undef)
    for i = 1:ND, j = 1:MD
        d[i, j] = derivativeat(p, x, ns[i, j, :])
    end
    return d
end

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

function derivative(p::PolynomialRnToR, ns::IntegerVec)
    a, e = _monomialsderivative(exponents(p), ns)
    return MPolynomial2(e, a .* coefficients(p))
end

function derivative(
    p::PolynomialRnToR{N,D,NT}, ns::StaticMatrix{ND,N,<:Integer}
) where {N,D,NT,ND}
    idx = 1
    e = MMatrix{N,NT * ND,Int}(undef)
    c = MMatrix{ND,NT * ND,eltype(coefficients(p))}(undef)
    c .= 0
    for i = 1:ND
        ai, ei = _monomialsderivative(exponents(p), ns[i, :])
        e[:, idx:(idx+NT-1)] = ei
        c[i, idx:(idx+NT-1)] = ai .* coefficients(p)
        idx += NT
    end
    return MPolynomial2(e, c)
end

function derivative(p::PolynomialRnToRm, ns::IntegerVec)
    a, e = _monomialsderivative(exponents(p), ns)
    return MPolynomial2(e, (a .* coefficients(p)')')
end


## Show

function Base.show(io::IO, p::PolynomialRnToR)
    for i = 1:nterms(p)
        c = coefficient(p, i)
        i > 1 && print(io, _signchar(c))
        print(io, abs(c), _prettymonomial(exponents(p, i)))
    end
end

function Base.show(io::IO, p::PolynomialRnToRm)
    for i = 1:nterms(p)
        i > 1 && print(io, " + ")
        print(io, _prettymonomial(exponents(p, i)), " ⋅ ", coefficient(p, i))
    end
end
