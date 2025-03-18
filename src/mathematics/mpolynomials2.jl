# -------------------------------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------------------------------

"""
    _factorialpower(m, n)

Return the product ``m (m - 1) ... (m - n)
"""
_factorialpower(m::Integer, n::Integer) = prod([m - i for i = 0:n-1])

function _monomialsat(exponents::SMatrix{N,M}, x::InRⁿ{N}) where {N,M}
    c = MVector{M,Float64}(undef)
    for i = 1:M
        v = 1
        for j = 1:N
            v *= Base.FastMath.pow_fast(x[j], exponents[j, i])
        end
        c[i] = v
    end
    return c
end

function _monomialsderivativeat(
    exponents::StaticMatrix{N,M}, x::InRⁿ{N}, ns::IntegerVec
) where {N,M}
    @assert N == length(ns)
    c = MVector{M,Float64}(undef)
    for i = 1:M
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

function _check(exponents, coefficients)
    @assert size(exponents, 2) == size(coefficients, ndims(coefficients))
    return true
end

function _lt(a::AbstractArray, b::AbstractArray)
    sa = sum(a)
    sb = sum(b)
    sa != sb && return sa < sb
    for (ea, eb) = zip(a, b)
        ea != eb && return ea < eb
    end
    return false
end

_arrange(::Type{Tuple{N}}, c::Array, p::AbstractVector) where {N} = c[p]
_arrange(::Type{Tuple{N,M}}, c::Array, p::AbstractVector) where {N,M} = c[:, p]


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

_dd(exponents) = InRⁿ{size(exponents, 1)}


# -------------------------------------------------------------------------------------------------
# Implementation
# -------------------------------------------------------------------------------------------------

struct MPolynomial2{NT,S,N,CT,D} <: AbstractMapping{InRⁿ{N},CT,D}

    exponents::SMatrix{N,NT,Int}
    coefficients::SArray{S}

    function MPolynomial2{NT,S,N,CT,D}(exponents, coefficients) where {NT,S,N,CT,D}
        p = sort!(collect(1:NT), lt=(i, j) -> _lt(exponents[:, i], exponents[:, j]), rev=true)
        return new{NT,S,N,CT,D}(exponents[:, p], _arrange(S, coefficients, p))
    end
end

function MPolynomial2(exponents::Matrix{Int}, coefficients::AbstractVector, D=_dd(exponents))
    @assert _check(exponents, coefficients)
    N = size(exponents, 1)
    NT = length(coefficients)
    S = Tuple{NT}
    return MPolynomial2{NT,S,N,InR,D}(exponents, coefficients)
end

function MPolynomial2(exponents::Matrix{Int}, coefficients::AbstractMatrix, D=_dd(exponents))
    @assert _check(exponents, coefficients)
    N = size(exponents, 1)
    M = size(coefficients, 1)
    NT = size(coefficients, 2)
    S = Tuple{M,NT}
    return MPolynomial2{NT,S,N,InRⁿ{M},D}(exponents, coefficients)
end


# Typedefs

const PolynomialRnToR{NT,N,D} = MPolynomial2{NT,Tuple{NT},N,InR,D}
const PolynomialRnToRm{NT,N,M,D} = MPolynomial2{NT,Tuple{M,NT},N,InRⁿ{M},D}


# Properties

exponents(p::MPolynomial2) = p.exponents
exponents(p::MPolynomial2, idx::Integer) = p.exponents[:, idx]
coefficients(p::MPolynomial2) = p.coefficients
coefficient(p::PolynomialRnToR, idx::Integer) = p.coefficients[idx]
coefficient(p::PolynomialRnToRm, idx::Integer) = p.coefficients[:, idx]
nterms(::Type{<:MPolynomial2{NT}}) where {NT} = NT
nterms(p::MPolynomial2) = nterms(typeof(p))


# Evaluation

_mulop(::PolynomialRnToR) = dot
_mulop(::PolynomialRnToRm) = *

valueat(p::MPolynomial2{NT,S,N}, x::InRⁿ{N}) where {NT,S,N} =
    _mulop(p)(coefficients(p), _monomialsat(exponents(p), x))

derivativeat(p::MPolynomial2{NT,S,N}, x::InRⁿ{N}, ns::IntegerVec) where {NT,S,N} =
    _mulop(p)(coefficients(p), _monomialsderivativeat(exponents(p), x, ns))

function derivativeat(
    p::MPolynomial2{NT,S,N}, x::InRⁿ{N}, ns::StaticMatrix{M,N,<:Integer}
) where {NT,S,N,M}
    d = MArray{S,Float64}(undef)
    println(S)
    # println(size(S))
    println(M)
    d
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
