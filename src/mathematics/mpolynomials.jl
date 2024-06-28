using Symbolics
using SymbolicUtils
using SymbolicUtils: Postwalk, Chain

import MMJMesh.Mathematics.FixedPolynomials as FP

# Multivariate polynomials
struct MPolynomial{N,D} <: FunctionRnToR{N,D}
    p::FP.Polynomial
end

function MPolynomial(exponents::Matrix{Int}, coefficients::AbstractVector{T}, d=nothing) where {T}
    n = size(exponents, 1)
    if isnothing(d)
        d = R^n
    end
    return MPolynomial{n,d}(FP.Polynomial(_simplify(exponents, coefficients)...))
end

valueat(p::MPolynomial{N}, x::SVector{N}) where {N} = p.p(x)

function derivative(f::MPolynomial{N,D}, ns::AbstractArray{<:Integer}) where {N,D}
    p = f.p
    for (idx, n) ∈ enumerate(ns)
        for _ = 1:n
            p = FP.differentiate(p, idx)
        end
    end
    return MPolynomial{N,D}(p)
end

function antiderivative(f::MPolynomial{N}, ns::AbstractArray{<:Integer}) where {N}
    e = copy(FP.exponents(f.p))
    c = Vector{promote_type(FP.eltype(f.p), Float64)}(FP.coefficients(f.p))
    for i = 1:N, _ = 1:ns[i], k = 1:FP.nterms(f.p)
        e[i, k] += 1
        c[k] = (1 // e[i, k]) * c[k]
    end
    return MPolynomial(e, c)
end

function Base.:(+)(p1::MPolynomial, p2::MPolynomial)
    _, e1, c1, _, e2, c2 = _extract(p1, p2)
    return MPolynomial(hcat(e1, e2), vcat(c1, c2), domain(p1))
end

function Base.:(*)(p1::MPolynomial{N}, p2::MPolynomial{N}) where {N}
    n1, e1, c1, n2, e2, c2 = _extract(p1, p2)
    n = n1 * n2
    e = zeros(Int, N, n)
    c = ones(promote_type(FP.eltype(p1.p), FP.eltype(p2.p)), n)
    cnt = 1
    for i = 1:n1
        for j = 1:n2
            e[:, cnt] = e1[:, i] + e2[:, j]
            c[cnt] = c1[i] * c2[j]
            cnt += 1
        end
    end
    return MPolynomial(e, c, domain(p1))
end

Base.:(*)(a::Real, p::MPolynomial{N,D}) where {N,D} =
    MPolynomial(FP.exponents(p.p), a * FP.coefficients(p.p), D)

Base.show(io::IO, p::MPolynomial) = print(io, p.p)
Base.:(==)(p1::MPolynomial{N,D}, p2::MPolynomial{N,D}) where {N,D} = p1.p == p2.p

"""
    mmonomials(n::Integer, p::Integer, dom=R^n, predicate=(...) -> true)

Generate multivariate momonials of `n` components up to degree `p`. Optionally, a domain
and a predicate can be specified
"""
function mmonomials(n::Integer, p::Integer, dom=R^n, predicate=(ps...) -> true; type=Float64)
    if n == 2
        return reshape(
            [
                MPolynomial([p1; p2;;], type.([1]), dom)
                for p1 in 0:p, p2 in 0:p
                if predicate(p1, p2)
            ],
            :
        )
    else
        error("Not implemented yet")
    end
end


# -------------------------------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------------------------------

function _simplify(e::Matrix{ET}, c::Vector{CT}) where {ET,CT}
    m = Dict{AbstractArray{Int},CT}()

    for (e, c) ∈ zip(eachcol(e), c)
        if !isequal(c, 0)
            if e ∈ keys(m)
                m[e] += c
            else
                m[e] = c
            end
        end
    end
    if !isempty(m)
        return stack(keys(m)), _integerize!(collect(values(m)))
    else
        return zeros(ET, size(e, 1), 1), zeros(CT, 1)
    end
end

function _extract(p1::MPolynomial, p2::MPolynomial)
    e1 = FP.exponents(p1.p)
    c1 = FP.coefficients(p1.p)
    n1 = FP.nterms(p1.p)
    e2 = FP.exponents(p2.p)
    c2 = FP.coefficients(p2.p)
    n2 = FP.nterms(p2.p)
    return n1, e1, c1, n2, e2, c2
end

_isintegervalue(x::T) where {T<:Integer} = true
_isintegervalue(x::Rational) = denominator(x) == 1
_isintegervalue(x::T) where {T<:AbstractFloat} = x == round(x, digits=0)
_isintegervalue(x) = false

function _integerize(expression::Num)
    @variables xone, xnull
    r = @rule ~x::_isintegervalue => xone * (Int(~x) + xnull)
    expression = simplify(expression)
    expression = simplify(expression, Postwalk(Chain([r])))
    expression = simplify(substitute(expression, Dict(xone => 1, xnull => 0)))
    return expression
end

function _integerize!(c::Vector{Num})
    for i = eachindex(c)
        c[i] = _integerize(c[i])
    end
    return c
end

_integerize!(c) = c

"""
    simplifyx(expression::Num) -> Num

Simplify expression and convert numbers to integers if possible.
"""
simplifyx(expression::Num) = _integerize(expression)
simplifyx(x::Float64) = x

