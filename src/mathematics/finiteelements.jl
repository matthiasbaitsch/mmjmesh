"""
	Space

General concept of a vector space over the real numbers.
"""
abstract type Space end
basis(::Space) = @abstractmethod
dimension(s::Space) = dimension(typeof(s))
Base.in(x, ::Type{<:Space}) = false
Base.in(x, s::Space) = in(x, typeof(s))

"""
	PolynomialSpace{N,K}

Space of possibly multivarite polynomials of N variables an maximum exponent K.
"""
abstract type PolynomialSpace{N,K} <: Space end
dimension(::Type{<:PolynomialSpace{1,K}}) where {K} = K + 1
Base.in(::AbstractArray{<:Integer}, ::Type{<:PolynomialSpace}) = @abstractmethod
Base.in(f::MappingFromR, s::Type{<:PolynomialSpace{1}}) = in([degree(f)], s)
Base.in(f::MPolynomial, s::Type{<:PolynomialSpace}) = all(in.(eachcol(f.p.exponents), s))

basis(::Type{<:PolynomialSpace{1,K}}, domain) where {K} = monomials(0:K, domain)
basis(s::Type{<:PolynomialSpace{N,K}}, domain) where {N,K} =
    mmonomials(N, K, domain, (k...) -> [k...] ∈ s)
basis(s::PolynomialSpace, domain) = basis(typeof(s), domain)

function _dimension(s::Type{<:PolynomialSpace{2,K}}) where {K}
    cnt = 0
    for i = 0:K, j = 0:K
        if [i, j] ∈ s
            cnt += 1
        end
    end
    return cnt
end

function printexponents(s::Type{<:PolynomialSpace{2,K}}) where {K}
    for sum = 0:2K
        k1 = sum
        k2 = 0
        b = false
        for _ = 0:sum
            if [k1, k2] ∈ s
                print("($k1,$k2) ")
                b = true
            end
            k1 -= 1
            k2 += 1
        end
        b && println()
    end
end


"""
	P{N,K}

Space of `n`-variate polynomials where the largest sum of exponents is not larger than `k`.
"""
struct P{N,K} <: PolynomialSpace{N,K} end
dimension(::Type{P{2,K}}) where {K} = Int((K + 1) * (K + 2) / 2)
Base.in(k::AbstractArray{<:Integer}, ::Type{P{N,K}}) where {N,K} = (length(k) == N && sum(k) <= K)

"""
	Q{N,K}

Space of `n`-variate polynomials where the largest exponent is not larger than `n`.
"""
struct Q{N,K} <: PolynomialSpace{N,K} end
dimension(::Type{Q{N,K}}) where {N,K} = (K + 1)^N
Base.in(k::AbstractArray{<:Integer}, ::Type{Q{N,K}}) where {N,K} = (length(k) == N && maximum(k) <= K)

"""
	S{N, K}

Space of serendipity polynomials. Only defined for n=2 and k=2,3.
"""
struct S{N,K} <: PolynomialSpace{N,K} end
dimension(s::Type{S{2,2}}) = 8
dimension(s::Type{S{2,3}}) = 12
Base.in(k::AbstractArray{<:Integer}, ::Type{S{2,K}}) where {K} = (k ∈ Q{2,K} && prod(k) <= K)
