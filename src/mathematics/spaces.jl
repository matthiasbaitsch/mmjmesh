"""
	FunctionSpace

General concept of a space containing functions defined on `D` which map type 
domain type `DT` to codomain type `CT`.
"""
abstract type FunctionSpace{DT,CT,D} end


"""
    basis(s::FunctionSpace{DT,CT}, D) -> Vector{<:AbstractMapping{DT,CT,D}}

Generate a basis for the specified function space defined on the domain `D`.
"""
basis(::Type{FunctionSpace}) = @abstractmethod
basis(s::FunctionSpace) = basis(typeof(s))


"""
    dimension
"""
dimension(::Type{FunctionSpace}) = @abstractmethod
dimension(s::FunctionSpace) = dimension(typeof(s))


""" Type of elements in the domain of the mapping. """
domaintype(::Type{<:FunctionSpace{DT}}) where {DT} = DT
domaintype(m::FunctionSpace) = domaintype(typeof(m))


""" Type of elements in the codomain of the mapping. """
codomaintype(::Type{<:FunctionSpace{DT,CT}}) where {DT,CT} = CT
codomaintype(m::FunctionSpace) = codomaintype(typeof(m))


""" Domain on which the functions are defined. """
domain(::Type{<:FunctionSpace{DT,CT,D}}) where {DT,CT,D} = D
domain(m::FunctionSpace) = domain(typeof(m))

Base.in(x, ::Type{<:FunctionSpace}) = false
Base.in(x, s::FunctionSpace) = in(x, typeof(s))


"""
	PolynomialSpace{N,K,D}

FunctionSpace of possibly multivarite polynomials of `N` variables an maximum exponent `K`.
"""
abstract type PolynomialSpace{N,K,D} <: FunctionSpace{InRⁿ{N},InR,D} end

dimension(::Type{<:PolynomialSpace{1,K}}) where {K} = K + 1

Base.in(::AbstractArray{<:Integer}, ::Type{<:PolynomialSpace}) = @abstractmethod
Base.in(f::MappingFromR, s::Type{<:PolynomialSpace{1}}) = in([degree(f)], s)
Base.in(f::MPolynomial, s::Type{<:PolynomialSpace}) = all(in.(eachcol(f.p.exponents), s))

basis(::Type{<:PolynomialSpace{1,K,D}}) where {K,D} = monomials(0:K, D)
basis(s::Type{<:PolynomialSpace{N,K,D}}) where {N,K,D} =
    mmonomials(N, K, D, (k...) -> [k...] ∈ s, type=Int)
basis(s::PolynomialSpace) = basis(typeof(s))

"""
    _dimension(s::Type{<:PolynomialSpace})

Fallback implemenentation to determine the dimension of space s. Useful for testing.
"""
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

FunctionSpace of `n`-variate polynomials where the largest sum of exponents is not larger than `k`.
"""
struct P{N,K,D} <: PolynomialSpace{N,K,D} end
dimension(::Type{<:P{2,K}}) where {K} = Int((K + 1) * (K + 2) / 2)
Base.in(k::AbstractArray{<:Integer}, ::Type{<:P{N,K}}) where {N,K} = (length(k) == N && sum(k) <= K)


"""
	Q{N,K,D}

FunctionSpace of `N`-variate polynomials where the largest exponent is smaller equal `K` defined on `D`.
"""
struct Q{N,K,D} <: PolynomialSpace{N,K,D} end
dimension(::Type{<:Q{N,K}}) where {N,K} = (K + 1)^N
Base.in(k::AbstractArray{<:Integer}, ::Type{<:Q{N,K}}) where {N,K} = (length(k) == N && maximum(k) <= K)


"""
	S{N, K}

FunctionSpace of serendipity polynomials. Only defined for n=2 and k=2,3.
"""
struct S{N,K,D} <: PolynomialSpace{N,K,D} end
dimension(::Type{<:S{2,K}}) where {K} = 4 * K
Base.in(k::AbstractArray{<:Integer}, ::Type{<:S{2,K}}) where {K} =
    (K <= 2 && k ∈ Q{2,K} && prod(k) <= K)


"""
    Q23R{K}

Reduced polynomial space `Q{2,3,K}` for nonconforming Kirchhoff rectangle.
"""
struct Q23R{K} <: PolynomialSpace{2,3,K} end
dimension(::Type{<:Q23R}) = 12
Base.in(k::AbstractArray{<:Integer}, ::Type{<:Q23R}) =
    maximum(k) <= 3 && sum(k) <= 4 && prod(k) < 4
