"""
	FunctionSpace

General concept of a space containing functions which map 
domain type `DT` to codomain type `CT`.
"""
abstract type FunctionSpace{DT,CT} end


"""
    basis(s::FunctionSpace{DT,CT}, d) -> Vector{<:AbstractMapping{DT,CT}}

Generate a basis for function space `s` defined on `d`.
"""
basis(::Type{FunctionSpace}, d) = @abstractmethod
basis(s::FunctionSpace, d) = basis(typeof(s), d)


"""
    dimension
"""
MMJMesh.dimension(::Type{FunctionSpace}) = @abstractmethod
MMJMesh.dimension(s::FunctionSpace) = dimension(typeof(s))


""" Type of elements in the domain of the mapping. """
domaintype(::Type{<:FunctionSpace{DT}}) where {DT} = DT
domaintype(m::FunctionSpace) = domaintype(typeof(m))


""" Type of elements in the codomain of the mapping. """
codomaintype(::Type{<:FunctionSpace{DT,CT}}) where {DT,CT} = CT
codomaintype(m::FunctionSpace) = codomaintype(typeof(m))


Base.in(x, ::Type{<:FunctionSpace}) = false
Base.in(x, s::FunctionSpace) = in(x, typeof(s))


"""
	PolynomialSpace{N,K}

FunctionSpace of possibly multivarite polynomials of `N` variables and maximum exponent `K`.
"""
abstract type PolynomialSpace{N,K} <: FunctionSpace{InRⁿ{N},InR} end

MMJMesh.dimension(::Type{<:PolynomialSpace{1,K}}) where {K} = K + 1

Base.in(::AbstractArray{<:Integer}, ::Type{<:PolynomialSpace}) = @abstractmethod
Base.in(f::MappingFromR, s::Type{<:PolynomialSpace{1}}) = in([degree(f)], s)
Base.in(f::MPolynomial2, s::Type{<:PolynomialSpace}) = all(in.(eachcol(f.p.exponents), s))

basis(::Type{<:PolynomialSpace{1,K}}, d) where {K} = monomials(0:K, d)
basis(s::Type{<:PolynomialSpace{N,K}}, d) where {N,K} =
    mmonomials2(N, K, d, (k...) -> [k...] ∈ s, type=Int)
basis(s::PolynomialSpace, d) = basis(typeof(s), d)

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

FunctionSpace of `N`-variate polynomials where the largest sum of exponents is not larger than `K`.
"""
struct P{N,K} <: PolynomialSpace{N,K} end
MMJMesh.dimension(::Type{<:P{2,K}}) where {K} = Int((K + 1) * (K + 2) / 2)
Base.in(k::AbstractArray{<:Integer}, ::Type{<:P{N,K}}) where {N,K} = (length(k) == N && sum(k) <= K)


"""
	Q{N,K}

FunctionSpace of `N`-variate polynomials where the largest exponent is smaller equal `K`.
"""
struct Q{N,K} <: PolynomialSpace{N,K} end
MMJMesh.dimension(::Type{<:Q{N,K}}) where {N,K} = (K + 1)^N
Base.in(k::AbstractArray{<:Integer}, ::Type{<:Q{N,K}}) where {N,K} = (length(k) == N && maximum(k) <= K)


"""
	S{N, K}

FunctionSpace of serendipity polynomials. Only defined for n=2 and k=2,3.
"""
struct S{N,K} <: PolynomialSpace{N,K} end
MMJMesh.dimension(::Type{<:S{2,K}}) where {K} = 4 * K
Base.in(k::AbstractArray{<:Integer}, ::Type{<:S{2,K}}) where {K} =
    (K <= 2 && k ∈ Q{2,K} && prod(k) <= K)


"""
    Q23R{K}

Reduced polynomial space `Q{2,3,K}` for nonconforming Kirchhoff rectangle.
"""
struct Q23R <: PolynomialSpace{2,3} end
MMJMesh.dimension(::Type{Q23R}) = 12
Base.in(k::AbstractArray{<:Integer}, ::Type{Q23R}) = maximum(k) <= 3 && sum(k) <= 4 && prod(k) < 4
