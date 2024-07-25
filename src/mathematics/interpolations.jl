# Helper type
struct TransposeMapping{DT,CT,D} <: AbstractMapping{DT,CT,D}
    u::AbstractMapping{DT}
    TransposeMapping(u::AbstractMapping{DT,InRⁿˣᵐ{N,M},D}) where {DT,N,M,D} =
        new{DT,InRⁿˣᵐ{M,N},D}(u)
end
valueat(u::TransposeMapping{DT}, x::DT) where {DT} = valueat(u.u, x)'
Base.:(==)(t1::TransposeMapping, t2::TransposeMapping) = (t1.u == t2.u)


"""
    Interpolation(h, C)

The interpolating mapping ``m: R^n \\to Y`` with

  ``m(x) = \\sum_{i=1}^m h_i(x) C_i``

with interpolation functions ``h_i : R^n \\to R``, coefficients 
``C_i \\in Y`` and ``Y \\in R^{k_1} \\times \\dots \\times R^{k_l}``.
"""
struct Interpolation{DT,CT,D} <: AbstractMapping{DT,CT,D}
    functions::AbstractMapping{DT,<:InRⁿ,D}
    coefficients::AbstractVecOrMat{<:Real}

    Interpolation(
        h::AbstractMapping{DT,<:AbstractVector,D}, C::AbstractVector{<:Real}
    ) where {DT,D} = new{DT,InR,D}(h, C)

    Interpolation(
        h::AbstractMapping{DT,<:AbstractVector,D}, C::AbstractMatrix{<:Real}
    ) where {DT,D} =
        new{DT,InRⁿ{size(C, 1)},D}(h, C)
end

valueat(m::Interpolation{DT}, x::DT) where {DT} = _mul(m.coefficients, m.functions(x))
Base.show(io::IO, m::Interpolation) =
    print(io, "Interpolation(", m.functions, ", ", m.coefficients, ")")
Base.:(==)(p1::Interpolation, p2::Interpolation) =
    (p1.functions == p2.functions) && (p1.coefficients == p2.coefficients)

# Derivatives generic case
derivative(m::Interpolation, n::Integer=1) =
    Interpolation(derivative(m.functions, n), m.coefficients)

derivativeat(m::Interpolation{DT}, x::DT, n::Integer=1) where {DT} =
    _mul(m.coefficients, derivativeat(m.functions, x, n))

# Derivatives for function Rn to R
function derivative(f::Interpolation{InRⁿ{N},InR}, n::Integer=1) where {N}
    @assert n == 1
    nc = length(f.coefficients)
    jh = jacobian(f.functions)
    colsjh = [MappingFromComponents([jh.components[i].components[j] for i = 1:nc]...) for j = 1:N]
    return MappingFromComponents([Interpolation(h, f.coefficients) for h = colsjh]...)
end

function derivativeat(m::Interpolation{InRⁿ{N},InR}, x::InRⁿ{N}, n::Int=1) where {N}
    return derivativeat(m.functions, x, n)' * m.coefficients
end

# Derivatives for mapping Rn to Rm
function derivative(u::Interpolation{InRⁿ{N},InRⁿ{M}}, n::Integer=1) where {N,M}
    @assert n == 1
    nc = size(u.coefficients, 2)
    jh = jacobian(u.functions)
    colsjh = [MappingFromComponents([jh.components[i].components[j] for i = 1:nc]...) for j = 1:N]
    return TransposeMapping(
        MappingFromComponents([Interpolation(h, u.coefficients) for h = colsjh]...)
    )
end

function derivativeat(u::Interpolation{InRⁿ{N},InRⁿ{M}}, x::InRⁿ{N}, n::Int=1) where {N,M}
    @assert n == 1
    ju = MMatrix{M,N,Float64}(undef)
    jh = jacobianat(u.functions, x)
    for j = 1:N
        ju[:, j] = u.coefficients * vec(jh[:, j])
    end
    return ju
end


# -------------------------------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------------------------------

_mul(C::AbstractVector, h::AbstractVector) = C ⋅ h
_mul(C::AbstractMatrix, h::AbstractVector) = C * h
_ncoefficients(C::AbstractVector) = length(C)
_ncoefficients(C::AbstractMatrix) = size(C, 2)

