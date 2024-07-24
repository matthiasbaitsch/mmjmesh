"""
    Interpolation(h, C)

The interpolating mapping ``m: R^n -> Y`` with

  ``m(x) = \\sum_{i=1}^m h_i(x) C_i``

with interpolation functions ``h_i : R^n -> R``, coefficients 
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
derivativeat(m::Interpolation{DT}, x::DT, n::Int=1) where {DT} =
    _mul(m.coefficients, derivativeat(m.functions, x, n))
derivative(m::Interpolation, n::Int=1) =
    Interpolation(derivative(m.functions, n), m.coefficients)
Base.show(io::IO, m::Interpolation) =
    print(io, "Interpolation(", m.functions, ", ", m.coefficients, ")")
Base.:(==)(p1::Interpolation, p2::Interpolation) =
    (p1.functions == p2.functions) && (p1.coefficients == p2.coefficients)


# -------------------------------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------------------------------

_mul(C::AbstractVector, h::AbstractVector) = C ⋅ h
_mul(C::AbstractMatrix, h::AbstractVector) = C * h
_ncoefficients(C::AbstractVector) = length(C)
_ncoefficients(C::AbstractMatrix) = size(C, 2)


