# -------------------------------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------------------------------

_mul(C::AbstractVector, h::AbstractVector) = C ⋅ h
_mul(C::AbstractMatrix, h::AbstractVector) = C * h
_ncoefficients(C::AbstractVector) = length(C)
_ncoefficients(C::AbstractMatrix) = size(C, 2)


# -------------------------------------------------------------------------------------------------
# Implementation
# -------------------------------------------------------------------------------------------------

"""
    Interpolation(h::MappingToRn, C)

The interpolating mapping ``m: R^n \\to Y``

``
    m(x) = \\sum_{i=1}^m h_i(x) C_i
``

with interpolation functions ``h_i : R^n \\to R`` and coefficients 
``C_i \\in Y`` where ``Y`` can be numbers or vectors.
"""
struct Interpolation{DT,CT,D} <: AbstractMapping{DT,CT,D}

    functions::MappingToRn
    coefficients::RealVecOrMat

    Interpolation(
        h::MappingToRn{DT,N,D}, C::AbstractVector{<:Real}
    ) where {DT,N,D} = new{DT,InR,D}(h, C)

    Interpolation(
        h::MappingToRn{DT,N,D}, C::AbstractMatrix{<:Real}
    ) where {DT,N,D} = new{DT,InRⁿ{size(C, 1)},D}(h, C)
end

Interpolation(
    h::AbstractVector{<:FunctionToR}, C::VecOrMat{<:Real}
) = Interpolation(MappingFromComponents(h...), C)

const InterpolationOnR{CT,D} = Interpolation{InR,CT,D}
const InterpolationRnToR{N,D} = Interpolation{InRⁿ{N},InR,D}
const InterpolationRnToRm{N,M,D} = Interpolation{InRⁿ{N},InRⁿ{M},D}


## Evaluation

valueat(m::Interpolation{DT}, x::DT) where {DT} = _mul(m.coefficients, m.functions(x))

derivative(m::InterpolationOnR, n::Integer=1) =
    Interpolation(derivative(m.functions, n), m.coefficients)

derivativeat(m::InterpolationOnR, x::InR, n::Integer=1) =
    _mul(m.coefficients, derivativeat(m.functions, x, n))

# Gradient
function derivative(m::InterpolationRnToR, n::Integer=1)
    n == 1 && return transpose(jacobian(m.functions)) * m.coefficients
    @notimplemented
end

function derivativeat(m::InterpolationRnToR{N}, x::InRⁿ{N}, n::Int=1) where {N}
    n == 1 && return jacobianat(m.functions, x)' * m.coefficients
    @notimplemented
end

# Jacobian
function derivative(m::InterpolationRnToRm, n::Integer=1)
    n == 1 && return m.coefficients * jacobian(m.functions)
    @notimplemented
end

function derivativeat(m::InterpolationRnToRm{N}, x::InRⁿ{N}, n::Int=1) where {N}
    n == 1 && return m.coefficients * jacobianat(m.functions, x)
    @notimplemented
end


## Show

Base.show(io::IO, m::Interpolation) =
    print(io, "Interpolation(", m.functions, ", ", m.coefficients, ")")
Base.:(==)(p1::Interpolation, p2::Interpolation) =
    (p1.functions == p2.functions) && (p1.coefficients == p2.coefficients)


