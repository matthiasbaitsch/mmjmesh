
# -------------------------------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------------------------------

mul(C::AbstractVector, h::AbstractVector) = C ⋅ h
mul(C::AbstractMatrix, h::AbstractVector) = C * h
ncoefficients(C::AbstractVector) = length(C)
ncoefficients(C::AbstractMatrix) = size(C, 2)


# -------------------------------------------------------------------------------------------------
# Lagrange shape functions
# -------------------------------------------------------------------------------------------------

"""
    lagrangebasis1d(p)

Constructs a 1D Lagrange basis of degree `p` on the interval `IHat`.
"""
lagrangebasis1d(p::Integer) =
    ParametricCurve(lagrangepolynomials(range(IHat, length=(p + 1)), IHat)...)


# -------------------------------------------------------------------------------------------------
# Interpolation
# -------------------------------------------------------------------------------------------------

"""
    Interpolation(c)

The interpolating mapping ``m : R^n -> Y`` with

  ``m(x) = \\sum_{i=1}^m h_i(x) C_i``

with interpolation functions ``h_i : R^n -> R``, coefficients 
``C_i \\in Y`` and ``Y \\in R^{k_1} \\times \\dots \\times R^{k_l}``.
"""
struct Interpolation{DT,CT,D} <: AbstractMapping{DT,CT,D}
    functions::AbstractMapping{DT,<:AbstractVector,D}
    coefficients::SArray
    Interpolation(h::AbstractMapping{DT,<:AbstractVector,D}, C::AbstractVector) where {DT,D} =
        new{DT,Float64,D}(h, SVector{length(C)}(C))
    Interpolation(h::AbstractMapping{DT,<:AbstractVector,D}, C::AbstractMatrix) where {DT,D} =
        new{DT,SVector{size(C, 1),Float64},D}(h, SMatrix{size(C, 1),size(C, 2)}(C))
end

function Interpolation(C::AbstractVecOrMat, d::Integer)
    nc = ncoefficients(C)
    if d == 1
        return Interpolation(lagrangebasis1d(nc - 1), C)
    else
        error("Not implemented")
    end
end

valueat(m::Interpolation, x) = mul(m.coefficients, m.functions(x))
derivativeat(m::Interpolation, x, n::Int=1) = mul(m.coefficients, derivativeat(m.functions, x, n))
derivative(m::Interpolation, n::Int=1) =
    Interpolation(derivative(m.functions, n), m.coefficients)
Base.show(io::IO, m::Interpolation) =
    print(io, "Interpolation(", m.functions, ", ", m.coefficients, ")")
Base.:(==)(p1::Interpolation, p2::Interpolation) =
    (p1.functions == p2.functions) && (p1.coefficients == p2.coefficients)




