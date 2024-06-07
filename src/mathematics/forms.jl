"""
    ValueAtLF(x)

Evaluates the function at `x`.
"""
struct ValueAtLF <: AbstractMapping{AbstractMapping,Real,Any}
    x
    ValueAtLF(x::Real) = new(x)
    ValueAtLF(x::AbstractArray{<:Real}) = new(Vector(x))
end
valueat(u::ValueAtLF, f::FunctionToR) = valueat(f, u.x)

"""
    DerivativeAtLF(x)

Evaluates the derivative at `x`.
"""
struct DerivativeAtLF <: AbstractMapping{AbstractMapping,Real,Any}
    x
    DerivativeAtLF(x::Real) = new(x)
    DerivativeAtLF(x::AbstractArray{<:Real}) = new(Vector(x))
end
valueat(u::DerivativeAtLF, f::FunctionRToR) = derivativeat(f, u.x)

"""
    PDerivativeAtLF(x, n)

Evaluates the partial derivative at `x` in directions `n` of a function 
``f : \\mathbb{R}^n \\to \\mathbb{R}``. The i-th entry of `n` specifies the
order of the partial derivative in direction `i`. For instance, `[1, 2]`
refers to the partial derivative ``f_{xyy}(\\mathbf{x})``.
"""
struct PDerivativeAtLF <: AbstractMapping{AbstractMapping,Real,Any}
    x
    n
    PDerivativeAtLF(x::AbstractArray{<:Real}, n::AbstractArray{<:Integer}) = new(Vector(x), Vector(n))
end
valueat(u::PDerivativeAtLF, f::FunctionToR) = derivativeat(f, u.x, u.n)

