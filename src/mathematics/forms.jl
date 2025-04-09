abstract type Form end
abstract type LinearForm <: AbstractMapping{AbstractMapping,Real,Any} end
abstract type BilinearForm <: AbstractMapping{SVector{2,AbstractMapping},Real,Any} end


"""
    ValueAtLF(x)

The linear form ``b`` with ``b(f) = f(x)`` for a specified position ``x``.
"""
struct ValueAtLF <: LinearForm
    x
    ValueAtLF(x::Real) = new(x)
    ValueAtLF(x::AbstractArray{<:Real}) = new(Vector(x))
end
valueat(u::ValueAtLF, f::FunctionToR) = valueat(f, u.x)

"""
    DerivativeAtLF(x)

The linear form ``b`` with ``b(f) = f'(x)`` for a specified position ``x``.
"""
struct DerivativeAtLF <: LinearForm
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
struct PDerivativeAtLF <: LinearForm
    x
    n
    PDerivativeAtLF(x::AbstractArray{<:Real}, n::AbstractArray{<:Integer}) =
        new(Vector(x), Vector(n))
end
valueat(u::PDerivativeAtLF, f::FunctionToR) = derivativeat(f, u.x, u.n)

∂xLF(x) = PDerivativeAtLF(x, [1, 0])
∂yLF(x) = PDerivativeAtLF(x, [0, 1])
∂xyLF(x) = PDerivativeAtLF(x, [1, 1])