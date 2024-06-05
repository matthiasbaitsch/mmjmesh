struct ValueAtLF <: AbstractMapping{AbstractMapping,Real,Any}
    x
    ValueAtLF(x::Real) = new(x)
    ValueAtLF(x::AbstractArray{<:Real}) = new(Vector(x))
end
valueat(u::ValueAtLF, f::FunctionToR) = f(u.x)

struct DerivativeAtLF <: AbstractMapping{AbstractMapping,Real,Any}
    x
    DerivativeAtLF(x::Real) = new(x)
    DerivativeAtLF(x::AbstractArray{<:Real}) = new(Vector(x))
end
valueat(u::DerivativeAtLF, f::FunctionRToR) = f'(u.x)

struct DDerivativeAtLF <: AbstractMapping{AbstractMapping,Real,Any}
    x
    d
    DDerivativeAtLF(x::AbstractArray{<:Real}, d::AbstractArray{<:Real}) = new(Vector(x), Vector(d))
end
valueat(u::DDerivativeAtLF, f::FunctionToR) = gradientat(f, u.x) ⋅ u.d

struct PDerivativeAtLF <: AbstractMapping{AbstractMapping,Real,Any}
    x
    n
    PDerivativeAtLF(x::AbstractArray{<:Real}, n::AbstractArray{<:Integer}) = new(Vector(x), Vector(n))
end
valueat(u::PDerivativeAtLF, f::FunctionToR) = derivativeat(f, u.x, u.n)

