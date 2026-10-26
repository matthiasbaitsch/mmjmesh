## Sine function

struct Sin{D} <: FunctionRToR{D}
    Sin(d=R) = new{d}()
end

Sin(::Identity{DT,D}) where {DT,D} = Sin(D)
Sin(f::T) where {T<:FunctionToR} = Sin() ∘ f
Base.sin(f::T) where {T<:FunctionToR} = Sin(f)

valueat(::Sin, x::InR) = sin(x)
derivativeat(::Sin, x::InR, n::Integer=1) = sin(x + n * π / 2)
derivative(::Sin{D}, n::Integer=1) where {D} = (-1.0)^(n ÷ 2) * (mod(n, 2) == 0 ? Sin(D) : Cos(D))
_antiderivative(::Sin{D}) where {D} = -Cos(D)
Base.show(io::IO, ::Sin) = print(io, "sin(x)")


## Cosine function

struct Cos{D} <: FunctionRToR{D}
    Cos(d=R) = new{d}()
end

Cos(::Identity{DT,D}) where {DT,D} = Cos(D)
Cos(f::T) where {T<:FunctionToR} = Cos() ∘ f
Base.cos(f::T) where {T<:FunctionToR} = Cos(f)

valueat(::Cos, x::InR) = cos(x)
derivativeat(::Cos, x::InR, n::Integer=1) = cos(x + n * π / 2)
derivative(::Cos{D}, n::Integer=1) where {D} =
    (-1.0)^((n + 1) ÷ 2) * (mod(n, 2) == 0 ? Cos(D) : Sin(D))
_antiderivative(::Cos{D}) where {D} = Sin(D)
Base.show(io::IO, ::Cos) = print(io, "cos(x)")


## Exp function

struct Exp{D} <: FunctionRToR{D}
    Exp(d=R) = new{d}()
end

Exp(::Identity{DT,D}) where {DT,D} = Exp(D)
Exp(f::T) where {T<:FunctionRToR} = Exp() ∘ f
Base.exp(f::T) where {T<:FunctionRToR} = Exp(f)

valueat(::Exp, x::InR) = exp(x)
derivativeat(::Exp, x::InR, ::Integer=1) = exp(x)
derivative(f::Exp, ::Integer=1) = f
antiderivative(f::Exp, ::Integer=1) = f
Base.show(io::IO, ::Exp) = print(io, "exp(x)")
