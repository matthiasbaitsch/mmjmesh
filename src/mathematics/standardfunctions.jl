# -------------------------------------------------------------------------------------------------
# Sine function
# -------------------------------------------------------------------------------------------------

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


# -------------------------------------------------------------------------------------------------
# Cosine function
# -------------------------------------------------------------------------------------------------

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


# -------------------------------------------------------------------------------------------------
# Exp function
# -------------------------------------------------------------------------------------------------

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


# -------------------------------------------------------------------------------------------------
# Hyperbolic sine function
# -------------------------------------------------------------------------------------------------

struct Sinh{D} <: FunctionRToR{D}
    Sinh(d=R) = new{d}()
end

Sinh(::Identity{DT,D}) where {DT,D} = Sinh(D)
Sinh(f::T) where {T<:FunctionToR} = Sinh() ∘ f
Base.sinh(f::T) where {T<:FunctionToR} = Sinh(f)

valueat(::Sinh, x::InR) = sinh(x)
derivativeat(::Sinh, x::InR, n::Integer=1) = mod(n, 2) == 0 ? sinh(x) : cosh(x)
derivative(::Sinh{D}, n::Integer=1) where {D} = mod(n, 2) == 0 ? Sinh(D) : Cosh(D)
_antiderivative(::Sinh{D}) where {D} = Cosh(D)
Base.show(io::IO, ::Sinh) = print(io, "sinh(x)")


# -------------------------------------------------------------------------------------------------
# Hyperbolic cosine function
# -------------------------------------------------------------------------------------------------

struct Cosh{D} <: FunctionRToR{D}
    Cosh(d=R) = new{d}()
end

Cosh(::Identity{DT,D}) where {DT,D} = Cosh(D)
Cosh(f::T) where {T<:FunctionToR} = Cosh() ∘ f
Base.cosh(f::T) where {T<:FunctionToR} = Cosh(f)

valueat(::Cosh, x::InR) = cosh(x)
derivativeat(::Cosh, x::InR, n::Integer=1) = mod(n, 2) == 0 ? cosh(x) : sinh(x)
derivative(::Cosh{D}, n::Integer=1) where {D} = mod(n, 2) == 0 ? Cosh(D) : Sinh(D)
_antiderivative(::Cosh{D}) where {D} = Sinh(D)
Base.show(io::IO, ::Cosh) = print(io, "cosh(x)")


# -------------------------------------------------------------------------------------------------
# Absolute value function
# -------------------------------------------------------------------------------------------------

struct Abs{D} <: FunctionRToR{D}
    Abs(d=R) = new{d}()
end

Abs(::Identity{DT,D}) where {DT,D} = Abs(D)
Abs(f::T) where {T<:FunctionToR} = Abs() ∘ f
Base.abs(f::T) where {T<:FunctionToR} = Abs(f)

valueat(::Abs, x::InR) = abs(x)
derivativeat(::Abs, x::InR, n::Integer=1) = n == 0 ? abs(x) : n == 1 ? sign(x) : 0.0
derivative(::Abs{D}, n::Integer=1) where {D} = n == 0 ? Abs(D) : n == 1 ? Sign(D) : Zero{InR,InR,D}()
_antiderivative(::Abs{D}) where {D} = 0.5 * Identity(D) * Abs(D)
pois(::Abs{D}) where {D} = 0 ∈ D ? [0.0] : Float64[]
Base.show(io::IO, ::Abs) = print(io, "abs(x)")


# -------------------------------------------------------------------------------------------------
# Sign function
# -------------------------------------------------------------------------------------------------

struct Sign{D} <: FunctionRToR{D}
    Sign(d=R) = new{d}()
end

Sign(::Identity{DT,D}) where {DT,D} = Sign(D)
Sign(f::T) where {T<:FunctionToR} = Sign() ∘ f
Base.sign(f::T) where {T<:FunctionToR} = Sign(f)

valueat(::Sign, x::InR) = sign(x)
derivativeat(::Sign, x::InR, n::Integer=1) = n == 0 ? sign(x) : 0.0
derivative(::Sign{D}, n::Integer=1) where {D} = n == 0 ? Sign(D) : Zero{InR,InR,D}()
_antiderivative(::Sign{D}) where {D} = Abs(D)
pois(::Sign{D}) where {D} = 0 ∈ D ? [0.0] : Float64[]
Base.show(io::IO, ::Sign) = print(io, "sign(x)")
