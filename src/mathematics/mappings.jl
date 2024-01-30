using IntervalSets
import Polynomials
import CairoMakie as cm

struct Everything end
Base.in(x, ::Everything) = true

# -------------------------------------------------------------------------------------------------
# General concept of mapping
# -------------------------------------------------------------------------------------------------
abstract type AbstractMapping{DT,CT,D} end
function valueat end
function derivativeat end
function derivative end
domain(::AbstractMapping{DT,CT,D}) where {DT,CT,D} = D
Base.similar(a::AbstractArray, ::Type{T}, dims::Base.OneTo...) where {T<:AbstractMapping} = similar(a, AbstractMapping, Base.to_shape(dims))

# Shorthand notations f' and f(x)
Base.adjoint(m::AbstractMapping) = derivative(m)
(m::AbstractMapping{DT,CT,D})(x::DT) where {DT,CT,D} = valueat(m, x)

# TODO: Do we need that?
struct Derivative{DT,CT,D} <: AbstractMapping{DT,CT,D}
    m::AbstractMapping
    Derivative(m::AbstractMapping{DT,CT,D}) where {DT,CT,D} = new{DT,CT,D}(m)
end
valueat(d::Derivative{DT,CT,D}, x::DT) where {DT,CT,D} = derivativeat(d.m, x)
derivative(m::AbstractMapping) = Derivative(m)
Base.show(io::IO, f::Derivative) = print(io, "($(f.m))'")


struct Zero{DT,CT,D} <: AbstractMapping{DT,CT,D} end
valueat(::Zero{DT,CT,D}, x::DT) where {DT,CT,D} = DT(0)
derivativeat(::Zero{DT,CT,D}, x::DT) where {DT,CT,D} = DT(0)
derivative(z::Zero) = z
Base.zero(::AbstractMapping{DT,CT,D}) where {DT,CT,D} = Zero{DT,CT,D}()
Base.show(io::IO, z::Zero) = print(io, "0")

struct Sum{DT,CT,D} <: AbstractMapping{DT,CT,D}
    m1::AbstractMapping
    m2::AbstractMapping
    Sum(m1::AbstractMapping{DT,CT,D1}, m2::AbstractMapping{DT,CT,D2}) where {DT,CT,D1,D2} = new{DT,CT,D1 ∩ D2}(m1, m2)
end
valueat(s::Sum, x) = valueat(s.m1, x) + valueat(s.m2, x)
derivative(f::Sum) = derivative(f.m1) + derivative(f.m2)
derivativeat(f::Sum{Real,Real,D}, x::Real) where {D} = derivativeat(f.m1, x) + derivativeat(f.m2, x)
Base.show(io::IO, s::Sum) = print(io, s.m1, " + ", s.m2)

struct ScaledMapping{DT,CT,D} <: AbstractMapping{DT,CT,D}
    a::Real
    m::AbstractMapping
    ScaledMapping(a::Real, m::AbstractMapping{DT,CT,D}) where {DT,CT,D} = new{DT,CT,D}(a, m)
end
valueat(s::ScaledMapping, x) = s.a * valueat(s.m, x)
derivative(f::ScaledMapping) = f.a * derivative(f.m)
derivativeat(f::ScaledMapping{Real,Real,D}, x::Real) where {D} = f.a * derivativeat(f.m, x)
Base.show(io::IO, s::ScaledMapping) = print(io, s.a, " * ", s.m)

# TODO: Concatenation

Base.:+(m1::AbstractMapping, m2::AbstractMapping) = Sum(m1, m2)
Base.:+(::Zero, m::AbstractMapping) = m
Base.:+(m::AbstractMapping, ::Zero) = m
Base.:+(m1::AbstractMapping, m2::ScaledMapping) = m1 === m2.m ? ScaledMapping(1 + m2.a, m2.m) : Sum(m1, m2)
Base.:+(m1::ScaledMapping, m2::AbstractMapping) = m2 + m1
Base.:+(::Zero, m::ScaledMapping) = m
Base.:*(a::Real, m::AbstractMapping) = a == 1 ? m : a == 0 ? zero(m) : ScaledMapping(a, m)
Base.:*(m::AbstractMapping, a::Real) = a * m
Base.:*(a::Real, m::ScaledMapping) = (a * m.a) * m.m
Base.:+(m1::ScaledMapping, m2::ScaledMapping) = m1.m === m2.m ? ScaledMapping(m1.a + m2.a, m1.m) : Sum(m1, m2)

# -------------------------------------------------------------------------------------------------
# Functions R → R
# -------------------------------------------------------------------------------------------------
const AbstractFunctionRToR{D} = AbstractMapping{Real,Real,D}

const R = -Inf .. Inf
const IHat = -1.0 .. 1.0

# TODO: Adaptive sampling
function sample(f::AbstractFunctionRToR, n=100)
    dom = domain(f)
    a, b = width(dom) != Inf ? endpoints(dom) : (0, 5)
    x = LinRange(a, b, n)
    y = f.(x)
    x, y
end

# TODO: Recipe
function plot(f::AbstractFunctionRToR, fs::AbstractFunctionRToR...)
    n = 100
    f = cm.lines(sample(f, n)...)
    for g ∈ fs
        cm.lines!(sample(g, n)...)
    end
    f
end

function antiderivative end
function integrate(f::AbstractFunctionRToR, I) 
    F = antiderivative(f)
    return F(rightendpoint(I)) - F(leftendpoint(I))
end

# -------------------------------------------------------------------------------------------------
# Special functions
# -------------------------------------------------------------------------------------------------

# Sine and Cosine
struct Sin{D} <: AbstractFunctionRToR{D}
    Sin(d=R) = new{d}()
end
valueat(::Sin, x::Real) = sin(x)
derivativeat(::Sin, x::Real) = cos(x)
Base.show(io::IO, ::Sin) = print(io, "sin(x)")

struct Cos{D} <: AbstractFunctionRToR{D}
    Cos(d=R) = new{d}()
end
valueat(::Cos, x::Real) = cos(x)
derivativeat(::Cos, x::Real) = -sin(x)
Base.show(io::IO, ::Cos) = print(io, "cos(x)")

derivative(::Sin{D}) where {D} = Cos(D)
derivative(::Cos{D}) where {D} = -1.0 * Sin(D)

# Polynomials
struct Polynomial{D} <: AbstractFunctionRToR{D}
    p::Polynomials.Polynomial
end
Polynomial(p::Polynomials.Polynomial, d=R) = Polynomial{d}(p)
Polynomial(c::AbstractArray, d=R) = Polynomial{d}(Polynomials.Polynomial(c))
valueat(p::Polynomial, x::Real) = p.p(x)
derivative(p::Polynomial{D}) where {D} = Polynomial(Polynomials.derivative(p.p), D)
antiderivative(p::Polynomial{D}) where {D} = Polynomial(Polynomials.integrate(p.p), D)

# TODO: derivativeat
Base.:+(p1::Polynomial{D1}, p2::Polynomial{D2}) where {D1,D2} = Polynomial(p1.p + p2.p, D1 ∩ D2)
Base.:*(p1::Polynomial{D1}, p2::Polynomial{D2}) where {D1,D2} = Polynomial(p1.p * p2.p, D1 ∩ D2)
Base.:*(a::Real, p::Polynomial{D}) where {D} = Polynomial(a * p.p, D)

# Monomials
function monomials(p::AbstractArray{Int}, d=R)
    function coeffs(pp::Int)
        c = zeros(pp + 1)
        c[pp+1] = 1
        return c
    end
    return [Polynomial(coeffs(i), d) for i in p]
end

fromroots(r::AbstractArray{<:Real}, d=R) = Polynomial(Polynomials.fromroots(r), d)

function lagrangepolynomials(c::AbstractArray{<:Real}, d=R)
    indices = 1:length(c)
    normalize(f, x) = 1 / f(x) * f
    return [
        normalize(
            fromroots(c[filter(j -> j != i, indices)], d), 
            c[i]
        ) for i in indices
    ]
end

