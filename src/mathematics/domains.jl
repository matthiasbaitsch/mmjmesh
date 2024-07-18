
""" The real numbers """
const R = -Inf .. Inf

""" Positive real numbers """
const RPlus = Interval{:open,:open}(0, Inf)
const R⁺ = RPlus

""" Non negative real numbers """
const R0Plus = Interval{:closed,:open}(0, Inf)
const R⁺₀ = R0Plus

""" Reference interval ``[-1, 1]`` """
const IHat = -1.0 .. 1.0
const ReferenceInterval = IHat

dimension(::AbstractInterval) = 1
Base.isfinite(I::Interval) = isfinite(width(I))

"""
    Base.intersect(s, x)
    Base.intersect(x, s)
    s ∩ x
    x ∩ s

Select all values from the vector `x` which are in the interval `s`.
"""
Base.intersect(s::AbstractInterval, a::AbstractVector{T}) where {T} = T[x for x in a if x ∈ s]
Base.intersect(a::AbstractVector, s::AbstractInterval) = intersect(s, a)

""" Real coordinate plane """
const R2 = R^2
const R² = R2

""" Real coordinate space """
const R3 = R^3
const R³ = R3

""" Reference square ``[-1, 1]^2`` """
const QHat = IHat × IHat
const ReferenceQuadrilateral = QHat

dimension(::Rectangle{<:SVector{N}}) where {N} = N
Base.isfinite(r::DomainSets.Rectangle) = all(isfinite.(r.a)) && all(isfinite.(r.b))

""" Element types """
const InR = Real
const InRⁿ{N} = SVector{N,<:Real}
const InRⁿˣᵐ{N,M} = SMatrix{N,M,<:Real}
const InR2 = InRⁿ{2}
const InR² = InR2
const InR3 = InRⁿ{3}
const InR³ = InR3
