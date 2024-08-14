# -------------------------------------------------------------------------------------------------
# Basic functionality
# -------------------------------------------------------------------------------------------------


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

""" Element of the real numbers """
const InR = Real
""" Element of the set of ``n``-tuples of real numbers """
const InRⁿ{N} = SVector{N,<:Real}
""" Element of the set of ``n \\times m``-matrices of real numbers """
const InRⁿˣᵐ{N,M} = SMatrix{N,M,<:Real}
""" Element of the set of pairs of real numbers """
const InR2 = InRⁿ{2}
""" Element of the set of pairs of real numbers """
const InR² = InR2
""" Element of the set of triples of real numbers """
const InR3 = InRⁿ{3}
""" Element of the set of triples of real numbers """
const InR³ = InR3


# -------------------------------------------------------------------------------------------------
# Points on domain
# -------------------------------------------------------------------------------------------------

""" 
    points(K, on, n=0)

Pick certain points from domain `K` at locations `on` which can take the following values

- `:corners` -- Returns points on corners of the domain

- `:sides` -- Gives `n` points equidistantly distributed on each side

- `:interior` -- Generates `n` by `n` points in the interior
"""
function points(K::AbstractInterval, on::Symbol, n::Integer=0)
    p1 = leftendpoint(K)
    p2 = rightendpoint(K)
    h = 1 // (n + 1)
    s = h:h:1-h

    if on == :corners
        return [p1, p2]
    elseif on == :interior
        return p1 .+ s' .* (p2 - p1) |> vec
    elseif on == :sides
        return []
    end

    error("Illegal value for on: $on")
end

function points(K::Rectangle, on::Symbol, n::Integer=0)
    p11, p12 = leftendpoint(K)
    p21, p22 = rightendpoint(K)
    c = [SA[p11, p12], SA[p21, p12], SA[p21, p22], SA[p11, p22]]
    h = 1 // (n + 1)
    s = h:h:1-h
    flatten(x) = reduce(vcat, x)

    if on == :corners
        return c
    elseif on == :sides
        return [
            [SVector{2,Float64}(c) for c in eachcol((c[i] .+ (s' .* (c[i%4+1] - c[i]))))]
            for i = 1:4
        ] |> flatten
    elseif on == :interior
        p1 = p11 .+ s' .* (p21 - p11)
        p2 = p12 .+ s' .* (p22 - p12)
        return reshape([SA[p1[i], p2[j]] for i = 1:n, j = 1:n], :, 1) |> vec
    end
    error("Illegal value for on: $on")
end

