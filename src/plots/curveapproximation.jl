# -------------------------------------------------------------------------------------------------
# Angle between vectors
# -------------------------------------------------------------------------------------------------

""" Clamp x to interval [-1, 1] """
clamp1(x::Real) = clamp(x, -1, 1)

""" Angle between vectors in degrees """
Base.angle(u::AbstractVector, v::AbstractVector) = u ⋅ v / (norm(u) * norm(v)) |>
                                                   clamp1 |> acos |> abs |> rad2deg


# -------------------------------------------------------------------------------------------------
# Parameter and corresponding point in 2D or 3D space
# -------------------------------------------------------------------------------------------------

struct PP{D}
    x::Float64
    point::SVector{D,Float64}
    PP(x::Real, y::Real) = new{2}(x, SA[x, y])
end


# -------------------------------------------------------------------------------------------------
# Curve approximation segment
# -------------------------------------------------------------------------------------------------

mutable struct Segment{D}
    level::Integer
    refine::Bool
    prev::Union{Segment,Nothing}
    pp1::PP{D}
    pp2::PP{D}
    next::Union{Segment,Nothing}
    Segment(p1::PP{D}, p2::PP{D}) where {D} = new{D}(0, false, nothing, p1, p2, nothing)
end

prev(s::Segment) = s.prev
next(s::Segment) = s.next
param1(s::Segment) = s.pp1.x
point1(s::Segment) = s.pp1.point
param2(s::Segment) = s.pp2.x
point2(s::Segment) = s.pp2.point
hasprev(s::Segment) = !isnothing(s.prev)
hasnext(s::Segment) = !isnothing(s.next)
Base.angle(s1::Nothing, s2::Segment) = 0
Base.angle(s1::Segment, s2::Segment) = angle(point2(s1) - point1(s1), point2(s2) - point1(s2))
Base.angle(s1::Segment, s2::Nothing) = 0
Base.isvalid(s::Segment) = isnan(s.pp1.point[end]) || isnan(s.pp2.point[end])
Base.show(io::IO, s::Segment) = Base.print(io, "Segment[$(s.level), $(s.refine), $(s.pp1.point), $(s.pp2.point)]")

function connect!(s1::Segment, s2::Segment)
    @assert s1.pp2 === s2.pp1
    s1.next = s2
    s2.prev = s1
    return s2
end

function connect!(s1::Segment, _::Nothing)
    s1.next = nothing
end

function refine!(s::Segment, pp::PP)
    nextsegment = next(s)
    newsegment = Segment(pp, s.pp2)
    s.pp2 = pp
    s.level += 1
    newsegment.level = s.level
    connect!(s, newsegment)
    connect!(newsegment, nextsegment)
    return newsegment
end


# -------------------------------------------------------------------------------------------------
# Curve approximation
# -------------------------------------------------------------------------------------------------

struct CurveApproximation
    head::Segment
    tail::Segment
end

function CurveApproximation(params::AbstractArray{<:Real}, f)
    xy = [PP(x, f(x)) for x ∈ params]
    head = tail = Segment(xy[1], xy[2])
    for i ∈ 2:length(params)-1
        tail = connect!(tail, Segment(xy[i], xy[i+1]))
    end
    return CurveApproximation(head, tail)
end

Base.iterate(p::CurveApproximation) = (p.head, next(p.head))
Base.iterate(_::CurveApproximation, state) = isnothing(state) ? nothing : (state, next(state))

function Base.length(p::CurveApproximation)
    cnt = 0
    for s ∈ p
        cnt += 1
    end
    return cnt
end

function Base.show(io::IO, p::CurveApproximation)
    for s ∈ p
        println(io, s)
        if hasnext(s)
            @printf(io, "%.1f°\n", angle(s, next(s)))
        end
    end
end

function mark!(p::CurveApproximation, maxangle::Real, maxrecursion::Integer)
    refine = false
    for s ∈ p
        s.refine = s.level < maxrecursion && (
            isvalid(s) ||
            angle(prev(s), s) > maxangle || angle(s, next(s)) > maxangle
        )
        refine |= s.refine
    end
    return refine
end

function refine!(p::CurveApproximation, f)
    for s ∈ p
        if s.refine
            xmid = 0.5 * (param1(s) + param2(s))
            refine!(s, PP(xmid, f(xmid)))
        end
    end
end

