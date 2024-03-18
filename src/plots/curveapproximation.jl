# -------------------------------------------------------------------------------------------------
# Angle between vectors in degrees
# -------------------------------------------------------------------------------------------------

clamp1(x::Real) = clamp(x, -1, 1)

Base.angle(u::AbstractVector, v::AbstractVector) = u ⋅ v / (norm(u) * norm(v)) |>
                                                   clamp1 |> acos |> abs |> rad2deg


# -------------------------------------------------------------------------------------------------
# Mapping parameter and corresponding point in 2D or 3D space
# -------------------------------------------------------------------------------------------------

struct PP{D}
    x::Float64
    point::SVector{D,Float64}
    PP(x, y) = new{2}(x, SA[x, y])
    PP(x, y::SVector{D,Float64}) where {D} = new{D}(x, y)
    PP(x, y::AbstractArray) = new{length(y)}(x, SVector{length(y),Float64}(y))
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
Base.isvalid(s::Segment) = !isfinite(s.pp1.point[end]) || !isfinite(s.pp2.point[end])
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

hasroot(s::Segment{2}) = point1(s)[end] * point2(s)[end] < -1e-12

function root(s::Segment{2})
    x1, y1 = point1(s)
    x2, y2 = point2(s)
    return abs(y2 - y1) > 1e-14 ? (x2 * y1 - x1 * y2) / (y1 - y2) : 0.5 * (x1 + x2)
end

function e1(s::Segment)
    u = point2(s) - point1(s)
    return u / norm(u)
end


# -------------------------------------------------------------------------------------------------
# Curve approximation
# -------------------------------------------------------------------------------------------------

struct CurveApproximation{D}
    head::Segment
    f::MappingFromR
end

function CurveApproximation(params::AbstractArray{<:Real}, f::MappingFromR)
    xy = [PP(x, f(x)) for x ∈ params]
    head = s = Segment(xy[1], xy[2])
    for i ∈ 2:length(params)-1
        s = connect!(s, Segment(xy[i], xy[i+1]))
    end
    return CurveApproximation{length(point1(head))}(head, f)
end

head(p::CurveApproximation) = p.head
Base.iterate(p::CurveApproximation) = (p.head, next(p.head))
Base.iterate(_::CurveApproximation, state) = isnothing(state) ? nothing : (state, next(state))

function Base.length(p::CurveApproximation)
    cnt = 0
    for _ ∈ p
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

function tail(p::CurveApproximation)
    s = p.head
    while hasnext(s)
        s = next(s)
    end
    return s
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

function refine!(p::CurveApproximation)
    for s ∈ p
        if s.refine
            xmid = 0.5 * (param1(s) + param2(s))
            refine!(s, PP(xmid, p.f(xmid)))
        end
    end
    return nothing
end

function handlenotfinite(x::SVector)
    if !isfinite(x[end])
        x = similar(x)
        x .= NaN
    end
    return x
end

function collectpoints(ca::CurveApproximation{D}, ir::Bool) where {D}

    # Segments containing a root
    segmentindex = 1
    segmentswithroot = Int[]
    for s ∈ ca
        ir && hasroot(s) && push!(segmentswithroot, segmentindex)
        segmentindex += 1
    end

    # Initialize
    rootindex = 1
    pointindex = 1
    segmentindex = 1
    segment = head(ca)
    np = length(ca) + length(segmentswithroot) + 1

    params = zeros(np)
    points = zeros(D, np)
    params[pointindex] = param1(segment)
    points[:, pointindex] = handlenotfinite(point1(segment))
    pointindex += 1

    # Collect
    for segment ∈ ca
        if rootindex <= length(segmentswithroot) && segmentindex == segmentswithroot[rootindex]
            p = root(segment)
            params[pointindex] = p
            points[1, pointindex] = p
            rootindex += 1
            pointindex += 1
        end
        params[pointindex] = param2(segment)
        points[:, pointindex] = handlenotfinite(point2(segment))
        pointindex += 1
        segmentindex += 1
    end

    return params, points
end
