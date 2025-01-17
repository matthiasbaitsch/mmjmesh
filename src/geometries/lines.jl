"""
    HLine(y)

A horizontal line in 2D which intercepts the y-axis at `y`.
"""
struct HLine <: GeometricObjectP{1,2}
    y::Real
end

parameterof(o::HLine, p::RealVec; atol::Real=1e-12) = abs(p[2] - o.y) <= atol ? p[1] : NaN


"""
    VLine(x)

A vertical line in 2D which intercepts the x-axis at `x`.
"""
struct VLine <: GeometricObjectP{1,2}
    x::Real
end

parameterof(o::VLine, p::RealVec; atol::Real=1e-12) = abs(p[1] - o.x) <= atol ? p[2] : NaN


"""
    Segment(p1, p2)

A line segment which connects the points `p1` and `p2`.
"""
struct Segment{DG} <: GeometricObjectP{1,DG}
    p1
    p2

    function Segment(p1, p2)
        @assert length(p1) == length(p2)
        return new{length(p1)}(p1, p2)
    end
end

function parameterof(s::Segment{2}, p::RealVec; atol::Real=1e-12)
    # Parameters
    u = s.p2 - s.p1
    n = normalize([-u[2], u[1]])

    # Intersection of lines
    # p1 + t1 * u = p + t2 * n
    # t1 * u - t2 * n = p - p1
    A = stack([u, -n], dims=2)
    b = p - s.p1
    t1, t2 = A \ b

    # Return condition
    if abs(t2) <= atol && -atol <= t1 && t1 <= 1 + atol
        return t1
    else
        return NaN
    end
end