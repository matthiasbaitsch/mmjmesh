"""
    sample1d(f, a, b; rp=false, maxrecursion=15, maxangle=2.5, npoints=5, yscale=1, ir=false) -> [params, ] points

Simple adaptive sampling of mapping `f` on the interval from `a` to `b`. The interval is initially
sampled at `npoints` and then refined recursively until either the angles are smaller than 
`maxangle` (in degrees) or the recursion depth reaches `maxrecursion`. The function returns the
sampled points in a matrix. If `rp=true`, the function returns parameters and corresponding points.

Parameters for real valued functions only: Angles are evaluated using the scaling factor `yscale`. 
Use `ir=true` to insert roots.

@see adapted_grid.jl in PlotUtils.jl
@see T. Bayer: Efficient plotting the functions with discontinuities
"""
function sample1d(
    f::MappingFromR, a::Real, b::Real;
    rp::Bool=false,
    maxrecursion::Integer=15, maxangle::Real=2.5, npoints::Integer=5, yscale::Real=1.0,
    ir::Bool=false
)

    # Check input
    @assert !ir || codomaintype(f) <: Real

    # Scaling factor
    f = yscale * SafeEval(f)

    # Helpers
    w = (b - a)
    h = w / (npoints - 1)

    # Initial sampling points
    x = wiggle!(collect(range(start=a, stop=b, length=npoints)), 0.5 * h)

    # Curve approximation and refinement
    ca = CurveApproximation(x, f)
    while mark!(ca, maxangle, maxrecursion)
        refine!(ca)
    end

    # Collect points with optional roots and scale
    params, points = collectpoints(ca, ir)
    points[2, :] /= yscale

    # Return
    if rp
        return params, points
    else
        return points
    end
end

sample1d(
    f::Function, a::Real, b::Real;
    maxrecursion::Integer=15, maxangle::Real=2.5, npoints::Integer=5, yscale::Real=1.0, ir::Bool=false
) = sample1d(
    MF(f), a, b,
    maxrecursion=maxrecursion, maxangle=maxangle, npoints=npoints, yscale=yscale, ir=ir
)


# -------------------------------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------------------------------

function wiggle!(x::Vector{<:Real}, eps::Real)
    n = length(x)
    x[2:n-1] += 0.5 * rand(Random.Random.Xoshiro(22421), -eps .. eps, n - 2)
    return x
end


struct MF <: FunctionRToR{Any}
    f::Function
end
MMJMesh.Mathematics.valueat(f::MF, x::Real) = f.f(x)


struct SafeEval{CT,D} <: AbstractMapping{Real,CT,D}
    f::MappingFromR{CT,D}
end

function MMJMesh.Mathematics.valueat(se::SafeEval{CT}, x::Real) where {CT}
    local y
    try
        y = se.f(x)
    catch err
        if err isa DomainError
            y = NaN
        else
            rethrow(err)
        end
    end
    return y
end
