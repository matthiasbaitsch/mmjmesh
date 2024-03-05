"""
    sampleadaptive2(f, a, b, maxrecursion, maxangle, npoints, yscale) -> Matrix{Real}

Simple adaptive sampling of mapping `f` on the interval from `a` to `b`. The interval is initially
sampled at `npoints` g and then refined recursively until either the angles are smaller than 
`maxangle` (in degrees) or the recursion depth reaches `maxrecursion`. For mappings into the real 
numbers, the angles are evaluated using the scaling factor `yscale`. The function returns the 
sampled points in a matrix.

@see adapted_grid.jl in PlotUtils.jl
@see T. Bayer: Efficient plotting the functions with discontinuities
"""
function sampleadaptive2(
    f::MappingFromR{CT}, a::Real, b::Real;
    maxrecursion::Integer=5, maxangle::Real=1, npoints::Integer=5, yscale::Real=1.0, ir::Bool=false
) where {CT}
    g = Vector{Pair{Real,CT}}()
    push!(g, Pair(a, f(a)))
    _sampleadaptive2!(g, f, a, b, maxrecursion, maxangle, npoints, yscale, Xoshiro(22421))
    push!(g, Pair(b, f(b)))
    return _tomatrix(g, ir)
end

sampleadaptive2(
    f::Function, a::Real, b::Real;
    maxrecursion::Integer=5, maxangle::Real=1, npoints::Integer=5, yscale::Real=1.0, ir::Bool=false
) = sampleadaptive2(
    MF(f), a, b,
    maxrecursion=maxrecursion, maxangle=maxangle, npoints=npoints, yscale=yscale, ir=ir
)


# -------------------------------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------------------------------


"""
    MF(f)

Construct `FunctionRToR` from Julia function.
"""
struct MF <: FunctionRToR{Any}
    f::Function
end
MMJMesh.Mathematics.valueat(f::MF, x::Real) = f.f(x)

""" Clamp x to interval [-1, 1] """
_clamp1(x::Real) = clamp(x, -1, 1)

""" Samples from fraph to points """
function _tomatrix(g::Vector{Pair{Real,Real}}, ir::Bool)
     return [p[i] for i ∈ 1:2, p ∈ g]
end

""" Angle between two vectors """
_angle(u::AbstractVector, v::AbstractVector) = u ⋅ v / (norm(u) * norm(v)) |>
                                               _clamp1 |> acos |> abs |> rad2deg

""" Angles between direction vectors """
function _angles(f::FunctionRToR, x::Vector{<:Real}, y::Vector{<:Real}, yscale::Real)
    n = length(x)
    d = diff(vcat(x', yscale * y'), dims=2)
    return [(i == 1 || i == n) ? 0 : _angle(d[:, i-1], d[:, i]) for i ∈ 1:n]
end

""" Move internal points of sampling points in `x` """
function _wiggle!(x::Vector{<:Real}, eps::Real, rng::AbstractRNG)
    n = length(x)
    x[2:n-1] += 0.5 * rand(rng, -eps .. eps, n - 2)
    return x
end

""" Implementation of adaptive sampling """
function _sampleadaptive2!(
    g::Vector{Pair{Real,CT}}, f::MappingFromR{CT}, a::Real, b::Real,
    d::Integer, maxangle::Real, npoints::Integer, yscale::Real, rng::AbstractRNG
) where {CT}

    # Helpers
    w = (b - a)
    h = w / (npoints - 1)
    nsegments = npoints - 1

    # Points and values (multiple evaluations of f at a, b for the sake of simplicity)
    x = _wiggle!(collect(range(start=a, stop=b, length=npoints)), 0.5 * h, rng)
    y = f.(x)

    # Compute angles
    φ = _angles(f, x, y, yscale)

    # Process intervals
    for i ∈ 1:nsegments
        if d > 0 && (φ[i] > maxangle || φ[i+1] > maxangle)
            _sampleadaptive2!(
                g, f, x[i], x[i+1], d - 1, maxangle, npoints, yscale, rng
            )
        end
        if i < nsegments
            push!(g, Pair(x[i+1], y[i+1]))
        end
    end
end

