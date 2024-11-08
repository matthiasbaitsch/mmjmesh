_midpoint(x, i) = (x[i] + x[i+1]) / 2

function _isascending(x)
    for i = 1:length(x)-1
        if x[i] >= x[i+1]
            return false
        end
    end
    return true
end

"""
    _makebreakpoints(pts, I)

Make breakpoints from Interval `I` and all points from `pts` which are in `I`.
"""
_makebreakpoints(pts, I) =
    unique!(vcat(leftendpoint(I), [p for p ∈ pts if p ∈ I], rightendpoint(I)))

"""
    PiecewiseFunction(breakpoints, functions)

Constructs a piecewise function.
"""
struct PiecewiseFunction{D} <: FunctionRToR{D}
    breakpoints::AbstractVector{<:Real}
    functions::AbstractVector{<:FunctionRToR}

    function PiecewiseFunction(
        breakpoints::AbstractVector{<:Real}, functions::AbstractVector{<:FunctionRToR}
    )
        @assert _isascending(breakpoints)
        @assert length(breakpoints) == length(functions) + 1
        return new{(breakpoints[1] .. breakpoints[end])}(breakpoints, functions)
    end
end

pois(f::PiecewiseFunction) = f.breakpoints[2:end-1]
npieces(f::PiecewiseFunction) = length(f.functions)
indexofpieceat(f::PiecewiseFunction, x::InR) =
    min(max(searchsortedfirst(f.breakpoints, x) - 1, 1), npieces(f))
functionat(f::FunctionRToR, ::InR) = f
functionat(f::PiecewiseFunction, x::InR) = f.functions[indexofpieceat(f, x)]
valueat(f::PiecewiseFunction, x::InR) = valueat(functionat(f, x), x)
derivativeat(f::PiecewiseFunction, x::InR, n::Integer=1) = derivativeat(functionat(f, x), x, n)

derivative(f::PiecewiseFunction, n::Integer=1) =
    PiecewiseFunction(f.breakpoints, derivative.(f.functions, n))

function integrate(f::PiecewiseFunction, I::Interval)
    s = 0.0
    pts = unique!(vcat([[leftendpoint(I)], [p for p ∈ pois(f) if p ∈ I], [rightendpoint(I)]]...))
    for i = 1:length(pts)-1
        sI = pts[i] .. pts[i+1]
        s += integrate(functionat(f, IntervalSets.mean(sI)), sI)
    end
    return s
end

function _applyoperator(f1::PiecewiseFunction, f2::FunctionRToR, op)
    breakpoints = _makebreakpoints(sort!(vcat(pois(f1), pois(f2))), domain(f1) ∩ domain(f2))

    functions = [
        begin
            mp = _midpoint(breakpoints, i)
            op(functionat(f1, mp), functionat(f2, mp))
        end
        for i = 1:length(breakpoints)-1
    ]

    return PiecewiseFunction(breakpoints, functions)
end

Base.:(+)(f1::FunctionRToR, f2::PiecewiseFunction) = _applyoperator(f2, f1, +)
Base.:(+)(f1::PiecewiseFunction, f2::FunctionRToR) = _applyoperator(f1, f2, +)
Base.:(+)(f1::PiecewiseFunction, f2::PiecewiseFunction) = _applyoperator(f1, f2, +)
Base.:(*)(f1::FunctionRToR, f2::PiecewiseFunction) = _applyoperator(f2, f1, *)
Base.:(*)(f1::PiecewiseFunction, f2::FunctionRToR) = _applyoperator(f1, f2, *)
Base.:(*)(f1::PiecewiseFunction, f2::PiecewiseFunction) = _applyoperator(f1, f2, *)
Base.:(*)(a::Real, f::PiecewiseFunction) = PiecewiseFunction(f.breakpoints, a * f.functions)

# TODO find better solution
Base.:(+)(f1::ScaledMapping{InR,InR}, f2::PiecewiseFunction) = _applyoperator(f2, f1, +)
Base.:(+)(f1::PiecewiseFunction, f2::ScaledMapping{InR,InR}) = _applyoperator(f1, f2, +)
Base.:(*)(f1::ScaledMapping{InR,InR}, f2::PiecewiseFunction) = _applyoperator(f2, f1, *)
Base.:(*)(f1::PiecewiseFunction, f2::ScaledMapping{InR,InR}) = _applyoperator(f1, f2, *)

"""
    interpolate(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}; order=1)

Construct a function which interpolates the points ``(x_1, y_1), (x_2, y_2), \\dots, (x_n, y_n)``. The entries in `x` have to be in ascending order and `x` and `y` need to have the same length.
"""
function interpolate(x::AbstractArray{<:Real}, y::AbstractArray{<:Real}; order=1)
    @assert order == 1
    @assert length(x) == length(y)
    @assert _isascending(x)

    return PiecewiseFunction(
        x, [affinefunction(x[i] .. x[i+1], y[i] .. y[i+1]) for i = 1:length(x)-1]
    )
end


