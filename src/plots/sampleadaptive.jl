"""
    sampleadaptive(f, a::Real, b::Real, ir::Bool, atol::Real, rtol::Real, level::Int)

Adaptive sampling of function `f` in the interval from `a` to `b`.

@see https://yacas.readthedocs.io/en/latest/book_of_algorithms/basic.html
"""
function sampleadaptive(f::Function, a::Real, b::Real; ir::Bool=false, tol::Real=1e-2, level::Int=10)
    r = valuerange(f, a, b, 20)
    return dosampleadaptive(f, a, b, ir, r * tol, level)
end

const X1 = 0:1/4:1
const L1 = lagrangepolynomials(X1, 0 .. 1)
const W1 = [integrate(L, 0 .. 1) for L in L1]
const X2 = [0, 0.253124, 0.4986745, 0.75834, 1]
const L2 = lagrangepolynomials(X2, 0 .. 1)
const W2 = [integrate(L, 0 .. 1) for L in L2]

function safeeval(f::Function, x::Real)
    y = 0.0
    try
        y = f(x)
    catch e
        if typeof(e) == DomainError
            y = NaN
        else
            rethrow(e)
        end
    end
    return y
end

function wiggle(x::LinRange, a::Real)
    y = collect(x)
    for i = 2:length(x)-1
        y[i] += a * (-0.5 + rand())
    end
    return y
end

function valuerange(f::Function, a::Real, b::Real, n::Int)
    x = wiggle(LinRange(a, b, n), 1e-5)
    x = x[isfinite.(x)]
    y = safeeval.(f, x) |> sort!
    o = Int(round(0.1 * n))
    return y[end-o] - y[1+o]
end

function approximationerror(h::Real, ix::AbstractVector{<:Real}, iw::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    m = (y[end] - y[1]) / (ix[end] - ix[1])
    z = y[1] .+ m * ix
    d2 = (y - z) .^ 2
    return sqrt(h * dot(iw, d2))
end

function dosampleadaptive(f::Function, a::Real, b::Real, ir::Bool, tol::Real, level::Int)
    X = X1
    W = W1

    w = b - a
    x = a .+ w * X
    y = safeeval.(f, x)
    ok = all(isfinite.(y))

    # Check if we have to refine
    if level > 0 && (!ok || approximationerror(0.5^level, X, W, y) > tol)
        return refine(f, a, b, ir, tol, level)
    end

    if ok
        if ir
            return insertroot(x, y)
        else
            return x, y
        end
    else
        return nothing, nothing
    end
end

function refine(f, a::Real, b::Real, ir::Bool, tol::Real, level::Int)
    m = (a + b) / 2
    x1, y1 = dosampleadaptive(f, a, m, ir, tol / 2, level - 1)
    x2, y2 = dosampleadaptive(f, m, b, ir, tol / 2, level - 1)
    if !isnothing(y1) && !isnothing(y2)
        return vcat(x1, x2[2:end]), vcat(y1, y2[2:end])
    elseif !isnothing(y1)
        return vcat(x1, NaN), vcat(y1, NaN)
    else
        return vcat(NaN, x2), vcat(NaN, y2)
    end
end

function insertroot(xx, yy)
    xn = similar(xx, 0)
    yn = similar(yy, 0)
    il = 1
    for i in 1:length(xx)-1
        y1, y2 = yy[i:i+1]
        if y1 * y2 < 0
            x1, x2 = xx[i:i+1]
            m = (y2 - y1) / (x2 - x1)
            x = x1 - y1 / m
            append!(xn, xx[il:i], [x])
            append!(yn, yy[il:i], [0])
            il = i
        end
    end
    if length(xn) > 0
        append!(xn, xx[il+1:end])
        append!(yn, yy[il+1:end])
    end
    if length(xn) > 0
        return xn, yn
    else
        return xx, yy
    end
end


