using Base.Iterators


"""
    linesegment(p1, p2)

Returns the line segment from `p1` to `p2` as parametric curve.
"""
function linesegment(p1::RealVec, p2::RealVec)
    Ns = MappingFromComponents(nodalbasis(makeelement(:lagrange, IHat, k=1))...)
    return Interpolation(Ns, tomatrix([p1, p2]))
end


"""
    polynomialinterpolation(points)
    polynomialinterpolation(p1, p2, p3, ...)
    polynomialinterpolation([p1, p2, p3, ...])

Constructs a polynomial parametric curve interpolating the specified points. 
The curve is defined on the unit interval [-1, 1].
"""
function polynomialinterpolation(pts::AbstractMatrix{<:Real})
    dps = [0, cumsum(norm.(eachcol(diff(pts, dims=2))))] |> flatten |> collect
    phi = affinefunction(0 .. dps[end], IHat)
    return Interpolation(lagrangepolynomials(phi.(dps), IHat), pts)
end
polynomialinterpolation(pts::AbstractVector) = polynomialinterpolation(tomatrix(pts))
polynomialinterpolation(pts::RealVec...) = polynomialinterpolation(collect(pts))
