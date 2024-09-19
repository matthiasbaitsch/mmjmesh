function linesegment(p1::RealVec, p2::RealVec)
    Ns = MappingFromComponents(nodalbasis(makeelement(:lagrange, IHat, k=1))...)
    return Interpolation(Ns, tomatrix([p1, p2]))
end
