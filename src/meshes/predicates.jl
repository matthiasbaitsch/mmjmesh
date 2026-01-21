using MMJMesh
using MMJMesh.Geometries


"""
    on(p::Point)
    on(s::GeometricObjectP)

Generates a predicate which tests if an object is on the specified geometric object. To be used in conjunction with the `entiy(m, dim, predicate)` function.
"""
on(p::Point; atol=1e-12) = n -> isapprox(p.coordinates, coordinates(n), atol=atol) ? 1 : NaN


_parameterof(s::GeometricObjectP, n::Node, atol) = parameterof(s, coordinates(n), atol=atol)

function _parameterof(s::GeometricObjectP, e::MeshEntity, atol)
    nn = 0
    pp = 0
    for c = eachcol(coordinates(e))
        p = parameterof(s, c, atol=atol)
        nn += 1
        pp += p
    end
    return pp / nn
end

on(s::GeometricObjectP; atol=1e-12) = e -> _parameterof(s, e, atol)

function hasnodes(nodes)
    nodeset = Set(nodes)

    function hasnodes(e)
        for n = nodeindices(e)
            if n ∉ nodeset
                return false
            end
        end
        return true
    end

    return e -> hasnodes(e) ? index(e) : NaN
end