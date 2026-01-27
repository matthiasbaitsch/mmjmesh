using MMJMesh
using MMJMesh.Geometries


"""
    on(o::Point; atol=1e-12)
    on(o::GeometricObjectP; atol=1e-12)

Generates a predicate which tests if an object is on the specified geometric object.
"""
on(o::Point; atol=1e-12) = node -> isapprox(o.coordinates, coordinates(node), atol=atol) ? 1 : NaN
on(o::GeometricObjectP; atol=1e-12) = entity -> _parameterof(o, entity, atol)

# Helpers

function _parameterof(o::GeometricObjectP, e::MeshEntity, atol)
    nn = 0
    pp = 0
    for c = eachcol(coordinates(e))
        p = parameterof(o, c, atol=atol)
        nn += 1
        pp += p
    end
    return pp / nn
end

_parameterof(o::GeometricObjectP, n::Node, atol) = parameterof(o, coordinates(n), atol=atol)


"""
    entities_in(pdim, idxs, select=all)

Generates a predicate which tests if indices of entities having parametric dimension `pdim` of an 
object are in the specified array vector `idxs`.
"""
function entities_in(pdim::Integer, idxs::IntegerVec; select)
    idxs = Set(idxs)
    return e -> select(i -> (i ∈ idxs), indices(e, pdim)) ? index(e) : NaN
end

