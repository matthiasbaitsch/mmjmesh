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
    entities_in(pdim, idxs; select)

Generates a predicate which tests if `pdim`-dimensional entities of an entity are in `idxs`. Use `select=all` to require all entities in `idxs`, use `select=any` if any entity in `idxs` is sufficient.
"""
function entities_in(pdim::Integer, idxs::IntegerVec; select)
    idxs = Set(idxs)
    return e -> select(i -> (i ∈ idxs), indices(e, pdim)) ? index(e) : NaN
end

