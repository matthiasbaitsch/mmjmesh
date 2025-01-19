using MMJMesh
using MMJMesh.Geometries


"""
    on(p::Point)
    on(s::GeometricObjectP)

Generates a predicate which tests if a node is on the specified geometric object. To be used
in conjunction with the `nodes(m, predicate)` function.
"""
on(s::GeometricObjectP; atol=1e-12) = n::Node -> parameterof(s, coordinates(n), atol=atol)
on(p::Point; atol=1e-12) = n -> isapprox(p.coordinates, coordinates(n), atol=atol) ? 1 : NaN