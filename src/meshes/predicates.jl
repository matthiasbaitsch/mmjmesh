using MMJMesh
using MMJMesh.Geometries

on(s::GeometricObjectP; atol=1e-12) = n::Node -> parameterof(s, coordinates(n), atol=atol)
on(p::Point; atol=1e-12) = n -> isapprox(p.coordinates, coordinates(n), atol=atol) ? 1 : NaN