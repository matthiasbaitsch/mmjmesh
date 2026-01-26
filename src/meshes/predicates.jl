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
    any_node_in(nodes::AbstractArray{<:Integer})
    all_nodes_in(nodes::AbstractArray{<:Integer})

Generates a predicate which tests if any or all nodes of an object are in the specified set of
nodes.
"""
all_nodes_in(nodes::AbstractArray{<:Integer}) = _nodes_in(nodes, select=all)
any_node_in(nodes::AbstractArray{<:Integer}) = _nodes_in(nodes, select=any)

# Helpers

function _nodes_in(nodes::AbstractArray{<:Integer}; select)
    nodeset = Set(nodes)
    return e -> select(n -> (n ∈ nodeset), nodeindices(e)) ? index(e) : NaN
end

