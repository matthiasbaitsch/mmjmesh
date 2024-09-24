# -------------------------------------------------------------------------------------------------
# Add nodes
# -------------------------------------------------------------------------------------------------

_lastnodeid(m::Mesh) = nnodes(m) == 0 ? 0 : MMJMesh.Topologies.entities(m.topology, 0)[end]

function addnode!(m::Mesh, x::RealVec)
    @assert length(x) == gdim(m)
    push!(m.geometry, Vector(x))
    push!(m.topology.entities[0], _lastnodeid(m) + 1)
    return nnodes(m)
end

function addnodes!(m::Mesh, x::RealVecVec)
    if length(x) == gdim(m)
        @assert all(length(x[1]) .== length.(x))
        return [addnode!(m, getindex.(x, i)) for i ∈ 1:length(x[1])]
    end
    if length(x[1]) == gdim(m)
        return [addnode!(m, xx) for xx ∈ x]
    end
end

addnodes!(m::Mesh, x::RealVec, y::RealVec) = addnodes!(m, [x, y])
addnodes!(m::Mesh, x::RealVec, y::RealVec, z::RealVec) = addnodes!(m, [x, y, z])

function addnodes!(m::Mesh, x::ParametricCurve{N}, n::Int) where {N}
    @assert n >= 2
    @assert N == gdim(m)
    @assert isfinite(domain(x))
    return [addnode!(m, x(t)) for t = range(domain(x), n)]
end

addnodes!(m::Mesh, x1::RealVec, x2::RealVec, n::Int) =
    addnodes!(m, linesegment(x1, x2), n)


# -------------------------------------------------------------------------------------------------
# Add elements
# -------------------------------------------------------------------------------------------------

function _addgroup!(m::Mesh, group::Symbol, indexes::IntegerVec)
    g = Group{MeshEntity{pdim(m)}}(indexes)
    if group ∉ names(m.groups)
        m.groups[group] = g
    else
        m.groups[group] = m.groups[group] ∪ g
    end
    return indexes
end

_addgroup!(::Mesh, ::Nothing) = indexes -> indexes
_addgroup!(m::Mesh, group::Symbol) = indexes -> _addgroup!(m, group, indexes)

addelement!(m::Mesh, nodes::Integer...; group::Union{Symbol,Nothing}=nothing) =
    addelement!(m, collect(nodes), group=group)

addelement!(m::Mesh, nodes::IntegerVec; group::Union{Symbol,Nothing}=nothing) =
    addlinks!(m.topology, pdim(m), 0, [nodes]) |> _addgroup!(m, group)

function addelements!(m::Mesh, nodes::IntegerVec...; n=-1, group::Union{Symbol,Nothing}=nothing)
    ne = round(Int, n)
    nn = [collect(n) for n = zip(nodes...)]

    ne != -1 && (nn = nn[1:ne])

    return addlinks!(m.topology, pdim(m), 0, nn) |> _addgroup!(m, group)
end
