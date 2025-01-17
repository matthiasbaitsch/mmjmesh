# -------------------------------------------------------------------------------------------------
# Add nodes
# -------------------------------------------------------------------------------------------------

_lastnodeid(m::Mesh) = nnodes(m) == 0 ? 0 : MMJMesh.Topologies.entities(m.topology, 0)[end]

addnodes!(m::Mesh, xs::RealVec, ys::RealVec) = addnodes!(m, [xs, ys])
addnodes!(m::Mesh, xs::RealVec, ys::RealVec, zs::RealVec) = addnodes!(m, [xs, ys, zs])
addnodes!(m::Mesh, p1::RealVec, p2::RealVec, n::Int) = addnodes!(m, linesegment(p1, p2), n)

function addnode!(m::Mesh, p::RealVec)
    @assert length(p) == gdim(m)
    push!(m.geometry, Vector(p))
    push!(m.topology.entities[0], _lastnodeid(m) + 1)
    return nnodes(m)
end

function addnodes!(m::Mesh, ps::RealVecVec)
    if length(ps) == gdim(m)
        @assert all(length(ps[1]) .== length.(ps))
        return [addnode!(m, getindex.(ps, i)) for i ∈ 1:length(ps[1])]
    end
    if length(ps[1]) == gdim(m)
        return [addnode!(m, xx) for xx ∈ ps]
    end
end

function addnodes!(m::Mesh{DT,DG}, x::ParametricCurve{DG}, n::Int) where {DT,DG}
    @assert n >= 2
    @assert isfinite(domain(x))
    return [addnode!(m, x(t)) for t = range(domain(x), n)]
end



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
