# -------------------------------------------------------------------------------------------------
# Type definitions
# -------------------------------------------------------------------------------------------------

const MeshEntityGroup = Group{<:MeshEntity}
const NodeGroup = MeshEntityGroup{Node}
const EdgeGroup = MeshEntityGroup{Edge}
const FaceGroup = MeshEntityGroup{Face}
const SolidGroup = MeshEntityGroup{Solid}

edim(::Type{<:Group{<:MeshEntity{DT}}}) where {DT} = DT
edim(g::Group{<:MeshEntity}) = edim(typeof(g))


# -------------------------------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------------------------------

function _collectboundary(m::Mesh{1})
    bn = Int[]
    for e ∈ nodes(m)
        if nedges(e) == 1
            append!(bn, nodeindices(e))
        end
    end
    return bn, Int[], Int[]
end


function _collectboundary(m::Mesh{2})
    bn = Int[]
    be = Int[]
    for e ∈ edges(m)
        if nfaces(e) == 1
            push!(be, e.index)
            append!(bn, nodeindices(e))
        end
    end
    return bn, be, Int[]
end


"""
    _idvector(values::AbstractVector) -> ids

Assigns an integer id to each unique entry in `values` and returns the `ids` vector
such that `ids[i]` is the id of the `i`-th value. Larger values have a larger id such 
that the resulting `ids` array satisfies  `cmp(values[i], values[j]) == cmp(ids[i], ids[j])`.
"""
function _idvector(values::AbstractVector)
    d = Dict([(val, idx) for (idx, val) in enumerate(unique!(sort(values)))])
    return [d[val] for val in values]
end


"""
    _registergrouprecipies(m::Mesh)

Registers recipes to generate groups for boundary entities and entities of a certain
parametric dimension. To be called once on `Mesh` creation.
"""
function _registergrouprecipies(m::Mesh)

    # On the boundary
    function boundaryrecipe()
        bn, be, bf = _collectboundary(m)
        m.groups.entries[:boundarynodes] = NodeGroup(bn)
        m.groups.entries[:boundaryedges] = EdgeGroup(be)
        m.groups.entries[:boundaryfaces] = FaceGroup(bf)
        return nothing
    end

    addrecipe!(m.groups, :boundarynodes, boundaryrecipe)
    addrecipe!(m.groups, :boundaryedges, boundaryrecipe)
    addrecipe!(m.groups, :boundaryfaces, boundaryrecipe)

    # By dimension
    addrecipe!(m.groups, :nodes, () -> NodeGroup(1:nnodes(m)))
    addrecipe!(m.groups, :edges, () -> EdgeGroup(1:nedges(m)))
    addrecipe!(m.groups, :faces, () -> FaceGroup(1:nfaces(m)))
    addrecipe!(m.groups, :solids, () -> SolidGroup(1:nsolids(m)))

    # Elements - TODO elegant solution
    if pdim(m) == 1
        addrecipe!(m.groups, :elements, () -> EdgeGroup(1:nedges(m)))
    elseif pdim(m) == 2
        addrecipe!(m.groups, :elements, () -> FaceGroup(1:nfaces(m)))
    elseif pdim(m) == 3
        addrecipe!(m.groups, :elements, () -> SolidGroup(1:nsolids(m)))
    end
end


"""
    _collectgroups(m::Mesh; d::Int, predefined::Bool=false)

Collect all groups the mesh entities of parametric dimension `d` belong to.
"""
function _collectgroups(m::Mesh; d::Int, predefined::Bool=false)
    names = [Symbol[] for _ in 1:nentities(m, d)]
    for name ∈ groupnames(m, d=d, predefined=predefined)
        for index ∈ m.groups[name]
            push!(names[index], name)
        end
    end
    return names
end


# -------------------------------------------------------------------------------------------------
# Interface
# -------------------------------------------------------------------------------------------------

"""
    groupnames(m::Mesh; d::Int=-1, predefined::Bool=false)

Returns all groups of dimension `d`. If `d` is not specified,
groups of all dimensions are included. If `predefined` is `false`, only 
user-defined groups are considered.
"""
function groupnames(m::Mesh; d::Int=-1, predefined::Bool=false)
    gc = m.groups
    isvalid(name) =
        (predefined || !ispredefined(gc, name)) &&
        (d == -1 || (!isempty(gc[name]) && edim(gc[name]) == d))
    return [name for name ∈ names(gc) if isvalid(name)]
end


"""
    groupids(m::Mesh; d::Int, predefined::Bool=false)

Unique numerical ids groups of specified dimension.
"""
groupids(m::Mesh; d::Int, predefined::Bool=false) =
    _idvector(_collectgroups(m, d=d, predefined=predefined))


"""
    ngroups(m::Mesh; d::Int=-1, predefined::Bool=false)

Number of groups.
"""
ngroups(m::Mesh; d::Int=-1, predefined::Bool=false) =
    length(groupnames(m, d=d, predefined=predefined))


"""
    hasgroups(m::Mesh; d::Int=-1, predefined::Bool=false)

Tests groups are defined.
"""
hasgroups(m::Mesh; d::Int=-1, predefined::Bool=false) =
    ngroups(m, d=d, predefined=predefined) > 0


"""
    groups(e::MeshEntity)

Names of groups entity e belongs to sorted by group size. That is, the most specific
group is returned first.
"""
groups(e::MeshEntity{DT}) where {DT} = return sort!(
    [
        group
        for group ∈ groupnames(e.mesh, d=DT, predefined=true)
        if e ∈ e.mesh.groups[group]
    ],
    lt=(n1, n2) -> length(e.mesh.groups[n1]) < length(e.mesh.groups[n2])
)


"""
    group(e, gc::GroupCollection)

Returns the smallest group containing `e``.
"""
group(e::MeshEntity) = groups(e)[1]
