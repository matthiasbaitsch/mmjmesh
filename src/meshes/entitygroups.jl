
# Short names
const NodeGroup = Group{Node}
const EdgeGroup = Group{Edge}
const FaceGroup = Group{Face}
const SolidGroup = Group{Solid}


# To my understanding, that should work like this
# edim(::Group{MeshEntity{DT}}) where {DT} = DT
# Why not?
edim(::NodeGroup) = 0
edim(::EdgeGroup) = 1
edim(::FaceGroup) = 2
edim(::SolidGroup) = 3


# Helper functions

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


function _idvector(s::AbstractVector)
    d = Dict([(j, i) for (i, j) in enumerate(Set(s))])
    return [d[k] for k in s]
end


# Mesh entity group functions

function populatepredfinedgroups!(m::Mesh)

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
end


function collectgroups(m::Mesh; d::Int=-1, predefined::Bool=false)
    namelists = [Symbol[] for _ in 1:nentities(m, d)]
    for name ∈ groupnames(m.groups, d=d, predefined=predefined)
        for index ∈ m.groups[name]
            push!(namelists[index], name)
        end
    end
    return namelists
end


"""
    groupnames(gc::GroupCollection; d::Int=-1, predefined::Bool=false)
    ngroups(gc::GroupCollection; d::Int=-1, predefined::Bool=false)
    hasgroups(gc::GroupCollection; d::Int=-1, predefined::Bool=false)
    groupids(m::Mesh; d::Int=-1, predefined::Bool=false)

Returns
    - names of groups
    - number of groups
    - whether there are groups
    - IDs for groups entities belong to

in the `GroupCollection` or `Mesh` `m` with dimension `d`. If `d` is not specified,
groups of all dimensions are included. If `predefined` is `false`, only user-defined
groups are considered.
"""

function groupnames(gc::GroupCollection; d::Int=-1, predefined::Bool=false)
    isvalid(name) =
        (predefined || !ispredefined(gc, name)) &&
        (d == -1 || (!isempty(gc[name]) && edim(gc[name]) == d))
    return [name for name ∈ names(gc) if isvalid(name)]
end

ngroups(gc::GroupCollection; d::Int=-1, predefined::Bool=false) =
    length(groupnames(gc, d=d, predefined=predefined))

hasgroups(gc::GroupCollection; d::Int=-1, predefined::Bool=false) =
    !isempty(groupnames(gc, d=d, predefined=predefined))

groupids(m::Mesh; d::Int=-1, predefined::Bool=false) =
    _idvector(collectgroups(m, d=d, predefined=predefined))


"""
    groupsof(e, gc::GroupCollection)

Names of groups entity e belongs to sorted by group size. That is, the most specific
group is returned first.
"""
groupsof(e::MeshEntity{DT}, gc::GroupCollection) where {DT} = return sort!(
    [name for name ∈ groupnames(gc, d=DT, predefined=true) if e ∈ gc[name]],
    lt=(n1, n2) -> length(gc[n1]) < length(gc[n2])
)

"""
    groupof(e, gc::GroupCollection)

Returns the smallest group containing e.
"""
groupof(e, gc::GroupCollection) = groupsof(e, gc)[1]
