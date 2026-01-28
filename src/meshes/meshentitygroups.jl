# -------------------------------------------------------------------------------------------------
# Type definitions
# -------------------------------------------------------------------------------------------------

"""
Group of mesh entities.
"""
const MeshEntityGroup = Group{<:MeshEntity}

"""
    edim(group)

Parametric dimension of entities in `group`.
"""
edim(::Group{<:MeshEntity{DT}}) where DT = DT


# -------------------------------------------------------------------------------------------------
# Interface
# -------------------------------------------------------------------------------------------------

"""
    name(g)

Name of group `g`.
"""
function name(g::Group)
    for n = groupnames(g.mesh, predefined=true)
        g.mesh.groups.entries[n] === g && return n
    end
end

"""
    definegroup!(name, entities)
    definegroup!(name, m, pdim, indices)

Define group `name`.
"""
definegroup!(name::Symbol, es::MeshEntityList{DT}) where {DT} =
    definegroup!(name, mesh(es), DT, indices(es))

function definegroup!(name::Symbol, m::Mesh, dim::Int, indices::IntegerVec)
    g = MeshEntityGroup{MeshEntity{dim}}(indices)
    g.mesh = m
    m.groups[name] = g
end

"""
    groupnames(m; d=-1, predefined=false)

Returns names of all non-empty groups in Mesh `m` of dimension `d`. If `d` is not specified,
groups of all dimensions are included. If `predefined` is `false`, only 
user defined groups are considered.
"""
function groupnames(m::Mesh; d::Int=-1, predefined::Bool=false)
    isvalid(name) =
        (predefined || !ispredefined(m.groups, name)) &&
        (d == -1 || (!isempty(m.groups[name]) && edim(m.groups[name]) == d))
    return sort([name for name ∈ names(m.groups) if isvalid(name)])
end

"""
    groupnames(e)

Names of groups entity e belongs to sorted by group size. That is, the most specific
group is returned first.
"""
groupnames(e::MeshEntity{DT}) where {DT} = return sort!(
    [
        groupname
        for groupname ∈ groupnames(e.mesh, d=DT, predefined=true)
        if e ∈ e.mesh.groups[groupname]
    ],
    lt=(name1, name2) -> length(e.mesh.groups[name1]) < length(e.mesh.groups[name2])
)


"""
    groupname(e)

Returns the smallest group containing `e``.
"""
groupname(e::MeshEntity) = groupnames(e) |> first

"""
    groupids(m::Mesh; d::Int, predefined::Bool=false)

Unique numerical ids groups of specified dimension.
"""
groupids(m::Mesh; d::Int, predefined::Bool=false) =
    _idvector(_entitygroupnames(m, d=d, predefined=predefined))


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


function group(m::Mesh, name::Symbol)
    # Hack to set group for recipes
    g = m.groups[name]
    g.mesh = m
    return g
end

# Entities from mesh via group
entities(m::Mesh, pdim::Integer, groupname::Symbol; select=all) =
    entities(m, pdim, indices(m, pdim, groupname; select=select))
entities(m::Mesh, groupname::Symbol) = entities(group(m, groupname))

# Entities from group
entities(g::MeshEntityGroup) = entities(g.mesh, edim(g), indices(g))
entities(g::MeshEntityGroup, pdim::Integer) = entities(entities(g), pdim)


# -------------------------------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------------------------------


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
    _entitygroupnames(m::Mesh; d::Int, predefined::Bool=false)

Returns a vector of the groups each entity of dimension `d` belongs to.
"""
function _entitygroupnames(m::Mesh; d::Int, predefined::Bool=false)
    groups = [Symbol[] for _ in 1:nentities(m, d)]
    for name ∈ groupnames(m, d=d, predefined=predefined)
        for index ∈ m.groups[name]
            push!(groups[index], name)
        end
    end
    sort!.(groups)
    return groups
end


function _collectboundary(m::Mesh{1})
    bn = Int[]
    for e ∈ nodes(m)
        if nedges(e) == 1
            append!(bn, indices(e, 0))
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
            append!(bn, indices(e, 0))
        end
    end
    return bn, be, Int[]
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
        m.groups.entries[:boundarynodes] = MeshEntityGroup{Node}(bn)
        m.groups.entries[:boundaryedges] = MeshEntityGroup{Edge}(be)
        m.groups.entries[:boundaryfaces] = MeshEntityGroup{Face}(bf)
        return nothing
    end

    addrecipe!(m.groups, :boundarynodes, boundaryrecipe)
    addrecipe!(m.groups, :boundaryedges, boundaryrecipe)
    addrecipe!(m.groups, :boundaryfaces, boundaryrecipe)

    # By dimension
    addrecipe!(m.groups, :nodes, () -> MeshEntityGroup{Node}(1:nnodes(m)))
    addrecipe!(m.groups, :edges, () -> MeshEntityGroup{Edge}(1:nedges(m)))
    addrecipe!(m.groups, :faces, () -> MeshEntityGroup{Face}(1:nfaces(m)))
    addrecipe!(m.groups, :solids, () -> MeshEntityGroup{Solid}(1:nsolids(m)))

    # Elements - TODO elegant solution
    if pdim(m) == 1
        addrecipe!(m.groups, :elements, () -> MeshEntityGroup{Edge}(1:nedges(m)))
    elseif pdim(m) == 2
        addrecipe!(m.groups, :elements, () -> MeshEntityGroup{Face}(1:nfaces(m)))
    elseif pdim(m) == 3
        addrecipe!(m.groups, :elements, () -> MeshEntityGroup{Solid}(1:nsolids(m)))
    end
end


