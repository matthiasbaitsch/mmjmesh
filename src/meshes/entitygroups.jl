const NodeGroup = EntityGroup{Node}
const EdgeGroup = EntityGroup{Edge}
const FaceGroup = EntityGroup{Face}
const SolidGroup = EntityGroup{Solid}


# To my understanding, that should work like this
# edim(::EntityGroup{MeshEntity{DT}}) where {DT} = DT
edim(::NodeGroup) = 0
edim(::EdgeGroup) = 1
edim(::FaceGroup) = 2
edim(::SolidGroup) = 3


function _collectboundary(m::Mesh{1})
    bn = Vector{Int}()
    for e ∈ nodes(m)
        if nedges(e) == 1
            append!(bn, nodeIdxs(e))
        end
    end
    return bn, Int[], Int[]
end


function _collectboundary(m::Mesh{2})
    bn = Vector{Int}()
    be = Vector{Int}()
    for e ∈ edges(m)
        if nfaces(e) == 1
            push!(be, e.index)
            append!(bn, nodeIdxs(e))
        end
    end
    return bn, be, Int[]
end


function _idvector(s::AbstractVector)
    d = Dict([(j, i) for (i, j) in enumerate(Set(s))])
    return [d[k] for k in s]
end


function populatepredfinedgroups!(m::Mesh)
    function recipe()
        bn, be, bf = _collectboundary(m)
        m.groups[:boundarynodes] = NodeGroup(bn)
        m.groups[:boundaryedges] = EdgeGroup(be)
        m.groups[:boundaryfaces] = FaceGroup(bf)
    end
    addrecipe!(m.groups, :boundarynodes, recipe)
    addrecipe!(m.groups, :boundaryedges, recipe)
    addrecipe!(m.groups, :boundaryfaces, recipe)
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


groupids(m::Mesh; d::Int=-1, predefined::Bool=false) = _idvector(collectgroups(m, d=d, predefined=predefined))


