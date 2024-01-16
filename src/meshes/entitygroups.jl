Base.in(target::MeshEntity, g::EntityGroup) = (pdim(target) == dimension(g) && target.index ∈ g.indexes)


function _collectboundary(m::Mesh{1,DG}) where {DG}
    bn = Vector{Int}()
    for e ∈ nodes(m)
        if nedges(e) == 1
            append!(bn, nodeIdxs(e))
        end
    end
    return bn, Int[], Int[]
end


function _collectboundary(m::Mesh{2,DG}) where {DG}
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
        m.groups[:boundarynodes] = EntityGroup(0, bn)
        m.groups[:boundaryedges] = EntityGroup(1, be)
        m.groups[:boundaryfaces] = EntityGroup(2, bf)
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


