import MMJMesh.Groups: addrecipe!, groupnames

module Detail
using MMJMesh.Meshes


function collectboundary(m::Mesh{1,DG}) where DG
    bn = Vector{Int}()
    for e ∈ nodes(m)
        if nedges(e) == 1
            append!(bn, nodeIdxs(e))
        end
    end
    return bn, Int[], Int[]
end


function collectboundary(m::Mesh{2,DG}) where DG
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

end


function populatepredfinedgroups!(m::Mesh)
    function recipe()
        bn, be, bf = Detail.collectboundary(m)
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


