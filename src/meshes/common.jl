nnodes(o) = nentities(o, 0)
nedges(o) = nentities(o, 1)
nfaces(o) = nentities(o, 2)
nsolids(o) = nentities(o, 3)
nodeIdxs(o) = indexes(o, 0)
edgeIdxs(o) = indexes(o, 1)
faceIdxs(o) = indexes(o, 2)
solidIdxs(o) = indexes(o, 3)
node(o, index::Int) = entity(o, 0, index)
edge(o, index::Int) = entity(o, 1, index)
face(o, index::Int) = entity(o, 2, index)
solid(o, index::Int) = entity(o, 3, index)
nodes(o) = entities(o, 0)
edges(o) = entities(o, 1)
faces(o) = entities(o, 2)
solids(o) = entities(o, 3)

# TODO: Use recipe for lazy creation, move to better place
function populatepredfinedgroups!(m::Mesh)
    # TODO 3D: Generalize
    be = Vector{Int}()
    bn = Vector{Int}()

    # TODO: Should be like this 
    # for e ∈ edges(m)
    # if nfaces(e) == 1
    # push!(be, e.index)

    l12 = links(m.topology, 1, 2)
    for (i, l) in enumerate(links(m.topology, 1, 0))
        if length(l12, i) == 1
            push!(be, i)
            append!(bn, l)
        end
    end

    # TODO
    m.groups.recipes[:boundarynodes] = sin
    m.groups[:boundarynodes] = EntityGroup(0, bn)
    m.groups.recipes[:boundaryedges] = sin
    m.groups[:boundaryedges] = EntityGroup(1, be)
end
