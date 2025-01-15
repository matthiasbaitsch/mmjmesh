# -------------------------------------------------------------------------------------------------
# Names for dimensions
# -------------------------------------------------------------------------------------------------

nnodes(o) = nentities(o, 0)
nedges(o) = nentities(o, 1)
nfaces(o) = nentities(o, 2)
nsolids(o) = nentities(o, 3)
node(o, index::Int) = entity(o, 0, index)
edge(o, index::Int) = entity(o, 1, index)
face(o, index::Int) = entity(o, 2, index)
solid(o, index::Int) = entity(o, 3, index)
nodes(o) = entities(o, 0)
edges(o) = entities(o, 1)
faces(o) = entities(o, 2)
solids(o) = entities(o, 3)
nodeindices(o) = indices(o, 0)
edgeindices(o) = indices(o, 1)
faceindices(o) = indices(o, 2)
solidindices(o) = indices(o, 3)


# -------------------------------------------------------------------------------------------------
# Iterate over all entities
# -------------------------------------------------------------------------------------------------

struct EntityIterator
    m
end

Base.iterate(it::EntityIterator) = node(it.m, 1), (0, 1)

function Base.iterate(it::EntityIterator, state)
    dim, index = state
    index += 1
    if index > nentities(it.m, dim)
        dim += 1
        index = 1
    end
    dim > pdim(it.m) && return nothing
    return entity(it.m, dim, index), (dim, index)
end

entities(m) = EntityIterator(m)
