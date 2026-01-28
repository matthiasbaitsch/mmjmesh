# -------------------------------------------------------------------------------------------------
# Names for dimensions
# -------------------------------------------------------------------------------------------------

nnodes(o) = nentities(o, 0)
nedges(o) = nentities(o, 1)
nfaces(o) = nentities(o, 2)
nsolids(o) = nentities(o, 3)

node(o, index) = entity(o, 0, index)
edge(o, index) = entity(o, 1, index)
face(o, index) = entity(o, 2, index)
solid(o, index) = entity(o, 3, index)

nodes(o, args...; kwargs...) = entities(o, 0, args...; kwargs...)
edges(o, args...; kwargs...) = entities(o, 1, args...; kwargs...)
faces(o, args...; kwargs...) = entities(o, 2, args...; kwargs...)
solids(o, args...; kwargs...) = entities(o, 3, args...; kwargs...)

nodeindices(o, args...; kwargs...) = indices(o, 0, args...; kwargs...)
edgeindices(o, args...; kwargs...) = indices(o, 1, args...; kwargs...)
faceindices(o, args...; kwargs...) = indices(o, 2, args...; kwargs...)
solidindices(o, args...; kwargs...) = indices(o, 3, args...; kwargs...)


# -------------------------------------------------------------------------------------------------
# Indices
# -------------------------------------------------------------------------------------------------

function indices(o::Union{MeshEntityGroup,MeshEntityList}, pdim::Int)
    if edim(o) == pdim
        return indices(o)
    else
        return Set(vcat(indices.(entities(o), pdim)...)) |> collect |> sort
    end
end


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

entities(m::Union{Mesh,MeshEntity}) = EntityIterator(m)