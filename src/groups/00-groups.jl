module Groups

import MMJMesh.MMJBase: SeqIntSet

struct EntityGroup{DT} <: AbstractVector{Int}
    indexes::SeqIntSet
end
EntityGroup(dt::Int, a::AbstractVector{Int}) = EntityGroup{dt}(SeqIntSet(a))

Base.length(g::EntityGroup) = length(g.indexes)
Base.size(g::EntityGroup) = size(g.indexes)
Base.isempty(g::EntityGroup) = isempty(g.indexes)
Base.in(target::Int, g::EntityGroup) = in(target, g)
Base.getindex(g::EntityGroup, i::Int) = g.indexes[i]
Base.show(io::IO, g::EntityGroup{DT}) where DT = (show(io, "EntityGroup{$DT}"); show(io, g.indexes))
Base.eltype(::EntityGroup) = Int
Base.iterate(g::EntityGroup) = iterate(g.indexes)
Base.iterate(g::EntityGroup, state) = iterate(g.indexes, state)

const NodeGroup = EntityGroup{0}
const EdgeGroup = EntityGroup{1}
const FaceGroup = EntityGroup{2}
const SlolidGroup = EntityGroup{3}

struct EntityGroupCollection
    groups::Dict{Symbol, EntityGroup}
end
EntityGroupCollection() = EntityGroupCollection(Dict{Symbol, EntityGroup}())

Base.getindex(gc::EntityGroupCollection, key::Symbol) = gc.groups[key]
Base.setindex!(gc::EntityGroupCollection, g::EntityGroup, key::Symbol) = gc.groups[key] = g

end


