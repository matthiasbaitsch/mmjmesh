module Groups

import MMJMesh.MMJBase: SeqIntSet

struct EntityGroup{DT} <: AbstractVector{Int}
    indexes::SeqIntSet
    EntityGroup{DT}(a::AbstractVector{Int}) where {DT} = new{DT}(SeqIntSet(a))
end
EntityGroup(dt::Int, a::AbstractVector{Int}) = EntityGroup{dt}(SeqIntSet(a))

Base.length(g::EntityGroup) = length(g.indexes)
Base.size(g::EntityGroup) = size(g.indexes)
Base.isempty(g::EntityGroup) = isempty(g.indexes)
Base.in(target::Int, g::EntityGroup) = in(target, g)
Base.getindex(g::EntityGroup, i::Int) = g.indexes[i]
Base.show(io::IO, g::EntityGroup{DT}) where {DT} = (print(io, "EntityGroup{$DT}"); show(io, g.indexes))
Base.eltype(::EntityGroup) = Int
Base.iterate(g::EntityGroup) = iterate(g.indexes)
Base.iterate(g::EntityGroup, state) = iterate(g.indexes, state)

const NodeGroup = EntityGroup{0}
const EdgeGroup = EntityGroup{1}
const FaceGroup = EntityGroup{2}
const SlolidGroup = EntityGroup{3}

struct EntityGroupCollection
    recipes::Dict{Symbol,Function}
    groups::Dict{Symbol,EntityGroup}
end
EntityGroupCollection() = EntityGroupCollection(Dict{Symbol,Function}(), Dict{Symbol,EntityGroup}())

ispredefined(gc::EntityGroupCollection, key::Symbol) = haskey(gc.recipes, key)
Base.getindex(gc::EntityGroupCollection, key::Symbol) = gc.groups[key]
Base.setindex!(gc::EntityGroupCollection, g::EntityGroup, key::Symbol) = gc.groups[key] = g

end


