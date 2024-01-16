module Groups

import MMJMesh.MMJBase: SeqIntSet

export NodeGroup, EdgeGroup, FaceGroup, SolidGroup, hasgroups, groupnames, ngroups

# -------------------------------------------------------------------------------------------------
# EntityGroup
# -------------------------------------------------------------------------------------------------

struct EntityGroup{DT} <: AbstractVector{Int}
    indexes::SeqIntSet
    EntityGroup{DT}(a::AbstractVector{Int}) where {DT} = new{DT}(SeqIntSet(a))
end

EntityGroup(dt::Int, a::AbstractVector{Int}) = EntityGroup{dt}(SeqIntSet(a))

# Short names
const NodeGroup = EntityGroup{0}
const EdgeGroup = EntityGroup{1}
const FaceGroup = EntityGroup{2}
const SolidGroup = EntityGroup{3}

# Functions from Base
Base.length(g::EntityGroup) = length(g.indexes)
Base.size(g::EntityGroup) = size(g.indexes)
Base.isempty(g::EntityGroup) = isempty(g.indexes)
Base.in(target::Int, g::EntityGroup) = in(target, g.indexes)
Base.getindex(g::EntityGroup, i::Int) = g.indexes[i]
Base.show(io::IO, g::EntityGroup{DT}) where {DT} = (print(io, "EntityGroup{$DT}"); show(io, g.indexes))
Base.eltype(g::EntityGroup) = eltype(g.indexes)
Base.iterate(g::EntityGroup) = iterate(g.indexes)
Base.iterate(g::EntityGroup, state) = iterate(g.indexes, state)
Base.union(g1::EntityGroup, g2::EntityGroup) = EntityGroup{dim(g1)}(union(g1.indexes, g2.indexes))
Base.intersect(g1::EntityGroup, g2::EntityGroup) = EntityGroup{dim(g1)}(intersect(g1.indexes, g2.indexes))
Base.setdiff(g1::EntityGroup, g2::EntityGroup) = EntityGroup{dim(g1)}(setdiff(g1.indexes, g2.indexes))

# Own functions
dim(::EntityGroup{DT}) where {DT} = DT

# -------------------------------------------------------------------------------------------------
# EntityGroupCollection
# -------------------------------------------------------------------------------------------------

struct EntityGroupCollection
    recipes::Dict{Symbol,Function}
    entries::Dict{Symbol,Union{Nothing,EntityGroup}}
end

EntityGroupCollection() = EntityGroupCollection(
    Dict{Symbol,Function}(),
    Dict{Symbol,Union{Nothing,EntityGroup}}()
)

# Functions from Base

function Base.getindex(gc::EntityGroupCollection, key::Symbol)
    if isnothing(gc.entries[key])
        gc.recipes[key]()
    end
    return gc.entries[key]
end

function Base.setindex!(gc::EntityGroupCollection, g::EntityGroup, key::Symbol)
    gc.entries[key] = g
end

function Base.show(io::IO, gc::EntityGroupCollection)
    println(io, "EntityGroupCollection")
    for key in sort(keys(gc.entries) |> collect)
        g = gc.entries[key]
        println(io, "  $key: $g")
    end 
end

# Own functions

function addrecipe!(gc::EntityGroupCollection, key::Symbol, recipe::Function)
    gc.recipes[key] = recipe
    gc.entries[key] = nothing
end

function groupnames(gc::EntityGroupCollection; d::Int=-1, predefined::Bool=false)
    names = Vector{Symbol}()
    for name in keys(gc.entries)
        if (predefined || !ispredefined(gc, name)) && (d == -1 || dim(gc[name]) == d)
            push!(names, name)
        end
    end
    return names
end

ngroups(gc::EntityGroupCollection; d::Int=-1, predefined::Bool=false) = length(groupnames(gc, d=d, predefined=predefined))

ispredefined(gc::EntityGroupCollection, key::Symbol) = haskey(gc.recipes, key)

hasgroups(gc::EntityGroupCollection, d::Int; predefined::Bool=false) = !isempty(groupnames(gc, d=d, predefined=predefined))


end


