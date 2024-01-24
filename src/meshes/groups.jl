# -------------------------------------------------------------------------------------------------
# EntityGroup
# -------------------------------------------------------------------------------------------------

struct EntityGroup{T} <: AbstractVector{Int}
    indexes::SeqIntSet
    EntityGroup{T}(a::AbstractVector{Int}) where {T} = new{T}(SeqIntSet(a))
end


# -------------------------------------------------------------------------------------------------
# Delegate methods
# -------------------------------------------------------------------------------------------------

# AbstractArray
Base.length(g::EntityGroup) = length(g.indexes)
Base.size(g::EntityGroup) = size(g.indexes)
Base.isempty(g::EntityGroup) = isempty(g.indexes)
Base.getindex(g::EntityGroup, i::Int) = g.indexes[i]

# Iterator
Base.eltype(g::EntityGroup) = eltype(g.indexes)
Base.iterate(g::EntityGroup) = iterate(g.indexes)
Base.iterate(g::EntityGroup, state) = iterate(g.indexes, state)

# Set operations
Base.union(g1::EntityGroup{T}, g2::EntityGroup{T}) where {T} =
    EntityGroup{T}(union(g1.indexes, g2.indexes))
Base.intersect(g1::EntityGroup{T}, g2::EntityGroup{T}) where {T} =
    EntityGroup{T}(intersect(g1.indexes, g2.indexes))
Base.setdiff(g1::EntityGroup{T}, g2::EntityGroup{T}) where {T} =
    EntityGroup{T}(setdiff(g1.indexes, g2.indexes))

# In 
Base.in(target::T1, g::EntityGroup{T2}) where {T1,T2} = T1 <: T2 && in(target.index, g.indexes)

# Show
Base.show(io::IO, g::EntityGroup{T}) where {T} = print(io, "$(T)Group$(g.indexes)")


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
        g = gc[name]
        if (predefined || !ispredefined(gc, name)) && (d == -1 || (!isempty(g) && edim(g) == d))
            push!(names, name)
        end
    end
    return names
end

ispredefined(gc::EntityGroupCollection, key::Symbol) = haskey(gc.recipes, key)
ngroups(gc::EntityGroupCollection; d::Int=-1, predefined::Bool=false) = length(groupnames(gc, d=d, predefined=predefined))
hasgroups(gc::EntityGroupCollection, d::Int; predefined::Bool=false) = !isempty(groupnames(gc, d=d, predefined=predefined))

"""
    groupsof(e, gc::EntityGroupCollection)

Names of groups entity e belongs to sorted by group size. That is, the most specific
group is returned first.
"""
groupsof(e, gc::EntityGroupCollection) = return sort!(
    [name for name ∈ groupnames(gc, d=pdim(e), predefined=true) if e ∈ gc[name]], 
    lt = (n1, n2) -> length(gc[n1]) < length(gc[n2])
)

groupof(e, gc::EntityGroupCollection) = groupsof(e, gc)[1]