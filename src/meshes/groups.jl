# -------------------------------------------------------------------------------------------------
# Group
# -------------------------------------------------------------------------------------------------

"""
    struct Group{T} <: AbstractVector{Int}

A set of indices of objects of type `T`. Requires a function `index(o::T) -> Int` in
this module.
"""
struct Group{T} <: AbstractVector{Int}
    indices::SeqIntSet
    Group{T}(a) where {T} = new{T}(SeqIntSet(a))
end

# Delegate methods
Base.length(g::Group) = length(g.indices)
Base.size(g::Group) = size(g.indices)
Base.isempty(g::Group) = isempty(g.indices)
Base.getindex(g::Group, i::Int) = g.indices[i]
Base.eltype(g::Group) = eltype(g.indices)
Base.iterate(g::Group) = iterate(g.indices)
Base.iterate(g::Group, state) = iterate(g.indices, state)

# Set operations
Base.union(g1::Group{T}, g2::Group{T}) where {T} = Group{T}(union(g1.indices, g2.indices))
Base.intersect(g1::Group{T}, g2::Group{T}) where {T} = Group{T}(intersect(g1.indices, g2.indices))
Base.setdiff(g1::Group{T}, g2::Group{T}) where {T} = Group{T}(setdiff(g1.indices, g2.indices))

# Others
Base.in(target::T1, g::Group{T2}) where {T1,T2} = T1 <: T2 && index(target) ∈ g.indices
Base.show(io::IO, g::Group{T}) where {T} = print(io, "$(T)Group$(g.indices)")


# -------------------------------------------------------------------------------------------------
# GroupCollection
# -------------------------------------------------------------------------------------------------

"""
    struct GroupCollection

A collection of `Group`s each identified by a `Symbol`. Predefined groups are created by
associating a `Symbol` with a recipe function. The recipe function is called when the group
is first accessed. The recipe function either returns a `Group` object which is then cached 
or `nothing`` in the case it manipulates the `GroupCollection` directly.
"""
struct GroupCollection
    recipes::Dict{Symbol,Function}
    entries::Dict{Symbol,Union{Nothing,Group}}
    GroupCollection() = new(Dict{Symbol,Function}(), Dict{Symbol,Union{Nothing,Group}}())
end


function Base.getindex(gc::GroupCollection, key::Symbol)
    if isnothing(gc.entries[key])
        r = gc.recipes[key]()
        if !isnothing(r)
            gc.entries[key] = r
        end
    end
    return gc.entries[key]
end

function Base.setindex!(gc::GroupCollection, g::Group, key::Symbol)
    @assert !ispredefined(gc, key) "Groupname '$key' is predefined"
    gc.entries[key] = g
end

function Base.show(io::IO, gc::GroupCollection)
    println(io, "GroupCollection")
    for key in sort(keys(gc.entries) |> collect)
        g = gc.entries[key]
        println(io, "  $key: $g")
    end
end

names(gc::GroupCollection) = keys(gc.entries)
ispredefined(gc::GroupCollection, key::Symbol) = haskey(gc.recipes, key)

function addrecipe!(gc::GroupCollection, key::Symbol, recipe::Function)
    gc.recipes[key] = recipe
    gc.entries[key] = nothing
end

function clearcache!(gc::GroupCollection)
    for key in keys(gc.recipes)
        gc.entries[key] = nothing
    end
end
