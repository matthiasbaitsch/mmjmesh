# -------------------------------------------------------------------------------------------------
# EntityData struct
# -------------------------------------------------------------------------------------------------

mutable struct EntityData{T}
    entity::Union{T,Nothing}
    EntityData{T}() where {T} = new{T}(nothing)
end

Base.getindex(ed::EntityData, name::Symbol) = ed.entity.mesh.data[name, ed.entity]


# -------------------------------------------------------------------------------------------------
# Data associated with groups of entities
# -------------------------------------------------------------------------------------------------

struct GroupData{T}
    values::Dict{Symbol,T}
    GroupData{T}() where {T} = new{T}(Dict{Symbol,T}())
end

function (gd::GroupData{T})(e::ET) where {T,ET}
    for g ∈ groups(e)
        haskey(gd.values, g) && return gd.values[g]
    end
end

"""
    m.data[:foo, :bar] = value

Associate `value` with the name `:foo` for entities in group `:bar`.
"""
function Base.setindex!(d::Data, value::T, name::Symbol, group::Symbol) where {T}
    if !haskey(d.mappings, name)
        d.mappings[name] = GroupData{T}()
    end
    d.mappings[name].values[group] = value
end
