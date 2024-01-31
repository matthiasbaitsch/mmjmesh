# -------------------------------------------------------------------------------------------------
# EntityData
# -------------------------------------------------------------------------------------------------

Base.getindex(ed::EntityData, name::Symbol) = ed.entity.mesh.data[name, ed.entity]


# -------------------------------------------------------------------------------------------------
# Data associated with entities
# -------------------------------------------------------------------------------------------------

struct _EntityMapping{T}
    mappings::Dict{Int,Vector{Union{T,Nothing}}}
    _EntityMapping{T}() where {T} = new{T}(Dict{Int,Vector{Union{T,Nothing}}}())
end

# Make _EntityMapping callable
function (gd::_EntityMapping{T})(e::MeshEntity{DT}) where {DT,T}
    return gd.mappings[DT][e.index]
end

# Add
function _addentitymapping!(mesh::Mesh, name::Symbol, dt::Int, t::Type)
    if !haskey(mesh.data.mappings, name)
        mesh.data.mappings[name] = _EntityMapping{t}()
    end
    em = mesh.data.mappings[name]
    if !haskey(em.mappings, dt)
        em.mappings[dt] = Vector{Union{t,Nothing}}(undef, nentities(mesh, dt))
        em.mappings[dt] .= nothing
    end
end

"""
    e.data[:foo] = value

Associate `value` with the name `:foo` for entity `e`.
"""
function Base.setindex!(ed::EntityData{MeshEntity{DT}}, value::T, name::Symbol) where {DT,T}
    _addentitymapping!(ed.entity.mesh, name, DT, T)
    ed.entity.mesh.data.mappings[name].mappings[DT][ed.entity.index] = value
end

function Base.setindex!(ed::EntityData{MeshEntity{DT}}, value::Function, name::Symbol) where {DT}
    _addentitymapping!(ed.entity.mesh, name, DT, Function)
    ed.entity.mesh.data.mappings[name].mappings[DT][ed.entity.index] = value
end

function Base.setindex!(d::MeshData, values::Vector{T}, name::Symbol, dt::Int) where {T}
    _addentitymapping!(d.mesh, name, dt, T)
    d.mappings[name].mappings[dt] .= values
end

# -------------------------------------------------------------------------------------------------
# Data associated with groups
# -------------------------------------------------------------------------------------------------

struct _GroupMapping{T}
    values::Dict{Symbol,T}
    _GroupMapping{T}() where {T} = new{T}(Dict{Symbol,T}())
end

# Make _GroupMapping callable
function (gd::_GroupMapping)(e::MeshEntity)
    for g ∈ groups(e)
        haskey(gd.values, g) && return gd.values[g]
    end
end

# Add
_addgroupmapping!(d::MeshData, name::Symbol, t::Type) =
    !haskey(d.mappings, name) && (d.mappings[name] = _GroupMapping{t}())

"""
    m.data[:foo, :bar] = value

Associate `value` with the name `:foo` for entities in group `:bar`.
"""
function Base.setindex!(d::MeshData, value::T, name::Symbol, group::Symbol) where {T}
    _addgroupmapping!(d, name, T)
    d.mappings[name].values[group] = value
end

function Base.setindex!(d::MeshData, f::Function, name::Symbol, group::Symbol)
    _addgroupmapping!(d, name, Function)
    d.mappings[name].values[group] = f
end
