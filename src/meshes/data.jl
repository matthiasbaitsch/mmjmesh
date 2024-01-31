# -------------------------------------------------------------------------------------------------
# MeshData
# -------------------------------------------------------------------------------------------------

mutable struct MeshData{T}
    mesh::T
    mappings::Dict{Symbol,Any}
    function MeshData{T}() where {T}
        md = new{T}()
        md.mappings = Dict{Symbol,Any}()
        return md
    end
end

"""
    m.data[:foo] = value

Associate `value` with the name `:foo` for everything in the mesh `m`.
"""
Base.setindex!(d::MeshData, value::Any, name::Symbol) = d.mappings[name] = (_...) -> value

"""
    v = m.data[:foo, e]

Get the value associated with the name `:foo` for entity `e` in the mesh `m`.
"""
Base.getindex(d::MeshData, name::Symbol, params...) = d.mappings[name](params...)


# -------------------------------------------------------------------------------------------------
# EntityData
# -------------------------------------------------------------------------------------------------

mutable struct EntityData{T}
    entity::T
    EntityData{T}() where {T} = new{T}()
end

