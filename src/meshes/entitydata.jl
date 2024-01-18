struct EntityData{DT}
    mesh::Mesh
    index::Int
end

Base.getindex(ed::EntityData{DT}, name::Symbol) where {DT} = ed.mesh.data[name, DT, ed.index]


"""
    m.data[:foo, :bar] = value

Associate `value` with the name `:foo` for entities in group `:bar`.
"""
function Base.setindex!(d::Data{Mesh{DT,DG}}, value, name::Symbol, group::Symbol) where {DT,DG}
    mesh = d.base
    group = mesh.groups[group]    
    function getvalue(m::Mesh, dt::Int, index::Int)
        if entity(m, dt, index) ∈ group
            return value
        end
        return nothing
    end
    addmapping!(d, name, getvalue)
end

