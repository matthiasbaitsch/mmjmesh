module Associations

export MeshData, MeshEntityData

mutable struct MeshData
    entries::Dict{Symbol, Any}

    MeshData() = new(Dict{Symbol, Any}())
end

Base.setindex!(md::MeshData, value, key::Symbol) = md.entries[key] = value
Base.getindex(md::MeshData, key::Symbol) = md.entries[key]

struct MeshEntityData end

end
