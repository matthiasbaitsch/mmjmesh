module Associations

export MeshData, addmapping!

mutable struct MeshData{T}
    source::Union{Nothing, T}
    mappings::Dict{Symbol, Function}
    MeshData(t::Type) = new{t}(nothing, Dict{Symbol, Any}())
end

addmapping!(md::MeshData, name::Symbol, m::Function) = md.mappings[name] = m

Base.setindex!(md::MeshData, value, name::Symbol) = addmapping!(md, name, (p...) -> value)

Base.getindex(md::MeshData, name::Symbol, params...) = md.mappings[name](md.source, params...)

end