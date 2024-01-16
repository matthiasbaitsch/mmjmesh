"""
    Data{T}

Data associates a `name` with a mapping  `(base, parameters...)` to `value`.
"""
mutable struct Data{T}
    base::Union{Nothing,T}
    mappings::Dict{Symbol,Function}
    Data(t::Type) = new{t}(nothing, Dict{Symbol,Any}())
end

addmapping!(md::Data, name::Symbol, m::Function) = md.mappings[name] = m
setbase!(md::Data{T}, b::T) where {T} = md.base = b
Base.setindex!(md::Data, value, name::Symbol) = addmapping!(md, name, (_...) -> value)
Base.getindex(md::Data, name::Symbol, params...) = md.mappings[name](md.base, params...)




