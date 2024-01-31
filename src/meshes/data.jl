"""
    Data{T}

Data associates a `name` with a mapping  `(base, params...) -> Union{Nothing, Any}`.
"""
mutable struct Data
    mappings::Dict{Symbol,Any}
    Data() = new(Dict{Symbol,Any}())
end

Base.setindex!(d::Data, value::Any, name::Symbol) = d.mappings[name] = (_...) -> value
Base.getindex(d::Data, name::Symbol, params...) = d.mappings[name](params...)
