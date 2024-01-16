"""
    Data{T}

Data associates a `name` with a mapping  `(base, parameters...)` to `value`.
"""
mutable struct Data{T}
    base::Union{Nothing,T}
    mappings::Dict{Symbol,Array{Function}}
    Data(t::Type) = new{t}(nothing, Dict{Symbol,Array{Function}}())
end


# Shortcut to associate everything with a value
Base.setindex!(md::Data, value, name::Symbol) = addmapping!(md, name, (_...) -> value)


function addmapping!(d::Data, name::Symbol, m::Function)
    !haskey(d.mappings, name) && (d.mappings[name] = Function[])
    push!(d.mappings[name], m)
end


function Base.getindex(d::Data, name::Symbol, params...)
    # TODO This is not efficient for many groups associated with the same name
    for m ∈ d.mappings[name]
        r = m(d.base, params...)
        !isnothing(r) && return r
    end
    return nothing
end




