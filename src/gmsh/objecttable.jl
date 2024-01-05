# -------------------------------------------------------------------------------------------------
# Pretty-print objects in an a table
# -------------------------------------------------------------------------------------------------

# Struct
struct ObjectTable{T}
    ids::Vector{Int}
    objects::Array{T}
end

# Consructors
function ObjectTable(a::Dict{Int, T}) where T 
    skeys = sort!([k for k ∈ keys(a)])
    ObjectTable(skeys, [a[k] for k in skeys])
end
ObjectTable(a::Array) = ObjectTable(collect(1:length(a)), a)

# Access
const ID_NAME = Symbol("#")
names(::ObjectTable{T}) where {T} = [ID_NAME; collect(fieldnames(T))]
indexof(ta::ObjectTable, nm::Symbol) = findfirst(x -> x == nm, names(ta))
Base.getindex(ta::ObjectTable, nm::Symbol, i::Int) = nm == ID_NAME ? ta.ids[i] : getfield(ta.objects[i], nm)

# Tables interface
Tables.istable(::Type{<:ObjectTable}) = true
Tables.columnaccess(::Type{<:ObjectTable}) = true
Tables.columns(ta::ObjectTable) = ta
Tables.columnnames(ta::ObjectTable) = names(ta)
#Tables.getcolumn(ta::ObjectTable, nm::Symbol) = Tables.getcolumn(ta, indexof(ta, nm))
Tables.getcolumn(ta::ObjectTable{T}, nm::Symbol) where {T} = [ta[nm, i] for i ∈ 1:length(ta.objects)]