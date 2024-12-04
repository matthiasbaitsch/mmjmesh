# -------------------------------------------------------------------------------------------------
# Data
# -------------------------------------------------------------------------------------------------

"""
    Data

`Data` makes it possible to associate data with mesh entities. This can happen on three levels:

1. `setdata!(m, :foo, 99)` sets property `:foo` to `99` for all mesh entities

1. `setdata!(group(m, :bar), :foo, 99)` sets property `:foo` for all mesh entities in group `:bar` to `99`

1. `setdata!(node(m, 1), :foo, 99)` sets property `:foo` for node 1 to `99`

If a property is specified on two or more levels, the most specific level is selected.
"""
mutable struct Data
    mesh
    meshdata::Dict{Symbol,Any}
    groupdata::Dict{Pair{Symbol,Symbol},Any}
    entitydata::Dict{Tuple{Symbol,Int,Int},Any}

    function Data()
        data = new()
        data.meshdata = Dict{Symbol,Any}()
        data.groupdata = Dict{Pair{Symbol,Symbol},Any}()
        data.entitydata = Dict{Tuple{Symbol,Int,Int},Any}()
        return data
    end
end
