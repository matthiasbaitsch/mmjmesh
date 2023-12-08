module MMJMesh

# Toplevel functions to be exported
"""
    entity(o, pdim::Int, idx::Int)

TODO
"""
function entity end

"""
    coordinates(o [, index::int])

TODO
"""
function coordinates end

export entity, coordinates

# Basic stuff
include("mmjbase/00-mmjbase.jl")

# Mesh consisits of a topology and a geometry
include("topologies/00-topologies.jl")
include("geometries/00-geometries.jl")
include("mathematics/00-mathematics.jl")
include("meshes/00-meshes.jl")

# Other things
include("plots/00-plots.jl")
include("utilities/00-utilities.jl")

end
