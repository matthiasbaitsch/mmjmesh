module MMJMesh

# Basic stuff
include("mmjbase/00-mmjbase.jl")

# Mesh consisits of a topology and a geometry
include("topologies/00-topologies.jl")
include("geometries/00-geometries.jl")
include("meshes/00-meshes.jl")

# Other things
include("plots/00-plots.jl")
include("utilities/00-utilities.jl")

end
