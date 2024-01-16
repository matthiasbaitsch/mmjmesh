module MMJMesh

# Basic stuff
include("mmjbase/00-mmjbase.jl")

# Mesh consists a topology, a geometry and data
include("topologies/00-topologies.jl")
include("geometries/00-geometries.jl")
include("mathematics/00-mathematics.jl")
include("groups/00-groups.jl")
include("meshes/00-meshes.jl")

# Other things
include("utilities/00-utilities.jl")
include("gmsh/00-gmsh.jl")
include("plots/00-plots.jl")

end
