module MMJMesh

# Basic stuff
include("mmjbase/00-mmjbase.jl")

# Mesh consisits of a topology and a geometry
include("topologies/00-topologies.jl")
include("geometries/00-geometries.jl")
include("mathematics/00-mathematics.jl")
include("associations/00-associations.jl")
include("meshes/00-meshes.jl")

# Other things
include("gmsh/00-gmsh.jl")
include("plots/00-plots.jl")
include("utilities/00-utilities.jl")

end
