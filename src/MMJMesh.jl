module MMJMesh

using Preferences
preference = @load_preference("preference", "default")

# Basic stuff
include("common.jl")
include("mmjbase/00-mmjbase.jl")
include("mathematics/00-mathematics.jl")

# Mesh consists a topology, a geometry and data
include("topologies/00-topologies.jl")
include("geometries/00-geometries.jl")
include("meshes/00-meshes.jl")

# Other things
include("utilities/00-utilities.jl")
include("gmsh/00-gmsh.jl")
include("plots/00-plots.jl")

# Export some
include("reexports.jl")

# Precompile
include("precompile.jl")

end
