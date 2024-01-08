module Utilities

# Modules needed by this module
using MMJMesh.Meshes
using MMJMesh.Topologies
using MMJMesh.Geometries

# Exports
## generatemeshes.jl
export Meshtype, QUADRANGLE, TRIANGLE
export makemeshoninterval, makemeshonrectangle

# Parts
include("generatemeshes.jl")

end
