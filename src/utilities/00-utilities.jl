module Utilities

# Exports

## generatemeshes.jl
export Meshtype, QUADRANGLE, TRIANGLE
export makemeshoninterval
export makemeshonrectangle

# Modules needed by this module
using MMJMesh.Meshes
using MMJMesh.Topologies
using MMJMesh.Geometries

# Parts
include("generatemeshes.jl")

end
